#include <stdio.h>
#include <queue>
#include <map>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>

#include "estimator.h"
#include "parameters.h"
#include "utility/visualization.h"


Estimator estimator;

std::condition_variable con;
double current_time = -1;
queue<sensor_msgs::ImuConstPtr> imu_buf;
queue<sensor_msgs::PointCloudConstPtr> feature_buf;
queue<sensor_msgs::PointCloudConstPtr> relo_buf;
int sum_of_wait = 0;

std::mutex m_buf;
std::mutex m_state;
std::mutex i_buf;
std::mutex m_estimator;

double latest_time;
Eigen::Vector3d tmp_P;
Eigen::Quaterniond tmp_Q;
Eigen::Vector3d tmp_V;
Eigen::Vector3d tmp_Ba;
Eigen::Vector3d tmp_Bg;
Eigen::Vector3d acc_0;
Eigen::Vector3d gyr_0;
bool init_feature = 0;
bool init_imu = 1;
double last_imu_t = 0;

/**
 * @brief 根据当前imu数据预测当前位姿
 * 
 * @param[in] imu_msg 
 */
// xc's todo: 这个函数是否仅仅用于可视化的逻辑，而并不会对算法的核心逻辑产生影响？
void predict(const sensor_msgs::ImuConstPtr &imu_msg)
{
    double t = imu_msg->header.stamp.toSec();
    if (init_imu)
    {
        latest_time = t;
        init_imu = 0;
        return;
    }
    double dt = t - latest_time;
    latest_time = t;

    // 得到加速度
    double dx = imu_msg->linear_acceleration.x;
    double dy = imu_msg->linear_acceleration.y;
    double dz = imu_msg->linear_acceleration.z;
    Eigen::Vector3d linear_acceleration{dx, dy, dz};
    // 得到角速度
    double rx = imu_msg->angular_velocity.x;
    double ry = imu_msg->angular_velocity.y;
    double rz = imu_msg->angular_velocity.z;
    Eigen::Vector3d angular_velocity{rx, ry, rz};
    // 上一时刻世界坐标系下加速度值
    /**
     * notes: 此处使用的是g=9.8，因此与加速度计建模推导中有所不同；如果给定的是-9.8，那么就与加速度计建模推导一致了
     */
    Eigen::Vector3d un_acc_0 = tmp_Q * (acc_0 - tmp_Ba) - estimator.g;

    // 中值陀螺仪的结果
    Eigen::Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - tmp_Bg;
    // 更新姿态
    /**
     * notes: 从理论上讲，tmp_Q并不是单位四元数，这是否会对后面产生影响
     */
    tmp_Q = tmp_Q * Utility::deltaQ(un_gyr * dt);
    // 当前时刻世界坐标系下的加速度值
    Eigen::Vector3d un_acc_1 = tmp_Q * (linear_acceleration - tmp_Ba) - estimator.g;
    // 加速度中值积分的值
    Eigen::Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
    // 经典物理中位置，速度更新方程
    tmp_P = tmp_P + dt * tmp_V + 0.5 * dt * dt * un_acc;
    tmp_V = tmp_V + dt * un_acc;

    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
}
// 用最新VIO结果更新最新imu对应的位姿
void update()
{
    TicToc t_predict;
    latest_time = current_time;
    tmp_P = estimator.Ps[WINDOW_SIZE];
    tmp_Q = estimator.Rs[WINDOW_SIZE];
    tmp_V = estimator.Vs[WINDOW_SIZE];
    tmp_Ba = estimator.Bas[WINDOW_SIZE];
    tmp_Bg = estimator.Bgs[WINDOW_SIZE];
    acc_0 = estimator.acc_0;
    gyr_0 = estimator.gyr_0;

    queue<sensor_msgs::ImuConstPtr> tmp_imu_buf = imu_buf;  // 遗留的imu的buffer，因为下面需要pop，所以copy了一份
    for (sensor_msgs::ImuConstPtr tmp_imu_msg; !tmp_imu_buf.empty(); tmp_imu_buf.pop())
        // 得到最新imu时刻的可靠的位姿
        predict(tmp_imu_buf.front());

}

// 获得匹配好的图像imu组
/**
 * notes: 直接一次性将数据对齐结束
 */
std::vector<std::pair<std::vector<sensor_msgs::ImuConstPtr>, sensor_msgs::PointCloudConstPtr>>
getMeasurements()
{
    std::vector<std::pair<std::vector<sensor_msgs::ImuConstPtr>, sensor_msgs::PointCloudConstPtr>> measurements;

    while (true)
    {
        if (imu_buf.empty() || feature_buf.empty())
            return measurements;
        // imu   *******
        // image          *****
        // 这就是image还没来
        if (!(imu_buf.back()->header.stamp.toSec() > feature_buf.front()->header.stamp.toSec() + estimator.td))
        {
            //ROS_WARN("wait for imu, only should happen at the beginning");
            sum_of_wait++;
            return measurements;
        }
        // imu        ****
        // image    ******
        // 这种只能扔掉一些image帧
        if (!(imu_buf.front()->header.stamp.toSec() < feature_buf.front()->header.stamp.toSec() + estimator.td))
        {
            ROS_WARN("throw img, only should happen at the beginning");
            feature_buf.pop();
            continue;
        }
        /**
         * notes: 为什么不判断imu_buf.front()->header.stamp.toSec() > feature_buf.back()->header.stamp.toSec() + estimator.td
         *      此时依然是没有重叠区域产生的，而且上面的if一定成立，所以会一直feature_buf.pop();
         */
        // 此时就保证了图像前一定有imu数据
        sensor_msgs::PointCloudConstPtr img_msg = feature_buf.front();
        feature_buf.pop();
        // 一般第一帧不会严格对齐，但是后面就都会对齐，当然第一帧也不会用到
        /**
         * IMU       ******************************
         * IMG            *     *     *     *
         *
         * notes: 可能第一帧图像之前可能有非常多的IMU数据，实际上第一帧的数据不会使用
         *
         * IMUs中存储的是img_msg前面的IMU信息以及其后面（也可能恰好相等）一个IMU信息
         */
        std::vector<sensor_msgs::ImuConstPtr> IMUs;
        while (imu_buf.front()->header.stamp.toSec() < img_msg->header.stamp.toSec() + estimator.td)
        {
            IMUs.emplace_back(imu_buf.front());
            imu_buf.pop();
        }
        // 保留图像时间戳后一个imu数据，但不会从buffer中扔掉
        // imu    *   *
        // image    *
        /**
         * notes: 此处主要是为了插值获取图像对应的IMU
         */
        IMUs.emplace_back(imu_buf.front());
        if (IMUs.empty())
            ROS_WARN("no imu between two image");
        measurements.emplace_back(IMUs, img_msg);

#if 0
        // 测试时间戳的对齐
        std::cout << std::fixed << std::setprecision(15);
        std::cout << "IMU的时间戳为：" << std::endl;
        for (int i = 0; i < IMUs.size(); i++) {
            if (i == IMUs.size() - 1 || i == IMUs.size() - 1)
                std::cout << "\t" << IMUs[i]->header.stamp.toSec() << std::endl;
        }
        std::cout << std::endl;
        std::cout << "image的时间戳为：\n\t" << img_msg->header.stamp.toSec() << std::endl;
        std::cout << std::endl;
#endif

    }
    return measurements;
}


/**
 * @brief imu消息存进buffer，同时按照imu频率预测位姿并发送，这样就可以提高里程计频率
 * 
 * @param[in] imu_msg 
 */
void imu_callback(const sensor_msgs::ImuConstPtr &imu_msg)
{
    if (imu_msg->header.stamp.toSec() <= last_imu_t)  // last_imu_t初值为0
    {
        ROS_WARN("imu message in disorder!");
        return;
    }
    last_imu_t = imu_msg->header.stamp.toSec();
    // 关于线程锁和条件变量：https://www.jianshu.com/p/c1dfa1d40f53
    m_buf.lock();
    imu_buf.push(imu_msg);
    m_buf.unlock();
    con.notify_one();  // 防止CPU占用过高

    // xc's todo: 为什么又进行赋值？
    last_imu_t = imu_msg->header.stamp.toSec();
    {
        std::lock_guard<std::mutex> lg(m_state);
        predict(imu_msg);
        std_msgs::Header header = imu_msg->header;
        header.frame_id = "world";
        // 只有初始化完成后才发送当前结果
        if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
            pubLatestOdometry(tmp_P, tmp_Q, tmp_V, header);
    }
}

/**
 * @brief 将前端信息送进buffer
 * 
 * @param[in] feature_msg 
 */
void feature_callback(const sensor_msgs::PointCloudConstPtr &feature_msg)
{
    if (!init_feature)  // 初值为false，第一帧忽略
    {
        //skip the first detected feature, which doesn't contain optical flow speed
        init_feature = 1;
        return;
    }
    m_buf.lock();
    feature_buf.push(feature_msg);
    m_buf.unlock();
    con.notify_one();
}

/**
 * @brief 将vins估计器复位
 * 
 * @param[in] restart_msg 
 */

void restart_callback(const std_msgs::BoolConstPtr &restart_msg)
{
    if (restart_msg->data == true)
    {
        ROS_WARN("restart the estimator!");
        m_buf.lock();
        while(!feature_buf.empty())
            feature_buf.pop();
        while(!imu_buf.empty())
            imu_buf.pop();
        m_buf.unlock();
        m_estimator.lock();
        estimator.clearState();
        estimator.setParameter();
        m_estimator.unlock();
        current_time = -1;
        last_imu_t = 0;
    }
    return;
}

// 将当前关键帧的闭环关键帧的结果保存到relo_buf之中
void relocalization_callback(const sensor_msgs::PointCloudConstPtr &points_msg)
{
    //printf("relocalization callback! \n");
    m_buf.lock();
    relo_buf.push(points_msg);
    m_buf.unlock();
}

// thread: visual-inertial odometry
void process()
{
    while (true)    // 这个线程是会一直循环下去
    {
        std::vector<std::pair<std::vector<sensor_msgs::ImuConstPtr>, sensor_msgs::PointCloudConstPtr>> measurements;
        std::unique_lock<std::mutex> lk(m_buf);
        con.wait(lk, [&]
                 {
            return (measurements = getMeasurements()).size() != 0;
                 });
        lk.unlock();    // 数据buffer的锁解锁，回调可以继续塞数据了
        m_estimator.lock(); // 进行后端求解，不能和复位重启冲突
        // 给予范围的for循环，这里就是遍历每组image imu组合
        for (auto &measurement : measurements)
        {
            auto img_msg = measurement.second;
            double dx = 0, dy = 0, dz = 0, rx = 0, ry = 0, rz = 0;
            // 遍历imu
            for (auto &imu_msg : measurement.first)
            {
                double t = imu_msg->header.stamp.toSec();
                double img_t = img_msg->header.stamp.toSec() + estimator.td;

                /**
                 * 这里仅仅只是对IMU集合的尾部做了插值，对头部没有做插值，也就是实际上没有与IMAGE对齐的IMU信息，从理论上构建的预积分不是两帧图像之间的预积分
                 * 如果单单从这里看的话，这里的逻辑是有问题的，但是在estimator.processIMU中，对第一帧IMU也做了预积分，通常来讲，第一帧是不可以做IMU预积分的
                 * 在IMU第一帧的时候做IMU预积分的策略是认为与IMAGE对齐的IMU也是IMU的第一帧，从而可以做IMU预积分；需要注意的是，只有第一次使用IMU的时候才有上面
                 * 的逻辑，一旦first_imu为true之后，acc_0和gyr_0已经是上一个预积分的最后一帧（当前预积分的第一帧）
                 *
                 * 这样写的坏处：
                 *      1. 代码可读性太差
                 *      2. 这里实际上认为起始的与IMAGE对齐的IMU信息与IMU集合的第一帧完全相同，这样做的精度弱于对头部也插值
                 */
                if (t <= img_t)
                { 
                    if (current_time < 0)
                        current_time = t;
                    double dt = t - current_time;
                    ROS_ASSERT(dt >= 0);
                    current_time = t;
                    dx = imu_msg->linear_acceleration.x;
                    dy = imu_msg->linear_acceleration.y;
                    dz = imu_msg->linear_acceleration.z;
                    rx = imu_msg->angular_velocity.x;
                    ry = imu_msg->angular_velocity.y;
                    rz = imu_msg->angular_velocity.z;
                    // 时间差和imu数据送进去
                    estimator.processIMU(dt, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz));
                    //printf("imu: dt:%f a: %f %f %f w: %f %f %f\n",dt, dx, dy, dz, rx, ry, rz);

                }
                else    // 这就是针对最后一个imu数据（最后一个IMU的时间戳也可能正好等于图像的时间戳，因此，也不用进行插值了），需要做一个简单的线性插值
                {
                    // current_time表示倒数第二个IMU时间戳，t表示倒数第一个IMU的时间戳
                    double dt_1 = img_t - current_time;
                    double dt_2 = t - img_t;
                    current_time = img_t;
                    ROS_ASSERT(dt_1 >= 0);
                    ROS_ASSERT(dt_2 >= 0);
                    ROS_ASSERT(dt_1 + dt_2 > 0);

                    double w1 = dt_2 / (dt_1 + dt_2);
                    double w2 = dt_1 / (dt_1 + dt_2);
                    dx = w1 * dx + w2 * imu_msg->linear_acceleration.x;
                    dy = w1 * dy + w2 * imu_msg->linear_acceleration.y;
                    dz = w1 * dz + w2 * imu_msg->linear_acceleration.z;
                    rx = w1 * rx + w2 * imu_msg->angular_velocity.x;
                    ry = w1 * ry + w2 * imu_msg->angular_velocity.y;
                    rz = w1 * rz + w2 * imu_msg->angular_velocity.z;
                    estimator.processIMU(dt_1, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz));
                    //printf("dimu: dt:%f a: %f %f %f w: %f %f %f\n",dt_1, dx, dy, dz, rx, ry, rz);
                }
            }
            // set relocalization frame
            // 回环相关部分
            sensor_msgs::PointCloudConstPtr relo_msg = NULL;
            while (!relo_buf.empty())   // 取出最新的回环帧，并不会都处理，性能有限
            {
                relo_msg = relo_buf.front();
                relo_buf.pop();
            }
            if (relo_msg != NULL)   // 有效回环信息
            {
                vector<Vector3d> match_points;
                double frame_stamp = relo_msg->header.stamp.toSec();    // 闭环关键帧对应的当前关键帧的时间戳
                // 闭环关键帧与其当前关键帧的匹配点，特征点是闭环关键帧的
                for (unsigned int i = 0; i < relo_msg->points.size(); i++)
                {
                    Vector3d u_v_id;
                    // 回环帧的归一化坐标和地图点idx
                    u_v_id.x() = relo_msg->points[i].x;
                    u_v_id.y() = relo_msg->points[i].y;
                    u_v_id.z() = relo_msg->points[i].z;  // 注意这里是特征点的id，不是1
                    match_points.push_back(u_v_id);
                }
                // 回环帧的位姿
                Vector3d relo_t(relo_msg->channels[0].values[0], relo_msg->channels[0].values[1], relo_msg->channels[0].values[2]);
                Quaterniond relo_q(relo_msg->channels[0].values[3], relo_msg->channels[0].values[4], relo_msg->channels[0].values[5], relo_msg->channels[0].values[6]);
                Matrix3d relo_r = relo_q.toRotationMatrix();
                int frame_index;
                frame_index = relo_msg->channels[0].values[7];  // 获取闭环的当前关键帧的索引
                estimator.setReloFrame(frame_stamp, frame_index, match_points, relo_t, relo_r);
            }

            ROS_DEBUG("processing vision data with stamp %f \n", img_msg->header.stamp.toSec());

            TicToc t_s;
            // 特征点id->特征点信息
            map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> image;
            for (unsigned int i = 0; i < img_msg->points.size(); i++)
            {
                // 从feature_tracker_node中的img_callback函数中的PUB_THIS_FRAME后的代码得知：
                int v = img_msg->channels[0].values[i] + 0.5;
                int feature_id = v / NUM_OF_CAM;
                int camera_id = v % NUM_OF_CAM;
                double x = img_msg->points[i].x;    // 去畸变后归一化相机系坐标
                double y = img_msg->points[i].y;
                double z = img_msg->points[i].z;
                double p_u = img_msg->channels[1].values[i];    // 特征点像素坐标
                double p_v = img_msg->channels[2].values[i];
                double velocity_x = img_msg->channels[3].values[i]; // 去畸变后的归一化相机系的特征点的运动速度
                double velocity_y = img_msg->channels[4].values[i];
                ROS_ASSERT(z == 1); // 检查是不是归一化
                Eigen::Matrix<double, 7, 1> xyz_uv_velocity;
                xyz_uv_velocity << x, y, z, p_u, p_v, velocity_x, velocity_y;
                /**
                 * notes: 对于双目而言，同一个特征点可能左边出现一次，右边出现一次
                 *
                 * 对某个特征点，其在哪个相机中被看到，看到的信息是什么
                 */
                image[feature_id].emplace_back(camera_id,  xyz_uv_velocity);
            }
            estimator.processImage(image, img_msg->header);

            // 一些打印以及topic的发送
            double whole_t = t_s.toc();  // 从定义这个变量到此时经过的时间
            printStatistics(estimator, whole_t);
            std_msgs::Header header = img_msg->header;
            header.frame_id = "world";

            pubOdometry(estimator, header);
            pubKeyPoses(estimator, header);
            pubCameraPose(estimator, header);
            pubPointCloud(estimator, header);
            pubTF(estimator, header);
            pubKeyframe(estimator);
            if (relo_msg != NULL)  // 将重定位的结果publish出去
                pubRelocalization(estimator);
            //ROS_ERROR("end: %f, at %f", img_msg->header.stamp.toSec(), ros::Time::now().toSec());
        }
        m_estimator.unlock();
        m_buf.lock();
        m_state.lock();
        if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
            update();
        m_state.unlock();
        m_buf.unlock();
    }
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "vins_estimator");
    ros::NodeHandle n("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);
    readParameters(n);
    estimator.setParameter();
#ifdef EIGEN_DONT_PARALLELIZE
    ROS_DEBUG("EIGEN_DONT_PARALLELIZE");
#endif
    ROS_WARN("waiting for image and imu...");

    // 注册一些publisher
    registerPub(n);
    // 接受imu消息
    ros::Subscriber sub_imu = n.subscribe(IMU_TOPIC, 2000, imu_callback, ros::TransportHints().tcpNoDelay());
    // 接受前端视觉光流结果
    ros::Subscriber sub_image = n.subscribe("/feature_tracker/feature", 2000, feature_callback);
    // 接受前端重启命令
    ros::Subscriber sub_restart = n.subscribe("/feature_tracker/restart", 2000, restart_callback);
    // 回环检测的fast relocalization响应，数据来源于对当前关键帧搜索到的闭环关键帧的地图点及其位姿
    ros::Subscriber sub_relo_points = n.subscribe("/pose_graph/match_points", 2000, relocalization_callback);

    // 核心处理线程
    std::thread measurement_process{process};
    ros::spin();

    return 0;
}
