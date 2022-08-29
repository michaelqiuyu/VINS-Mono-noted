#include "feature_manager.h"

int FeaturePerId::endFrame()
{
    return start_frame + feature_per_frame.size() - 1;
}

FeatureManager::FeatureManager(Matrix3d _Rs[])
    : Rs(_Rs)
{
    for (int i = 0; i < NUM_OF_CAM; i++)
        ric[i].setIdentity();
}

void FeatureManager::setRic(Matrix3d _ric[])
{
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        ric[i] = _ric[i];
    }
}

void FeatureManager::clearState()
{
    feature.clear();
}

/**
 * @brief 得到有效的地图点的数目
 * 
 * @return int 
 */
int FeatureManager::getFeatureCount()
{
    int cnt = 0;
    for (auto &it : feature)
    {

        it.used_num = it.feature_per_frame.size();
        // 两帧以上看到就可以三角化地图点；不要是最后几帧才看到的；
        if (it.used_num >= 2 && it.start_frame < WINDOW_SIZE - 2)
        {
            cnt++;
        }
    }
    return cnt;
}

/**
 * @brief 增加特征点信息，同时检查上一帧是否时关键帧
 * 
 * @param[in] frame_count 
 * @param[in] image 
 * @param[in] td 
 * @return true 
 * @return false 
 */
// 仅仅由Estimator::processImage调用
// xc's todo: 检查的是当前帧还是当前帧的上一帧？：检查的是当前帧跟踪的好不好，如果不好，返回true，将其作为独立的一帧；如果好的话，返回false，将其与倒数第二帧合并
bool FeatureManager::addFeatureCheckParallax(int frame_count, const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &image, double td)
{
    ROS_DEBUG("input feature: %d", (int)image.size());
    ROS_DEBUG("num of feature: %d", getFeatureCount());
    double parallax_sum = 0;
    int parallax_num = 0;
    last_track_num = 0;
    // 遍历每个特征点
    for (auto &id_pts : image)
    {
        // 用特征点信息构造一个对象
        // id_pts.second[0]表示vector<pair<int, Eigen::Matrix<double, 7, 1>>>中的第一个元素；由于本工程是单目，实际上也只有一个元素
        FeaturePerFrame f_per_fra(id_pts.second[0].second, td);

        int feature_id = id_pts.first;
        // 在已有的id中寻找是否是有相同的特征点
        auto it = find_if(feature.begin(), feature.end(), [feature_id](const FeaturePerId &it)
                          {
            return it.feature_id == feature_id;
                          });
        // 这是一个新的特征点
        if (it == feature.end())
        {
            // 在特征点管理器中，新创建一个特征点id，这里的frame_count就是该特征点在滑窗中的当前位置，作为这个特征点的起始位置
            feature.push_back(FeaturePerId(feature_id, frame_count));
            feature.back().feature_per_frame.push_back(f_per_fra);
        }
        // 如果这是一个已有的特征点，就在对应的“组织”下增加一个帧属性
        else if (it->feature_id == feature_id)
        {
            // xc's todo: 并没有记录帧号，后面是如何知道这个特征点信息是在那一帧下的？：
            // 从start_frame开始递增即可，这里跟ORB-SLAM3的特征点法不同，一旦某个特征点在某一帧不可见，那么在后面的帧都不可见，可见的帧一定是连续的
            it->feature_per_frame.push_back(f_per_fra);
            last_track_num++;   // 追踪到上一帧的特征点数目
        }
    }
#if 0
    std::cout << "last_track_num = " << last_track_num << std::endl;
#endif
    // 前两帧都设置为KF，追踪过少也认为是KF
    if (frame_count < 2 || last_track_num < 20)
        return true;

#if 0
    std::cout << "feature.size = " << feature.size() << std::endl;  // 并不会一直膨胀
#endif

    for (auto &it_per_id : feature)
    {
        // 计算的实际上是frame_count-1,也就是前一帧是否为关键帧
        // 因此起始帧至少得是frame_count - 2,同时至少覆盖到frame_count - 1帧（也就是被它看到）
        /**
         * notes: 实际上就是某个特征点至少能够同时被倒数第三帧和倒数第二帧看到
         */
        // xc's todo: 随着系统的运行，feature会不断的膨胀，那么还是会一直执行这个判断吗？有必要对所有的feature进行判断吗？:经过测试，后面会对feature处理，从而使得feature不会一直膨胀
        if (it_per_id.start_frame <= frame_count - 2 &&
            it_per_id.start_frame + int(it_per_id.feature_per_frame.size()) - 1 >= frame_count - 1)
        {
            parallax_sum += compensatedParallax2(it_per_id, frame_count);
            parallax_num++;
        }
    }
    // 这个和上一帧没有相同的特征点
    if (parallax_num == 0)
    {
        // xc's todo: 没有跟踪到的特征点该怎么办？系统跟踪失败重定位？
        return true;
    }
    else
    {
        ROS_DEBUG("parallax_sum: %lf, parallax_num: %d", parallax_sum, parallax_num);
        ROS_DEBUG("current parallax: %lf", parallax_sum / parallax_num * FOCAL_LENGTH);
        // 看看平均视差是否超过一个阈值
#if 0
        std::cout << "parallax_sum / parallax_num * FOCAL_LENGTH = " << parallax_sum / parallax_num * FOCAL_LENGTH << std::endl;
        std::cout << "parallax_sum / parallax_num = " << parallax_sum / parallax_num << std::endl;
#endif
        // 经过测试，这里的判断大部分时候都是true
        return parallax_sum / parallax_num >= MIN_PARALLAX;  // MIN_PARALLAX已经使用虚拟焦距处理过了
    }
}

void FeatureManager::debugShow()
{
    ROS_DEBUG("debug show");
    for (auto &it : feature)
    {
        ROS_ASSERT(it.feature_per_frame.size() != 0);
        ROS_ASSERT(it.start_frame >= 0);
        ROS_ASSERT(it.used_num >= 0);

        ROS_DEBUG("%d,%d,%d ", it.feature_id, it.used_num, it.start_frame);
        int sum = 0;
        for (auto &j : it.feature_per_frame)
        {
            ROS_DEBUG("%d,", int(j.is_used));
            sum += j.is_used;
            printf("(%lf,%lf) ",j.point(0), j.point(1));
        }
        ROS_ASSERT(it.used_num == sum);
    }
}

/**
 * @brief 得到同时被frame_count_l frame_count_r帧看到的特征点在各自的坐标
 * 
 * @param[in] frame_count_l 
 * @param[in] frame_count_r 
 * @return vector<pair<Vector3d, Vector3d>> 
 */

vector<pair<Vector3d, Vector3d>> FeatureManager::getCorresponding(int frame_count_l, int frame_count_r)
{
    vector<pair<Vector3d, Vector3d>> corres;
    for (auto &it : feature)
    {
        // 保证需要的特征点被这两帧都观察到
        if (it.start_frame <= frame_count_l && it.endFrame() >= frame_count_r)
        {
            Vector3d a = Vector3d::Zero(), b = Vector3d::Zero();
            // 获得在feature_per_frame中的索引
            int idx_l = frame_count_l - it.start_frame;
            int idx_r = frame_count_r - it.start_frame;

            a = it.feature_per_frame[idx_l].point;  // 在frame_count_l中的去畸变的归一化相机系坐标

            b = it.feature_per_frame[idx_r].point;  // 在frame_count_r中的去畸变的归一化相机系坐标
            
            corres.push_back(make_pair(a, b));  // 返回相机坐标系下的坐标对
        }
    }
    return corres;
}

void FeatureManager::setDepth(const VectorXd &x)
{
    int feature_index = -1;
    for (auto &it_per_id : feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        // xc's todo: 经过测试，有大量观测大于等于2，但是起始帧大于等于window_size - 2的，这些特征点不处理吗？
        //  在初始化环节，地图点存储在sfm_f中，并没有赋值给feature，因此初始化的时候，所有的地图点的深度都是-1
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;

        it_per_id.estimated_depth = 1.0 / x(++feature_index);
        //ROS_INFO("feature id %d , start_frame %d, depth %f ", it_per_id->feature_id, it_per_id-> start_frame, it_per_id->estimated_depth);
        if (it_per_id.estimated_depth < 0)
        {
            it_per_id.solve_flag = 2;
        }
        else
            it_per_id.solve_flag = 1;
    }
}

/**
 * @brief 移除一些不能被三角化的点
 * 
 */
void FeatureManager::removeFailures()
{
    for (auto it = feature.begin(), it_next = feature.begin();
         it != feature.end(); it = it_next)
    {
        it_next++;
        if (it->solve_flag == 2)
            feature.erase(it);
    }
}

/**
 * @brief 把给定的深度赋值给各个特征点作为逆深度
 * 
 * @param[in] x 
 */
void FeatureManager::clearDepth(const VectorXd &x)
{
    int feature_index = -1;
    for (auto &it_per_id : feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        // xc's todo: 经过测试，有大量观测大于等于2，但是起始帧大于等于window_size - 2的，这些特征点不处理吗？
        //  在初始化环节，地图点存储在sfm_f中，并没有赋值给feature，因此初始化的时候，所有的地图点的深度都是-1
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        it_per_id.estimated_depth = 1.0 / x(++feature_index);
    }
}

/**
 * @brief 得到特征点的逆深度
 * 
 * @return VectorXd 
 */
VectorXd FeatureManager::getDepthVector()
{
    VectorXd dep_vec(getFeatureCount());  // 获得有效的地图点的数目，并且这些地图点不是从倒数第三帧才开始看到的
    int feature_index = -1;
    for (auto &it_per_id : feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
#if 0
//        std::cout << "it_per_id.used_num = " << it_per_id.used_num << std::endl;
//        std::cout << "it_per_id.start_frame = " << it_per_id.start_frame << std::endl;
        if (it_per_id.used_num >= 2 && it_per_id.start_frame >= WINDOW_SIZE - 2) {
            std::cout << "这个正常的地图点也不处理吗" << std::endl;
            std::cout << "这个特征点的深度为：" << it_per_id.estimated_depth << std::endl;
        }

#endif
        // xc's todo: 经过测试，有大量观测大于等于2，但是起始帧大于等于window_size - 2的，这些特征点不处理吗？
        //  在初始化环节，地图点存储在sfm_f中，并没有赋值给feature，因此初始化的时候，所有的地图点的深度都是-1
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
#if 1
        dep_vec(++feature_index) = 1. / it_per_id.estimated_depth;  // 使用的是逆深度
#else
        dep_vec(++feature_index) = it_per_id->estimated_depth;
#endif
    }
    return dep_vec;
}

/**
 * @brief 利用观测到该特征点的所有位姿来三角化特征点
 * 
 * @param[in] Ps 
 * @param[in] tic 
 * @param[in] ric 
 */
void FeatureManager::triangulate(Vector3d Ps[], Vector3d tic[], Matrix3d ric[])
{
    /**
     * 多帧三角化，利用所有的观测来构建约束，得到最小二乘解
     * 注意这里的计算是在start_frame下进行的，得到的深度也是start_frame坐标系下的深度
     * 特征点管理器中的特征点并没有在sfm中重建出来
     */
    // 遍历每一个特征点
    for (auto &it_per_id : feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        // 如果是最近才看到的特征点也不会处理，尽管其观测大于2
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;

        if (it_per_id.estimated_depth > 0)  // 代表已经三角化过了
            continue;
        int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;

        ROS_ASSERT(NUM_OF_CAM == 1);
        // 被一帧看到就可以构建2个方程
        Eigen::MatrixXd svd_A(2 * it_per_id.feature_per_frame.size(), 4);
        int svd_idx = 0;

        Eigen::Matrix<double, 3, 4> P0;
        // Twi -> Twc,第一个观察到这个特征点的KF的位姿
        /**
         * 在初始化里面，world是枢纽帧，tic[0]是零向量，因为这个时候的Ps并没有表示IMU的信息，还是表示枢纽帧下的相机平移
         * 因此t0 = Ps[imu_i] + Rs[imu_i] * tic[0] = Ps[imu_i]实际上就是枢纽帧（也就是当前的世界系）下的相机平移，正好跟R0对应上
         *
         * 这里的代码写的非常不好，可读性太差了，很容易让人感到疑惑，代码最好还是解耦一些比较好
         */
        Eigen::Vector3d t0 = Ps[imu_i] + Rs[imu_i] * tic[0];
        Eigen::Matrix3d R0 = Rs[imu_i] * ric[0];
        P0.leftCols<3>() = Eigen::Matrix3d::Identity();
        P0.rightCols<1>() = Eigen::Vector3d::Zero();
        // 遍历所有看到这个特征点的KF
        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            // 得到该KF的相机坐标系位姿
            Eigen::Vector3d t1 = Ps[imu_j] + Rs[imu_j] * tic[0];
            Eigen::Matrix3d R1 = Rs[imu_j] * ric[0];
            // T_w_cj -> T_c0_cj
            Eigen::Vector3d t = R0.transpose() * (t1 - t0);
            Eigen::Matrix3d R = R0.transpose() * R1;
            Eigen::Matrix<double, 3, 4> P;
            // T_c0_cj -> T_cj_c0相当于把c0当作世界系
            P.leftCols<3>() = R.transpose();
            P.rightCols<1>() = -R.transpose() * t;
            // 这一步实际上没有必要
            Eigen::Vector3d f = it_per_frame.point.normalized();
            // 构建超定方程的其中两个方程：每一帧提供两个方程
            // 使用的形式是：(u * p3.t - p1.t) * p = 0, (v * p3.t - p2.t) * p = 0;
            // 注意这里使用的是归一化相机系坐标，不是图像坐标，但是推导过程是一样的
            svd_A.row(svd_idx++) = f[0] * P.row(2) - f[2] * P.row(0);
            svd_A.row(svd_idx++) = f[1] * P.row(2) - f[2] * P.row(1);

            if (imu_i == imu_j)
                continue;
        }
        ROS_ASSERT(svd_idx == svd_A.rows());
        Eigen::Vector4d svd_V = Eigen::JacobiSVD<Eigen::MatrixXd>(svd_A, Eigen::ComputeThinV).matrixV().rightCols<1>();
        // 求解齐次坐标下的深度
        double svd_method = svd_V[2] / svd_V[3];
        //it_per_id->estimated_depth = -b / A;
        //it_per_id->estimated_depth = svd_V[2] / svd_V[3];
        // 得到的深度值实际上就是第一个观察到这个特征点的相机坐标系下的深度值
        // xc's todo: 仅仅使用了计算的深度，在后面的优化中直接按照在start_frame下的归一化相机系坐标乘以深度得到相机系的三维坐标：
        /**
         * 这里隐含了假设：地图点在start_frame下的重投影误差为0，也就是在start_frame的射线上寻找一个深度而已
         * 但是其实际上的计算又不是按照这个假设求解的，还是按照正常的三角化求解，只获取其深度直接赋值，这不是前后矛盾吗？
         * 如果按照假设的话，应该只是求解start_frame下的深度s，也就是pw = s * R1.t * u1 - R1.t * t1
         *
         * 对于任何一个观测有，u2.hat * (R2 * pw + t2) = 0，可以化简得到sp1 = p2，构建∑||sp1 - p2||^2来求解
         *
         * 在initial_sfm.cpp中的triangulateTwoFrames，使用两视图三角化进行测试，发现：
         *      1. 使用我提出的方法与DLT的方法计算的世界点的坐标大部分时候都相差在2cm以内
         *      2. 还是有时候相差达到10cm以上
         * 
         * 能够运行的原因：
         *      这里的深度值仅仅只是一个初值，并不会是最终的值，因此只要初值不是非常偏离真值，对最后的结果的影响可能是可以忽略不计的
         */
        // xc's todo: 后续可以完整实现仅仅计算深度值的方式构建最小二乘，跑通slam，与原始版本进行对比
        it_per_id.estimated_depth = svd_method;
        //it_per_id->estimated_depth = INIT_DEPTH;

        if (it_per_id.estimated_depth < 0.1)
        {
            it_per_id.estimated_depth = INIT_DEPTH; // 距离太近就设置成默认值
        }

    }
}

void FeatureManager::removeOutlier()
{
    ROS_BREAK();
    int i = -1;
    for (auto it = feature.begin(), it_next = feature.begin();
         it != feature.end(); it = it_next)
    {
        it_next++;
        i += it->used_num != 0;
        if (it->used_num != 0 && it->is_outlier == true)
        {
            feature.erase(it);
        }
    }
}

/**
 * @brief 
 * 
 * @param[in] marg_R  被移除的位姿
 * @param[in] marg_P 
 * @param[in] new_R    转接地图点的位姿
 * @param[in] new_P 
 */

void FeatureManager::removeBackShiftDepth(Eigen::Matrix3d marg_R, Eigen::Vector3d marg_P, Eigen::Matrix3d new_R, Eigen::Vector3d new_P)
{
    for (auto it = feature.begin(), it_next = feature.begin();
         it != feature.end(); it = it_next)
    {
        it_next++;

        if (it->start_frame != 0)   // 如果不是被移除的帧看到，那么该地图点对应的起始帧id减一
            it->start_frame--;
        else
        {
            Eigen::Vector3d uv_i = it->feature_per_frame[0].point;    // 取出归一化相机坐标系坐标
            it->feature_per_frame.erase(it->feature_per_frame.begin()); // 该点不再被原来的第一帧看到，因此从中移除
            if (it->feature_per_frame.size() < 2)   // 如果这个地图点没有至少被两帧看到：注意执行这个判断的时候已经删除了一个观测了
            {
                feature.erase(it);  // 那他就没有存在的价值了
                continue;
            }
            else    // 进行管辖权的转交
            {
                Eigen::Vector3d pts_i = uv_i * it->estimated_depth; // 起始相机坐标系下的坐标
                Eigen::Vector3d w_pts_i = marg_R * pts_i + marg_P;  // 转到世界坐标系下
                Eigen::Vector3d pts_j = new_R.transpose() * (w_pts_i - new_P);  // 转到新的最老帧的相机坐标系下，由于起始帧被移除了，因此需要重新计算起始帧下的坐标
                double dep_j = pts_j(2);
                if (dep_j > 0)  // 看看深度是否有效
                    it->estimated_depth = dep_j;    // 有效的话就得到在现在最老帧下的深度值
                else
                    it->estimated_depth = INIT_DEPTH;   // 无效就设置默认值
            }
        }
        // remove tracking-lost feature after marginalize
        /*
        if (it->endFrame() < WINDOW_SIZE - 1)
        {
            feature.erase(it);
        }
        */
    }
}
/**
 * @brief 这个还没初始化结束，因此相比刚才，不进行地图点新的深度的换算，因为此时还有进行视觉惯性对齐
 * 
 */
void FeatureManager::removeBack()
{
    for (auto it = feature.begin(), it_next = feature.begin();
         it != feature.end(); it = it_next)
    {
        it_next++;

        if (it->start_frame != 0)
            it->start_frame--;
        else
        {
            it->feature_per_frame.erase(it->feature_per_frame.begin());
            if (it->feature_per_frame.size() == 0)
                feature.erase(it);
        }
    }
}

// 对margin倒数第二帧进行处理
void FeatureManager::removeFront(int frame_count)
{
    for (auto it = feature.begin(), it_next = feature.begin(); it != feature.end(); it = it_next)
    {
        it_next++;

        if (it->start_frame == frame_count) // 如果地图点被最后一帧看到，由于滑窗，他的起始帧减1
        {
            it->start_frame--;
        }
        else
        {
            int j = WINDOW_SIZE - 1 - it->start_frame;  // 倒数第二帧在这个地图点对应KF vector的idx
            if (it->endFrame() < frame_count - 1)   // 如果该地图点不能被倒数第二帧看到，那没什么好做的
                continue;
            it->feature_per_frame.erase(it->feature_per_frame.begin() + j); // 能被倒数第二帧看到，erase掉这个索引，也就是这个观测
            if (it->feature_per_frame.size() == 0)  // 如果这个地图点没有别的观测了
                feature.erase(it);  // 就没有存在的价值了
        }
    }
}

// xc's todo: 此处使用的是归一化的相机坐标系的坐标，如果按照针孔模型的话，要达到10的间隔，说明像素间隔有几千了，是否在前面做了什么处理？：阈值MIN_PARALLAX已经使用虚拟焦距处理过了
// 仅仅由FeatureManager::addFeatureCheckParallax调用
double FeatureManager::compensatedParallax2(const FeaturePerId &it_per_id, int frame_count)
{
    //check the second last frame is keyframe or not
    //parallax between second last frame and third last frame
    // 找到相邻两帧
    const FeaturePerFrame &frame_i = it_per_id.feature_per_frame[frame_count - 2 - it_per_id.start_frame];  // 倒数第三帧
    const FeaturePerFrame &frame_j = it_per_id.feature_per_frame[frame_count - 1 - it_per_id.start_frame];  // 倒数第二帧

    double ans = 0;
    Vector3d p_j = frame_j.point;

    double u_j = p_j(0);
    double v_j = p_j(1);

    Vector3d p_i = frame_i.point;
    Vector3d p_i_comp;

    //int r_i = frame_count - 2;
    //int r_j = frame_count - 1;
    //p_i_comp = ric[camera_id_j].transpose() * Rs[r_j].transpose() * Rs[r_i] * ric[camera_id_i] * p_i;
    p_i_comp = p_i;
    double dep_i = p_i(2);
    double u_i = p_i(0) / dep_i;
    double v_i = p_i(1) / dep_i;
    double du = u_i - u_j, dv = v_i - v_j;  // 归一化相机坐标系的坐标差
    // 当都是归一化坐标系时，他们两个都是一样的
    double dep_i_comp = p_i_comp(2);
    double u_i_comp = p_i_comp(0) / dep_i_comp;
    double v_i_comp = p_i_comp(1) / dep_i_comp;
    double du_comp = u_i_comp - u_j, dv_comp = v_i_comp - v_j;

    /**
     * 对去畸变的归一化相机系而言：
     * u = fx * x + cx, v = fy * y + cy
     * (u2 - u1)^2 + (v2 - v1)^2 = fx^2 * (x2 - x1)^2 + fy^2 * (y2 - y1)^2
     *
     * 本系统中大量使用虚拟焦距来做统一的阈值处理，但是如果fx!=fy，那么使用一个虚拟焦距来处理的话并不一定合适，尤其是二者相差较大的时候，最好使用两个虚拟焦距
     * 比如在这里，使用统一的虚拟焦距f，那么有[(u2 - u1)^2 + (v2 - v1)^2] / f^2 = (fx/f)^2 * [(x2 - x1)^2 + (fy/fx)^2 * (y2 - y1)^2]，很显然，后面中括号里面的表达式并不是归一化相机系下的距离
     * 而在这里直接使用了归一化相机系下的距离来与阈值判断，这显然是不合适的
     *
     * 合理的设想：设置两个虚拟焦距f1和f2，且f1 / fx = f2 / fy，但是实际上在fx和fy位置的情况下，不可能设置出f1和f2；因此，此工程实际上限制的fx和fy的相对大小，也就是此工程要求fx和fy尽量接近
     */

    // 从函数的调用逻辑看，这里传入的point的坐标一定是去畸变的归一化相机系坐标，来源是estimator_node.cpp：xyz_uv_velocity << x, y, z, p_u, p_v, velocity_x, velocity_y;
    ans = max(ans, sqrt(min(du * du + dv * dv, du_comp * du_comp + dv_comp * dv_comp)));

    return ans;
}