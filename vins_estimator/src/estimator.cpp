#include "estimator.h"

int test_index = 0;
int ini_num = 0;

Estimator::Estimator(): f_manager{Rs}
{
    ROS_INFO("init begins");
    clearState();
}

/**
 * @brief 外参，重投影置信度，延时设置
 * 
 */
void Estimator::setParameter()
{
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = TIC[i];
        ric[i] = RIC[i];
    }
    f_manager.setRic(ric);
    // 这里可以看到虚拟相机的用法
    /**
     * notes: 重投影误差的阈值为1.5个像素
     */
    ProjectionFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    ProjectionTdFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    td = TD;
}

// 所有状态全部重置
void Estimator::clearState()
{
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        Rs[i].setIdentity();
        Ps[i].setZero();
        Vs[i].setZero();
        Bas[i].setZero();
        Bgs[i].setZero();
        dt_buf[i].clear();
        linear_acceleration_buf[i].clear();
        angular_velocity_buf[i].clear();

        if (pre_integrations[i] != nullptr)
            delete pre_integrations[i];
        pre_integrations[i] = nullptr;
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = Vector3d::Zero();
        ric[i] = Matrix3d::Identity();
    }

    for (auto &it : all_image_frame)
    {
        if (it.second.pre_integration != nullptr)
        {
            delete it.second.pre_integration;
            it.second.pre_integration = nullptr;
        }
    }

    solver_flag = INITIAL;
    first_imu = false,
    sum_of_back = 0;
    sum_of_front = 0;
    frame_count = 0;
    solver_flag = INITIAL;
    initial_timestamp = 0;
    all_image_frame.clear();
    td = TD;


    if (tmp_pre_integration != nullptr)
        delete tmp_pre_integration;
    if (last_marginalization_info != nullptr)
        delete last_marginalization_info;

    tmp_pre_integration = nullptr;
    last_marginalization_info = nullptr;
    last_marginalization_parameter_blocks.clear();

    f_manager.clearState();

    failure_occur = 0;
    relocalization_info = 0;

    drift_correct_r = Matrix3d::Identity();
    drift_correct_t = Vector3d::Zero();
}

/**
 * @brief 对imu数据进行处理，包括更新预积分量，和提供优化状态量的初始值
 * 
 * @param[in] dt 
 * @param[in] linear_acceleration 
 * @param[in] angular_velocity 
 */
void Estimator::processIMU(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity)
{
    if (!first_imu)
    {
        // notes: 注意只有第一次使用IMU的时候才会执行这个逻辑，这个逻辑是为了在第一次构建IMU预积分的时候与第一帧图像对齐；
        // 后面没有这个逻辑了，因为acc_0和gyr_0已经更新并与初始图像帧对齐了
        // 这里的处理并不好，逻辑耦合太强了，可读性太差了
        first_imu = true;  // 初始化时候为false
        // acc_0和gyr_0并不是表示初始的IMU信息，而是上一帧的信息
        acc_0 = linear_acceleration;
        gyr_0 = angular_velocity;
#if 0
        std::cout << std::fixed << std::setprecision(15);
        std::cout << "acc_0 = " << acc_0 << ", gry_0 = " << gyr_0 << std::endl;
#endif
    }
#if 0
    if (pre_integrations[frame_count]) {
        std::cout << "帧号：" << frame_count << "对应的预积分已经new出来" << std::endl;
    } else {
        std::cout << "帧号：" << frame_count << "对应的预积分还没有new出来" << std::endl;
    }
#endif
    // 滑窗中保留11帧，frame_count表示现在处理第几帧，一般处理到第11帧时就保持不变了
    // 由于预积分是帧间约束，因此第1个预积分量实际上是用不到的
    if (!pre_integrations[frame_count])
    {
        // 初始化的时候都是空，还有这里new出来后并没有return，而是继续往后面运行
        pre_integrations[frame_count] = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};
    }
    // 所以只有大于0才处理
    if (frame_count != 0)  // 第一帧不处理
    {
        // 保存dt、acc、gyr，并计算雅克比和协方差
        // 注意：即使pre_integrations[frame_count]刚刚new出来，也会做IMU预积分，这也是为什么IMU头部没有插值的原因
        pre_integrations[frame_count]->push_back(dt, linear_acceleration, angular_velocity);
        //if(solver_flag != NON_LINEAR)
            // 这个量用来做初始化用的，注意当frame_count=0的时候是个空指针
            tmp_pre_integration->push_back(dt, linear_acceleration, angular_velocity);  // 在estimator.processImage中new
        // 保存传感器数据
        // dt_buf保存的是所有的IMU信息，当然是分vector存储每一个区间的IMU信息；而IntegrationBase中仅仅存放当前区间中的IMU信息
        dt_buf[frame_count].push_back(dt);
        linear_acceleration_buf[frame_count].push_back(linear_acceleration);
        angular_velocity_buf[frame_count].push_back(angular_velocity);
        // 又是一个中值积分，更新滑窗中状态量，本质是给非线性优化提供可信的初始值
        int j = frame_count;         
        Vector3d un_acc_0 = Rs[j] * (acc_0 - Bas[j]) - g;  // 此时Rs[j]还没有更新
        Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - Bgs[j];
        Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
        Vector3d un_acc_1 = Rs[j] * (linear_acceleration - Bas[j]) - g;
        Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
        Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc;
        Vs[j] += dt * un_acc;
    }
    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
}

void Estimator::processImage(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &image, const std_msgs::Header &header)
{
    ROS_DEBUG("new image coming ------------------------------------------");
    ROS_DEBUG("Adding feature points %lu", image.size());
    // Step 1 将特征点信息加到f_manager这个特征点管理器中，同时进行是否关键帧的检查
    if (f_manager.addFeatureCheckParallax(frame_count, image, td))
        // 如果上一帧是关键帧，则滑窗中最老的帧就要被移出滑窗
        marginalization_flag = MARGIN_OLD;
    else
        // 否则移除上一帧
        marginalization_flag = MARGIN_SECOND_NEW;
    /**
     * 通过下面的测试代码可以看到：
     *      1. 当为old的时候就会生成关键帧，当为new的时候只会有里程计的信息，不会生成关键帧
     *      2. 从后面可以看到，如果为new，那么只是会把新的一帧与最后一帧的预积分进行合并，并删除最后一帧看到的地图点，当然这一帧的位姿也就没有了
     *      3. 这样的删除操作并没有充分使用数据，而且也会使得关键帧的生成非常频繁
     */
#if 0
    getchar();
    std::cout << "marginalization_flag = ";
    if (marginalization_flag == MARGIN_OLD) {
        std::cout << "MARGIN_OLD" << std::endl;
    } else {
        std::cout << "MARGIN_SECOND_NEW" << std::endl;
    }
#endif

    ROS_DEBUG("this frame is--------------------%s", marginalization_flag ? "reject" : "accept");
    ROS_DEBUG("%s", marginalization_flag ? "Non-keyframe" : "Keyframe");
    ROS_DEBUG("Solving %d", frame_count);
    ROS_DEBUG("number of feature: %d", f_manager.getFeatureCount());
    Headers[frame_count] = header;  // 对应的是当前预积分的结束时刻
#if 0
    std::cout << "frame_count = " << frame_count << std::endl;
    std::cout << "Headers[frame_count]" << Headers[frame_count].stamp.toSec() << std::endl;
#endif

    // all_image_frame用来做初始化相关操作，他保留滑窗起始到当前的所有帧
    // 有一些帧会因为不是KF，被MARGIN_SECOND_NEW，但是及时较新的帧被margin，他也会保留在这个容器中，因为初始化要求使用所有的帧，而非只要KF
    ImageFrame imageframe(image, header.stamp.toSec());
    // 对于图像的第一帧，其tmp_pre_intergration为nullptr，也对应了图像第一帧之前的IMU不做预积分，不使用
    imageframe.pre_integration = tmp_pre_integration;  // 注意这里使用的是tmp_pre_intergration，而不是pre_integrations
    // 这里就是简单的把图像和预积分绑定在一起，这里预积分就是两帧之间的，滑窗中实际上是两个KF之间的
    // 实际上是准备用来初始化的相关数据
    all_image_frame.insert(make_pair(header.stamp.toSec(), imageframe));
    // tmp_pre_integration是对每一个普通帧都做预积分
    tmp_pre_integration = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};

    // 没有外参初值
    // Step 2： 外参初始化
    if(ESTIMATE_EXTRINSIC == 2)
    {
        ROS_INFO("calibrating extrinsic param, rotation movement is needed");
        if (frame_count != 0)
        {
            // 这里标定imu和相机的旋转外参的初值
            // 因为预积分是相邻帧的约束，因为这里得到的图像关联也是相邻的
            vector<pair<Vector3d, Vector3d>> corres = f_manager.getCorresponding(frame_count - 1, frame_count);
            Matrix3d calib_ric;
            // pre_integrations[frame_count]代表着从frame_count-1到frame_count之间的预积分
            if (initial_ex_rotation.CalibrationExRotation(corres, pre_integrations[frame_count]->delta_q, calib_ric)) {
                ROS_WARN("initial extrinsic rotation calib success");
                ROS_WARN_STREAM("initial extrinsic rotation: " << endl << calib_ric);
                ric[0] = calib_ric;  // ric表示多个相机与IMU的外参，这里只标定了第一个
                RIC[0] = calib_ric;
                // 标志位设置成可信的外参初值
                ESTIMATE_EXTRINSIC = 1;
            }
        }
    }

#if 0
    std::cout << "frame_count = " << frame_count << std::endl;
#endif

    if (solver_flag == INITIAL)
    {
        if (frame_count == WINDOW_SIZE) // 有足够的帧数
        {
            bool result = false;
            // 要有可信的外参值，同时距离上次初始化不成功至少相邻0.1s: 不成功的附近依然不会成功
            // Step 3： VIO初始化
            if( ESTIMATE_EXTRINSIC != 2 && (header.stamp.toSec() - initial_timestamp) > 0.1)  // initial_timestamp初始化为0
            {
               result = initialStructure();
               initial_timestamp = header.stamp.toSec();
            }
            if(result)
            {
                solver_flag = NON_LINEAR;
                // Step 4： 非线性优化求解VIO
                solveOdometry();
                // Step 5： 滑动窗口
                slideWindow();
                // Step 6： 移除无效地图点
                f_manager.removeFailures(); // 移除无效地图点
                ROS_INFO("Initialization finish!");
                last_R = Rs[WINDOW_SIZE];   // 滑窗里最新的位姿
                last_P = Ps[WINDOW_SIZE];
                last_R0 = Rs[0];    // 滑窗里最老的位姿
                last_P0 = Ps[0];
                
            }
            else
                slideWindow();
        }
        else
            // xc's todo: frame_count只在这里递增，如果上面初始化失败的话，frame_count就不变了吗？：frame_count到达window_size之后就不变了
            frame_count++;  //
    }
    else
    {
        TicToc t_solve;
        solveOdometry();
        ROS_DEBUG("solver costs: %fms", t_solve.toc());
        // 检测VIO是否正常
        if (failureDetection())
        {
            ROS_WARN("failure detection!");
            failure_occur = 1;
            // 如果异常，重启VIO
            clearState();
            setParameter();
            ROS_WARN("system reboot!");
            return;
        }

        TicToc t_margin;
        slideWindow();
        f_manager.removeFailures();
        ROS_DEBUG("marginalization costs: %fms", t_margin.toc());
        // prepare output of VINS
        // 给可视化用的
        key_poses.clear();
        for (int i = 0; i <= WINDOW_SIZE; i++)
            key_poses.push_back(Ps[i]);

        last_R = Rs[WINDOW_SIZE];
        last_P = Ps[WINDOW_SIZE];
        last_R0 = Rs[0];
        last_P0 = Ps[0];
    }

#if 0
    std::cout << "frame_count = " << frame_count << std::endl;
#endif
}

/**
 * @brief VIO初始化，将滑窗中的P V Q恢复到第0帧并且和重力对齐
 * 
 * @return true 
 * @return false 
 */
bool Estimator::initialStructure()
{
    TicToc t_sfm;
    // Step 1 check imu observibility：检测运动速度的方差
    {
        // 此操作不会产生时机影响，即使IMU激励不够，也会往后面执行

        /**
         * 此处使用预积分的V，而不是直接使用加速度计测量的比力来计算的原因：
         *      1. 预积分的V只需要用到IMU的测量值就可以计算
         *      2. 如果要使用比力，那么a_w = Rwb * (f_b - ba - na) + g，此时需要知道Rwb，也就是需要知道IMU的位姿，这在初始化里面是未知的
         *      3. 因此，这里使用delta_v来计算三个轴加速度方差的总和
         *      4. beita / delta_t = R * (a - gw)，由于gw是固定值，因此beita / delta_t的波动也就代表了a的波动
         *
         * 做这个操作的原因：
         *      1. IMU预积分是基于牛顿第二定律的，如果body不动或者匀速运动，那么aw=0，那么就没有办法积分出速度和加速度了，IMU预积分也无法计算了
         *      2. 因此，这里的操作就是为了避免body静止或者匀速运动导致aw接近0
         */
        // xc's todo: 为什么使用加速度的标准差来判断激励，即使方差为0，也只能说明是匀加速运动而非匀速运动呀
        map<double, ImageFrame>::iterator frame_it;
        Vector3d sum_g;
#if 0
        Vector3d sum_v;
#endif
        // 从第二帧开始检查imu
        for (frame_it = all_image_frame.begin(), frame_it++; frame_it != all_image_frame.end(); frame_it++)
        {
            double dt = frame_it->second.pre_integration->sum_dt;
            Vector3d tmp_g = frame_it->second.pre_integration->delta_v / dt;
            // 累加每一帧带重力加速度的deltav：预积分中的delta_v
            sum_g += tmp_g;
            /**
             * 1. 使用预积分的位移除以时间来做速度的判断是不对的，这个可以从预积分的位移公式看出，因为预积分的位姿并不代表两帧之间的位移
             *      alpha / delta_t = R * (v - v_bk + 0.5 * g * delta_t)；由于delta_t和v_bk都是变量，因此alpha / delta_t的波动不能代表v的波动
             * 2. 直接使用预积分的速度是否接近某个常数来判断也是不行的，因为beita = R * (delta_v - g * delta_t)，由于delta_t是变量，因此beita为常数也不能说明delta_v等于0
             * 3. 视觉初始化之后，使用枢纽帧下的平移参数来计算也是不行的，因为此时还没有尺度（尺度还是依赖预积分得到的，互为条件了），阈值无法给出
             *
             * 在没有其他传感器的情况下，真实尺度的速度是拿不到的，因此没办法做速度的方差校验；这里适应加速度来判断，相当于将条件变严格了；
             * 如果加速度的条件都满足了，速度的条件一定能满足
             */
#if 0
            Vector3d tmp_v = frame_it->second.pre_integration->delta_p / dt;
            std::cout << "frame_it->second.pre_integration->delta_p = " << frame_it->second.pre_integration->delta_p << std::endl;
            sum_v += tmp_v;
#endif
        }
        Vector3d aver_g;
        aver_g = sum_g * 1.0 / ((int)all_image_frame.size() - 1);
        double var = 0;
#if 0
        Vector3d aver_v;
        aver_v = sum_v * 1.0 / ((int)all_image_frame.size() - 1);
        std::cout << "aver_v = " << aver_v << std::endl;
        double var_v = 0;
#endif
        for (frame_it = all_image_frame.begin(), frame_it++; frame_it != all_image_frame.end(); frame_it++)
        {
            double dt = frame_it->second.pre_integration->sum_dt;
            Vector3d tmp_g = frame_it->second.pre_integration->delta_v / dt;
            // 求方差: 三个轴上的方差一起计算
            var += (tmp_g - aver_g).transpose() * (tmp_g - aver_g);
            //cout << "frame g " << tmp_g.transpose() << endl;

#if 0
            Vector3d tmp_v = frame_it->second.pre_integration->delta_p / dt;
            var_v += (tmp_v - aver_v).transpose() * (tmp_v - aver_v);
#endif

        }
        // 得到的标准差：使用的是无偏估计的方差
        var = sqrt(var / ((int)all_image_frame.size() - 1));
#if 0
        var_v = sqrt(var_v / ((int)all_image_frame.size() - 1));
        std::cout << "var = " << var << std::endl;
        std::cout << "var_v = " << var_v << std::endl;
#endif
        //ROS_WARN("IMU variation %f!", var);
        // 实际上检查结果并没有用
        if(var < 0.25)
        {
            ROS_INFO("IMU excitation not enouth!");
            //return false;
        }
    }
    // Step 2 global sfm
    // 做一个纯视觉slam
    Quaterniond Q[frame_count + 1];
    Vector3d T[frame_count + 1];
    map<int, Vector3d> sfm_tracked_points;
    vector<SFMFeature> sfm_f;   // 保存每个特征点的信息
    // 遍历所有的特征点
    for (auto &it_per_id : f_manager.feature)
    {
        int imu_j = it_per_id.start_frame - 1;  // 这个跟imu无关，就是存储观测特征点的帧的索引
        SFMFeature tmp_feature; // 用来后续做sfm
        tmp_feature.state = false;
        tmp_feature.id = it_per_id.feature_id;
        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;  // 注意这里递增了，这样也解释了上面为什么减1
            Vector3d pts_j = it_per_frame.point;
            // 帧号以及去畸变的归一化相机系坐标
            tmp_feature.observation.push_back(make_pair(imu_j, Eigen::Vector2d{pts_j.x(), pts_j.y()}));
        }
        sfm_f.push_back(tmp_feature);
    } 
    Matrix3d relative_R;
    Vector3d relative_T;
    int l;
    // 如果以枢纽帧为世界系，那么得到的位姿的表达应该是Twc，world代表枢纽帧，camera代表最后一帧
    if (!relativePose(relative_R, relative_T, l))
    {
        ROS_INFO("Not enough features or parallax; Move device around");
        return false;
    }
#if 0
    std::cout << "枢纽帧 = " << l << std::endl;
    std::cout << "ini_num = " << ini_num++ << std::endl;
#endif
    GlobalSFM sfm;
    // 进行sfm的求解：需要注意的是，relative_R和relative_T仅仅通过E得到，并没有经过优化环节，后续应该会有一个统一的GBA
    // 这里的sfm求解的是位姿，尽管内部也求解了地图点（sfm_f），但是后面没有使用，而是利用位姿重新三角化了
    if(!sfm.construct(frame_count + 1, Q, T, l,
              relative_R, relative_T,
              sfm_f, sfm_tracked_points))
    {
#if 0
        std::cout << "求解失败1" << std::endl;
#endif
        ROS_DEBUG("global SFM failed!");
        marginalization_flag = MARGIN_OLD;
        return false;
    }
#if 0
    std::cout << "求解失败2" << std::endl;
    if (test_index < 1) {
        map<double, ImageFrame>::iterator frame_it_;
        map<int, Vector3d>::iterator it_;
        frame_it_ = all_image_frame.begin();
        std::cout << std::fixed << std::setprecision(15);
        for (int i = 0; frame_it_ != all_image_frame.end( ); frame_it_++) {
            cv::Mat r, rvec, t, D, tmp_r;
            std::cout << "i = " << i << std::endl;
            std::cout << "frame_it->first = " << frame_it_->first << std::endl;
            std::cout << "Headers[i].stamp.toSec() = " << Headers[i].stamp.toSec() << std::endl;
            i++;
        }
        test_index++;
        return false;
    }

    std::cout << "求解失败3" << std::endl;
    if (test_index < 2) {
        std::cout << "second" << std::endl;
        map<double, ImageFrame>::iterator frame_it_;
        map<int, Vector3d>::iterator it_;
        frame_it_ = all_image_frame.begin();
        std::cout << std::fixed << std::setprecision(15);
        for (int i = 0; frame_it_ != all_image_frame.end( ); frame_it_++) {
            cv::Mat r, rvec, t, D, tmp_r;
            std::cout << "i = " << i << std::endl;
            std::cout << "frame_it->first = " << frame_it_->first << std::endl;
            std::cout << "Headers[i].stamp.toSec() = " << Headers[i].stamp.toSec() << std::endl;
            i++;
        }
    }
#endif

#if 0
    // 测试all_image_frame与Headers的大小
    std::cout << "all_image_frame.size = " << all_image_frame.size() << std::endl;
    std::cout << "Headers[frame_count].stamp.toSec() = " << Headers[frame_count].stamp.toSec() << std::endl;
#endif

    // Step 3 solve pnp for all frame
    // step2只是针对KF进行sfm，初始化需要all_image_frame中的所有元素，因此下面通过KF来求解其他的非KF的位姿
    map<double, ImageFrame>::iterator frame_it;
    map<int, Vector3d>::iterator it;
    frame_it = all_image_frame.begin( );
    // i代表跟这个帧最近的KF的索引，用来提供pnp的初始值
    for (int i = 0; frame_it != all_image_frame.end( ); frame_it++)
    {
        // provide initial guess
        cv::Mat r, rvec, t, D, tmp_r;
        // 这一帧本身就是KF，因此可以直接得到位姿
        // 对于正常的初始化，初始化的那些帧都会被认为是关键帧
        if((frame_it->first) == Headers[i].stamp.toSec())
        {
#if 0
            std::cout << "关键帧：" << i << std::endl;
#endif
            frame_it->second.is_key_frame = true;
            // notes: 这里的R和T表示的对象是不同的，因此后面的某些操作是比较具有迷惑性的，但只要知道这一点就能很清晰的理解
            frame_it->second.R = Q[i].toRotationMatrix() * RIC[0].transpose();  // 得到Rwi，world对应着枢纽帧
            frame_it->second.T = T[i];  // 初始化不估计平移外参，实际上表示的是枢纽帧下的相机的平移，与R表示的IMU信息是不同的
            i++;
            continue;
        }
#if 0
        std::cout << "third" << std::endl;
        std::cout << std::fixed << std::setprecision(15);
        std::cout << "i = " << i << std::endl;
        std::cout << "frame_it->first = " << frame_it->first << std::endl;
        std::cout << "Headers[i].stamp.toSec() = " << Headers[i].stamp.toSec() << std::endl;
#endif
        /**
         * notes:
         *      1. 试想有这样一种场景：frame_it对应的时间戳总是大于Headers对应的时间戳，因此二者都会执行frame_it++和i++，因此二者的时间戳总是没有办法对齐，也就不会有关键帧了
         *      2. 在测试中，如果我们margin了倒数第二帧，那么就会产生时间戳不对应的情况发生，但此时由于这个判断不通过，i就会滞后frame_it，后面二者就对齐了，因此后面的图像就被认为是关键帧了
         */
        if((frame_it->first) > Headers[i].stamp.toSec())
        {
#if 0
            std::cout << "i递增了1" << std::endl;
#endif
            i++;
        }
        // 最近的KF提供一个初始值，Twc -> Tcw
        Matrix3d R_inital = (Q[i].inverse()).toRotationMatrix();
        Vector3d P_inital = - R_inital * T[i];
        // eigen -> cv
        cv::eigen2cv(R_inital, tmp_r);
        cv::Rodrigues(tmp_r, rvec);
        cv::eigen2cv(P_inital, t);

        frame_it->second.is_key_frame = false;
        vector<cv::Point3f> pts_3_vector;
        vector<cv::Point2f> pts_2_vector;
        // 遍历这一帧对应的特征点，获得其对应的3维点（如果有的话），用于后面的pnp求解
        for (auto &id_pts : frame_it->second.points)
        {
            int feature_id = id_pts.first;
            // 由于是单目，这里id_pts.second大小就是1
            for (auto &i_p : id_pts.second)
            {
                it = sfm_tracked_points.find(feature_id);
                if(it != sfm_tracked_points.end())  // 有对应的三角化出来的3d点
                {
                    Vector3d world_pts = it->second;    // 地图点的世界坐标
                    cv::Point3f pts_3(world_pts(0), world_pts(1), world_pts(2));
                    pts_3_vector.push_back(pts_3);
                    // 去畸变的归一化相机系坐标
                    Vector2d img_pts = i_p.second.head<2>();
                    cv::Point2f pts_2(img_pts(0), img_pts(1));
                    pts_2_vector.push_back(pts_2);
                }
            }
        }
        cv::Mat K = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);     
        if(pts_3_vector.size() < 6)  // 匹配的数量少于6
        {
            cout << "pts_3_vector size " << pts_3_vector.size() << endl;
            ROS_DEBUG("Not enough points for solve pnp !");
            return false;
        }
        // 依然是调用opencv求解pnp接口：使用优化的方式，并使用初始的外参
        if (! cv::solvePnP(pts_3_vector, pts_2_vector, K, D, rvec, t, 1))
        {
            ROS_DEBUG("solve pnp fail!");
            return false;
        }
        // cv -> eigen,同时Tcw -> Twc
        cv::Rodrigues(rvec, r);
        MatrixXd R_pnp,tmp_R_pnp;
        cv::cv2eigen(r, tmp_R_pnp);
        R_pnp = tmp_R_pnp.transpose();
        MatrixXd T_pnp;
        cv::cv2eigen(t, T_pnp);
        T_pnp = R_pnp * (-T_pnp);
        // Twc -> Twi
        // 由于尺度未恢复，因此平移暂时不转到imu系
        frame_it->second.R = R_pnp * RIC[0].transpose();
        frame_it->second.T = T_pnp;
    }
    // 到此就求解出用来做视觉惯性对齐的所有视觉帧的位姿
    // Step 4 视觉惯性对齐
    if (visualInitialAlign())
        return true;
    else
    {
        ROS_INFO("misalign visual structure with IMU");
        return false;
    }

}

bool Estimator::visualInitialAlign()
{
    TicToc t_g;
    VectorXd x;
    //solve scale
    bool result = VisualIMUAlignment(all_image_frame, Bgs, g, x);
    if(!result)
    {
        ROS_DEBUG("solve g failed!");
        return false;
    }

    // change state
    // 首先把对齐后KF的位姿附给滑窗中的值，Rwi twi，到目前位置，world对应着枢纽帧
    for (int i = 0; i <= frame_count; i++)
    {
        Matrix3d Ri = all_image_frame[Headers[i].stamp.toSec()].R;  // Rc0_bi，表示枢纽帧下的IMU位姿
        Vector3d Pi = all_image_frame[Headers[i].stamp.toSec()].T;  // 这里的T还是枢纽帧下面的相机平移，跟R是在枢纽帧下的IMU旋转是不同的
        Ps[i] = Pi;
        Rs[i] = Ri;
        // 在initialStructure中已经将其设置为关键帧了
        all_image_frame[Headers[i].stamp.toSec()].is_key_frame = true;
    }

#if 0
    for (auto &it_per_id : f_manager.feature) {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
//        std::cout << "it_per_id.used_num = " << it_per_id.used_num << std::endl;
//        std::cout << "it_per_id.start_frame = " << it_per_id.start_frame << ", window_size - 2 = " << WINDOW_SIZE -2 << std::endl;
        if (std::abs(it_per_id.estimated_depth + 1) > 1e-6)
            std::cout << "这个特征点的深度为：" << it_per_id.estimated_depth << std::endl;
        std::cout << std::endl;
    }
#endif

    // 以下操作的目的是将所有观测数目在2以上的特征点的深度全部赋值为-1，后面会做多帧的三角化
    // 经测试，此时特征点管理器中的特征点的深度值都是-1，还没有被赋值
    VectorXd dep = f_manager.getDepthVector();  // 根据有效特征点数初始化这个动态向量
    for (int i = 0; i < dep.size(); i++)
        dep[i] = -1;    // 深度预设都是-1
    f_manager.clearDepth(dep);  // 特征管理器把所有的特征点逆深度也设置为-1

    //triangulate on cam pose , no tic
    Vector3d TIC_TMP[NUM_OF_CAM];
    for(int i = 0; i < NUM_OF_CAM; i++)
        TIC_TMP[i].setZero();
    ric[0] = RIC[0];
    f_manager.setRic(ric);
    // 多约束三角化所有的特征点，注意，仍带是尺度模糊的
#if 0
    std::cout << "tic[0] = " << TIC_TMP[0] << std::endl;
#endif
    // xc's todo: 在这里传入的平移为零向量，三角化的结果准确吗?：这里的Ps实际上还是枢纽帧下的相机平移，因此实际上不应该传入平移，也可以说应该传入零向量
    // 在单目sfm中确实已经两帧三角化点了，但是并没有赋值给Ps
    f_manager.triangulate(Ps, &(TIC_TMP[0]), &(RIC[0]));

    double s = (x.tail<1>())(0);
#if 0
    std::cout << "***************************尺度s = " << s << std::endl;
#endif

    // 将滑窗中的预积分重新计算
    // 需要注意的是，在陀螺仪零偏初始值计算中，对all_image_frame中的预积分（tmp_pre_integration）重新计算了，这里是对pre_integrations中的预积分重新计算
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        pre_integrations[i]->repropagate(Vector3d::Zero(), Bgs[i]);
    }
    // 下面开始把所有的状态对齐到第0帧的imu坐标系
    for (int i = frame_count; i >= 0; i--)
        /**
         * s * Pc0_ci - s * Pc0_bi = Rc0_b0 * Pb_c，见论文公式6，c0代表枢纽帧的相机系
         *
         * s * Ps[i] - Rs[i] * TIC[0] = s * Pc0_bi
         * s * Ps[0] - Rs[0] * TIC[0] = s * Pc0_b0
         *
         * Ps[i]代表枢纽帧下b0到bi的向量，并且尺度已经恢复了
         */
        Ps[i] = s * Ps[i] - Rs[i] * TIC[0] - (s * Ps[0] - Rs[0] * TIC[0]);
    int kv = -1;
    map<double, ImageFrame>::iterator frame_i;
    // 把求解出来KF的速度赋给滑窗中
    for (frame_i = all_image_frame.begin(); frame_i != all_image_frame.end(); frame_i++)
    {
        if(frame_i->second.is_key_frame)
        {
            kv++;
            // 当时求得速度是第k帧IMU坐标系下面的速度，现在转到world系
            // second.R目前为止表示的是Twi，world是枢纽帧的相机系
            // Rw_bk * Vbk_bk = Rw_bk，注意这里的速度仍然是枢纽帧下面的，还没有转换到第0帧IMU坐标系下
            // xc's todo: Vs保存的是哪个坐标系下面的速度，他跟上面的Ps的坐标系是不同的？：后面会进一步更新
            Vs[kv] = frame_i->second.R * x.segment<3>(kv * 3);
        }
    }
    // 把尺度模糊的3d点恢复到真实尺度下
    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        // s是枢纽帧坐标系下的尺度因子，也就是相机系到真实世界的尺度变换；上面三角化得到的深度也是基于start_frame计算的，二者必须对应
        // 枢纽帧坐标系是相机系，实际上对于任意的相机系，其尺度因子与枢纽帧的尺度因子是一样的
        it_per_id.estimated_depth *= s;
    }
    // 所有的P V Q全部对齐到第0帧的，同时和对齐到重力方向
    Matrix3d R0 = Utility::g2R(g);  // g是枢纽帧下的重力方向，得到枢纽帧到某个坐标系的变换（实际上是一个Z轴与重力重合的坐标系，并且补偿了yaw）
    // Rw_c0 * Rc0_b0 = Rw_b0，也就是第0帧IMU坐标系到世界系，这个世界系也就是某个Z轴与重力重合的坐标系（并且补偿yaw角）
    double yaw = Utility::R2ypr(R0 * Rs[0]).x();
#if 0
    std::cout << "Utility::R2ypr(R0 * Rs[0]).x() = " << Utility::R2ypr(R0 * Rs[0]).x() << std::endl;
    std::cout << "Utility::R2ypr(R0 * Rs[0]).y() = " << Utility::R2ypr(R0 * Rs[0]).y() << std::endl;
    std::cout << "Utility::R2ypr(R0 * Rs[0]).z() = " << Utility::R2ypr(R0 * Rs[0]).z() << std::endl;
    std::cout << "R0 * g = " << R0 * g << std::endl;
#endif
    /**
     * notes:
     *      1. 首先，bi并不是重力与Z轴重合的坐标系，应该说每一个imu都对应着一个重力与Z轴重合的坐标系
     *      2. 本身，从world→b0需要经过欧拉角-yaw, -pitch, -roll
     *      3. Eigen::Vector3d{-yaw, 0, 0} * R0表示将枢纽帧变换到world，然后在经过欧拉角-yaw得到的坐标系
     *      4. world是一个重力与Z轴重合的坐标系，仅仅变动yaw角度，并不是以b0为世界系，这里的新的世界系遵循上述操作产生的，其Z轴与重力方向重合
     */
    // xc's todo: 为什么要补偿yaw角？只要是一个重力与Z轴重合的坐标系不就可以了吗
    R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;  // R0表示RnewWorld_c0
    g = R0 * g;
#if 0
    std::cout << "g = " << g << std::endl;
#endif
    //Matrix3d rot_diff = R0 * Rs[0].transpose();
    // 在这里构建的世界系的原点在b0的原点，因此只有旋转而没有平移
    Matrix3d rot_diff = R0;
    for (int i = 0; i <= frame_count; i++)
    {
        // 在此之前，Ps、Rs和Vs均是枢纽帧坐标系下的物理量，将其变换到上面计算得到的世界系下面
        Ps[i] = rot_diff * Ps[i];  // 这里仅仅只是旋转了向量，而没有平移，说明新的世界系的原点与b0重合
        Rs[i] = rot_diff * Rs[i];   // 全部对齐到重力下，同时yaw角对齐到第一帧
        Vs[i] = rot_diff * Vs[i];
    }
    // g变成了标准重力，Rs[0]的位姿并不是单位帧，一定要牢记，担其yaw角为0
    ROS_DEBUG_STREAM("g0     " << g.transpose());
    ROS_DEBUG_STREAM("my R0  " << Utility::R2ypr(Rs[0]).transpose());
#if 0
    // yaw角为0
    std::cout << "Utility::R2ypr(Rs[0]).transpose() = " << Utility::R2ypr(Rs[0]).transpose() << std::endl;
#endif

    return true;
}

/**
 * @brief 寻找滑窗内一个帧作为枢纽帧，要求和最后一帧既有足够的共视也要有足够的视差
 *        这样其他帧都对齐到这个枢纽帧上
 *        得到T_l_last
 * @param[in] relative_R 
 * @param[in] relative_T 
 * @param[in] l 
 * @return true 
 * @return false 
 */

bool Estimator::relativePose(Matrix3d &relative_R, Vector3d &relative_T, int &l)
{
    // find previous frame which contians enough correspondance and parallex with newest frame
    // 优先从最前面开始
    /**
     * notes:
     *      1. 我们希望两帧离的足够远，使得有足够的平移
     *      2. 我们又不希望离的过远，避免没有足够的匹配点
     */
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        vector<pair<Vector3d, Vector3d>> corres;
        corres = f_manager.getCorresponding(i, WINDOW_SIZE);  // 最后一帧是固定的，选择的是开始的一帧
        // 要求共视的特征点足够多
        if (corres.size() > 20)  // 这里的条件是希望足够接近，从而有足够的匹配点
        {
            double sum_parallax = 0;
            double average_parallax;
            for (int j = 0; j < int(corres.size()); j++)
            {
                Vector2d pts_0(corres[j].first(0), corres[j].first(1));
                Vector2d pts_1(corres[j].second(0), corres[j].second(1));
                double parallax = (pts_0 - pts_1).norm();   // 计算了视差
                sum_parallax = sum_parallax + parallax;

            }
            // 计算每个特征点的平均视差：注意这个视差是去畸变的归一化相机系下，需要变换到图像坐标系下
            average_parallax = 1.0 * sum_parallax / int(corres.size());
            // 有足够的视差在通过本质矩阵恢复第i帧和最后一帧之间的 R t T_i_last
            if(average_parallax * 460 > 30 && m_estimator.solveRelativeRT(corres, relative_R, relative_T))
            {
                l = i;
                ROS_DEBUG("average_parallax %f choose l %d and newest frame to triangulate the whole structure", average_parallax * 460, l);
                return true;
            }
        }
    }
    return false;
}

void Estimator::solveOdometry()
{
    // 保证滑窗中帧数满了
    if (frame_count < WINDOW_SIZE)
        return;
    // 其次要求初始化完成
    if (solver_flag == NON_LINEAR)
    {
        TicToc t_tri;
        // 先把应该三角化但是没有三角化的特征点三角化
        // xc's todo: 为什么还有应该初始化而没有初始化的特征点：初始化结束后执行这个函数并没有需要三角化的地图点了
        // 初始化结束后进入这个函数会有需要三角化点地图点
        f_manager.triangulate(Ps, tic, ric);
        ROS_DEBUG("triangulation costs %f", t_tri.toc());
        optimization();
    }
}

/**
 * @brief 由于ceres的参数块都是double数组，因此这里把参数块从eigen的表示转成double数组
 * 
 */
void Estimator::vector2double()  // 赋初值
{
    /**
     * 同一种类型的状态量放在一个二维double中，每一个二维double的元素都是一个待优化的状态量
     */
    // KF的位姿
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        para_Pose[i][0] = Ps[i].x();
        para_Pose[i][1] = Ps[i].y();
        para_Pose[i][2] = Ps[i].z();
        Quaterniond q{Rs[i]};
        para_Pose[i][3] = q.x();
        para_Pose[i][4] = q.y();
        para_Pose[i][5] = q.z();
        para_Pose[i][6] = q.w();

        para_SpeedBias[i][0] = Vs[i].x();
        para_SpeedBias[i][1] = Vs[i].y();
        para_SpeedBias[i][2] = Vs[i].z();

        para_SpeedBias[i][3] = Bas[i].x();
        para_SpeedBias[i][4] = Bas[i].y();
        para_SpeedBias[i][5] = Bas[i].z();

        para_SpeedBias[i][6] = Bgs[i].x();
        para_SpeedBias[i][7] = Bgs[i].y();
        para_SpeedBias[i][8] = Bgs[i].z();
    }
    // 外参
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        para_Ex_Pose[i][0] = tic[i].x();
        para_Ex_Pose[i][1] = tic[i].y();
        para_Ex_Pose[i][2] = tic[i].z();
        Quaterniond q{ric[i]};
        para_Ex_Pose[i][3] = q.x();
        para_Ex_Pose[i][4] = q.y();
        para_Ex_Pose[i][5] = q.z();
        para_Ex_Pose[i][6] = q.w();
    }
    /**
     * 优化的是特征点的start_frame下的逆深度，这样做的理由是：
     *      1. 从理论上，三维点投影到相机应该是完全与像素重合，也就是没有重投影误差，但实际上是不可能的
     *      2. 这里认为在start_frame下是没有重投影误差的，也就是在地图点一定在初次观测到的start_frame的射线上，这在理论上也是合理的，这个特征点在start_frame上通过特征提取获得的，不是跟踪得到的
     *      3. 后面对这个特征点的观测是通过光流跟踪得到的，准确度就下降了，这个意思就是说，光流也是有误差的，地图点投影到这些帧上就不一定与光流跟踪的像素点重合了
     *      4. 所以优化start_frame下的逆深度是合理的，而且这种优化方式限制的地图点的运动范围，使得优化的地图点不会大幅偏离真值
     *
     * 这样做的好处是：
     *      1. 从优化的变量是三维变成了一维，有利于提高实时性
     *      2. 优化的是逆深度，对于非常远的点也会有比较好的鲁棒性，不至于出现数值非常大导致的数值精度问题；而一般不会有距离相机非常近的点
     *
     * 拓展：
     *      1. 这种思路同样可以用于ORB-SLAM3等特征点法，因为在特征点匹配过程中也会存在误匹配，使用这种方式仅仅优化逆深度也能够根据重投影误差删除错误的匹配
     */
    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < f_manager.getFeatureCount(); i++)
        para_Feature[i][0] = dep(i);
    // 传感器时间同步
    if (ESTIMATE_TD)
        para_Td[0][0] = td;
}

/**
 * @brief double -> eigen 同时fix第一帧的yaw和平移，固定了四自由度的零空间
 * 
 */
void Estimator::double2vector()
{
    // xc's todo: 窗口第一帧的yaw角和位移为什么不应该更新？
    /**
     * 以下阐述yaw角的波动和平移的增减对系统没有影响，也就是损失函数不发生变化
     * 
     * 对于视觉重投影的残差有：e = pcj / pcj(2) - obs
     *                     pcj = Rbc.t * (Rwbj.t * (Rwbi * (Rbc * pci / λ + tbc) + twbi) - Rwbj.t * twbj) - Rbc.t * tbc
     *      如果对于R，P(这里实际上对应twbi和twbj)都右乘一个delta_R，那么残差e不会发生改变；如果对于P增减一个delta_P，残差e也不会发生改变
     *      其实就是说，整个视觉系统可以整体旋转和平移，系统自身的残差不会发生变化
     * 
     * 对于IMU预积分的残差有：
     *      1. e1 = log(Ri.t * Rj)
     *      2. e2 = Ri.t * (Vj - Vi + gw * delta_t)
     *      3. e3 = Ri.t * (Pj - Pi - Vi * delta_t + 0.5 * gw * delta_t^2)
     *
     *      如果对于R，P，V都右乘一个R(delta_yaw)，那么残差e1、e2和e3都不会发生改变；如果对于P增减一个delta_P，残差e1、e2和e3也都不会发生改变
     *      对于e1，上面的阐述很容易理解，对于e2和e3，只要证明R(-delta_yaw) * gw = gw即可，而实际上这是成立的
     *
     * 对于边缘化的残差有：边缘化的残差实际上来源于视觉重投影残差和IMU预积分的残差；注意这里的边缘化的残差和jacobian并没有严格的按照FEJ来推导
     *      边缘化的目的实际上就是为了限制窗口中关键帧yaw和position的任意变动，比如在ORB-SLAM中限定第一帧不动，也就是将不可观的pose变的可观了（提供了一个世界系）；
     *      如果我们按照严格的FEJ来推导残差及其雅克比，那么理论上这里并不需要进行补偿；然而由于没有按照严格的FEJ，不可观的yaw和position依旧未知，
     *      并且由于变量在边缘化约束中按照FEJ的方式计算jacobian，而在视觉重投影误差和IMU预积分误差中没有按照FEJ的方式计算jacobian，因此yaw会发生
     *      较大的变动（添加边缘化约束与不添加边缘化约束相比，不添加边缘化的约束很快就收敛了，状态量没有较大的变化）
     *
     * 这里可以这样理解：边缘化的残差本来是用来约束不可观的yaw和position的（当然也会对pitch和roll产生约束），但由于没有严格按照FEJ，导致我们实际上没有限制住yaw和position，只能将其补偿回优化前的状态，
     * 以此来保证轨迹的连续性，这只是一种折中的方式，这种方式也就意味着窗口第一帧的yaw和position是不变的，这也是不是很准确的；
     * 所以，实际上，这里的边缘化的约束实际上只是约束了pitch和roll，对yaw和position的约束会被补偿回去
     * 
     * 如果不添加边缘化的约束，整个系统有4自由度不可观，并且也会丧失窗口外的观测对窗口内的状态的约束，这样做是不可取的；由于采用的是迭代优化，因此，
     * 实际上，不添加边缘化的约束的话，系统是比较快的就收敛了，并且各个状态量的变动比较小，当然，这也取决于状态量的初值
     *
     * 在视觉惯性系统中，不可观测的自由度为4，由于重力的存在，使得pitch和roll变的可以观测了；系统的不可观测状态量是4个---yaw, x, y, z，改变这4个状态量，残差不会发生改变
     *
     * https://github.com/HKUST-Aerial-Robotics/VINS-Mono/issues/18
     */

    // 取出优化前的第一帧的位姿
    // 初始化之前，其yaw角为0，这也可以从视觉惯性对齐中获取世界系的逻辑中得知；初始化之后，窗口的第一帧不再是系统的第一帧了，因此yaw不为0
    Vector3d origin_R0 = Utility::R2ypr(Rs[0]);
    Vector3d origin_P0 = Ps[0];
#if 0
    std::cout << "优化更新前第一帧的yaw为：" << origin_R0.x() << std::endl;
    std::cout << "优化更新前第一帧的平移为：" << origin_P0.transpose() << std::endl;
#endif

    /**
     * notes: 以下操作的目的是为了维持窗口第一帧的yaw角不变，第一帧的位移不变
     *
     * 值得注意的是，即使在初始化之前，也没有固定住某一帧的位姿，也就是说位姿可以任意漂移和旋转；然后再通过yaw角补偿变换一下
     */

    if (failure_occur)  // 初始值为false，如果failureDetection检测到失败，就会变成true
    {
        origin_R0 = Utility::R2ypr(last_R0);
        origin_P0 = last_P0;
        failure_occur = 0;
    }
    // 优化后的第一帧的位姿
    Vector3d origin_R00 = Utility::R2ypr(Quaterniond(para_Pose[0][6],
                                                      para_Pose[0][3],
                                                      para_Pose[0][4],
                                                      para_Pose[0][5]).toRotationMatrix());
    // yaw角差
    double y_diff = origin_R0.x() - origin_R00.x();
#if 0
    std::cout << "origin: x = " << origin_R00.x() << ", y = " << origin_R00.y() << ", z = " << origin_R00.z() << std::endl;
    std::cout << "opt: x = " << origin_R0.x() << ", y = " << origin_R0.y() << ", z = " << origin_R0.z() << std::endl;
#endif
    /**
     * Rs[0] = Ry2Rp2Rr2, Ropt = Ry1Rp1Rr1
     * y_diff = y2 - y1
     * rot_diff = R(y2 - y1)
     *
     * rot_diff * Ropt = R(y2 - y1) * Ry1Rp1Rr1 = Ry2Rp1Rr1，因此窗口第一帧的yaw角度并没有发生变化，其他帧的位姿做相同的变换
     */
    Matrix3d rot_diff = Utility::ypr2R(Vector3d(y_diff, 0, 0));
    // 接近万象节死锁的问题 https://blog.csdn.net/AndrewFan/article/details/60981437
    if (abs(abs(origin_R0.y()) - 90) < 1.0 || abs(abs(origin_R00.y()) - 90) < 1.0)
    {
        // xc's todo: 万向节死锁的时候不信任优化结果
        ROS_DEBUG("euler singular point!");
        rot_diff = Rs[0] * Quaterniond(para_Pose[0][6],
                                       para_Pose[0][3],
                                       para_Pose[0][4],
                                       para_Pose[0][5]).toRotationMatrix().transpose();
        /**
         * rot_diff = Rs[0] * Ropt.t
         * rot_diff * Ropt = Rs[0] * Ropt.t * Ropt = Rs[0]
         * 在万向节死锁的时候，直接将窗口第一帧的位姿恢复，其他帧的位姿做相同的变换
         */
    }

    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        // 保持第1帧的yaw不变
        Rs[i] = rot_diff * Quaterniond(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]).normalized().toRotationMatrix();
        // 保持第1帧的位移不变：pi - p0 = R * (pi - p0)，对向量执行相同的变换
        Ps[i] = rot_diff * Vector3d(para_Pose[i][0] - para_Pose[0][0],
                                para_Pose[i][1] - para_Pose[0][1],
                                para_Pose[i][2] - para_Pose[0][2]) + origin_P0;

        Vs[i] = rot_diff * Vector3d(para_SpeedBias[i][0],
                                    para_SpeedBias[i][1],
                                    para_SpeedBias[i][2]);

        Bas[i] = Vector3d(para_SpeedBias[i][3],
                          para_SpeedBias[i][4],
                          para_SpeedBias[i][5]);

        Bgs[i] = Vector3d(para_SpeedBias[i][6],
                          para_SpeedBias[i][7],
                          para_SpeedBias[i][8]);
    }
#if 0
    Vector3d R0_temp = Utility::R2ypr(Rs[0]);
    std::cout << "优化更新后的第一帧yaw为：" << R0_temp.x() << std::endl;
    std::cout << "优化更新前第一帧的平移为：" << Ps[0].transpose() << std::endl;
#endif
    // notes: 经过测试后，优化更新前后的第一帧的yaw是一样的，优化更新前后的第一帧的位移是一样的，这与理论是一样的；
    // ORB-SLAM中使用共视图，将局部地图之外的关键帧的位姿固定从而限定了局部地图中关键帧的优化，否则的话，每次优化都可以整体任意偏移和旋转
    // xc's todo: 在这里是如何限制窗口内的帧的位姿的变动的？在这里使用了FEJ来保证对窗口内关键帧的约束，但是又没有严格按照FEJ来执行，因此对yaw和position进行补偿，保证轨迹的平滑
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = Vector3d(para_Ex_Pose[i][0],
                          para_Ex_Pose[i][1],
                          para_Ex_Pose[i][2]);
        ric[i] = Quaterniond(para_Ex_Pose[i][6],
                             para_Ex_Pose[i][3],
                             para_Ex_Pose[i][4],
                             para_Ex_Pose[i][5]).toRotationMatrix();
    }
    // 重新设置各个特征点的逆深度
    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < f_manager.getFeatureCount(); i++)
        dep(i) = para_Feature[i][0];
    f_manager.setDepth(dep);
    if (ESTIMATE_TD)
        td = para_Td[0][0];

    // relative info between two loop frame
#if 0
    std::cout << "relocalization_info = " << relocalization_info << std::endl;
#endif
    if(relocalization_info) // 类似进行一个调整
    {
        /**
         * 在vins_estimator中，我们并没有使用闭环帧的任何信息，只是根据其与窗口中的某一关键帧的共视来构建视觉重投影误差，进而优化得到闭环帧在VIO窗口中的位姿
         * 既然它可以被认为是窗口中的某一帧，那么其yaw和position也应该像窗口中的其他关键帧一样被补偿
         * 这个位姿当然不会被直接使用
         */
        Matrix3d relo_r;
        Vector3d relo_t;
        relo_r = rot_diff * Quaterniond(relo_Pose[6], relo_Pose[3], relo_Pose[4], relo_Pose[5]).normalized().toRotationMatrix();
        relo_t = rot_diff * Vector3d(relo_Pose[0] - para_Pose[0][0],
                                     relo_Pose[1] - para_Pose[0][1],
                                     relo_Pose[2] - para_Pose[0][2]) + origin_P0;
        double drift_correct_yaw;
        drift_correct_yaw = Utility::R2ypr(prev_relo_r).x() - Utility::R2ypr(relo_r).x();
#if 0
        double drift_correct_pitch = Utility::R2ypr(prev_relo_r).y() - Utility::R2ypr(relo_r).y();
        double drift_correct_roll = Utility::R2ypr(prev_relo_r).z() - Utility::R2ypr(relo_r).z();
        std::cout << "yaw = " << drift_correct_yaw << ", pitch = " << drift_correct_pitch << ", roll = " << drift_correct_roll << std::endl;

        Eigen::Matrix3d drift_R = prev_relo_r * relo_r.inverse();
        double drift_correct_yaw1 = Utility::R2ypr(drift_R).x();
        double drift_correct_pitch1 = Utility::R2ypr(drift_R).y();
        double drift_correct_roll1 = Utility::R2ypr(drift_R).z();
        std::cout << "yaw1 = " << drift_correct_yaw1 << ", pitch1 = " << drift_correct_pitch1 << ", roll1 = " << drift_correct_roll1 << std::endl;
#endif
        /**
         * 这里的drift_correct_yaw实际上描述的是当前VIO坐标系到闭环帧坐标系的旋转，这个量描述了系统旋转的漂移程度
         *
         * relo_r = Ry1 * Rp1 * Rr1, prev_relo_r = Ry2 * Rp2 * Rr2
         * 由于一个VIO系统是yaw和position不可观，因此yaw存在漂移（就像尺度存在漂移一样），这里操作的目的就是修正VIO系统的yaw角度和position的漂移
         *
         * drift_correct_r = R(y2 - y1)，那么drift_correct_r * relo_r = Ry2 * Rp1 * Rr1
         *
         * (drift_correct_r    drift_correct_t)  *  (relo_t)  = prev_relo_t；实际上就是将当前VIO的位姿进行修正
         * (       0                   1      )     (   1  )
         */
        // 注意即使到了这里，我们并没有直接修改当前VIO或者闭环帧的位姿，只是在求解一些相对位姿
        drift_correct_r = Utility::ypr2R(Vector3d(drift_correct_yaw, 0, 0));
        drift_correct_t = prev_relo_t - drift_correct_r * relo_t;   
        // T_loop_w * T_w_cur = T_loop_cur
        relo_relative_t = relo_r.transpose() * (Ps[relo_frame_local_index] - relo_t);  // 当前关键帧到闭环帧的平移
        relo_relative_q = relo_r.transpose() * Rs[relo_frame_local_index];  // 当前关键帧到闭环帧的旋转
        // 返回的角度在-0.5pi~0.5pi
        relo_relative_yaw = Utility::normalizeAngle(Utility::R2ypr(Rs[relo_frame_local_index]).x() - Utility::R2ypr(relo_r).x());
        //cout << "vins relo " << endl;
        //cout << "vins relative_t " << relo_relative_t.transpose() << endl;
        //cout << "vins relative_yaw " <<relo_relative_yaw << endl;
        relocalization_info = 0;  // 复位为0，等待下一次回环

    }
}

bool Estimator::failureDetection()
{
    if (f_manager.last_track_num < 2)   // 跟踪上一帧的地图点的数目
    {
        ROS_INFO(" little feature %d", f_manager.last_track_num);
        //return true;  // 注意没有返回true，有可能矫正回来
    }
    if (Bas[WINDOW_SIZE].norm() > 2.5)  // 加速度零偏是否正常
    {
        ROS_INFO(" big IMU acc bias estimation %f", Bas[WINDOW_SIZE].norm());
        return true;
    }
    if (Bgs[WINDOW_SIZE].norm() > 1.0)  // 陀螺仪零偏是否正常
    {
        ROS_INFO(" big IMU gyr bias estimation %f", Bgs[WINDOW_SIZE].norm());
        return true;
    }
    /*
    if (tic(0) > 1)
    {
        ROS_INFO(" big extri param estimation %d", tic(0) > 1);
        return true;
    }
    */
    Vector3d tmp_P = Ps[WINDOW_SIZE];
    if ((tmp_P - last_P).norm() > 5)    // 两帧之间运动是否过大
    {
        ROS_INFO(" big translation");
        return true;
    }
    if (abs(tmp_P.z() - last_P.z()) > 1)    // 重力方向运动是否过大，通常情况下不会发生
    {
        ROS_INFO(" big z translation");
        return true; 
    }
    Matrix3d tmp_R = Rs[WINDOW_SIZE];
    Matrix3d delta_R = tmp_R.transpose() * last_R;
    Quaterniond delta_Q(delta_R);
    double delta_angle;
    delta_angle = acos(delta_Q.w()) * 2.0 / 3.14 * 180.0;  // 使用的是相对旋转对应的轴角来判断
    if (delta_angle > 50)   // 两帧姿态变化是否过大
    {
        ROS_INFO(" big delta_angle ");
        //return true;  // 注意没有返回true，有可能矫正回来
    }
    return false;
}

/**
 * @brief 进行非线性优化
 * 
 */
void Estimator::optimization()
{
    /**
     * 优化函数包括的残差类型有：
     *      1. 视觉重投影误差（分为优化时间同步和不优化时间同步）
     *      2. 视觉与IMU预积分构建的残差
     *      3. 边缘化获取的约束
     *      4. 回环检测的约束
     */
    // 借助ceres进行非线性优化
    ceres::Problem problem;
    ceres::LossFunction *loss_function;
    //loss_function = new ceres::HuberLoss(1.0);
    // Cauchy: ρ(s) = log(1 + s)   →   ρ'(s) = (1 + s)^(-1)   →   ρ''(s) = -(1 + s)^(-2)
    loss_function = new ceres::CauchyLoss(1.0);
    // Step 1 定义待优化的参数块，类似g2o的顶点
    // 参数块 1： 滑窗中位姿包括位置和姿态，共11帧
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        // 由于姿态不满足正常的加法，也就是李群上没有加法，因此需要自己定义他的加法
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Pose[i], SIZE_POSE, local_parameterization);
        problem.AddParameterBlock(para_SpeedBias[i], SIZE_SPEEDBIAS);  // 默认的加法
    }
    // 参数块 2： 相机imu间的外参
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Ex_Pose[i], SIZE_POSE, local_parameterization);
        if (!ESTIMATE_EXTRINSIC)
        {
            ROS_DEBUG("fix extinsic param");
            // 如果不需要优化外参就设置为fix
            problem.SetParameterBlockConstant(para_Ex_Pose[i]);
        }
        else
            ROS_DEBUG("estimate extinsic param");
    }
    // 传感器的时间同步
    if (ESTIMATE_TD)
    {
        problem.AddParameterBlock(para_Td[0], 1);
        //problem.SetParameterBlockConstant(para_Td[0]);
    }
    // 实际上还有地图点，其实平凡的（不需要定义“加法”）参数块不需要调用AddParameterBlock，增加残差块接口时会自动绑定
    TicToc t_whole, t_prepare;
    // eigen -> double
    vector2double();  // 给予优化的初值
    // Step 2 通过残差约束来添加残差块，类似g2o的边
    // 上一次的边缘化结果作为这一次的先验
    /**
     * 1. 实际上，当我们将边缘化的约束注释掉的时候，优化前后的R变动非常小，这个也是非常容易理解的，每次滑窗去除最老帧，加入最新帧，整体优化幅度其实很小（系统变化很小）
     *
     * 2. 但是当我们加入边缘化的约束后，yaw角的波动变大，pitch和roll的波动依然很小；这是由于边缘化的约束添加后，系统的优化幅度会大幅增加（系统变动比较大，需要大幅波动才能再次稳定），
     * 而yaw角的波动不会影响到系统的残差，因此在yaw角上的优化幅度会非常大，反之，由于pitch和roll相比于yaw会显著影响到系统的残差，而且当前值已经使得窗口中的约束近似局部最优（也就是1中所说的优化前后变动很小），
     * 而边缘化的约束也不会大幅波动pitch和roll，因此pitch和roll的波动相比于yaw而言，小了非常多
     */
    if (last_marginalization_info)
    {
        // construct new marginlization_factor
        MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
        problem.AddResidualBlock(marginalization_factor, NULL,
                                 last_marginalization_parameter_blocks);
    }
    // imu预积分的约束
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        int j = i + 1;
        // 时间过长这个约束就不可信了
        if (pre_integrations[j]->sum_dt > 10.0)
            continue;
        IMUFactor* imu_factor = new IMUFactor(pre_integrations[j]);  // costFunction，注意j是从1开始的
        problem.AddResidualBlock(imu_factor, NULL, para_Pose[i], para_SpeedBias[i], para_Pose[j], para_SpeedBias[j]);
    }
    int f_m_cnt = 0;
    int feature_index = -1;
    // 视觉重投影的约束
    // 遍历每一个特征点
    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        // 进行特征点有效性的检查
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
 
        ++feature_index;
        // 第一个观测到这个特征点的帧idx
        int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
        // 特征点在第一个帧下的归一化相机系坐标
        Vector3d pts_i = it_per_id.feature_per_frame[0].point;
        // 遍历看到这个特征点的所有KF
        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            if (imu_i == imu_j) // 自己跟自己不能形成重投影
            {
                continue;
            }
#if 0
            std::cout << "it_per_id.feature_per_frame[0].cur_td = " << it_per_id.feature_per_frame[0].cur_td << std::endl;
            std::cout << "it_per_frame.cur_td = " << it_per_frame.cur_td << std::endl;
#endif
            // 取出另一帧的归一化相机坐标
            Vector3d pts_j = it_per_frame.point;
            // 带有时间延时的是另一种形式
            if (ESTIMATE_TD)
            {
                    // 这里的td都是从配置文件中读取的，实际上都是一个td
                    ProjectionTdFactor *f_td = new ProjectionTdFactor(pts_i, pts_j, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocity,
                                                                     it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td,
                                                                     it_per_id.feature_per_frame[0].uv.y(), it_per_frame.uv.y());
                    problem.AddResidualBlock(f_td, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index], para_Td[0]);
                    /*
                    double **para = new double *[5];
                    para[0] = para_Pose[imu_i];
                    para[1] = para_Pose[imu_j];
                    para[2] = para_Ex_Pose[0];
                    para[3] = para_Feature[feature_index];
                    para[4] = para_Td[0];
                    f_td->check(para);
                    */
            }
            else
            {
                ProjectionFactor *f = new ProjectionFactor(pts_i, pts_j);   // 构造函数就是同一个特征点在不同帧的观测
                // 约束的变量是该特征点的第一个观测帧以及其他一个观测帧，加上外参和特征点逆深度
                problem.AddResidualBlock(f, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index]);
            }
            f_m_cnt++;
        }
    }

    ROS_DEBUG("visual measurement count: %d", f_m_cnt);
    ROS_DEBUG("prepare for ceres: %f", t_prepare.toc());
    // xc's todo: 为什么优化闭环关键帧，其应该没有漂移，应该是后面的关键帧才会有漂移发生？这个优化发生在4自由度位姿图优化前还是后？
    // xc's todo: 如果优化回环关键帧，那么那些与回环关键帧构成共视的关键帧的位姿难道不应该同步优化吗？
    // 回环检测相关的约束
    if(relocalization_info)  // 找到了有效的回环信息
    {
        //printf("set relocalization factor! \n");
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(relo_Pose, SIZE_POSE, local_parameterization);    // 需要优化的回环帧位姿
        int retrive_feature_index = 0;
        int feature_index = -1;
        // 遍历现有地图点
        for (auto &it_per_id : f_manager.feature)
        {
            it_per_id.used_num = it_per_id.feature_per_frame.size();
            if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
                continue;
            ++feature_index;
            int start = it_per_id.start_frame;
            // 闭环关键帧与窗口中的relo_frame_local_index对应的关键帧构成回环，如果relo_frame_local_index对应的关键帧看不到，那么闭环关键帧也一定看不到
            if(start <= relo_frame_local_index)   
            {   
                // 寻找回环帧能看到的地图点
                while((int)match_points[retrive_feature_index].z() < it_per_id.feature_id)
                {
                    retrive_feature_index++;
                }
                // 这个地图点也能被回环帧看到
                if((int)match_points[retrive_feature_index].z() == it_per_id.feature_id)
                {
                    // 构建一个重投影约束，这个地图点的起始帧和该回环帧之间
                    Vector3d pts_j = Vector3d(match_points[retrive_feature_index].x(), match_points[retrive_feature_index].y(), 1.0);
                    Vector3d pts_i = it_per_id.feature_per_frame[0].point;
                    
                    ProjectionFactor *f = new ProjectionFactor(pts_i, pts_j);
                    // 某个特征点既可以被回环关键帧看到，又可以被窗口中的某一个关键帧看到，从而构建视觉重投影误差
                    // 误匹配的可能性还是比较大的，使用核函数
                    problem.AddResidualBlock(f, loss_function, para_Pose[start], relo_Pose, para_Ex_Pose[0], para_Feature[feature_index]);
                    retrive_feature_index++;
                }     
            }
        }

    }
    // Step 3 ceres优化求解
    ceres::Solver::Options options;

    options.linear_solver_type = ceres::DENSE_SCHUR;
    //options.num_threads = 2;
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = NUM_ITERATIONS;
    //options.use_explicit_schur_complement = true;
    //options.minimizer_progress_to_stdout = true;
    //options.use_nonmonotonic_steps = true;
    if (marginalization_flag == MARGIN_OLD)
        // 下面的边缘化老的操作比较多，因此给他优化时间就少一些
        options.max_solver_time_in_seconds = SOLVER_TIME * 4.0 / 5.0;
    else
        options.max_solver_time_in_seconds = SOLVER_TIME;
    TicToc t_solver;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);  // ceres优化求解
    //cout << summary.BriefReport() << endl;
    ROS_DEBUG("Iterations : %d", static_cast<int>(summary.iterations.size()));
    ROS_DEBUG("solver costs: %f", t_solver.toc());
    // 把优化后double -> eigen
    double2vector();  // 将优化后的状态量保存，某些状态量的优化结果会被进一步更新
    // Step 4 边缘化
    TicToc t_whole_marginalization;
    /**
     * 需要注意的是：上一次的边缘化形成的约束一定对本次待边缘化的变量形成约束，才会将上一次边缘化形成的约束加入到整个marginalization_info中，否则是不会添加的
     * 这里所说的对应着边缘化的时候的B矩阵不应该有全0行，如果有的话，直接删除即可，不需要这个约束
     *
     * 对整个边缘化过程的剖析：注意一定要出现old，否则一直没有办法形成last_marginalization_info
     *      1. 为了产生last_marginalization_info，从某一帧开始一定要出现old；如果刚开始是new，那么实际上窗口的第一帧就是整个系统的第一帧，系统一直在合并新来的帧与窗口最后一帧；
     *      2. 一旦某一帧开始出现old，last_marginalization_info被创建出来，边缘化后对窗口内的某些参数块形成了约束，在后面会使用到
     *      3. 如果是从old到old，那么需要添加的约束有窗口第一帧与第二帧的IMU预积分约束，第一帧就看到的地图点的视觉重投影约束，以及上一次边缘化形成的约束对本次第一帧的位姿、速度和零偏有约束
     *      4. 如果是old到new，那么只需要添加上一次边缘化形成的约束对本次窗口最后一帧的位姿有约束
     *          4.1 为什么本次不添加被最后一帧看到的重投影约束：预积分会合并，但是窗口最后一帧的视觉信息（观测到的地图点）是直接删除的，因此最后一帧并没有与窗口中的其他帧产生视觉重投影的约束了，因此也就没有必要构建与最后一帧相关的重投影约束了
     *              而在边缘化第1帧和第2帧的时候，第1帧的观测依然存在，其依然应该对窗口内的帧产生约束，故而添加重投影约束；从因子图的角度理解，我们要边缘化的约束是待边缘化变量与非边缘化变量有边连接，在这里如果将最后一帧的地图点删除，连接也就断了
     *              比较奇怪的是，后面保留了前面的预积分约束中与最后一帧有关的约束，这里可以这样理解：直接删除位姿和地图点的话，就不需要做边缘化构建约束了，但是前面边缘化的结果依然在，如果可以的话，依然需要构建约束
     *          4.2 为什么不添加窗口倒数第二帧与窗口最后一帧的预积分约束：最后一帧与新进来的帧合并成为最后一帧，倒数第二帧与窗口最后一帧的预积分约束包含在合并的预积分中了，并没有被删除；
     *              不像第1帧和第二帧直接的预积分约束，由于第1帧不在窗口中而消失了
     *
     *      5. 如果是new到new，那么依然只需要添加上一次边缘化形成的约束对本次窗口最后一帧的位姿有约束
     *      6. 如果是new到old，那么需要添加的约束有窗口第一帧与第二帧的IMU预积分约束，第一帧就看到的地图点的视觉重投影约束，以及上一次边缘化形成的约束对本次第一帧的位姿、速度和零偏有约束
     */
    if (marginalization_flag == MARGIN_OLD)
    {
        // 一个用来边缘化操作的对象
        MarginalizationInfo *marginalization_info = new MarginalizationInfo();
        // 这里类似手写高斯牛顿，因此也需要都转成double数组
        vector2double();  // 将经过double2vector后的结果重新赋值给待优化的状态量作为初值
        // 关于边缘化有几点注意的地方
        // 1、找到需要边缘化的参数块，这里是地图点，第0帧位姿，第0帧速度零偏
        // 2、找到构造高斯牛顿下降时跟这些待边缘化相关的参数块有关的残差约束，那就是预积分约束，重投影约束，以及上一次边缘化约束
        // 3、这些约束连接的参数块中，不需要被边缘化的参数块，就是被提供先验约束的部分，也就是滑窗中剩下的位姿和速度零偏

        // 上一次的边缘化结果
        if (last_marginalization_info)
        {
            vector<int> drop_set;
            // last_marginalization_parameter_blocks是上一次边缘化对哪些当前参数块有约束
            for (int i = 0; i < static_cast<int>(last_marginalization_parameter_blocks.size()); i++)
            {
                // 上一次边缘化中留下来的参数块包含了本次滑窗起始的位姿、零偏和速度，这个约束也需要考虑进去
                if (last_marginalization_parameter_blocks[i] == para_Pose[0] ||
                    last_marginalization_parameter_blocks[i] == para_SpeedBias[0])
                    drop_set.push_back(i);
            }
            // 处理方式和其他残差块相同
            // construct new marginlization_factor
            MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(marginalization_factor, NULL,
                                                                           last_marginalization_parameter_blocks,
                                                                           drop_set);

            marginalization_info->addResidualBlockInfo(residual_block_info);
        }
        /**
         * notes:
         *      1. drop set：表示parameter block中下标为drop set的元素是待边缘化的
         */
        // 只有第1个预积分和待边缘化参数块相连
        {
            if (pre_integrations[1]->sum_dt < 10.0)  // 时间太长，预积分的结果就不可信了
            {
                // 跟构建ceres约束问题一样，这里也需要得到残差和雅克比
                IMUFactor* imu_factor = new IMUFactor(pre_integrations[1]);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(imu_factor, NULL,
                                                                           vector<double *>{para_Pose[0], para_SpeedBias[0], para_Pose[1], para_SpeedBias[1]},
                                                                           vector<int>{0, 1});  // 这里就是第0和1个参数块是需要被边缘化的，也就是para_Pose[0], para_SpeedBias[0]
                marginalization_info->addResidualBlockInfo(residual_block_info);
            }
        }
        // 遍历视觉重投影的约束
        {
            int feature_index = -1;
            for (auto &it_per_id : f_manager.feature)
            {
                it_per_id.used_num = it_per_id.feature_per_frame.size();
                if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
                    continue;

                ++feature_index;

                int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
                // notes：只找能被第0帧看到的特征点，这些点才是要被边缘化的，因为第0帧是要被边缘化掉的
                if (imu_i != 0)
                    continue;

                Vector3d pts_i = it_per_id.feature_per_frame[0].point;
                // 遍历看到这个特征点的所有KF，通过这个特征点，建立和第0帧的约束
                for (auto &it_per_frame : it_per_id.feature_per_frame)
                {
                    imu_j++;
                    if (imu_i == imu_j)
                        continue;

                    Vector3d pts_j = it_per_frame.point;
                    // 根据是否约束延时确定残差阵
                    if (ESTIMATE_TD)
                    {
                        /**
                         * 在视觉重投影的边缘化的过程中，如果我们只是边缘化位姿，那么会导致fill-in现象，使得整个矩阵不再是稀疏矩阵
                         * 在这里，我们将位姿和地图点都边缘化了，这样就可以保证整个矩阵的稀疏性了，详细的原理性推导可以参考SLAM14讲第10章的部分
                         */
                        ProjectionTdFactor *f_td = new ProjectionTdFactor(pts_i, pts_j, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocity,
                                                                          it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td,
                                                                          it_per_id.feature_per_frame[0].uv.y(), it_per_frame.uv.y());
                        ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f_td, loss_function,
                                                                                        vector<double *>{para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index], para_Td[0]},
                                                                                        vector<int>{0, 3});  // 这里第0个和第3个参数块是需要被边缘化的，也就是para_Pose[imu_i]和para_Feature[feature_index]
                        marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                    else
                    {
                        ProjectionFactor *f = new ProjectionFactor(pts_i, pts_j);
                        ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                                                                                       vector<double *>{para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index]},
                                                                                       vector<int>{0, 3});  // 这里第0个和第3个参数块是需要被边缘化的，也就是para_Pose[imu_i]和para_Feature[feature_index]
                        marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                }
            }
        }
        // 所有的残差块都收集好了
        TicToc t_pre_margin;
        // 进行预处理：计算每个因子的残差和jacobian，根据Robust Kernel对残差也jacobian进行缩放，对所有的参数块进行备份
        marginalization_info->preMarginalize();
        ROS_DEBUG("pre marginalization %f ms", t_pre_margin.toc());
        
        TicToc t_margin;
        // 边缘化操作
        marginalization_info->marginalize();
        ROS_DEBUG("marginalization %f ms", t_margin.toc());
        // 即将滑窗，因此记录新地址对应的老地址
        std::unordered_map<long, double *> addr_shift;  // 随滑窗而变动的参数块的地址才需要变动
        for (int i = 1; i <= WINDOW_SIZE; i++)  // 注意这里是从1开始递增
        {
            // 位姿和速度都要滑窗移动
            addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];
            addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
        }
        // 外参和时间延时不变
        for (int i = 0; i < NUM_OF_CAM; i++)
            addr_shift[reinterpret_cast<long>(para_Ex_Pose[i])] = para_Ex_Pose[i];
        if (ESTIMATE_TD)
        {
            addr_shift[reinterpret_cast<long>(para_Td[0])] = para_Td[0];
        }
        // parameter_blocks实际上就是addr_shift中的没有边缘化的参数块的地址
        vector<double *> parameter_blocks = marginalization_info->getParameterBlocks(addr_shift);

        if (last_marginalization_info)
            delete last_marginalization_info;
        last_marginalization_info = marginalization_info;   // 本次边缘化的所有信息
        last_marginalization_parameter_blocks = parameter_blocks;   // 代表该次边缘化对没有边缘化的参数块形成约束，这些参数块在滑窗之后的地址
        
    }
    else    // 边缘化倒数第二帧
    {
        /**
         * 如果一开始就是margine old并且一直都是margin old，那么实际上边缘化部分的代码没有起到任何作用，一直都没有last_marginalization_info创建成功
         */

        // 要求有上一次边缘化的结果同时，即将被margin掉的(窗口最后一帧的信息)在上一次边缘化后的约束的参数块中
        // 预积分结果合并，因此只有位姿margin掉
        if (last_marginalization_info &&
            std::count(std::begin(last_marginalization_parameter_blocks), std::end(last_marginalization_parameter_blocks), para_Pose[WINDOW_SIZE - 1]))
        {

            MarginalizationInfo *marginalization_info = new MarginalizationInfo();
            vector2double();  // 将经过double2vector后的结果重新赋值给待优化的状态量作为初值
            // xc's todo: 重复判断？
            if (last_marginalization_info)
            {
                vector<int> drop_set;
                for (int i = 0; i < static_cast<int>(last_marginalization_parameter_blocks.size()); i++)
                {
                    // 速度零偏只会margin第1个，不可能出现倒数第二个；para_SpeedBias出现在预积分的约束中
                    ROS_ASSERT(last_marginalization_parameter_blocks[i] != para_SpeedBias[WINDOW_SIZE - 1]);
                    // 这种case只会margin掉倒数第二个位姿，应该可以在找到之后直接break的，最多只会找到1个
                    if (last_marginalization_parameter_blocks[i] == para_Pose[WINDOW_SIZE - 1])
                        drop_set.push_back(i);
                }
                // construct new marginlization_factor
                // 这里只会更新一下margin factor
                MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(marginalization_factor, NULL,
                                                                               last_marginalization_parameter_blocks,
                                                                               drop_set);

                marginalization_info->addResidualBlockInfo(residual_block_info);
            }
            // 这里的操作如出一辙
            TicToc t_pre_margin;
            ROS_DEBUG("begin marginalization");
            marginalization_info->preMarginalize();
            ROS_DEBUG("end pre marginalization, %f ms", t_pre_margin.toc());

            TicToc t_margin;
            ROS_DEBUG("begin marginalization");
            marginalization_info->marginalize();
            ROS_DEBUG("end marginalization, %f ms", t_margin.toc());
            
            std::unordered_map<long, double *> addr_shift;
            for (int i = 0; i <= WINDOW_SIZE; i++)
            {
                if (i == WINDOW_SIZE - 1)
                    continue;
                else if (i == WINDOW_SIZE)  // 滑窗，最新帧成为次新帧
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];
                    addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
                }
                else    // 其他不变
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i];
                    addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i];
                }
            }
            for (int i = 0; i < NUM_OF_CAM; i++)
                addr_shift[reinterpret_cast<long>(para_Ex_Pose[i])] = para_Ex_Pose[i];
            if (ESTIMATE_TD)
            {
                addr_shift[reinterpret_cast<long>(para_Td[0])] = para_Td[0];
            }
            
            vector<double *> parameter_blocks = marginalization_info->getParameterBlocks(addr_shift);
            if (last_marginalization_info)
                delete last_marginalization_info;
            last_marginalization_info = marginalization_info;
            last_marginalization_parameter_blocks = parameter_blocks;
            
        }
    }
    ROS_DEBUG("whole marginalization costs: %f", t_whole_marginalization.toc());
    
    ROS_DEBUG("whole time for ceres: %f", t_whole.toc());
}

// 滑动窗口 
void Estimator::slideWindow()
{
    TicToc t_margin;
    // 根据边缘化种类的不同，进行滑窗的方式也不同
#if 0
    std::cout << "fist" << std::endl;
    std::cout << "marginalization_flag = " << marginalization_flag << std::endl;
    std::cout << std::fixed << std::setprecision(15);
    for (int i = 0; i < WINDOW_SIZE; i++) {
        std::cout << "headers[i].t = " << Headers[i].stamp.toSec() << std::endl;
    }
#endif
    if (marginalization_flag == MARGIN_OLD)
    {
        double t_0 = Headers[0].stamp.toSec();
        back_R0 = Rs[0];
        back_P0 = Ps[0];
        // 必须是填满了滑窗才可以
        if (frame_count == WINDOW_SIZE)
        {
            // 一帧一帧交换过去：注意这里仅仅只是更新了预积分相关的信息，all_image_frame的信息在后面更新
            for (int i = 0; i < WINDOW_SIZE; i++)
            {
                Rs[i].swap(Rs[i + 1]);

                std::swap(pre_integrations[i], pre_integrations[i + 1]);

                dt_buf[i].swap(dt_buf[i + 1]);
                linear_acceleration_buf[i].swap(linear_acceleration_buf[i + 1]);
                angular_velocity_buf[i].swap(angular_velocity_buf[i + 1]);

                Headers[i] = Headers[i + 1];
                Ps[i].swap(Ps[i + 1]);
                Vs[i].swap(Vs[i + 1]);
                Bas[i].swap(Bas[i + 1]);
                Bgs[i].swap(Bgs[i + 1]);
            }
#if 0
            std::cout << std::fixed << std::setprecision(15);
            std::cout << "Ps[WINDOW_SIZE]" << Ps[WINDOW_SIZE] << ", Ps[WINDOW_SIZE - 1] = " << Ps[WINDOW_SIZE - 1] << std::endl;
            std::cout << "Bas[WINDOW_SIZE] = " << Bas[WINDOW_SIZE] << ", Bas[WINDOW_SIZE - 1] = " << Bas[WINDOW_SIZE - 1] << std::endl;
#endif
            // 最后一帧的状态量赋上当前值，最为初始值
            // 注意上面是swap，而不是直接赋值覆盖，因此才有下面的操作，上面的for循环之后，类似Ps[WINDOW_SIZE]已经发生了变化
            Headers[WINDOW_SIZE] = Headers[WINDOW_SIZE - 1];
            Ps[WINDOW_SIZE] = Ps[WINDOW_SIZE - 1];
            Vs[WINDOW_SIZE] = Vs[WINDOW_SIZE - 1];
            Rs[WINDOW_SIZE] = Rs[WINDOW_SIZE - 1];
            Bas[WINDOW_SIZE] = Bas[WINDOW_SIZE - 1];
            Bgs[WINDOW_SIZE] = Bgs[WINDOW_SIZE - 1];
            // 预积分量就得置零
            delete pre_integrations[WINDOW_SIZE];
#if 0
            std::cout << std::fixed << std::setprecision(15);
            std::cout << "acc_0 = " << acc_0 << ", gyr_0 = " << gyr_0 << std::endl;
#endif
            pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};
            // buffer清空，等待新的数据来填
            dt_buf[WINDOW_SIZE].clear();
            linear_acceleration_buf[WINDOW_SIZE].clear();
            angular_velocity_buf[WINDOW_SIZE].clear();
            // 清空all_image_frame最老帧之前的状态
            if (true || solver_flag == INITIAL)  // 始终都会执行
            {
                // 预积分量是堆上的空间，因此需要手动释放
                map<double, ImageFrame>::iterator it_0;
                it_0 = all_image_frame.find(t_0);
                delete it_0->second.pre_integration;
                it_0->second.pre_integration = nullptr;
#if 0
                std::cout << "开始执行判断" << std::endl;
                std::cout << "all_image_frame.size = " << all_image_frame.size() << std::endl;
                for (map<double, ImageFrame>::iterator it = all_image_frame.begin(); it != all_image_frame.end(); ++it) {
                    std::cout << "it = " << it->first << std::endl;
                    std::cout << "it_0 = " << it_0->first << std::endl;
                    if (it != it_0) {
                        std::cout << "发生不相等了" << std::endl;
                    }
                }
#endif
                // 在初始化阶段，it_0对应的一般是all_image_frame的第一帧，因此这个for循环一般进不去就直接结束了
                for (map<double, ImageFrame>::iterator it = all_image_frame.begin(); it != it_0; ++it)
                {
#if 0
                    std::cout << "释放预积分" << std::endl;
#endif
                    if (it->second.pre_integration)
                        delete it->second.pre_integration;
                    it->second.pre_integration = NULL;
                }
                // 释放完空间之后再erase
                all_image_frame.erase(all_image_frame.begin(), it_0);
                all_image_frame.erase(t_0);

            }
            slideWindowOld();
        }
    }
    else
    {
        if (frame_count == WINDOW_SIZE)
        {
            // 将最后两个预积分观测合并成一个
            // 需要注意的是：此时frame_count已经等于window_size了，也就是已经不在窗口内了，是窗口右边的第一帧
            for (unsigned int i = 0; i < dt_buf[frame_count].size(); i++)
            {
                double tmp_dt = dt_buf[frame_count][i];
                Vector3d tmp_linear_acceleration = linear_acceleration_buf[frame_count][i];
                Vector3d tmp_angular_velocity = angular_velocity_buf[frame_count][i];

                pre_integrations[frame_count - 1]->push_back(tmp_dt, tmp_linear_acceleration, tmp_angular_velocity);

                dt_buf[frame_count - 1].push_back(tmp_dt);
                linear_acceleration_buf[frame_count - 1].push_back(tmp_linear_acceleration);
                angular_velocity_buf[frame_count - 1].push_back(tmp_angular_velocity);
            }
            // 简单的滑窗交换：将变量与合并后的预积分相对应
            Headers[frame_count - 1] = Headers[frame_count];  // 预积分已经合并了，因此其结束时刻也要更新，下面的变量也是一样
            Ps[frame_count - 1] = Ps[frame_count];
            Vs[frame_count - 1] = Vs[frame_count];
            Rs[frame_count - 1] = Rs[frame_count];
            Bas[frame_count - 1] = Bas[frame_count];
            Bgs[frame_count - 1] = Bgs[frame_count];
            // reset最新预积分量
            delete pre_integrations[WINDOW_SIZE];
#if 0
            std::cout << std::fixed << std::setprecision(15);
            std::cout << "acc_0 = " << acc_0 << ", gyr_0 = " << gyr_0 << std::endl;
#endif
            // 注意，这里的acc_0和gyr_0并不是当前预积分的第一个读数，而是已经更新了，是当前预积分的最后一个读数，是下一个预积分的第一个读数
            // 但从初始化的逻辑看，这里不new也没关系，在图像处理函数中也会new出来
            pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};
            // clear相关buffer
            dt_buf[WINDOW_SIZE].clear();
            linear_acceleration_buf[WINDOW_SIZE].clear();
            angular_velocity_buf[WINDOW_SIZE].clear();

            slideWindowNew();
        }
    }
#if 0
    std::cout << "second" << std::endl;
    std::cout << "marginalization_flag = " << marginalization_flag << std::endl;
    std::cout << std::fixed << std::setprecision(15);
    for (int i = 0; i < WINDOW_SIZE; i++) {
        std::cout << "headers[i].t = " << Headers[i].stamp.toSec() << std::endl;
    }
#endif
}

// real marginalization is removed in solve_ceres()
// 对被移除的倒数第二帧的地图点进行处理
void Estimator::slideWindowNew()
{
    sum_of_front++;
    f_manager.removeFront(frame_count);
}
// real marginalization is removed in solve_ceres()
// 由于地图点是绑定在第一个看见它的位姿上的，因此需要对被移除的帧看见的地图点进行解绑，以及每个地图点的首个观测帧id减1
void Estimator::slideWindowOld()
{
    sum_of_back++;

    bool shift_depth = solver_flag == NON_LINEAR ? true : false;
    if (shift_depth)    // 如果初始化过了
    {
        Matrix3d R0, R1;
        Vector3d P0, P1;
        // back_R0 back_P0是被移除的帧的位姿
        R0 = back_R0 * ric[0];  // 被移除的相机的姿态
        R1 = Rs[0] * ric[0];    // 当前最老的相机姿态
        P0 = back_P0 + back_R0 * tic[0];    // 被移除的相机的位置
        P1 = Ps[0] + Rs[0] * tic[0];    // 当前最老的相机位置
        // 下面要做的事情把被移除帧看见地图点的管理权交给当前的最老帧
        f_manager.removeBackShiftDepth(R0, P0, R1, P1);
    }
    else
        f_manager.removeBack();
}

/**
 * @brief 接受回环帧的消息
 * 
 * @param[in] _frame_stamp 闭环关键帧对应的当前关键帧的时间戳
 * @param[in] _frame_index 获取闭环的当前关键帧的索引
 * @param[in] _match_points 闭环关键帧的与当前关键帧匹配的归一化相机系的特征点
 * @param[in] _relo_t 闭环关键帧的位置
 * @param[in] _relo_r 闭环关键帧的旋转
 */
void Estimator::setReloFrame(double _frame_stamp, int _frame_index, vector<Vector3d> &_match_points, Vector3d _relo_t, Matrix3d _relo_r)
{
    relo_frame_stamp = _frame_stamp;
    relo_frame_index = _frame_index;
    match_points.clear();
    match_points = _match_points;
    prev_relo_t = _relo_t;
    prev_relo_r = _relo_r;
    // 在滑窗中寻找当前关键帧，因为VIO送给回环结点的是倒数第三帧，因此，很有可能这个当前帧还在滑窗里
    for(int i = 0; i < WINDOW_SIZE; i++)
    {
        if(relo_frame_stamp == Headers[i].stamp.toSec())
        {
            relo_frame_local_index = i; // 对应滑窗中的第i帧
            relocalization_info = 1;    // 这是一个有效的回环信息
            for (int j = 0; j < SIZE_POSE; j++)
                // 注意这里直接赋值，对应的是同一个pose，因此必须逐个元素赋值
                relo_Pose[j] = para_Pose[i][j]; // 借助VIO优化回环帧位姿，初值先设为当前帧位姿
        }
    }
}

