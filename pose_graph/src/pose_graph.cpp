#include "pose_graph.h"

PoseGraph::PoseGraph()
{
    // 用于rviz的可视化
    posegraph_visualization = new CameraPoseVisualization(1.0, 0.0, 1.0, 1.0);
    posegraph_visualization->setScale(0.1);
    posegraph_visualization->setLineWidth(0.01);
    // 生成一个线程，该线程用于进行4自由度全局优化
	t_optimization = std::thread(&PoseGraph::optimize4DoF, this);
    // 初始化一些变量
    earliest_loop_index = -1;
    t_drift = Eigen::Vector3d(0, 0, 0);
    yaw_drift = 0;
    r_drift = Eigen::Matrix3d::Identity();
    w_t_vio = Eigen::Vector3d(0, 0, 0);
    w_r_vio = Eigen::Matrix3d::Identity();
    global_index = 0;
    sequence_cnt = 0;
    sequence_loop.push_back(0);
    base_sequence = 1;

}

PoseGraph::~PoseGraph()
{
	t_optimization.join();
}
/**
 * @brief 注册一些发布publisher
 * 
 * @param[in] n 
 */
void PoseGraph::registerPub(ros::NodeHandle &n)
{
    pub_pg_path = n.advertise<nav_msgs::Path>("pose_graph_path", 1000);
    pub_base_path = n.advertise<nav_msgs::Path>("base_path", 1000);
    pub_pose_graph = n.advertise<visualization_msgs::MarkerArray>("pose_graph", 1000);
    for (int i = 1; i < 10; i++)
        pub_path[i] = n.advertise<nav_msgs::Path>("path_" + to_string(i), 1000);
}

// 加载二进制词袋库
void PoseGraph::loadVocabulary(std::string voc_path)
{
    voc = new BriefVocabulary(voc_path);
    db.setVocabulary(*voc, false, 0);
}

// flag_detect_loop传入的始终是1
void PoseGraph::addKeyFrame(KeyFrame* cur_kf, bool flag_detect_loop)
{
    //shift to base frame
    Vector3d vio_P_cur;
    Matrix3d vio_R_cur;
    // 发生了sequence的跳变，可以理解成不是一个地图，这里的sequence_cnt递增之后，就到了一个地图了
    if (sequence_cnt != cur_kf->sequence)
    {
        sequence_cnt++;
        // 复位一些变量
        sequence_loop.push_back(0);
        w_t_vio = Eigen::Vector3d(0, 0, 0);
        w_r_vio = Eigen::Matrix3d::Identity();
        m_drift.lock();
        t_drift = Eigen::Vector3d(0, 0, 0);
        r_drift = Eigen::Matrix3d::Identity();
        m_drift.unlock();
    }
    // 更新一下VIO位姿
    cur_kf->getVioPose(vio_P_cur, vio_R_cur);   // 得到VIO节点的位姿
    // 注意这里的w_r_vio和w_t_vio只有在地图合并的时候才会赋值，这里实际上是将当前关键帧的位姿变化到上一次闭环的世界系下，也就是统一到一个坐标系下面
    // 如果只是闭环，也就是本来就是在一个坐标系下，那么这个w_r_vio始终为单位阵，w_t_vio始终为0向量
    vio_P_cur = w_r_vio * vio_P_cur + w_t_vio;
    vio_R_cur = w_r_vio *  vio_R_cur;
    cur_kf->updateVioPose(vio_P_cur, vio_R_cur);    // 更新VIO位姿
    cur_kf->index = global_index;   // 赋值索引
    global_index++;
	int loop_index = -1;
    if (flag_detect_loop)  // 始终检测
    {
        TicToc tmp_t;
        loop_index = detectLoop(cur_kf, cur_kf->index);
    }
    else
    {
        addKeyFrameIntoVoc(cur_kf);  // 将新提的FAST角点的信息加入到DBOW的数据库中
    }
	if (loop_index != -1)   // 代表找到了有效的回环帧
	{
        //printf(" %d detect loop with %d \n", cur_kf->index, loop_index);
        KeyFrame* old_kf = getKeyFrame(loop_index); // 得到回环帧的指针

        if (cur_kf->findConnection(old_kf)) // 如果确定两者回环
        {
            // 更新最早回环帧，用来确定全局优化的范围
            if (earliest_loop_index > loop_index || earliest_loop_index == -1)
                earliest_loop_index = loop_index;

            Vector3d w_P_old, w_P_cur, vio_P_cur;
            Matrix3d w_R_old, w_R_cur, vio_R_cur;
            old_kf->getVioPose(w_P_old, w_R_old);
            cur_kf->getVioPose(vio_P_cur, vio_R_cur);

            Vector3d relative_t;
            Quaterniond relative_q;
            // T_old_cur
            relative_t = cur_kf->getLoopRelativeT();
            relative_q = (cur_kf->getLoopRelativeQ()).toRotationMatrix();
            // T_w_old * T_old_cur = T_w_cur,这就是回环矫正后当前帧的位姿
            /**
             * relative_pose中，闭环帧的位姿是通过当前关键帧的地图点计算得到，因此就得到了没有漂移的当前关键帧到闭环关键帧时间的位姿变化
             * 再通过没有漂移的闭环关键帧的位姿就可以计算得到当前关键帧的位姿
             * Tw_old * Told_cur即为所得
             */
            w_P_cur = w_R_old * relative_t + w_P_old;
            w_R_cur = w_R_old * relative_q;
            double shift_yaw;
            Matrix3d shift_r;
            Vector3d shift_t; 
            // 回环矫正前的位姿认为是T_w'_cur
            /**
             * T_w_cur * T_cur_w' = T_w_w'，计算的是矫正前的位姿到矫正后的位姿的变换
             *
             * 这里的shift_yaw只是在两个地图合并的时候才会使用，两个地图的Z轴都是重力方向，因此两个地图的世界系实际上只会有yaw角的变动
             */
            // xc's todo: 为什么这里不使用相对位姿的yaw角度，而是分解求解yaw角度来做减法；实际上两者是不一样的，除非两个坐标系真的只有yaw角变化，而pitch和roll完全没有变化，此时两者才有可能是一样的
            // 在这里已经假设两个世界系只有yaw角的变动了，因此就直接这样计算了，如果不满足这个假设也没办法了，毕竟系统运行到这里，也没办法继续优化世界系的Z轴方向了
            // 实际上选择相对旋转的yaw角与这里的shift_yaw一般差距很小，因此误差应该也不会大
            shift_yaw = Utility::R2ypr(w_R_cur).x() - Utility::R2ypr(vio_R_cur).x();
            shift_r = Utility::ypr2R(Vector3d(shift_yaw, 0, 0));
            shift_t = w_P_cur - w_R_cur * vio_R_cur.transpose() * vio_P_cur; 
            // shift vio pose of whole sequence to the world frame
            // 如果这两个不是同一个sequence，并且当前sequence没有跟之前合并，多地图合并
            if (old_kf->sequence != cur_kf->sequence && sequence_loop[cur_kf->sequence] == 0)
            {
                /**
                 * 从打印的结果看，两个世界系的pitch和roll并不相等，这是不正常的，这说明系统的重力方向并没有优化到准确方向，当然这个差值一般相对于yaw很小
                 * 有没有可能根据两个世界系的pitch和roll不同来进一步优化系统的世界系Z轴方向
                 */
#if 0
                double shift_pitch = Utility::R2ypr(w_R_cur).y() - Utility::R2ypr(vio_R_cur).y();
                double shift_roll = shift_yaw = Utility::R2ypr(w_R_cur).z() - Utility::R2ypr(vio_R_cur).z();
                std::cout << "shift_yaw = " << shift_yaw << ", shift_pitch = " << shift_pitch << ", shift_roll = " << shift_roll << std::endl;

                Eigen::Matrix3d Rw_w = w_R_cur * vio_R_cur.inverse();
                double shift_yaw1 = Utility::R2ypr(Rw_w).x();
                double shift_pitch1 = Utility::R2ypr(Rw_w).y();
                double shift_roll1 = Utility::R2ypr(Rw_w).z();
                std::cout << "shift_yaw1 = " << shift_yaw1 << ", shift_pitch = " << shift_pitch1 << ", shift_roll = " << shift_roll1 << std::endl;
#endif
                // 以下实际上是将当前序列的位姿转换到闭环关键帧所在序列的世界系下
                w_r_vio = shift_r;
                w_t_vio = shift_t;
                // T_w_w' * T_w'_cur
                vio_P_cur = w_r_vio * vio_P_cur + w_t_vio;
                vio_R_cur = w_r_vio *  vio_R_cur;
                cur_kf->updateVioPose(vio_P_cur, vio_R_cur);    // 更新当前帧位姿
                list<KeyFrame*>::iterator it = keyframelist.begin();
                // notes: 没有将这个序列的地图点也更新一下
                // 同时把这个序列当前帧之间的位姿都更新过来
                for (; it != keyframelist.end(); it++)   
                {
                    if((*it)->sequence == cur_kf->sequence)
                    {
                        Vector3d vio_P_cur;
                        Matrix3d vio_R_cur;
                        (*it)->getVioPose(vio_P_cur, vio_R_cur);
                        vio_P_cur = w_r_vio * vio_P_cur + w_t_vio;
                        vio_R_cur = w_r_vio *  vio_R_cur;
                        (*it)->updateVioPose(vio_P_cur, vio_R_cur);
                    }
                }
                // 代表这个序列已经跟之前的序列合并过了，这里实现的也就是一个序列的合并
                sequence_loop[cur_kf->sequence] = 1;
            }
            // notes：合并之后还是会做4自由度位姿图优化，这个跟ORB-SLAM3是不一样的，在ORB-SLAM3中，地图合并之后就结束了
            m_optimize_buf.lock();
            optimize_buf.push(cur_kf->index);   // 相当于通知4dof优化线程开始干活
            m_optimize_buf.unlock();
        }
	}
	m_keyframelist.lock();
    Vector3d P;
    Matrix3d R;
    cur_kf->getVioPose(P, R);
    // 根据全局优化后进行位姿调整，注意这里是r_drift和t_drift
#if 0
    if (loop_index != -1)
        std::cout << "在addKeyFrame中更新当前关键帧的位姿" << std::endl;
#endif
    // 4DOF优化耗时比较大，因此r_drift与t_drift并不总是在其中计算，更有可能在updateKeyFrameLoop得到计算；测试表明，更新经常性发生在这里
    P = r_drift * P + t_drift;
    R = r_drift * R;
    // 注意在这里更新的是世界系的位姿，不是VIO的位姿，后面如果要做4DOF优化的时候，在其中使用的是VIO的位姿
    cur_kf->updatePose(P, R);
    // 下面是可视化部分
    Quaterniond Q{R};
    geometry_msgs::PoseStamped pose_stamped;
    pose_stamped.header.stamp = ros::Time(cur_kf->time_stamp);
    pose_stamped.header.frame_id = "world";
    pose_stamped.pose.position.x = P.x() + VISUALIZATION_SHIFT_X;
    pose_stamped.pose.position.y = P.y() + VISUALIZATION_SHIFT_Y;
    pose_stamped.pose.position.z = P.z();
    pose_stamped.pose.orientation.x = Q.x();
    pose_stamped.pose.orientation.y = Q.y();
    pose_stamped.pose.orientation.z = Q.z();
    pose_stamped.pose.orientation.w = Q.w();
    path[sequence_cnt].poses.push_back(pose_stamped);
    path[sequence_cnt].header = pose_stamped.header;

    if (SAVE_LOOP_PATH) {
        ofstream loop_path_file(VINS_RESULT_PATH, ios::app);
        loop_path_file.setf(ios::fixed, ios::floatfield);
        loop_path_file.precision(0);
        loop_path_file << cur_kf->time_stamp * 1e9 << ",";
        loop_path_file.precision(5);
        loop_path_file  << P.x() << ","
              << P.y() << ","
              << P.z() << ","
              << Q.w() << ","
              << Q.x() << ","
              << Q.y() << ","
              << Q.z() << ","
              << endl;
        loop_path_file.close();
    }
    //draw local connection
    if (SHOW_S_EDGE)
    {
        list<KeyFrame*>::reverse_iterator rit = keyframelist.rbegin();
        for (int i = 0; i < 4; i++)
        {
            if (rit == keyframelist.rend())
                break;
            Vector3d conncected_P;
            Matrix3d connected_R;
            if((*rit)->sequence == cur_kf->sequence)
            {
                (*rit)->getPose(conncected_P, connected_R);
                posegraph_visualization->add_edge(P, conncected_P);
            }
            rit++;
        }
    }
    if (SHOW_L_EDGE)
    {
        if (cur_kf->has_loop)
        {
            //printf("has loop \n");
            KeyFrame* connected_KF = getKeyFrame(cur_kf->loop_index);
            Vector3d connected_P,P0;
            Matrix3d connected_R,R0;
            connected_KF->getPose(connected_P, connected_R);
            //cur_kf->getVioPose(P0, R0);
            cur_kf->getPose(P0, R0);
            if(cur_kf->sequence > 0)
            {
                //printf("add loop into visual \n");
                posegraph_visualization->add_loopedge(P0, connected_P + Vector3d(VISUALIZATION_SHIFT_X, VISUALIZATION_SHIFT_Y, 0));
            }
            
        }
    }
    //posegraph_visualization->add_pose(P + Vector3d(VISUALIZATION_SHIFT_X, VISUALIZATION_SHIFT_Y, 0), Q);

	keyframelist.push_back(cur_kf); // 当前帧送进KF容器中
    publish();
	m_keyframelist.unlock();
}

/**
 * @brief 加载KF
 * 
 * @param[in] cur_kf 
 * @param[in] flag_detect_loop 是否做回环检测
 */
void PoseGraph::loadKeyFrame(KeyFrame* cur_kf, bool flag_detect_loop)
{
    cur_kf->index = global_index;
    global_index++;
    int loop_index = -1;
    if (flag_detect_loop)
        // 进行回环检测
       loop_index = detectLoop(cur_kf, cur_kf->index);
    else
    {
        // 把当前帧的信息加载进数据库中，用来后续进行回环检测
        addKeyFrameIntoVoc(cur_kf);
    }
    if (loop_index != -1)
    {
        printf(" %d detect loop with %d \n", cur_kf->index, loop_index);
        KeyFrame* old_kf = getKeyFrame(loop_index);
        if (cur_kf->findConnection(old_kf))
        {
            if (earliest_loop_index > loop_index || earliest_loop_index == -1)
                earliest_loop_index = loop_index;
            m_optimize_buf.lock();
            optimize_buf.push(cur_kf->index);
            m_optimize_buf.unlock();
        }
    }
    m_keyframelist.lock();
    Vector3d P;
    Matrix3d R;
    cur_kf->getPose(P, R);
    Quaterniond Q{R};
    geometry_msgs::PoseStamped pose_stamped;
    pose_stamped.header.stamp = ros::Time(cur_kf->time_stamp);
    pose_stamped.header.frame_id = "world";
    pose_stamped.pose.position.x = P.x() + VISUALIZATION_SHIFT_X;
    pose_stamped.pose.position.y = P.y() + VISUALIZATION_SHIFT_Y;
    pose_stamped.pose.position.z = P.z();
    pose_stamped.pose.orientation.x = Q.x();
    pose_stamped.pose.orientation.y = Q.y();
    pose_stamped.pose.orientation.z = Q.z();
    pose_stamped.pose.orientation.w = Q.w();
    base_path.poses.push_back(pose_stamped);
    base_path.header = pose_stamped.header;

    //draw local connection
    if (SHOW_S_EDGE)
    {
        list<KeyFrame*>::reverse_iterator rit = keyframelist.rbegin();
        for (int i = 0; i < 1; i++)
        {
            if (rit == keyframelist.rend())
                break;
            Vector3d conncected_P;
            Matrix3d connected_R;
            if((*rit)->sequence == cur_kf->sequence)
            {
                (*rit)->getPose(conncected_P, connected_R);
                posegraph_visualization->add_edge(P, conncected_P);
            }
            rit++;
        }
    }
    /*
    if (cur_kf->has_loop)
    {
        KeyFrame* connected_KF = getKeyFrame(cur_kf->loop_index);
        Vector3d connected_P;
        Matrix3d connected_R;
        connected_KF->getPose(connected_P,  connected_R);
        posegraph_visualization->add_loopedge(P, connected_P, SHIFT);
    }
    */

    keyframelist.push_back(cur_kf);
    //publish();
    m_keyframelist.unlock();
}

/**
 * @brief 根据index找到对应的KF对象的指针
 * 
 * @param[in] index 
 * @return KeyFrame* 
 */
KeyFrame* PoseGraph::getKeyFrame(int index)
{
//    unique_lock<mutex> lock(m_keyframelist);
    list<KeyFrame*>::iterator it = keyframelist.begin();
    for (; it != keyframelist.end(); it++)   
    {
        if((*it)->index == index)
            break;
    }
    if (it != keyframelist.end())
        return *it;
    else
        return NULL;
}

// 进行回环检测，寻找候选的回环帧
int PoseGraph::detectLoop(KeyFrame* keyframe, int frame_index)
{
    // put image into image_pool; for visualization
    cv::Mat compressed_image;
    if (DEBUG_IMAGE)
    {
        int feature_num = keyframe->keypoints.size();
        cv::resize(keyframe->image, compressed_image, cv::Size(376, 240));
        putText(compressed_image, "feature_num:" + to_string(feature_num), cv::Point2f(10, 10), CV_FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(255));
        image_pool[frame_index] = compressed_image;
    }
    TicToc tmp_t;
    //first query; then add this frame into database!
    QueryResults ret;
    TicToc t_query;
    // 调用词袋查询接口，查询结果是ret，最多返回4个备选KF，查找距离当前至少50帧的KF
    db.query(keyframe->brief_descriptors, ret, 4, frame_index - 50);
    //printf("query time: %f", t_query.toc());
    //cout << "Searching for Image " << frame_index << ". " << ret << endl;

    TicToc t_add;
    // 当然也会把当前帧送进数据库中，便于后续帧的查询
    db.add(keyframe->brief_descriptors);
    //printf("add feature time: %f", t_add.toc());
    // ret[0] is the nearest neighbour's score. threshold change with neighour score
    bool find_loop = false;
    cv::Mat loop_result;
    if (DEBUG_IMAGE)
    {
        loop_result = compressed_image.clone();  // 当前关键帧对应的图像
        if (ret.size() > 0)
            putText(loop_result, "neighbour score:" + to_string(ret[0].Score), cv::Point2f(10, 50), CV_FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255));
    }
    // visual loop result 
    if (DEBUG_IMAGE)
    {
        for (unsigned int i = 0; i < ret.size(); i++)  // 回环帧对应的图像
        {
            int tmp_index = ret[i].Id;
            auto it = image_pool.find(tmp_index);
            cv::Mat tmp_image = (it->second).clone();
            putText(tmp_image, "index:  " + to_string(tmp_index) + "loop score:" + to_string(ret[i].Score), cv::Point2f(10, 50), CV_FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255));
            cv::hconcat(loop_result, tmp_image, loop_result);  // 将当前关键帧及闭环帧的图像全部水平合并
        }
    }
    // notes: 最高分要大于0.05，次高分要大于0.015，就认为找到了回环
    // a good match with its nerghbour
    // ret按照得分大小降序排列的，这里确保返回的候选KF数目至少一个且得分满足要求
    if (ret.size() >= 1 &&ret[0].Score > 0.05)  // 最高分大于0.05
        // 开始遍历其他候选帧
        for (unsigned int i = 1; i < ret.size(); i++)
        {
            //if (ret[i].Score > ret[0].Score * 0.3)
            // 如果有大于一个候选帧且得分满足要求
            if (ret[i].Score > 0.015)
            {          
                find_loop = true;   // 就认为找到回环了
                int tmp_index = ret[i].Id;
                if (DEBUG_IMAGE && 0)
                {
                    auto it = image_pool.find(tmp_index);
                    cv::Mat tmp_image = (it->second).clone();
                    putText(tmp_image, "loop score:" + to_string(ret[i].Score), cv::Point2f(10, 50), CV_FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(255));
                    cv::hconcat(loop_result, tmp_image, loop_result);
                }
                // notes：添加一个break应该更符合逻辑
            }

        }
/*
    if (DEBUG_IMAGE)
    {
        cv::imshow("loop_result", loop_result);
        cv::waitKey(20);
    }
*/
    if (find_loop && frame_index > 50)  // 认为找到了回环并且当前帧是第50帧以后，也就是前面50帧不会做回环的操作
    {
        int min_index = -1;
        // 寻找得分大于0.015的idx最小的那个帧，这样可以尽可能调整更多帧的位姿
        for (unsigned int i = 0; i < ret.size(); i++)  // 注意这里认为0.015以上的得分都是有效的闭环帧，并不是直接找得分最高的候选闭环帧
        {
            if (min_index == -1 || (ret[i].Id < min_index && ret[i].Score > 0.015))
                min_index = ret[i].Id;
        }
        return min_index;   // 当前帧和min_index形成了回环
    }
    else
        return -1;  // -1代表回环失败

}

void PoseGraph::addKeyFrameIntoVoc(KeyFrame* keyframe)
{
    // put image into image_pool; for visualization
    cv::Mat compressed_image;
    if (DEBUG_IMAGE)
    {
        int feature_num = keyframe->keypoints.size();
        cv::resize(keyframe->image, compressed_image, cv::Size(376, 240));
        putText(compressed_image, "feature_num:" + to_string(feature_num), cv::Point2f(10, 10), CV_FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(255));
        image_pool[keyframe->index] = compressed_image;
    }

    db.add(keyframe->brief_descriptors);  // 仅仅添加了新提的FAST角点的描述子，并没有窗口中Harris角点的描述子
}

/**
 * @brief 如果找到回环，该线程就执行位姿优化
 * 
 */
void PoseGraph::optimize4DoF()
{
    while(true)
    {
        int cur_index = -1;
        int first_looped_index = -1;
        m_optimize_buf.lock();
        // 取出最新的形成回环的当前帧
        while(!optimize_buf.empty())
        {
            cur_index = optimize_buf.front();  // 获取的是闭环的当前关键帧的索引
            first_looped_index = earliest_loop_index;   // 找到当前关键帧的闭环关键帧的索引
            optimize_buf.pop();
        }
        m_optimize_buf.unlock();
        if (cur_index != -1)
        {
            printf("optimize pose graph \n");
            TicToc tmp_t;
            m_keyframelist.lock();
            KeyFrame* cur_kf = getKeyFrame(cur_index);  // 取出当前帧对应的KF指针

            int max_length = cur_index + 1; // 预设最大长度，总之优化帧数不可能超过这么多

            // w^t_i   w^q_i
            double t_array[max_length][3];
            Quaterniond q_array[max_length];
            double euler_array[max_length][3];
            double sequence_array[max_length];
            // 定义一个ceres优化问题，这里只优化位移和yaw角
            ceres::Problem problem;
            ceres::Solver::Options options;
            options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
            //options.minimizer_progress_to_stdout = true;
            //options.max_solver_time_in_seconds = SOLVER_TIME * 3;
            options.max_num_iterations = 5;
            ceres::Solver::Summary summary;
            ceres::LossFunction *loss_function;
            loss_function = new ceres::HuberLoss(0.1);
            //loss_function = new ceres::CauchyLoss(1.0);
            // 由于yaw不满足简单的增量更新方式，因此需要自定义其local param形式
            ceres::LocalParameterization* angle_local_parameterization =
                AngleLocalParameterization::Create();

            list<KeyFrame*>::iterator it;

            int i = 0;
            // 遍历KF的list
            for (it = keyframelist.begin(); it != keyframelist.end(); it++)
            {
                if ((*it)->index < first_looped_index)  // idx小于最早回环帧就算了，不在回环调整的范围内
                    continue;
                (*it)->local_index = i; // 这个是在本次优化中的idx
                Quaterniond tmp_q;
                Matrix3d tmp_r;
                Vector3d tmp_t;
                //  获取VIO下的位姿，注意，即使是构成闭环的当前关键帧依然选择的是VIO的位姿，注意在addKeyFrame中可能更新过世界系的位姿，但是VIO的位姿没有更新过
                (*it)->getVioPose(tmp_t, tmp_r);
                tmp_q = tmp_r;
                t_array[i][0] = tmp_t(0);
                t_array[i][1] = tmp_t(1);
                t_array[i][2] = tmp_t(2);
                q_array[i] = tmp_q;
                // 四元数转欧拉角
                Vector3d euler_angle = Utility::R2ypr(tmp_q.toRotationMatrix());
                euler_array[i][0] = euler_angle.x();
                euler_array[i][1] = euler_angle.y();
                euler_array[i][2] = euler_angle.z();

                sequence_array[i] = (*it)->sequence;
                // 只有yaw角参与优化，成为参数块
                problem.AddParameterBlock(euler_array[i], 1, angle_local_parameterization);
                problem.AddParameterBlock(t_array[i], 3);
                // 条件1：闭环关键帧的位姿不优化；条件2：加载的地图的关键帧不优化
                if ((*it)->index == first_looped_index || (*it)->sequence == 0)
                {   
                    problem.SetParameterBlockConstant(euler_array[i]);
                    problem.SetParameterBlockConstant(t_array[i]);
                }

                //add edge
                // 建立约束，每一帧和之前4帧建立约束关系，找到这一帧前面4帧，并且需要他们是一个序列中的KF
                /**
                 * 以下的约束实际上是“稳态的”，在没有其他约束的情况下，系统实际上并不需要优化就已经达到最优了；
                 * 后面会根据当前关键帧及其闭环关键帧构建约束，促使系统发生变化，开始优化，知道收敛
                 *
                 * 第一种约束：每一个关键帧与其前面4个关键帧的相对平移和yaw角度的插值，当前状态即可使得约束达到最优
                 * 第二种约束：当前关键帧与其闭环关键帧构成的相对平移和yaw角度的差值，当前状态并不能使得优化最优，相反，此时的残差应该是比较大的
                 *
                 * 需要注意的是，第二种约束的观测是在当前关键帧下面的地图点计算闭环关键帧的位姿并经过VIO中的滑窗优化得到的，是一个正确的、不含有漂移的相对约束
                 */
                for (int j = 1; j < 5; j++)
                {
                    // i从0开始递增，i- j < 0显然取不到
                  if (i - j >= 0 && sequence_array[i] == sequence_array[i-j])
                  {
#if 0
                      std::cout << "first_looped_index" << first_looped_index << std::endl;
                      std::cout << "i = " << i << ", j = " << j << ", i - j = " << i - j << std::endl;
#endif
                    // 计算T_i-j_w * T_w_i = T_i-j_i
                    Vector3d euler_conncected = Utility::R2ypr(q_array[i-j].toRotationMatrix());
                    Vector3d relative_t(t_array[i][0] - t_array[i-j][0], t_array[i][1] - t_array[i-j][1], t_array[i][2] - t_array[i-j][2]);
                    relative_t = q_array[i-j].inverse() * relative_t;  // 统一到第i帧，对应着Ti-j_i的平移部分
                    double relative_yaw = euler_array[i][0] - euler_array[i-j][0];  // 两帧之间的位姿的yaw角度的差值
                    ceres::CostFunction* cost_function = FourDOFError::Create( relative_t.x(), relative_t.y(), relative_t.z(),
                                                   relative_yaw, euler_conncected.y(), euler_conncected.z());
                    // 对i-j帧和第i帧都成约束
                    problem.AddResidualBlock(cost_function, NULL, euler_array[i-j], 
                                            t_array[i-j], 
                                            euler_array[i], 
                                            t_array[i]);
                  }
                }
                /**
                 * 我们将第一帧设为i，最后一帧设为j，并且第j帧也就是构成闭环的当前关键帧，那么
                 *  根据上面的约束有：
                 *      yaw_j - yaw_i = (yaw_j - yaw_j-1) + (yaw_j-1 - yaw_j-2) + ... + (yaw_i+1 - yaw_i)
                 *  根据下面的约束有：
                 *      yaw_j - yaw_i = (*it)->getLoopRelativeYaw()
                 *
                 *  对平移也有类似的式子，同样可以构成约束
                 *
                 * 这两个等式的坐标是相同的，右边是不同的，因此构成了图优化，推动了优化的进行
                 */

                //add loop edge
                // 如果这一帧有回环帧
                if((*it)->has_loop)
                {
                    assert((*it)->loop_index >= first_looped_index);
                    int connected_index = getKeyFrame((*it)->loop_index)->local_index;  // 得到回环帧在这次优化中的idx
                    Vector3d euler_conncected = Utility::R2ypr(q_array[connected_index].toRotationMatrix());
                    Vector3d relative_t;
                    relative_t = (*it)->getLoopRelativeT(); // 得到当前帧和回环帧的相对位姿
                    double relative_yaw = (*it)->getLoopRelativeYaw();
                    ceres::CostFunction* cost_function = FourDOFWeightError::Create( relative_t.x(), relative_t.y(), relative_t.z(),
                                                                               relative_yaw, euler_conncected.y(), euler_conncected.z());
                    problem.AddResidualBlock(cost_function, loss_function, euler_array[connected_index], 
                                                                  t_array[connected_index], 
                                                                  euler_array[i], 
                                                                  t_array[i]);
                }
                
                if ((*it)->index == cur_index)  // 到当前帧了，不会再有添加了，结束
                    break;
                i++;
            }
            m_keyframelist.unlock();

            ceres::Solve(options, &problem, &summary);
            //std::cout << summary.BriefReport() << "\n";
            
            //printf("pose optimization time: %f \n", tmp_t.toc());
            /*
            for (int j = 0 ; j < i; j++)
            {
                printf("optimize i: %d p: %f, %f, %f\n", j, t_array[j][0], t_array[j][1], t_array[j][2] );
            }
            */
#if 0
            std::cout << "在optimize4DOF中更新当前关键帧的位姿" << std::endl;
#endif
            m_keyframelist.lock();
            i = 0;
            // 将优化后的位姿恢复
            for (it = keyframelist.begin(); it != keyframelist.end(); it++)
            {
                if ((*it)->index < first_looped_index)
                    continue;
                Quaterniond tmp_q;
                tmp_q = Utility::ypr2R(Vector3d(euler_array[i][0], euler_array[i][1], euler_array[i][2]));
                Vector3d tmp_t = Vector3d(t_array[i][0], t_array[i][1], t_array[i][2]);
                Matrix3d tmp_r = tmp_q.toRotationMatrix();
                (*it)-> updatePose(tmp_t, tmp_r);   // 更新世界系的位姿，VIO的位姿没有发生变动

                if ((*it)->index == cur_index)
                    break;
                i++;
            }

            Vector3d cur_t, vio_t;
            Matrix3d cur_r, vio_r;
            cur_kf->getPose(cur_t, cur_r);  // 最新优化后的位姿
            cur_kf->getVioPose(vio_t, vio_r);   // VIO的位姿
            m_drift.lock();
            // 计算当前帧的VIO位姿和优化后位姿差
            /**
             * cur_r = Ry2 * Rp2 * Rr2, vio_r = Ry1 * Rp1 * Rr1
             * shift_r = R(y2 - y1)
             * R(y2 - y1) * vio_r = Ry2 * Rp1 * Rr1
             *
             * 如果令vio对应的世界系是w'，那么cur_r * vio_r表示w'到w的旋转；Tw_cur * Tcur_w'
             * 使用r_shift描述两个世界系之间的相对旋转，由于VIO系统，只有yaw角漂移，因此仅仅使用yaw
             * 使用t_shift描述两个世界系之间的相对平移，注意这里的计算与updateKeyFrameLoop不同，直接使用的r_shift，而不再是使用原始结果了；
             *
             * 从优化的原理上看，这里并没有优化pitch和roll，因此优化前后的pitch和roll实际上没有发生变化（这一点可以从测试上得到验证）
             *
             * 既然pitch和roll没有发生变动，因此使用原始数据计算和直接使用r_shift是一样的，这是区别于updateKeyFrameLoop的地方
             */
#if 0
            Eigen::Matrix3d drift_R = cur_r * vio_r.inverse();
            double drift_yaw = Utility::R2ypr(drift_R).x();
            double drift_pitch = Utility::R2ypr(drift_R).y();
            double drift_roll = Utility::R2ypr(drift_R).z();
            std::cout << "yaw = " << drift_yaw << ", pitch = " << drift_pitch << ", roll = " << drift_roll << std::endl;
#endif
            yaw_drift = Utility::R2ypr(cur_r).x() - Utility::R2ypr(vio_r).x();
            r_drift = Utility::ypr2R(Vector3d(yaw_drift, 0, 0));
            t_drift = cur_t - r_drift * vio_t;
            m_drift.unlock();
            //cout << "t_drift " << t_drift.transpose() << endl;
            //cout << "r_drift " << Utility::R2ypr(r_drift).transpose() << endl;
            //cout << "yaw drift " << yaw_drift << endl;

            it++;
            // 遍历当前帧之后的所有位姿，根据算得的位姿差进行调整
            /**
             * notes:
             *      1. 这里主要是多线程运行，其他线程还在产生关键帧
             *      2. 需要留意的是，其他线程产生的关键帧不能太多（距离当前关键帧不能太久），否则依然会有漂移
             */
            for (; it != keyframelist.end(); it++)
            {
                // 将构成闭环的当前关键帧之后的关键帧的位姿根据当前关键帧的漂移进行修正，这里认为他们的漂移程度与当前关键帧的漂移程度相同
                Vector3d P;
                Matrix3d R;
                (*it)->getVioPose(P, R);
                P = r_drift * P + t_drift;
                R = r_drift * R;
                (*it)->updatePose(P, R);  // 将关键帧的位姿更新到世界系，注意VIO下的位姿没有发生改变
            }
            m_keyframelist.unlock();
            // 可视化部分
            updatePath();
        }

        std::chrono::milliseconds dura(2000);
        std::this_thread::sleep_for(dura);  // 间隔两秒，不会频繁的触发
    }
}

void PoseGraph::updatePath()
{
    m_keyframelist.lock();
    list<KeyFrame*>::iterator it;
    for (int i = 1; i <= sequence_cnt; i++)
    {
        path[i].poses.clear();
    }
    base_path.poses.clear();
    posegraph_visualization->reset();

    if (SAVE_LOOP_PATH)
    {
        ofstream loop_path_file_tmp(VINS_RESULT_PATH, ios::out);
        loop_path_file_tmp.close();
    }

    for (it = keyframelist.begin(); it != keyframelist.end(); it++)
    {
        Vector3d P;
        Matrix3d R;
        (*it)->getPose(P, R);
        Quaterniond Q;
        Q = R;
//        printf("path p: %f, %f, %f\n",  P.x(),  P.z(),  P.y() );

        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header.stamp = ros::Time((*it)->time_stamp);
        pose_stamped.header.frame_id = "world";
        pose_stamped.pose.position.x = P.x() + VISUALIZATION_SHIFT_X;
        pose_stamped.pose.position.y = P.y() + VISUALIZATION_SHIFT_Y;
        pose_stamped.pose.position.z = P.z();
        pose_stamped.pose.orientation.x = Q.x();
        pose_stamped.pose.orientation.y = Q.y();
        pose_stamped.pose.orientation.z = Q.z();
        pose_stamped.pose.orientation.w = Q.w();
        if((*it)->sequence == 0)
        {
            base_path.poses.push_back(pose_stamped);
            base_path.header = pose_stamped.header;
        }
        else
        {
            path[(*it)->sequence].poses.push_back(pose_stamped);
            path[(*it)->sequence].header = pose_stamped.header;
        }

        if (SAVE_LOOP_PATH)
        {
            ofstream loop_path_file(VINS_RESULT_PATH, ios::app);
            loop_path_file.setf(ios::fixed, ios::floatfield);
            loop_path_file.precision(0);
            loop_path_file << (*it)->time_stamp * 1e9 << ",";
            loop_path_file.precision(5);
            loop_path_file  << P.x() << ","
                  << P.y() << ","
                  << P.z() << ","
                  << Q.w() << ","
                  << Q.x() << ","
                  << Q.y() << ","
                  << Q.z() << ","
                  << endl;
            loop_path_file.close();
        }
        //draw local connection
        if (SHOW_S_EDGE)
        {
            list<KeyFrame*>::reverse_iterator rit = keyframelist.rbegin();
            list<KeyFrame*>::reverse_iterator lrit;
            for (; rit != keyframelist.rend(); rit++)  
            {  
                if ((*rit)->index == (*it)->index)
                {
                    lrit = rit;
                    lrit++;
                    for (int i = 0; i < 4; i++)
                    {
                        if (lrit == keyframelist.rend())
                            break;
                        if((*lrit)->sequence == (*it)->sequence)
                        {
                            Vector3d conncected_P;
                            Matrix3d connected_R;
                            (*lrit)->getPose(conncected_P, connected_R);
                            posegraph_visualization->add_edge(P, conncected_P);
                        }
                        lrit++;
                    }
                    break;
                }
            } 
        }
        if (SHOW_L_EDGE)
        {
            if ((*it)->has_loop && (*it)->sequence == sequence_cnt)
            {
                
                KeyFrame* connected_KF = getKeyFrame((*it)->loop_index);
                Vector3d connected_P;
                Matrix3d connected_R;
                connected_KF->getPose(connected_P, connected_R);
                //(*it)->getVioPose(P, R);
                (*it)->getPose(P, R);
                if((*it)->sequence > 0)
                {
                    posegraph_visualization->add_loopedge(P, connected_P + Vector3d(VISUALIZATION_SHIFT_X, VISUALIZATION_SHIFT_Y, 0));
                }
            }
        }

    }
    publish();
    m_keyframelist.unlock();
}

/**
 * @brief 保存视觉地图
 * 
 */
void PoseGraph::savePoseGraph()
{
    m_keyframelist.lock();
    TicToc tmp_t;
    FILE *pFile;
    printf("pose graph path: %s\n",POSE_GRAPH_SAVE_PATH.c_str());
    printf("pose graph saving... \n");
    string file_path = POSE_GRAPH_SAVE_PATH + "pose_graph.txt";
    pFile = fopen (file_path.c_str(),"w");
    //fprintf(pFile, "index time_stamp Tx Ty Tz Qw Qx Qy Qz loop_index loop_info\n");
    list<KeyFrame*>::iterator it;
    // 遍历所有的KF
    for (it = keyframelist.begin(); it != keyframelist.end(); it++)
    {
        std::string image_path, descriptor_path, brief_path, keypoints_path;
        if (DEBUG_IMAGE)
        {
            image_path = POSE_GRAPH_SAVE_PATH + to_string((*it)->index) + "_image.png";
            imwrite(image_path.c_str(), (*it)->image);
        }
        // 分别存储VIO位姿和全局优化后位姿以及对应的回环信息
        Quaterniond VIO_tmp_Q{(*it)->vio_R_w_i};
        Quaterniond PG_tmp_Q{(*it)->R_w_i};
        Vector3d VIO_tmp_T = (*it)->vio_T_w_i;
        Vector3d PG_tmp_T = (*it)->T_w_i;

#if 0
        std::cout << "VIO_tmp_Q = " << VIO_tmp_Q.coeffs().transpose() << std::endl;
        std::cout << "VIO_tmp_T = " << VIO_tmp_T.transpose() << std::endl;
        std::cout << "PG_tmp_Q = " << PG_tmp_Q.coeffs().transpose() << std::endl;
        std::cout << "PG_tmp_T = " << PG_tmp_T.transpose() << std::endl;
#endif

        fprintf (pFile, " %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f %f %f %f %f %f %f %d\n",(*it)->index, (*it)->time_stamp, 
                                    VIO_tmp_T.x(), VIO_tmp_T.y(), VIO_tmp_T.z(), 
                                    PG_tmp_T.x(), PG_tmp_T.y(), PG_tmp_T.z(), 
                                    VIO_tmp_Q.w(), VIO_tmp_Q.x(), VIO_tmp_Q.y(), VIO_tmp_Q.z(), 
                                    PG_tmp_Q.w(), PG_tmp_Q.x(), PG_tmp_Q.y(), PG_tmp_Q.z(), 
                                    (*it)->loop_index, 
                                    (*it)->loop_info(0), (*it)->loop_info(1), (*it)->loop_info(2), (*it)->loop_info(3),
                                    (*it)->loop_info(4), (*it)->loop_info(5), (*it)->loop_info(6), (*it)->loop_info(7),
                                    (int)(*it)->keypoints.size());
        // 存储这一帧的特征点和描述子
        // write keypoints, brief_descriptors   vector<cv::KeyPoint> keypoints vector<BRIEF::bitset> brief_descriptors;
        assert((*it)->keypoints.size() == (*it)->brief_descriptors.size());
        brief_path = POSE_GRAPH_SAVE_PATH + to_string((*it)->index) + "_briefdes.dat";
        std::ofstream brief_file(brief_path, std::ios::binary);
        keypoints_path = POSE_GRAPH_SAVE_PATH + to_string((*it)->index) + "_keypoints.txt";
        FILE *keypoints_file;
        keypoints_file = fopen(keypoints_path.c_str(), "w");
        // 存储描述子，像素坐标和归一化相机坐标
        for (int i = 0; i < (int)(*it)->keypoints.size(); i++)
        {
            brief_file << (*it)->brief_descriptors[i] << endl;
            fprintf(keypoints_file, "%f %f %f %f\n", (*it)->keypoints[i].pt.x, (*it)->keypoints[i].pt.y, 
                                                     (*it)->keypoints_norm[i].pt.x, (*it)->keypoints_norm[i].pt.y);
        }
        brief_file.close();
        fclose(keypoints_file);
    }
    fclose(pFile);

    printf("save pose graph time: %f s\n", tmp_t.toc() / 1000);
    m_keyframelist.unlock();
}

/**
 * @brief 加载已有的视觉地图
 * 
 */
void PoseGraph::loadPoseGraph()
{
    TicToc tmp_t;
    FILE * pFile;
    string file_path = POSE_GRAPH_SAVE_PATH + "pose_graph.txt";
    printf("lode pose graph from: %s \n", file_path.c_str());
    printf("pose graph loading...\n");
    pFile = fopen (file_path.c_str(),"r");
    if (pFile == NULL)
    {
        printf("lode previous pose graph error: wrong previous pose graph path or no previous pose graph \n the system will start with new pose graph \n");
        return;
    }
    int index;
    double time_stamp;
    double VIO_Tx, VIO_Ty, VIO_Tz;
    double PG_Tx, PG_Ty, PG_Tz;
    double VIO_Qw, VIO_Qx, VIO_Qy, VIO_Qz;
    double PG_Qw, PG_Qx, PG_Qy, PG_Qz;
    double loop_info_0, loop_info_1, loop_info_2, loop_info_3;
    double loop_info_4, loop_info_5, loop_info_6, loop_info_7;
    int loop_index;
    int keypoints_num;
    Eigen::Matrix<double, 8, 1 > loop_info;
    int cnt = 0;
    // 按照存储的视觉地图的格式，加载视觉地图
    while (fscanf(pFile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %d", &index, &time_stamp, 
                                    &VIO_Tx, &VIO_Ty, &VIO_Tz, 
                                    &PG_Tx, &PG_Ty, &PG_Tz, 
                                    &VIO_Qw, &VIO_Qx, &VIO_Qy, &VIO_Qz, 
                                    &PG_Qw, &PG_Qx, &PG_Qy, &PG_Qz, 
                                    &loop_index,
                                    &loop_info_0, &loop_info_1, &loop_info_2, &loop_info_3, 
                                    &loop_info_4, &loop_info_5, &loop_info_6, &loop_info_7,
                                    &keypoints_num) != EOF) 
    {
        /*
        printf("I read: %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %d\n", index, time_stamp, 
                                    VIO_Tx, VIO_Ty, VIO_Tz, 
                                    PG_Tx, PG_Ty, PG_Tz, 
                                    VIO_Qw, VIO_Qx, VIO_Qy, VIO_Qz, 
                                    PG_Qw, PG_Qx, PG_Qy, PG_Qz, 
                                    loop_index,
                                    loop_info_0, loop_info_1, loop_info_2, loop_info_3, 
                                    loop_info_4, loop_info_5, loop_info_6, loop_info_7,
                                    keypoints_num);
        */
        cv::Mat image;
        std::string image_path, descriptor_path;
        if (DEBUG_IMAGE)
        {
            image_path = POSE_GRAPH_SAVE_PATH + to_string(index) + "_image.png";
            image = cv::imread(image_path.c_str(), 0);
        }
        // 同样加载VIO位姿和PG位姿
        Vector3d VIO_T(VIO_Tx, VIO_Ty, VIO_Tz);
        Vector3d PG_T(PG_Tx, PG_Ty, PG_Tz);
        Quaterniond VIO_Q;
        VIO_Q.w() = VIO_Qw;
        VIO_Q.x() = VIO_Qx;
        VIO_Q.y() = VIO_Qy;
        VIO_Q.z() = VIO_Qz;
        Quaterniond PG_Q;
        PG_Q.w() = PG_Qw;
        PG_Q.x() = PG_Qx;
        PG_Q.y() = PG_Qy;
        PG_Q.z() = PG_Qz;
        Matrix3d VIO_R, PG_R;
        VIO_R = VIO_Q.toRotationMatrix();
        PG_R = PG_Q.toRotationMatrix();
        // 记载回环信息
        Eigen::Matrix<double, 8, 1 > loop_info;
        loop_info << loop_info_0, loop_info_1, loop_info_2, loop_info_3, loop_info_4, loop_info_5, loop_info_6, loop_info_7;

        if (loop_index != -1)
            if (earliest_loop_index > loop_index || earliest_loop_index == -1)
            {
                earliest_loop_index = loop_index;
            }

        // load keypoints, brief_descriptors   
        string brief_path = POSE_GRAPH_SAVE_PATH + to_string(index) + "_briefdes.dat";
        std::ifstream brief_file(brief_path, std::ios::binary);
        string keypoints_path = POSE_GRAPH_SAVE_PATH + to_string(index) + "_keypoints.txt";
        FILE *keypoints_file;
        keypoints_file = fopen(keypoints_path.c_str(), "r");
        vector<cv::KeyPoint> keypoints;
        vector<cv::KeyPoint> keypoints_norm;
        // 记载特征点的像素坐标，归一化相机坐标和描述子
        vector<BRIEF::bitset> brief_descriptors;
        for (int i = 0; i < keypoints_num; i++)
        {
            BRIEF::bitset tmp_des;
            brief_file >> tmp_des;
            brief_descriptors.push_back(tmp_des);
            cv::KeyPoint tmp_keypoint;
            cv::KeyPoint tmp_keypoint_norm;
            double p_x, p_y, p_x_norm, p_y_norm;
            if(!fscanf(keypoints_file,"%lf %lf %lf %lf", &p_x, &p_y, &p_x_norm, &p_y_norm))
                printf(" fail to load pose graph \n");
            tmp_keypoint.pt.x = p_x;
            tmp_keypoint.pt.y = p_y;
            tmp_keypoint_norm.pt.x = p_x_norm;
            tmp_keypoint_norm.pt.y = p_y_norm;
            keypoints.push_back(tmp_keypoint);
            keypoints_norm.push_back(tmp_keypoint_norm);
        }
        brief_file.close();
        fclose(keypoints_file);
        // 根据加载的信息生成KF
        KeyFrame* keyframe = new KeyFrame(time_stamp, index, VIO_T, VIO_R, PG_T, PG_R, image, loop_index, loop_info, keypoints, keypoints_norm, brief_descriptors);
        loadKeyFrame(keyframe, 0);
        if (cnt % 20 == 0)
        {
            publish();
        }
        cnt++;
    }
    fclose (pFile);
    printf("load pose graph time: %f s\n", tmp_t.toc()/1000);
    base_sequence = 0;
}

void PoseGraph::publish()
{
    for (int i = 1; i <= sequence_cnt; i++)
    {
        //if (sequence_loop[i] == true || i == base_sequence)
        if (1 || i == base_sequence)
        {
            pub_pg_path.publish(path[i]);
            pub_path[i].publish(path[i]);
            posegraph_visualization->publish_by(pub_pose_graph, path[sequence_cnt].header);
        }
    }
    base_path.header.frame_id = "world";
    pub_base_path.publish(base_path);
    //posegraph_visualization->publish_by(pub_pose_graph, path[sequence_cnt].header);
}

// 注意这个函数并没有update任何关键帧的位姿
void PoseGraph::updateKeyFrameLoop(int index, Eigen::Matrix<double, 8, 1 > &_loop_info)
{
    KeyFrame* kf = getKeyFrame(index);  // 根据获取闭环的当前关键帧的索引获取关键帧
    // 虽然其在获得闭环的过程中使用pnp已经计算过了，但是这里的相对位姿会更准确，因此的将其更新到获取闭环的当前关键帧
    kf->updateLoop(_loop_info);
    // 如果yaw角度比比较大或者距离比较远，那么看到的可能就不是一个场景了，此时计算得到的闭环帧也就不那么可信了
    if (abs(_loop_info(7)) < 30.0 && Vector3d(_loop_info(0), _loop_info(1), _loop_info(2)).norm() < 20.0)
    {
        if (FAST_RELOCALIZATION)    // 肯定只有这种情况下触发
        {
            KeyFrame* old_kf = getKeyFrame(kf->loop_index); // 得到回环帧信息
            Vector3d w_P_old, w_P_cur, vio_P_cur;
            Matrix3d w_R_old, w_R_cur, vio_R_cur;
            // 注意这里闭环关键帧获取的是pose，其可能在前面矫正过，位姿是准确的
            old_kf->getPose(w_P_old, w_R_old);
            // 注意这里的当前关键帧获取的VIO位姿，这个位姿存在漂移
            kf->getVioPose(vio_P_cur, vio_R_cur);

            Vector3d relative_t;
            Quaterniond relative_q;
            // 得到T_loop_cur
            relative_t = kf->getLoopRelativeT();
            relative_q = (kf->getLoopRelativeQ()).toRotationMatrix();
            // T_w_loop * T_loop_cur = T_w_cur
            w_P_cur = w_R_old * relative_t + w_P_old;
            w_R_cur = w_R_old * relative_q;
            double shift_yaw;
            Matrix3d shift_r;
            Vector3d shift_t; 
            // 更新VIO位姿和修正位姿的delta pose
            /**
             * w_R_cur = Ry2 * Rp2 * Rr2, vio_R_cur = Ry1 * Rp1 * Rr1
             * shift_r = R(y2 - y1)
             * R(y2 - y1) * vio_R_cur = Ry2 * Rp1 * Rr1
             *
             * 如果令vio对应的世界系是w'，那么w_R_cur * vio_R_cur表示w'到w的旋转；Tw_cur * Tcur_w'
             * 使用shift_r描述两个世界系之间的相对旋转，由于VIO系统，只有yaw角漂移，因此仅仅使用yaw
             * 使用shift_t描述两个世界系之间的相对平移，注意相对平移的计算与相对旋转不同，相对平移的计算并没有使用shift_r，而是使用原始的w_R_cur * vio_R_cur
             * 这里的逻辑其实不太能自洽，主要是因为重力方向并没真是的重力方向，因此实际上并不仅仅只有yaw角度上有漂移，这是系统设计上的问题；系统并没有完全解决这个问题，而是采用了折中的方式
             */
            shift_yaw = Utility::R2ypr(w_R_cur).x() - Utility::R2ypr(vio_R_cur).x();
            shift_r = Utility::ypr2R(Vector3d(shift_yaw, 0, 0));
            shift_t = w_P_cur - w_R_cur * vio_R_cur.transpose() * vio_P_cur; 

            m_drift.lock();
            yaw_drift = shift_yaw;
            r_drift = shift_r;
            t_drift = shift_t;
            m_drift.unlock();
        }
    }
}