#include "feature_tracker.h"

int FeatureTracker::n_id = 0;

bool inBorder(const cv::Point2f &pt)
{
    const int BORDER_SIZE = 1;
    int img_x = cvRound(pt.x);
    int img_y = cvRound(pt.y);
    return BORDER_SIZE <= img_x && img_x < COL - BORDER_SIZE && BORDER_SIZE <= img_y && img_y < ROW - BORDER_SIZE;
}

// 根据状态位，进行“瘦身”
void reduceVector(vector<cv::Point2f> &v, vector<uchar> status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i])
            v[j++] = v[i];
    v.resize(j);
}

void reduceVector(vector<int> &v, vector<uchar> status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i])
            v[j++] = v[i];
    v.resize(j);
}


FeatureTracker::FeatureTracker()
{
}

// 给现有的特征点设置mask，目的为了特征点的均匀化：第一帧不做这个操作
void FeatureTracker::setMask()
{
    if(FISHEYE)
        mask = fisheye_mask.clone();
    else
        mask = cv::Mat(ROW, COL, CV_8UC1, cv::Scalar(255));

    // prefer to keep features that are tracked for long time
    vector<pair<int, pair<cv::Point2f, int>>> cnt_pts_id;

    // 执行光流跟踪之后，track_cnt在reduceVector中被resize了，在还没有赋值的时候，使用的是默认值0
    for (unsigned int i = 0; i < forw_pts.size(); i++)  // 第一帧的时候，forw_pts为空
        cnt_pts_id.push_back(make_pair(track_cnt[i], make_pair(forw_pts[i], ids[i])));
    // 利用光流特点，追踪多的稳定性好，排前面
    sort(cnt_pts_id.begin(), cnt_pts_id.end(), [](const pair<int, pair<cv::Point2f, int>> &a, const pair<int, pair<cv::Point2f, int>> &b)
         {
            return a.first > b.first;
         });

    forw_pts.clear();
    ids.clear();
    track_cnt.clear();

    /*
     * 这里的均匀化比较简单，先对特征点的质量进行排序，然后画圈进行均匀化
     */
    for (auto &it : cnt_pts_id)
    {
        if (mask.at<uchar>(it.second.first) == 255)
        {
            // 把挑选剩下的特征点重新放进容器
            forw_pts.push_back(it.second.first);
            ids.push_back(it.second.second);
            track_cnt.push_back(it.first);
            // opencv函数，把周围一个圆内全部置0,这个区域不允许别的特征点存在，避免特征点过于集中
            /*
             * author: xiongchao
             * 一旦画圆之后，这个院内的mask的值就变了，就不会通过这个if判断了
             */
            cv::circle(mask, it.second.first, MIN_DIST, 0, -1);
        }
    }
}

// 把新的点加入容器，id给-1作为区分
void FeatureTracker::addPoints()
{
    for (auto &p : n_pts)
    {
        forw_pts.push_back(p);
        ids.push_back(-1);
        track_cnt.push_back(1);
    }
}

/**
 * @brief 
 * 
 * @param[in] _img 输入图像
 * @param[in] _cur_time 图像的时间戳
 * 1、图像均衡化预处理
 * 2、光流追踪
 * 3、提取新的特征点（如果发布）
 * 4、所有特征点去畸变，计算速度
 */
void FeatureTracker::readImage(const cv::Mat &_img, double _cur_time)
{
    cv::Mat img;
    TicToc t_r;
    cur_time = _cur_time;  // 令人比较迷惑的是，cur_img和cur_pts表示上一帧的图像和特征点，而cur_time表示当前时间

    if (EQUALIZE)
    {
        // 图像太暗或者太亮，提特征点比较难，所以均衡化一下
        // ! opencv 函数看一下
        cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE(3.0, cv::Size(8, 8));
        TicToc t_c;
        clahe->apply(_img, img);
        ROS_DEBUG("CLAHE costs: %fms", t_c.toc());
    }
    else
        img = _img;

    // 这里forw表示当前，cur表示上一帧
    if (forw_img.empty())   // 第一次输入图像，prev_img这个没用
    {
        prev_img = cur_img = forw_img = img;
    }
    else
    {
        forw_img = img;
    }

    forw_pts.clear();

    if (cur_pts.size() > 0) // 上一帧有特征点，就可以进行光流追踪了
    {
        TicToc t_o;
        vector<uchar> status;
        vector<float> err;
        // 调用opencv函数进行光流追踪
        // Step 1 通过opencv光流追踪给的状态位剔除outlier
        cv::calcOpticalFlowPyrLK(cur_img, forw_img, cur_pts, forw_pts, status, err, cv::Size(21, 21), 3);

        for (int i = 0; i < int(forw_pts.size()); i++)
            // Step 2 通过图像边界剔除outlier
            /*
             * author: xiongchao
             * 这里的范围是1 <= x <= row - 1; 1 <= y <= col - 1
             */
            if (status[i] && !inBorder(forw_pts[i]))    // 追踪状态好检查在不在图像范围
                status[i] = 0;
        reduceVector(prev_pts, status); // 没用到
        // notes：在这里会将cur_pts和forw_pts的size处理成一样，后面会进行对极几何约束
        reduceVector(cur_pts, status);
        reduceVector(forw_pts, status);
        reduceVector(ids, status);  // 特征点的id
        reduceVector(cur_un_pts, status);   // 去畸变后的坐标
        reduceVector(track_cnt, status);    // 追踪次数
        ROS_DEBUG("temporal optical flow costs: %fms", t_o.toc());
    }
    // 被追踪到的是上一帧就存在的，因此追踪数+1
    for (auto &n : track_cnt)
        n++;

    /**
     * author: xiongchao
     * 如果不发布的话，我们就仅仅只是做光流跟踪和去畸变操作
     *
     */
    if (PUB_THIS_FRAME)
    {
        // Step 3 通过对级约束来剔除outlier
        rejectWithF();  // 第二次运行这个函数的时候，track_cnt依旧为空
        ROS_DEBUG("set mask begins");
        TicToc t_m;
        setMask();  // 均匀化；第二次运行这个函数的时候，track_cnt才开始计数
        ROS_DEBUG("set mask costs %fms", t_m.toc());

        ROS_DEBUG("detect feature begins");
        TicToc t_t;
        /*
         * author: xiongchao
         * 新提的特征点的数目
         */
        int n_max_cnt = MAX_CNT - static_cast<int>(forw_pts.size());
        if (n_max_cnt > 0)
        {
            if(mask.empty())
                cout << "mask is empty " << endl;
            if (mask.type() != CV_8UC1)
                cout << "mask type wrong " << endl;
            if (mask.size() != forw_img.size())
                cout << "wrong size " << endl;
            // 只有发布才可以提取更多特征点，同时避免提的点进mask
            // 会不会这些点集中？会，不过没关系，他们下一次作为老将就得接受均匀化的洗礼
            /**
             * author: xiongchao
             * 第一帧的特征提取就是这里，第一次提取的特征没有经过均匀化，但对后续跟踪到的特征点会进行均匀化
             *
             * MIN_DIST用于提取的特征点之间的最小间距（像素单位）
             *
             * 想法：这里的逻辑是：先提特征点→光流跟踪→均匀化，需要注意的是提特征点的接口已经有均匀化的操作了，不然的话逻辑就会有问题，因为如果提的特征点太密集，可能在均匀化的时候删掉了大量的特征点
             */
            cv::goodFeaturesToTrack(forw_img, n_pts, MAX_CNT - forw_pts.size(), 0.01, MIN_DIST, mask);
        }
        else
            n_pts.clear();
        ROS_DEBUG("detect feature costs: %fms", t_t.toc());

        ROS_DEBUG("add feature begins");
        TicToc t_a;
        addPoints();
        ROS_DEBUG("selectFeature costs: %fms", t_a.toc());
    }
    prev_img = cur_img;
    prev_pts = cur_pts;
    prev_un_pts = cur_un_pts;   // 以上三个量无用
    cur_img = forw_img; // 实际上是上一帧的图像
    cur_pts = forw_pts; // 上一帧的特征点
    // 在做对极约束的时候，已经对特征点做了去畸变，在这里实际上只需要对新增的特征点去畸变即可；第一帧的特征点全部来自新增
    undistortedPoints();
    prev_time = cur_time;
}

/**
 * @brief 
 * 
 */
void FeatureTracker::rejectWithF()
{
    // 当前被追踪到的光流至少8个点
    if (forw_pts.size() >= 8)  // 第一帧不做这个操作
    {
        ROS_DEBUG("FM ransac begins");
        TicToc t_f;
        // 在做完光流跟踪后，会将cur_pts和forw_pts的size处理成一样
        vector<cv::Point2f> un_cur_pts(cur_pts.size()), un_forw_pts(forw_pts.size());
        for (unsigned int i = 0; i < cur_pts.size(); i++)
        {
            Eigen::Vector3d tmp_p;
            // 得到相机坐标系的值：2D点反向投影到去畸变的3D坐标，针孔相机返回的是归一化相机系坐标，鱼眼返回的是一个单位球面坐标
            // 这里的操作里面含有x/z，实际上就是在计算归一化相机系坐标了
            m_camera->liftProjective(Eigen::Vector2d(cur_pts[i].x, cur_pts[i].y), tmp_p);
            // 这里用一个虚拟相机，原因同样参考https://github.com/HKUST-Aerial-Robotics/VINS-Mono/issues/48
            // 这里有个好处就是对F_THRESHOLD和相机无关
            // 投影到虚拟相机的像素坐标系
            /**
             * 使用虚拟相机的原因是：对极几何的阈值受到fx和fy的影响(例如fx = alpha * f，alpha表示像素/米，也就是说fx代表了焦距等价于多少个像素)
             * 但是在对极几何中，需要给定点到极线的距离的阈值约束，对于不同的fx，应该给与不同的阈值（这一点与ORB-SLAM3是不同的，ORB3认为在像素平面上是1自由度的卡方分布阈值，与fx和fy的值是无关的）
             * 在这里有一个假设是：FOCAL_LENGTH对应于F_THRESHOLD；二者根据真实大小同等比例缩放，例如fx=2*FOCAL_LENGTH，那么对极几何的阈值为2*F_THRESHOLD
             *
             * 从以下可以推导：
             *      x / z = (u - cx) / fx, y / z = (v - cy) / fy
             *      u' = f0 * x / z + cx / 2; v' = f0 * y / z + cy / 2
             *
             *      u'   f0/fx    0    (1-f0/fx)*cx         u
             *      v' =   0    f0/fy  (1-f0/fy)*cy    *    v     →      p' = A * p
             *      1      0      0          1              1
             *
             *      对于对极几何约束而言，p2.t * k2.inv.t * t.hat * R * k1.inv * p1 = 0 →
             *          p2'.t * A.int.t * K2.inv * t.hat * R * k1.inv * A.inv * p1' = 0 →
             *
             *                   f0/fx    0    (1-f0/fx)*cx         fx    0    cx       f0   0   cx
             *      其中A * K1 =   0     f0/fy  (1-f0/fy)*cy    *    0     fy   cy   =   0    f0  cy
             *                    0       0         1               0     0    1        0    0   1
             *      通过虚拟相机，使得对极几何使用的内参矩阵实际上变成了一个设定好的、固定的矩阵，这样就可以对任何fx和fy不用每次都变动阈值，
             *      而是使用一个设定好的F_THRESHOLD（我们当然也可以每次动态的计算阈值）
             *
             */
            tmp_p.x() = FOCAL_LENGTH * tmp_p.x() / tmp_p.z() + COL / 2.0;
            tmp_p.y() = FOCAL_LENGTH * tmp_p.y() / tmp_p.z() + ROW / 2.0;
            un_cur_pts[i] = cv::Point2f(tmp_p.x(), tmp_p.y());

            m_camera->liftProjective(Eigen::Vector2d(forw_pts[i].x, forw_pts[i].y), tmp_p);
            tmp_p.x() = FOCAL_LENGTH * tmp_p.x() / tmp_p.z() + COL / 2.0;
            tmp_p.y() = FOCAL_LENGTH * tmp_p.y() / tmp_p.z() + ROW / 2.0;
            un_forw_pts[i] = cv::Point2f(tmp_p.x(), tmp_p.y());
        }

        vector<uchar> status;
        // opencv接口计算本质矩阵，某种意义也是一种对级约束的outlier剔除
        /**
         * author: xiongchao
         *
         * F_THRESHOLD：点到极限的距离，这里的阈值的设定是否应该参考ORB-SLAM3中依赖卡方分布检验的结果
         */
        cv::findFundamentalMat(un_cur_pts, un_forw_pts, cv::FM_RANSAC, F_THRESHOLD, 0.99, status);
        int size_a = cur_pts.size();
        reduceVector(prev_pts, status);
        reduceVector(cur_pts, status);
        reduceVector(forw_pts, status);
        reduceVector(cur_un_pts, status);
        reduceVector(ids, status);
        reduceVector(track_cnt, status);
        ROS_DEBUG("FM ransac: %d -> %lu: %f", size_a, forw_pts.size(), 1.0 * forw_pts.size() / size_a);
        ROS_DEBUG("FM ransac costs: %fms", t_f.toc());
    }
}

/**
 * @brief 
 * 
 * @param[in] i 
 * @return true 
 * @return false 
 *  给新的特征点赋上id,越界就返回false
 */
bool FeatureTracker::updateID(unsigned int i)
{
    if (i < ids.size())
    {
        if (ids[i] == -1)
            ids[i] = n_id++;  // n_id静态变量，从零开始递增
        return true;
    }
    else
        return false;
}

void FeatureTracker::readIntrinsicParameter(const string &calib_file)
{
    ROS_INFO("reading paramerter of camera %s", calib_file.c_str());
    // 读到的相机内参赋给m_camera
    m_camera = CameraFactory::instance()->generateCameraFromYamlFile(calib_file);
}

void FeatureTracker::showUndistortion(const string &name)
{
    cv::Mat undistortedImg(ROW + 600, COL + 600, CV_8UC1, cv::Scalar(0));
    vector<Eigen::Vector2d> distortedp, undistortedp;
    for (int i = 0; i < COL; i++)
        for (int j = 0; j < ROW; j++)
        {
            Eigen::Vector2d a(i, j);
            Eigen::Vector3d b;
            m_camera->liftProjective(a, b);
            distortedp.push_back(a);
            undistortedp.push_back(Eigen::Vector2d(b.x() / b.z(), b.y() / b.z()));
            //printf("%f,%f->%f,%f,%f\n)\n", a.x(), a.y(), b.x(), b.y(), b.z());
        }
    for (int i = 0; i < int(undistortedp.size()); i++)
    {
        cv::Mat pp(3, 1, CV_32FC1);
        pp.at<float>(0, 0) = undistortedp[i].x() * FOCAL_LENGTH + COL / 2;
        pp.at<float>(1, 0) = undistortedp[i].y() * FOCAL_LENGTH + ROW / 2;
        pp.at<float>(2, 0) = 1.0;
        //cout << trackerData[0].K << endl;
        //printf("%lf %lf\n", p.at<float>(1, 0), p.at<float>(0, 0));
        //printf("%lf %lf\n", pp.at<float>(1, 0), pp.at<float>(0, 0));
        if (pp.at<float>(1, 0) + 300 >= 0 && pp.at<float>(1, 0) + 300 < ROW + 600 && pp.at<float>(0, 0) + 300 >= 0 && pp.at<float>(0, 0) + 300 < COL + 600)
        {
            undistortedImg.at<uchar>(pp.at<float>(1, 0) + 300, pp.at<float>(0, 0) + 300) = cur_img.at<uchar>(distortedp[i].y(), distortedp[i].x());
        }
        else
        {
            //ROS_ERROR("(%f %f) -> (%f %f)", distortedp[i].y, distortedp[i].x, pp.at<float>(1, 0), pp.at<float>(0, 0));
        }
    }
    cv::imshow(name, undistortedImg);
    cv::waitKey(0);
}

// 当前帧所有点统一去畸变，同时计算特征点速度，用来后续时间戳标定
/**
 * author: xiongchao
 * 需要注意的是，在桶形畸变的情形下，如果更新后的点与原始点不在同一个1/4图像内，这个去畸变算法就会发散，但是这种情形不会出现
 */
void FeatureTracker::undistortedPoints()
{
    cur_un_pts.clear();
    cur_un_pts_map.clear();
    //cv::undistortPoints(cur_pts, un_pts, K, cv::Mat());
    // notes: 在这里对所有的特征点都做了去畸变的操作，但是在对极几何中已经对特征点做了去畸变了，因此在这里可以选择对没有做对极几何的点去畸变
    for (unsigned int i = 0; i < cur_pts.size(); i++)
    {
        // 有的之前去过畸变了，这里连同新的特征点重新做一次，有些费时；ids[i]!=-1表示非新增，在对极几何中已经做了，不需要重新做了，此处逻辑可以进行修改，加快运行速度
        Eigen::Vector2d a(cur_pts[i].x, cur_pts[i].y);
        Eigen::Vector3d b;
        m_camera->liftProjective(a, b);
        // 获取的是归一化相机系坐标
        cur_un_pts.push_back(cv::Point2f(b.x() / b.z(), b.y() / b.z()));
        // id->坐标的map
        cur_un_pts_map.insert(make_pair(ids[i], cv::Point2f(b.x() / b.z(), b.y() / b.z())));
        //printf("cur pts id %d %f %f", ids[i], cur_un_pts[i].x, cur_un_pts[i].y);
    }
    // caculate points velocity
    if (!prev_un_pts_map.empty())
    {
        double dt = cur_time - prev_time;
        pts_velocity.clear();
        for (unsigned int i = 0; i < cur_un_pts.size(); i++)
        {
            if (ids[i] != -1)
            {
                std::map<int, cv::Point2f>::iterator it;
                it = prev_un_pts_map.find(ids[i]);
                // 找到同一个特征点
                /**
                 * author: xiongchao
                 * ids保存的是从第一帧开始能够跟踪到的特征点的id以及后面新加特征点的id，可能有的id从第一帧开始一直被跟踪到，
                 * 导致id一直没变，但是prev_un_pts_map在一直更新，因此此处的逻辑是没有问题的
                 *
                 * 思考：ids[i] != -1是否意味着这个特征点被跟踪到了，因此it != prev_un_pts_map.end()一定成立
                 */
                if (it != prev_un_pts_map.end())
                {
                    double v_x = (cur_un_pts[i].x - it->second.x) / dt;
                    double v_y = (cur_un_pts[i].y - it->second.y) / dt;
                    // 得到在归一化平面的速度
                    pts_velocity.push_back(cv::Point2f(v_x, v_y));
                }
                else
                    pts_velocity.push_back(cv::Point2f(0, 0));
            }
            else

            {
                pts_velocity.push_back(cv::Point2f(0, 0));
            }
        }
    }
    else
    {
        // 第一帧的情况
        for (unsigned int i = 0; i < cur_pts.size(); i++)
        {
            pts_velocity.push_back(cv::Point2f(0, 0));
        }
    }
    prev_un_pts_map = cur_un_pts_map;
}
