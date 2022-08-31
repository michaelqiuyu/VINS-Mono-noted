#include "initial_ex_rotation.h"

InitialEXRotation::InitialEXRotation(){
    frame_count = 0;
    Rc.push_back(Matrix3d::Identity());
    Rc_g.push_back(Matrix3d::Identity());
    Rimu.push_back(Matrix3d::Identity());
    ric = Matrix3d::Identity();
}

// 标定imu和相机之间的旋转外参，通过imu和图像计算的旋转使用手眼标定计算获得
bool InitialEXRotation::CalibrationExRotation(vector<pair<Vector3d, Vector3d>> corres, Quaterniond delta_q_imu, Matrix3d &calib_ric_result)
{
    /**
     * 对Ax = 0，且||x|| = 1的求解，通过最小二乘求得x，实际上-x同样满足要求
     *      1. 在此处，q与-q等价，对应同一个旋转
     *      2. 在DLT三角化处，由于需要除以第四个维度，因此x与-x获得的三维点等同
     */
    frame_count++;
    // 根据特征关联求解两个连续帧相机的旋转R12，2→1
    Rc.push_back(solveRelativeR(corres));
    // IMU预积分的结果是j → i的，也就是Rj = Ri * delta_R；Rj * pj = Ri * detla_R * pj，因此detla_R = Rij
    Rimu.push_back(delta_q_imu.toRotationMatrix());
    // 通过外参把imu的旋转转移到相机坐标系
    /**
     * author: xiongchao
     *
     * notes:
     *      1. q_ckbk * q_bkbk+1 * q_bk+1ck+1 = q_ckck+1
     *      2. 当frame_count < WINDOW_SIZE的时候，会一直执行这个函数，也就是一直求解外参；
     *      3. 之所以不断求解，而不是最后求解一次的原因是仅仅使用两帧求解的参数不一定具有普遍性（泛化能力）
     *      4. ric.inverse() * delta_q_imu * ric = Rc1i1 * Ri1i2 * Ri2c2 = Rc1c2，应该与Rc[i]接近
     *      5. 使用ric来评估本组的delta_q_imu与Rc之间是否匹配（意思是能否根据ric直接建立关系），如果匹配性不够，就降低权重
     */
    Rc_g.push_back(ric.inverse() * delta_q_imu * ric);  // ric是上一次求解得到的外参，初始化为单位阵，Rc_g用来计算核函数的

    /**
     * 从A的由来的角度解释标定的时候为什么需要三个轴都有激励：
     *      1. 假设我们始终都在激励一个轴并且运动速度十分均匀，那么对任意的i，i与i-1构建的方程都几乎相同，也就是无论有多少帧，实际上只有0和1构建的方程是有效方程
     *      2. 从数学上讲，就是A的前四行可以表示A的任何一行，Ax=0仅仅有前四个有效方程，后面的都是冗余的方程
     *      3. 我们在这里并不是求解Ax=0，而是求解最小二乘解，如果仅仅只有4个方程，那么求得的解不一定适合真实的场景（一旦运动不是匀速，也可以求解，两组解差异很大），没有体现最小二乘取“平均”（对每一个方程的误差都不会很大）的思想
     *      4. 只有在不同的时刻对不同的轴实行激励，那么A中的有效方程才会更多，求解的外参才会更加接近正确结果
     *      5. 如果body静止或者没有转动，那么就不可以根据陀螺仪的读数计算相对旋转了，此时也会产生不好的求解结果（注意：我们并不需要a，因此也不会使用到牛顿第二定律）
     *
     * 使用euroc_config的数据进行测试，如果我们每次仅仅使用一帧来求解，得到的结果是:
     *      1. 每次的结果都不一样，并且差距非常大
     *      2. 通过与正确的外参进行对比，发现优化的结果经常偏离正确值，偶尔才会有比较正确的结果
     */
    Eigen::MatrixXd A(frame_count * 4, 4);  // 每一帧提供4个方程
#if 0
    Eigen::MatrixXd A(4, 4);  // 每一帧提供4个方程
#endif
    A.setZero();
    int sum_ok = 0;
    /**
     * author: xiongchao
     *
     * notes:
     *      1. 从1开始因为至少两帧才有对极几何
     *      2. q_cb * q_bkbk+1 = q_ckck+1 * q_cb = q_ckbk+1
     *      3. 注意到，不管frame_count等于多少都会执行这个遍历，也就是不止使用当前帧及其前一帧来计算
     */
    for (int i = 1; i <= frame_count; i++)
    {
#if 0
        if (i != frame_count)
            continue;
#endif

        Quaterniond r1(Rc[i]);  // 通过视觉的对极几何得到
        Quaterniond r2(Rc_g[i]);  // 通过IMU和外参得到

        /**
         * author: xiongchao
         *
         * notes:
         *      1. q1 * q2.conj：对于单位四元数而言，逆等于共轭，也就是求二者之间的相对旋转
         *      2. 虚部 / 实部 = 单位向量n * tan(theta / 2)：反向求解theta
         *      3. 这里的angular_distance实际上求的是相对旋转的轴角
         *      4. 如果angular_distance偏大，说明这一组的delta_q_imu与Rc不太匹配，如果直接使用他们来计算外参，误差会比较大，所以权重应该降低
         *      5. 使用这种方式也保证了每次计算的Ric的波动不会大，避免每次计算得到的Ric不停地跳动，那样是不正常的
         *      6. 当然，这样也就要求第一次得到的ric是比较准确的，后面的计算才会在一个正确的方向上
         */
        double angular_distance = 180 / M_PI * r1.angularDistance(r2);

        /**
         * notes:
         *      1. 实际上，两帧图像之间的相对旋转非常小，也就是q接近(1, 0)；两帧IMU之间的相对旋转也非常小，也就是q接近(1, 0)；
         *         因此后面的L-R应该比较接近0矩阵（元素数值上），因此A应该比较接近0矩阵（元素数值上）
         *      2. 为了避免L和R非常接近而产生的数值稳定性的问题，可以使用第i帧和第j帧来构建方程，而不是使用连续帧来构建方程；后续可以尝试实现验证效果
         */
#if 0
        std::cout << "r1 = " << r1.w() << " " << r1.vec().transpose() << std::endl;
        std::cout << "r2 = " << r2.w() << " " << r2.vec().transpose() << std::endl;
        std::cout << "delra_angle = " << r1.angularDistance(r2) << std::endl;
#endif
        ROS_DEBUG(
            "%d %f", i, angular_distance);
        // 一个简单的核函数
        /**
         * author: xiongchao
         *
         * notes:
         *      1. kernel function: 防止某些误差占据了过大的权重，导致计算错误；常用于存在异常值的时候
         *      2. 当第一帧和第二帧进行外参标定的时候，由于此时ric为单位阵，因此误差较大
         *      3. 这是使用核函数的可能原因是：前面的IMU运动激励不够，因此标定的不太准确，因此降低他们的权重
         */
        double huber = angular_distance > 5.0 ? 5.0 / angular_distance : 1.0;
        ++sum_ok;
        Matrix4d L, R;

        double w = Quaterniond(Rc[i]).w();
        Vector3d q = Quaterniond(Rc[i]).vec();
        L.block<3, 3>(0, 0) = w * Matrix3d::Identity() + Utility::skewSymmetric(q);
        L.block<3, 1>(0, 3) = q;
        L.block<1, 3>(3, 0) = -q.transpose();
        L(3, 3) = w;

        Quaterniond R_ij(Rimu[i]);
        w = R_ij.w();
        q = R_ij.vec();
        R.block<3, 3>(0, 0) = w * Matrix3d::Identity() - Utility::skewSymmetric(q);
        R.block<3, 1>(0, 3) = q;
        R.block<1, 3>(3, 0) = -q.transpose();
        R(3, 3) = w;

        // A.block<4, 4>((i - 1) * 4, 0) = huber * (L - R);    // 作用在残差上面
        A.block<4, 4>(0, 0) = huber * (L - R);    // 作用在残差上面
#if 0
        std::cout << "L = " << L << std::endl;
        std::cout << "R = " << R << std::endl;
        std::cout << "L - R" << L - R << std::endl;
#endif
    }

#if 0
    std::cout << "A.t() * A = " << (A.transpose() * A).determinant() << std::endl;
#endif

    JacobiSVD<MatrixXd> svd(A, ComputeFullU | ComputeFullV);
    Matrix<double, 4, 1> x = svd.matrixV().col(3);
    Quaterniond estimated_R(x);  // IMU→camera
    ric = estimated_R.toRotationMatrix().inverse();  // camera→IMU
    //cout << svd.singularValues().transpose() << endl;
    //cout << ric << endl;
#if 0
    if (frame_count < WINDOW_SIZE) {
        std::cout << "frame_count = " << frame_count << std::endl;
        std::cout << "ric = " << ric << std::endl;
        std::cout << "奇异值为：" << svd.singularValues().transpose() << std::endl;
        std::cout << "svd.singularValues().tail<3>() = " << svd.singularValues().tail<3>() << std::endl;

        Eigen::Matrix4d mat = A.transpose() * A;
        std::cout << "mat = " << mat << std::endl;

    }

#endif
    Vector3d ric_cov;
    ric_cov = svd.singularValues().tail<3>();  // 最后三个
    // 倒数第二个奇异值，因为旋转是3个自由度，因此检查一下第三小的奇异值是否足够大，通常需要足够的运动激励才能保证得到没有奇异的解
    /**
     * author: xiongchao
     *
     * notes:对于Ax = 0且||x|| = 1
     *      1. 特征值 = 奇异值 * 奇异值
     *      2. 如果奇异值特别小，说明特征值更小，接近0，说明A的秩接近2
     *      3. 当A的秩为2的时候，x = k1v1 + k2v2，且k1 * v1_30 + k2 * v2_30 = 1，保证齐次性；此时有无穷个解，不符合要求
     *      4. 实际上在三角化那里也有同样的问题，但是那里并没有这个判断，是因为三角化的点会经过重投影误差的检验，而这里并没有类似的检验
     *
     * A.t * A的特征值之和等于其对角线元素之和；随着frame_count的增大，A.t * A的对角线上的元素必然递增，也就意味着特征值必然递增
     * 又其行列式的值等于特征值的乘积，也就意味着，行列式的值也在递增；通过使用euroc_config进行测试，结果与上面的结论相同
     */
    if (frame_count >= WINDOW_SIZE && ric_cov(1) > 0.25)
    {
        calib_ric_result = ric;
        return true;
    }
    else
        return false;
}

Matrix3d InitialEXRotation::solveRelativeR(const vector<pair<Vector3d, Vector3d>> &corres)
{
    if (corres.size() >= 9)
    {
        vector<cv::Point2f> ll, rr;
        for (int i = 0; i < int(corres.size()); i++)
        {
            ll.push_back(cv::Point2f(corres[i].first(0), corres[i].first(1)));
            rr.push_back(cv::Point2f(corres[i].second(0), corres[i].second(1)));
        }
        // 这里用的是相机坐标系，因此这个函数得到的也就是E矩阵
        cv::Mat E = cv::findFundamentalMat(ll, rr);
        cv::Mat_<double> R1, R2, t1, t2;
        decomposeE(E, R1, R2, t1, t2);

        // 旋转矩阵的行列式应该是1,这里如果是-1就取一下反
        /**
         * author: xiongchao
         *
         * notes: 此处的检查是否会有浮点数运算导致出现异常，如det(R1) = -1 + 1e-9，本应该取反，却没有 取反
         *      直接判断det(R1) < 0不是更加直接吗；计算需要使用绝对值的操作，使用1e-6会跟好
         */
        if (determinant(R1) + 1.0 < 1e-09)
        {
            E = -E;
            decomposeE(E, R1, R2, t1, t2);
        }
        double ratio1 = max(testTriangulation(ll, rr, R1, t1), testTriangulation(ll, rr, R1, t2));
        double ratio2 = max(testTriangulation(ll, rr, R2, t1), testTriangulation(ll, rr, R2, t2));
        // 获得内点比例最多的那一组解
        cv::Mat_<double> ans_R_cv = ratio1 > ratio2 ? R1 : R2;

        // 解出来的是R21，也就是1→2

        /**
         * author: xiongchao
         *
         * notes: 直接使用转置并将opencv转化为eigen即可
         */
        Matrix3d ans_R_eigen;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                ans_R_eigen(j, i) = ans_R_cv(i, j); // 这里转换成R12
        return ans_R_eigen;
    }
    return Matrix3d::Identity();
}

/**
 * @brief 通过三角化来检查R t是否合理
 * 
 * @param[in] l l相机的观测
 * @param[in] r r相机的观测
 * @param[in] R 旋转矩阵
 * @param[in] t 位移
 * @return double 
 */
double InitialEXRotation::testTriangulation(const vector<cv::Point2f> &l,
                                          const vector<cv::Point2f> &r,
                                          cv::Mat_<double> R, cv::Mat_<double> t)
{
    cv::Mat pointcloud;
    // 其中一帧设置为单位阵
    cv::Matx34f P = cv::Matx34f(1, 0, 0, 0,
                                0, 1, 0, 0,
                                0, 0, 1, 0);
    // 第二帧就设置为R t对应的位姿
    cv::Matx34f P1 = cv::Matx34f(R(0, 0), R(0, 1), R(0, 2), t(0),
                                 R(1, 0), R(1, 1), R(1, 2), t(1),
                                 R(2, 0), R(2, 1), R(2, 2), t(2));
    // 看一下opencv的定义：在triangulate.cpp中使用DLT算法
    cv::triangulatePoints(P, P1, l, r, pointcloud);
    int front_count = 0;
    for (int i = 0; i < pointcloud.cols; i++)
    {
        // 因为是齐次的，所以要求最后一维等于1
        double normal_factor = pointcloud.col(i).at<float>(3);
        // 得到在各自相机坐标系下的3d坐标，等价于Rp + t
        cv::Mat_<double> p_3d_l = cv::Mat(P) * (pointcloud.col(i) / normal_factor);
        cv::Mat_<double> p_3d_r = cv::Mat(P1) * (pointcloud.col(i) / normal_factor);
        // 通过深度是否大于0来判断是否合理
        if (p_3d_l(2) > 0 && p_3d_r(2) > 0)
            front_count++;
    }
    ROS_DEBUG("MotionEstimator: %f", 1.0 * front_count / pointcloud.cols);
    return 1.0 * front_count / pointcloud.cols;
}

// 具体解法参考多视图几何
void InitialEXRotation::decomposeE(cv::Mat E,
                                 cv::Mat_<double> &R1, cv::Mat_<double> &R2,
                                 cv::Mat_<double> &t1, cv::Mat_<double> &t2)
{
    cv::SVD svd(E, cv::SVD::MODIFY_A);
    cv::Matx33d W(0, -1, 0,
                  1, 0, 0,
                  0, 0, 1);
    cv::Matx33d Wt(0, 1, 0,
                   -1, 0, 0,
                   0, 0, 1);
    R1 = svd.u * cv::Mat(W) * svd.vt;
    R2 = svd.u * cv::Mat(Wt) * svd.vt;
    t1 = svd.u.col(2);
    t2 = -svd.u.col(2);
}
