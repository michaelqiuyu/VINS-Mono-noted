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
    frame_count++;
    // 根据特征关联求解两个连续帧相机的旋转R12，2→1
    Rc.push_back(solveRelativeR(corres));
    Rimu.push_back(delta_q_imu.toRotationMatrix());  // 用来计算核函数的
    // 通过外参把imu的旋转转移到相机坐标系
    /**
     * author: xiongchao
     *
     * notes:
     *      1. q_ckbk * q_bkbk+1 * q_bk+1ck+1 = q_ckck+1
     *      2. 当frame_count < WINDOW_SIZE的时候，会一直执行这个函数，也就是一直求解外参；
     *      3. 之所以不断求解，而不是最后求解一次的原因是仅仅使用两帧求解的参数不一定具有普遍性（泛化能力）
     *      4. ric.inverse() * delta_q_imu * ric = Rc1i1 * Ri1i2 * Ri2c2 = Rc1c2，应该与Rc[i]接近
     */
    Rc_g.push_back(ric.inverse() * delta_q_imu * ric);  // ric是上一次求解得到的外参

    Eigen::MatrixXd A(frame_count * 4, 4);
    A.setZero();
    int sum_ok = 0;
    /**
     * author: xiongchao
     *
     * notes:
     *      1. 从1开始因为至少两帧才有对极几何
     *      2. q_cb * q_bkbk+1 = q_ckck+1 * q_cb = q_ckbk+1
     */
    for (int i = 1; i <= frame_count; i++)
    {
        Quaterniond r1(Rc[i]);  // 通过视觉的对极几何得到
        Quaterniond r2(Rc_g[i]);  // 通过IMU和外参得到

        /**
         * author: xiongchao
         *
         * notes:
         *      1. q1 * q2.conj：对于单位四元数而言，逆等于共轭，也就是求二者之间的相对旋转
         *      2. 虚部 / 实部 = 单位向量n * tan(theta / 2)：反向求解theta
         *      3. 这里的angular_distance实际上求的是相对旋转的轴角
         */
        double angular_distance = 180 / M_PI * r1.angularDistance(r2);
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

        A.block<4, 4>((i - 1) * 4, 0) = huber * (L - R);    // 作用在残差上面
    }

    JacobiSVD<MatrixXd> svd(A, ComputeFullU | ComputeFullV);
    Matrix<double, 4, 1> x = svd.matrixV().col(3);
    Quaterniond estimated_R(x);  // IMU→camera
    ric = estimated_R.toRotationMatrix().inverse();  // camera→IMU
    //cout << svd.singularValues().transpose() << endl;
    //cout << ric << endl;
    Vector3d ric_cov;
    ric_cov = svd.singularValues().tail<3>();
    // 倒数第二个奇异值，因为旋转是3个自由度，因此检查一下第三小的奇异值是否足够大，通常需要足够的运动激励才能保证得到没有奇异的解
    /**
     * author: xiongchao
     *
     * notes:对于Ax = 0且||x|| = 1
     *      1. 特征值 = 奇异值 * 奇异值
     *      2. 如果奇异值特别小，说明特征值更小，接近0，说明A的秩接近2
     *      3. 当A的秩为2的时候，x = k1v1 + k2v2，且k1 * v1_30 + k2 * v2_30 = 1，保证齐次性；此时有无穷个解，不符合要求
     *      4. 实际上在三角化那里也有同样的问题，但是那里并没有这个判断，是因为三角化的点会经过重投影误差的检验，而这里并没有类似的检验
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
         * notes: 此处的检查是否会有浮点数运算导致出现异常，如det(R1) = -0.99999998，本应该取反，却没有 取反
         *      直接判断det(R1) < 0不是更加直接吗
         */
        if (determinant(R1) + 1.0 < 1e-09)
        {
            E = -E;
            decomposeE(E, R1, R2, t1, t2);
        }
        double ratio1 = max(testTriangulation(ll, rr, R1, t1), testTriangulation(ll, rr, R1, t2));
        double ratio2 = max(testTriangulation(ll, rr, R2, t1), testTriangulation(ll, rr, R2, t2));
        cv::Mat_<double> ans_R_cv = ratio1 > ratio2 ? R1 : R2;

        // 解出来的是R21，也就是1→2

        /**
         * author: xiongchao
         *
         * notes: 直接使用转置不行吗？
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
    // 看一下opencv的定义
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

// 具体解法参考多视角几何
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
