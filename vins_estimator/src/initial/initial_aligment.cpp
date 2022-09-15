#include "initial_alignment.h"

/**
 * @brief 求解陀螺仪零偏，同时利用求出来的零偏重新进行预积分
 * 
 * @param[in] all_image_frame 
 * @param[in] Bgs 
 */
void solveGyroscopeBias(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs)
{
    /**
     * 思路：
     *      1. 通过相机位姿和相机与IMU的外参计算两帧IMU之间的相对旋转（这个不含有陀螺仪零偏）；预积分也能得到两帧之间的相对旋转；二者理论上应该相等
     *      2. 由于陀螺仪零偏初始为0，因此不会严格相等，因此可以构建最小二乘进行陀螺仪零偏的求解
     *      3. 二者不等的原因当然不仅仅是陀螺仪零偏，还有视觉估计的误差，预积分的误差（中值法本身自带误差）等；因此这里实际上是调整陀螺仪零偏使其适应视觉的结果
     *
     * 计算过程：
     *      1. qbk_c0 = qbk_ck * qck_c0，这里的旋转都是已知的，因此可以计算qbk_c0
     *      2. qbk_bk+1 = qbk_c0 * qc0_bk+1，这里的旋转都是已知的
     *      3. qbk_bk+1.inv * rbk_bk+1的虚部应该是零向量，据此即可构建最小二乘问题，当然rbk_bk+1 = rbk_bk+1_观测 * (1, 0.5 * Jbw * delta_bw)
     *      4. 执行变换qbk_bk+1.inv * rbk_bk+1_观测 * (1, 0.5 * Jbw * delta_bw) = I(q) → (1, 0.5 * Jbw * delta_bw) = rbk_bk+1_观测.inv * qbk_bk+1
     *      5. 均取虚部即可构建一个Ax = b的模型
     *      6. 由于k的取值不唯一，因此模型转变为∑||Ai * x - bi||^2 → x.t() * ∑Ai.t() * Ai * x - ∑(2 * bi.t() * Ai) * x + C
     *      7. 令A =  ∑Ai.t() * Ai，B = ∑(bi.t() * Ai)，则问题转化为min(x.t() * A * x - 2 * B * x)
     *      8. 上述问题为一个凸优化问题，有闭式解，即求解c
     */
    Matrix3d A;
    Vector3d b;
    Vector3d delta_bg;
    A.setZero();
    b.setZero();
    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++)
    {
        frame_j = next(frame_i);
        MatrixXd tmp_A(3, 3);
        tmp_A.setZero();
        VectorXd tmp_b(3);
        tmp_b.setZero();
        // frame_i->second.R是q_cobk = q_cock * q_ckbk
        Eigen::Quaterniond q_ij(frame_i->second.R.transpose() * frame_j->second.R);  // qbk_bk+1
        tmp_A = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_R, O_BG);  // Jbw
        tmp_b = 2 * (frame_j->second.pre_integration->delta_q.inverse() * q_ij).vec();  // b
        A += tmp_A.transpose() * tmp_A;  // Ai.t() * Ai
        b += tmp_A.transpose() * tmp_b;  // bi.t() * Ai
    }
    delta_bg = A.ldlt().solve(b);  // Ax = B
    ROS_WARN_STREAM("gyroscope bias initial calibration " << delta_bg.transpose());
    // 滑窗中的零偏设置为求解出来的零偏
    for (int i = 0; i <= WINDOW_SIZE; i++)
        Bgs[i] += delta_bg;  // 窗口内的零偏都相同
    // 对all_image_frame中预积分量根据当前零偏重新积分
    // all_image_frame的begin没有重新预积分，是因为第一帧不做预积分
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
    {
        frame_j = next(frame_i);
        // 只在初始化的时候重新传播
        frame_j->second.pre_integration->repropagate(Vector3d::Zero(), Bgs[0]);
    }
}


MatrixXd TangentBasis(Vector3d &g0)
{
    Vector3d b, c;
    Vector3d a = g0.normalized();
    Vector3d tmp(0, 0, 1);
    if(a == tmp)
        tmp << 1, 0, 0;
    b = (tmp - a * (a.transpose() * tmp)).normalized();
    c = a.cross(b);
    MatrixXd bc(3, 2);
    bc.block<3, 1>(0, 0) = b;
    bc.block<3, 1>(0, 1) = c;
    return bc;
}

/**
 * @brief 得到了一个初始的重力向量，引入重力大小作为先验，再进行几次迭代优化，求解最终的变量
 * 
 * @param[in] all_image_frame 
 * @param[in] g 
 * @param[in] x 
 */
void RefineGravity(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    /**
     * 在这里再次求解了x，那么还有必要执行LinearAlignment吗？有必要
     *      1. 在本函数中，需要提供一个初始的重力方向g，这个g是通过LinearAlignment来得到初值的，如果直接使用一个随意的初值，结果可能不一定准确
     *      2. g的求解不仅仅与g相关，实际上与所有的待求解的参数都有关系，因为g影响了矩阵A和b的构建，因此对最终所有的参数都有影响
     *      3. 在这里迭代求解g的过程中，实际上速度和尺度也在发生变化，而我们并没有使用到LinearAlignment中计算得到的速度和尺度
     *      4. LinearAlignment仅仅只是为了给本函数提供一个重力方向的初值
     */
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
    /**
     * 求解的变量有3 * all_frame_count + 3，能够构建的方程的个数有6 * (all_frame_count - 1)
     * 如果想要正常求解，那么6 * (all_frame_count - 1) >= 3 * all_frame_count + 3，也就是all_frame_count >= 2
     */
    int n_state = all_frame_count * 3 + 2 + 1;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);  // 注意g0在不停地更新
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);

            MatrixXd tmp_A(6, 9);
            tmp_A.setZero();
            VectorXd tmp_b(6);
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;


            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 2>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;
            tmp_A.block<3, 1>(0, 8) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;     
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0;

            tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 2>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0;


            Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
            cov_inv.setIdentity();

            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

            A.bottomRightCorner<3, 3>() += r_A.bottomRightCorner<3, 3>();
            b.tail<3>() += r_b.tail<3>();

            A.block<6, 3>(i * 3, n_state - 3) += r_A.topRightCorner<6, 3>();
            A.block<3, 6>(n_state - 3, i * 3) += r_A.bottomLeftCorner<3, 6>();
        }
            A = A * 1000.0;
            b = b * 1000.0;
            x = A.ldlt().solve(b);
            VectorXd dg = x.segment<2>(n_state - 3);
            // 始终保证重力的模长为G的模长：g0的方向在不停的更新；并不能保证求得的结果直接就在球面上，而是需要不停地归一化
            g0 = (g0 + lxly * dg).normalized() * G.norm();
            //double s = x(n_state - 1);
    }   
    g = g0;
}

/**
 * @brief 求解各帧的速度，枢纽帧的重力方向，以及尺度
 * 
 * @param[in] all_image_frame 
 * @param[in] g 
 * @param[in] x 
 * @return true 
 * @return false 
 */
bool LinearAlignment(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    // 这一部分内容对照论文进行理解
    int all_frame_count = all_image_frame.size();
    /**
     * 求解的变量有3 * all_frame_count + 4，能够构建的方程的个数有6 * (all_frame_count - 1)
     * 如果想要正常求解，那么6 * (all_frame_count - 1) >= 3 * all_frame_count + 4，也就是all_frame_count > 2
     */
    int n_state = all_frame_count * 3 + 3 + 1;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int i = 0;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
    {
        // 此处的代码比较复杂，可读性很差，其实直接构建A矩阵即可，没有必要先通过tmp_A.t() * tmp_A然后赋值来构建
        // 直接按照求解的整个向量的形式来构建这个很大的A，其中A的维度是(6 * window_size) * n_state
        // 每相邻两帧可以提供6个方程，其中平移和速度的预积分各三个，具体见论文公式10~12
        frame_j = next(frame_i);

        MatrixXd tmp_A(6, 10);
        tmp_A.setZero();
        VectorXd tmp_b(6);
        tmp_b.setZero();

        double dt = frame_j->second.pre_integration->sum_dt;

        tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
        tmp_A.block<3, 3>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
        // 对尺度部分的元素除了100，也就是求解从s变成了100s
        tmp_A.block<3, 1>(0, 9) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;     
        tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0];
        //cout << "delta_p   " << frame_j->second.pre_integration->delta_p.transpose() << endl;
        tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
        tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
        tmp_A.block<3, 3>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
        tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v;
        //cout << "delta_v   " << frame_j->second.pre_integration->delta_v.transpose() << endl;

        Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
        //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
        //MatrixXd cov_inv = cov.inverse();
        cov_inv.setIdentity();

        MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
        VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

        A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
        b.segment<6>(i * 3) += r_b.head<6>();

        A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
        b.tail<4>() += r_b.tail<4>();

        A.block<6, 4>(i * 3, n_state - 4) += r_A.topRightCorner<6, 4>();
        A.block<4, 6>(n_state - 4, i * 3) += r_A.bottomLeftCorner<4, 6>();
    }
    // 增强数值稳定性
    A = A * 1000.0;
    b = b * 1000.0;
    x = A.ldlt().solve(b);
    // 由于求解的是100s，因此这里需要除以100
    double s = x(n_state - 1) / 100.0;
    ROS_DEBUG("estimated scale: %f", s);
    g = x.segment<3>(n_state - 4);
    ROS_DEBUG_STREAM(" result g     " << g.norm() << " " << g.transpose());
    // 做一些检查
    if(fabs(g.norm() - G.norm()) > 1.0 || s < 0)
    {
        return false;
    }
    // 重力修复
    RefineGravity(all_image_frame, g, x);
    // 得到真实尺度
    s = (x.tail<1>())(0) / 100.0;
    (x.tail<1>())(0) = s;
    ROS_DEBUG_STREAM(" refine     " << g.norm() << " " << g.transpose());
    if(s < 0.0 )
        return false;   
    else
        return true;
}

/**
 * @brief 
 * 
 * @param[in] all_image_frame 每帧的位姿和对应的预积分量
 * @param[out] Bgs 陀螺仪零偏
 * @param[out] g 重力向量
 * @param[out] x 其他状态量
 * @return true 
 * @return false 
 */

bool VisualIMUAlignment(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs, Vector3d &g, VectorXd &x)
{
    solveGyroscopeBias(all_image_frame, Bgs);

    if(LinearAlignment(all_image_frame, g, x))
        return true;
    else 
        return false;
}
