#pragma once
#include <ros/assert.h>
#include <iostream>
#include <eigen3/Eigen/Dense>

#include "../utility/utility.h"
#include "../parameters.h"
#include "integration_base.h"

#include <ceres/ceres.h>

class IMUFactor : public ceres::SizedCostFunction<15, 7, 9, 7, 9>
{
  public:
    IMUFactor() = delete;  // 使用“=delete”来将不想要的函数定义成删除的函数，我们将不能以任何方式使用它们
    IMUFactor(IntegrationBase* _pre_integration):pre_integration(_pre_integration)
    {
    }
    /**
     * @brief  使用ceres解析求导，必须重载这个函数
     * 
     * @param[in] parameters 这是一个二维数组，每个参数块都是一个double数组，而一个观测会对多个参数块形成约束
     * @param[in] residuals 残差的计算结果，是一个一维数组，残差就是该观测量和约束的状态量通过某种关系形成残差
     * @param[in] jacobians 残差对参数块的雅克比矩阵，这也是一个二维数组，对任意一个参数块的雅克比矩阵都是一个一维数组
     * @return true 
     * @return false 
     */
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {
        // 便于后续计算，把参数块都转换成eigen
        // imu预积分的约束的参数是相邻两帧的位姿 速度和零偏
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Vector3d Bai(parameters[1][3], parameters[1][4], parameters[1][5]);
        Eigen::Vector3d Bgi(parameters[1][6], parameters[1][7], parameters[1][8]);

        Eigen::Vector3d Pj(parameters[2][0], parameters[2][1], parameters[2][2]);
        Eigen::Quaterniond Qj(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

        Eigen::Vector3d Vj(parameters[3][0], parameters[3][1], parameters[3][2]);
        Eigen::Vector3d Baj(parameters[3][3], parameters[3][4], parameters[3][5]);
        Eigen::Vector3d Bgj(parameters[3][6], parameters[3][7], parameters[3][8]);

//Eigen::Matrix<double, 15, 15> Fd;
//Eigen::Matrix<double, 15, 12> Gd;

//Eigen::Vector3d pPj = Pi + Vi * sum_t - 0.5 * g * sum_t * sum_t + corrected_delta_p;
//Eigen::Quaterniond pQj = Qi * delta_q;
//Eigen::Vector3d pVj = Vi - g * sum_t + corrected_delta_v;
//Eigen::Vector3d pBaj = Bai;
//Eigen::Vector3d pBgj = Bgi;

//Vi + Qi * delta_v - g * sum_dt = Vj;
//Qi * delta_q = Qj;

//delta_p = Qi.inverse() * (0.5 * g * sum_dt * sum_dt + Pj - Pi);
//delta_v = Qi.inverse() * (g * sum_dt + Vj - Vi);
//delta_q = Qi.inverse() * Qj;

#if 0
        if ((Bai - pre_integration->linearized_ba).norm() > 0.10 ||
            (Bgi - pre_integration->linearized_bg).norm() > 0.01)
        {
            pre_integration->repropagate(Bai, Bgi);
        }
#endif
        /**
         * notes:
         *      1. 注意这里传入的位姿表达是PoseLocalParameterization，其前3维是平移，后4维是旋转，并且旋转使用的是四元数；四元数的前3维是虚部，最后一维是实部；
         *         在其成员函数plus中规定了其更新方式是右乘，因此此处的求导也都是根据右乘推导的；由于四元数是过参数的表达形式，因此这里可以看到，作者构建的是一个15 * 7
         *         的矩阵，但是最后一列是0向量，在这里并没有按照标准的ceres方式对q（global_size）进行jacobian求解，而是直接对theta（local_size）进行求解，这也是为什么
         *         PoseLocalParameterization中的ComputeJacobian是(I, 0).t的形式了，也就是e对(theta, 0)求导，(theta, 0)对theta求导
         *         由于jacobian是(I, 0).t，所以global_jacobian的最后一列实际上没有必要求解，直接赋值为0即可
         *
         */

        Eigen::Map<Eigen::Matrix<double, 15, 1>> residual(residuals);
        // notes: 注意残差的构建方式是Visual - IMU，因此对零偏的导数有负号
        residual = pre_integration->evaluate(Pi, Qi, Vi, Bai, Bgi,
                                            Pj, Qj, Vj, Baj, Bgj);
        // 因为ceres没有g2o设置信息矩阵的接口，因此置信度直接乘在残差上，这里通过LLT分解，相当于将信息矩阵开根号
        // 执行LU分解中的乔累斯基分解（LLT分解），并获取LT矩阵（matrixL为获取L矩阵）
        Eigen::Matrix<double, 15, 15> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 15, 15>>(pre_integration->covariance.inverse()).matrixL().transpose();
        //sqrt_info.setIdentity();
        // 这就是带有信息矩阵的残差
        /**
         * H = LLT
         * (LT * e)T * (LT * e) = eT * L * LT * e = eT * H * e
         */
        residual = sqrt_info * residual;
        // 关于雅克比的计算手动推导一下
        if (jacobians)
        {
            double sum_dt = pre_integration->sum_dt;
            Eigen::Matrix3d dp_dba = pre_integration->jacobian.template block<3, 3>(O_P, O_BA);
            Eigen::Matrix3d dp_dbg = pre_integration->jacobian.template block<3, 3>(O_P, O_BG);

            Eigen::Matrix3d dq_dbg = pre_integration->jacobian.template block<3, 3>(O_R, O_BG);

            Eigen::Matrix3d dv_dba = pre_integration->jacobian.template block<3, 3>(O_V, O_BA);
            Eigen::Matrix3d dv_dbg = pre_integration->jacobian.template block<3, 3>(O_V, O_BG);

            // 对雅克比矩阵中最大值和最小值的检验：数值过大或者过小都不好
            if (pre_integration->jacobian.maxCoeff() > 1e8 || pre_integration->jacobian.minCoeff() < -1e8)
            {
                ROS_WARN("numerical unstable in preintegration");
                //std::cout << pre_integration->jacobian << std::endl;
                // ROS_BREAK();
            }

            /**
             * notes: 以下公式的推导的详情请见：崔华坤文档的公式20~24
             */

            // e对四个参数求导
            if (jacobians[0])  // 对i时刻的平移和旋转的jacobian：崔华坤文档的公式20和21
            {
                Eigen::Map<Eigen::Matrix<double, 15, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
                jacobian_pose_i.setZero();

                /**
                 * notes:
                 *      1. 对四元数扰动的形式（右乘）
                 *      2. 不用过分专注于求导的形式，因为只需要求导和更新一致即可
                 *      3. 求导过程中可能使用了近似，需要注意近似的精度，如果精度太差，求导也就失去了意义；至少应该在分子趋近于0的时候，分母也逼近0
                 *      4. 对四元数求导与对旋转矩阵求导的结果是等价的，只要他们的扰动方向一致即可
                 */

                // 对位移的预积分对i时刻的平移和旋转的jacobian
                jacobian_pose_i.block<3, 3>(O_P, O_P) = -Qi.inverse().toRotationMatrix();  // 对i时刻的平移
                jacobian_pose_i.block<3, 3>(O_P, O_R) = Utility::skewSymmetric(Qi.inverse() * (0.5 * G * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt));  // 对i时刻的旋转

#if 0
            jacobian_pose_i.block<3, 3>(O_R, O_R) = -(Qj.inverse() * Qi).toRotationMatrix();
#else
                // 对旋转的预积分对位姿的jacobian：对i时刻平移的jacobian为0
                Eigen::Quaterniond corrected_delta_q = pre_integration->delta_q * Utility::deltaQ(dq_dbg * (Bgi - pre_integration->linearized_bg));
                jacobian_pose_i.block<3, 3>(O_R, O_R) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(corrected_delta_q)).bottomRightCorner<3, 3>();  // 对i时刻的旋转
#endif
                // 对速度的预积分对位姿的jacobian：对i时刻的平移为0
                jacobian_pose_i.block<3, 3>(O_V, O_R) = Utility::skewSymmetric(Qi.inverse() * (G * sum_dt + Vj - Vi));  // 对i时刻的旋转

                // 陀螺仪零偏的残差的加速度零偏的残差：对i时刻的旋转和平移的jacobian都是0

                // 链式求导，sqrt_info * residual对residual求导，residual对优化变量求导
                jacobian_pose_i = sqrt_info * jacobian_pose_i;

                if (jacobian_pose_i.maxCoeff() > 1e8 || jacobian_pose_i.minCoeff() < -1e8)
                {
                    ROS_WARN("numerical unstable in preintegration");
                    //std::cout << sqrt_info << std::endl;
                    //ROS_BREAK();
                }
            }
            if (jacobians[1])  // 对i时刻的速度、加速度零偏和陀螺仪零偏的jacobian：崔华坤文档的公式20和22
            {
                Eigen::Map<Eigen::Matrix<double, 15, 9, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[1]);
                jacobian_speedbias_i.setZero();
                // 对位移的预积分对i时刻的速度、加速度零偏和陀螺仪零偏的jacobian
                jacobian_speedbias_i.block<3, 3>(O_P, O_V - O_V) = -Qi.inverse().toRotationMatrix() * sum_dt;  // 对i时刻的速度
                jacobian_speedbias_i.block<3, 3>(O_P, O_BA - O_V) = -dp_dba;  // 对i时刻的加速度的零偏
                jacobian_speedbias_i.block<3, 3>(O_P, O_BG - O_V) = -dp_dbg;  // 对i时刻的陀螺仪的零偏

#if 0
            jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -dq_dbg;
#else
                //Eigen::Quaterniond corrected_delta_q = pre_integration->delta_q * Utility::deltaQ(dq_dbg * (Bgi - pre_integration->linearized_bg));
                //jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -Utility::Qleft(Qj.inverse() * Qi * corrected_delta_q).bottomRightCorner<3, 3>() * dq_dbg;
                // 对旋转的预积分对陀螺仪零偏的jacobian
                /**
                 * 此处的求解虽然与四元数类似，但是还是有区别的，在对四元数扰动求解中，我们直接操作的就是四元数的扰动项delta_theta
                 * 在此处，对旋转的预积分是零偏的函数，但是其扰动并不是(1, 0.5 * delta_bg)，而依然是(1, 0.5 * delta_theta)，且delta_theta = JRg * delta_bg
                 */
                 // 对旋转的预积分对i时刻的速度、加速度零偏和陀螺仪零偏的jacobian：对i时刻的速度和加速度零偏的jacobian为0
                jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -Utility::Qleft(Qj.inverse() * Qi * pre_integration->delta_q).bottomRightCorner<3, 3>() * dq_dbg;  // 对i时刻的陀螺仪的零偏
#endif

                // 对速度的预积分对i时刻的速度、加速度零偏和陀螺仪零偏的jacobian：
                jacobian_speedbias_i.block<3, 3>(O_V, O_V - O_V) = -Qi.inverse().toRotationMatrix();  // 对i时刻的速度
                jacobian_speedbias_i.block<3, 3>(O_V, O_BA - O_V) = -dv_dba;  // 对i时刻的加速度零偏
                jacobian_speedbias_i.block<3, 3>(O_V, O_BG - O_V) = -dv_dbg;  // 对i时刻的陀螺仪零偏

                // 加速度零偏的残差：对陀螺仪零偏的jacobian为0
                jacobian_speedbias_i.block<3, 3>(O_BA, O_BA - O_V) = -Eigen::Matrix3d::Identity();  // 对i时刻的加速度零偏
                // 陀螺仪零偏的残差：对加速度零偏的jacobian为0
                jacobian_speedbias_i.block<3, 3>(O_BG, O_BG - O_V) = -Eigen::Matrix3d::Identity();  // 对i时刻的陀螺仪零偏

                jacobian_speedbias_i = sqrt_info * jacobian_speedbias_i;

                //ROS_ASSERT(fabs(jacobian_speedbias_i.maxCoeff()) < 1e8);
                //ROS_ASSERT(fabs(jacobian_speedbias_i.minCoeff()) < 1e8);
            }
            if (jacobians[2])  // 对j时刻的平移和旋转的jacobian：崔华坤文档的公式20和23
            {
                Eigen::Map<Eigen::Matrix<double, 15, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[2]);
                jacobian_pose_j.setZero();
                // 对位移的预积分对j时刻的平移和旋转的jacobian：对j时刻的旋转的jacobian为0
                jacobian_pose_j.block<3, 3>(O_P, O_P) = Qi.inverse().toRotationMatrix();  // 对j时刻的平移

#if 0
            jacobian_pose_j.block<3, 3>(O_R, O_R) = Eigen::Matrix3d::Identity();
#else
                // 对旋转的预积分对j时刻的平移和旋转的jacobian：对j时刻的平移的jacobian为0
                Eigen::Quaterniond corrected_delta_q = pre_integration->delta_q * Utility::deltaQ(dq_dbg * (Bgi - pre_integration->linearized_bg));
                jacobian_pose_j.block<3, 3>(O_R, O_R) = Utility::Qleft(corrected_delta_q.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();  // 对j时刻的旋转

                // 对速度的预积分对j时刻的平移和旋转的jacobian都是0
                // 加速度零偏的残差对j时刻的平移和旋转的jacobian都是0
                // 陀螺仪零偏的残差对j时刻的平移和旋转的jacobian都是0
#endif

                jacobian_pose_j = sqrt_info * jacobian_pose_j;

                //ROS_ASSERT(fabs(jacobian_pose_j.maxCoeff()) < 1e8);
                //ROS_ASSERT(fabs(jacobian_pose_j.minCoeff()) < 1e8);
            }
            if (jacobians[3])  // 对j时刻的速度、加速度零偏和陀螺仪零偏的jacobian：崔华坤文档的公式20和24
            {
                Eigen::Map<Eigen::Matrix<double, 15, 9, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[3]);
                jacobian_speedbias_j.setZero();

                // 对位移的预积分对j时刻的速度、加速度零偏和陀螺仪零偏的jacobian都是0
                // 对旋转的预积分对j时刻的速度、加速度零偏和陀螺仪零偏的jacobian都是0

                // 对速度的预积分对j时刻的速度、加速度零偏和陀螺仪零偏的jacobian：对j时刻的加速度零偏和陀螺仪零偏的jacobian都是0
                jacobian_speedbias_j.block<3, 3>(O_V, O_V - O_V) = Qi.inverse().toRotationMatrix();

                // 加速度的残差对j时刻的速度、加速度零偏和陀螺仪零偏的jacobian：对j时刻速度和陀螺仪零偏的jacobian为0
                jacobian_speedbias_j.block<3, 3>(O_BA, O_BA - O_V) = Eigen::Matrix3d::Identity();  // 对j时刻的加速度零偏

                // 陀螺仪的残差对j时刻的速度、加速度零偏和陀螺仪零偏的jacobian：对j时刻速度和加速度零偏的jacobian为0
                jacobian_speedbias_j.block<3, 3>(O_BG, O_BG - O_V) = Eigen::Matrix3d::Identity();  // 对j时刻的陀螺仪零偏

                jacobian_speedbias_j = sqrt_info * jacobian_speedbias_j;

                //ROS_ASSERT(fabs(jacobian_speedbias_j.maxCoeff()) < 1e8);
                //ROS_ASSERT(fabs(jacobian_speedbias_j.minCoeff()) < 1e8);
            }
        }

        return true;
    }

    //bool Evaluate_Direct(double const *const *parameters, Eigen::Matrix<double, 15, 1> &residuals, Eigen::Matrix<double, 15, 30> &jacobians);

    //void checkCorrection();
    //void checkTransition();
    //void checkJacobian(double **parameters);
    IntegrationBase* pre_integration;

};

