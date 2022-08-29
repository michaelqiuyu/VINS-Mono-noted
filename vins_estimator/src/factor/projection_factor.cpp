#include "projection_factor.h"

Eigen::Matrix2d ProjectionFactor::sqrt_info;  // sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity()
double ProjectionFactor::sum_t;

ProjectionFactor::ProjectionFactor(const Eigen::Vector3d &_pts_i, const Eigen::Vector3d &_pts_j) : pts_i(_pts_i), pts_j(_pts_j)
{
#ifdef UNIT_SPHERE_ERROR
    Eigen::Vector3d b1, b2;
    Eigen::Vector3d a = pts_j.normalized();
    Eigen::Vector3d tmp(0, 0, 1);
    if(a == tmp)
        tmp << 1, 0, 0;
    b1 = (tmp - a * (a.transpose() * tmp)).normalized();
    b2 = a.cross(b1);
    tangent_base.block<1, 3>(0, 0) = b1.transpose();
    tangent_base.block<1, 3>(1, 0) = b2.transpose();
#endif
};

bool ProjectionFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    TicToc tic_toc;
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

    double inv_dep_i = parameters[3][0];
    // 地图点在i帧相机坐标系下坐标
    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    // 转成第i帧imu坐标系
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    // 转成世界坐标系
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    // 转到第j帧imu坐标系
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    // 转到第j帧相机坐标系
    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    Eigen::Map<Eigen::Vector2d> residual(residuals);

    /**
     * u = fx * x + cx, v = fy * y + cy
     *
     * (u - u0)^2 + (v - v0)^2 = fx^2 * (x - x0)^2 + fy^2 * (y - y0)^2 ！= (x - x0)^2 + (y - y0)^2
     * 如果fx和fy并不接近，那么优化归一化相机系的残差相当于认为给重投影误差加了权重，这样做并不符合重投影误差的定义，其与重投影误差并不等效
     */

#ifdef UNIT_SPHERE_ERROR 
    residual =  tangent_base * (pts_camera_j.normalized() - pts_j.normalized());
#else
    double dep_j = pts_camera_j.z();    // 第j帧相机系下深度
    residual = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();  // 归一化相机系下的重投影误差
#endif

    /**
     * 如果fx = fy，那么(x - x0)^2 + (y - y0)^2 = [(u - u0)^2 + (v - v0)^2] / fx^2
     * residual^2 = sqrt_info^2 * residual^2 = (FOCAL_LENGTH / 1.5)^2 * (x - x0)^2 + (y - y0)^2 = [(u - u0)^2 + (v - v0)^2] * (FOCAL_LENGTH / (1.5 * fx))^2
     * 令fx = k * FOCAL_LENGTH, e1 = sqrt((u - u0)^2 + (v - v0)^2), e2 = residual, 那么
     * e2 = e1 / (1.5 * k)，其中e1是重投影误差；因为这里不止视觉残差，这里的1.5有降低视觉权重的逻辑在里面，这里的参数应该是都是调参的结果
     */
    residual = sqrt_info * residual;    // 误差乘上信息矩阵，1.5个像素

    /**
     * 总体投影过程：相机系i帧下的地图点→body系i帧→世界系→body系j帧→相机系j帧→相机系j帧归一化
     * 总体投影方程：e = pcj / pcj(2) - obs
     *            pcj = Rbc.t * (Rwbj.t * (Rwbi * (Rbc * pci / λ + tbc) + twbi) - Rwbj.t * twbj) - Rbc.t * tbc
     *
     * 对旋转一律使用右扰动，对旋转矩阵求导和对四元数求导等价，在这里都是对轴角进行求导，因此在四元数的情况下，最后一列总是0，相应的ComputeJacobian函数是(I, 0).t
     *
     * e = (e1) = (x / z - u)，其中pcj = (x, y, z).t
     *     (e2) = (y / z - v)
     * e对pcj的导数为reduce，如下所示
     */

    if (jacobians)
    {
        Eigen::Matrix3d Ri = Qi.toRotationMatrix();
        Eigen::Matrix3d Rj = Qj.toRotationMatrix();
        Eigen::Matrix3d ric = qic.toRotationMatrix();
        Eigen::Matrix<double, 2, 3> reduce(2, 3);
#ifdef UNIT_SPHERE_ERROR
        double norm = pts_camera_j.norm();
        Eigen::Matrix3d norm_jaco;
        double x1, x2, x3;
        x1 = pts_camera_j(0);
        x2 = pts_camera_j(1);
        x3 = pts_camera_j(2);
        norm_jaco << 1.0 / norm - x1 * x1 / pow(norm, 3), - x1 * x2 / pow(norm, 3),            - x1 * x3 / pow(norm, 3),
                     - x1 * x2 / pow(norm, 3),            1.0 / norm - x2 * x2 / pow(norm, 3), - x2 * x3 / pow(norm, 3),
                     - x1 * x3 / pow(norm, 3),            - x2 * x3 / pow(norm, 3),            1.0 / norm - x3 * x3 / pow(norm, 3);
        reduce = tangent_base * norm_jaco;
#else
        // e对pcj的jacobian
        reduce << 1. / dep_j, 0, -pts_camera_j(0) / (dep_j * dep_j),    // 重投影误差对j帧相机坐标系下坐标求导
            0, 1. / dep_j, -pts_camera_j(1) / (dep_j * dep_j);
#endif
        reduce = sqrt_info * reduce;  // 链式求导，添加一个系数阵

        if (jacobians[0])  // 视觉重投影误差对i帧的旋转和平移的jacobian
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix<double, 3, 6> jaco_i;
            jaco_i.leftCols<3>() = ric.transpose() * Rj.transpose();  // pcj对i帧平移的jacobian
            jaco_i.rightCols<3>() = ric.transpose() * Rj.transpose() * Ri * -Utility::skewSymmetric(pts_imu_i);  // pcj对i帧旋转的jacobian

            jacobian_pose_i.leftCols<6>() = reduce * jaco_i;
            jacobian_pose_i.rightCols<1>().setZero();
        }

        if (jacobians[1])  // 视觉重投影误差对j帧的旋转和平移的jacobian
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);

            Eigen::Matrix<double, 3, 6> jaco_j;
            jaco_j.leftCols<3>() = ric.transpose() * -Rj.transpose();  // pcj对j帧平移的jacobian
            jaco_j.rightCols<3>() = ric.transpose() * Utility::skewSymmetric(pts_imu_j);  // pcj对j帧旋转的jacobian

            jacobian_pose_j.leftCols<6>() = reduce * jaco_j;
            jacobian_pose_j.rightCols<1>().setZero();
        }
        if (jacobians[2])  // 视觉重投影误差对外参的旋转和平移的jacobian
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_ex_pose(jacobians[2]);
            Eigen::Matrix<double, 3, 6> jaco_ex;
            jaco_ex.leftCols<3>() = ric.transpose() * (Rj.transpose() * Ri - Eigen::Matrix3d::Identity());  // pcj对外参平移的jacobian
            Eigen::Matrix3d tmp_r = ric.transpose() * Rj.transpose() * Ri * ric;
            jaco_ex.rightCols<3>() = -tmp_r * Utility::skewSymmetric(pts_camera_i) + Utility::skewSymmetric(tmp_r * pts_camera_i) +
                                     Utility::skewSymmetric(ric.transpose() * (Rj.transpose() * (Ri * tic + Pi - Pj) - tic));  // pcj对外参旋转的jacobian
            jacobian_ex_pose.leftCols<6>() = reduce * jaco_ex;
            jacobian_ex_pose.rightCols<1>().setZero();
        }
        if (jacobians[3])  // 视觉重投影误差对尺度的jacobian
        {
            Eigen::Map<Eigen::Vector2d> jacobian_feature(jacobians[3]);
#if 1
            jacobian_feature = reduce * ric.transpose() * Rj.transpose() * Ri * ric * pts_i * -1.0 / (inv_dep_i * inv_dep_i);  // pcj对尺度的jacobian
#else
            jacobian_feature = reduce * ric.transpose() * Rj.transpose() * Ri * ric * pts_i;
#endif
        }
    }
    sum_t += tic_toc.toc();

    return true;
}

// 并未使用
void ProjectionFactor::check(double **parameters)
{
    double *res = new double[15];
    double **jaco = new double *[4];
    jaco[0] = new double[2 * 7];
    jaco[1] = new double[2 * 7];
    jaco[2] = new double[2 * 7];
    jaco[3] = new double[2 * 1];
    Evaluate(parameters, res, jaco);
    puts("check begins");

    puts("my");

    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[2]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Vector2d>(jaco[3]) << std::endl
              << std::endl;

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    double inv_dep_i = parameters[3][0];

    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);


    Eigen::Vector2d residual;
#ifdef UNIT_SPHERE_ERROR 
    residual =  tangent_base * (pts_camera_j.normalized() - pts_j.normalized());
#else
    double dep_j = pts_camera_j.z();
    residual = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();
#endif
    residual = sqrt_info * residual;

    puts("num");
    std::cout << residual.transpose() << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 2, 19> num_jacobian;
    for (int k = 0; k < 19; k++)
    {
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

        Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
        Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
        double inv_dep_i = parameters[3][0];

        int a = k / 3, b = k % 3;
        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

        if (a == 0)
            Pi += delta;
        else if (a == 1)
            Qi = Qi * Utility::deltaQ(delta);
        else if (a == 2)
            Pj += delta;
        else if (a == 3)
            Qj = Qj * Utility::deltaQ(delta);
        else if (a == 4)
            tic += delta;
        else if (a == 5)
            qic = qic * Utility::deltaQ(delta);
        else if (a == 6)
            inv_dep_i += delta.x();

        Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
        Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
        Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
        Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
        Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

        Eigen::Vector2d tmp_residual;
#ifdef UNIT_SPHERE_ERROR 
        tmp_residual =  tangent_base * (pts_camera_j.normalized() - pts_j.normalized());
#else
        double dep_j = pts_camera_j.z();
        tmp_residual = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();
#endif
        tmp_residual = sqrt_info * tmp_residual;
        num_jacobian.col(k) = (tmp_residual - residual) / eps;
    }
    std::cout << num_jacobian << std::endl;
}
