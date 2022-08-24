#include "utility.h"

Eigen::Matrix3d Utility::g2R(const Eigen::Vector3d &g)
{
    Eigen::Matrix3d R0;
    Eigen::Vector3d ng1 = g.normalized();
    // 注意：重力并不总是(0, 0, 1).t，还有可能是(0, 0, -1)，取决于预积分推导中使用的形式
    // 在这里建议使用配置参数避免不同的地方使用了不同的重力定义
    Eigen::Vector3d ng2{0, 0, 1.0};
    // ng2 = R0 * ng1，也就是说R0表示从ng1对应的坐标系变换到世界系
    R0 = Eigen::Quaterniond::FromTwoVectors(ng1, ng2).toRotationMatrix();
    double yaw = Utility::R2ypr(R0).x();
    // yaw不可观，因此将yaw角补偿掉
    // R0 = Ry * Rp * Rr → R(-y) * R0 = Rp * Rr
    R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;
    // R0 = Utility::ypr2R(Eigen::Vector3d{-90, 0, 0}) * R0;
    return R0;
}
