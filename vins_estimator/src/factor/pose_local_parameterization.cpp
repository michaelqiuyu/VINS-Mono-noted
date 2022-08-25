#include "pose_local_parameterization.h"

bool PoseLocalParameterization::Plus(const double *x, const double *delta, double *x_plus_delta) const
{
    // Map类用于经过C++中普通的连续指针或者数组 （raw C/C++ arrays）来构造Eigen里的Matrix类，这就比如Eigen里的Matrix类的数据和raw C++array 共享了一片地址，也就是引用
    Eigen::Map<const Eigen::Vector3d> _p(x);
    Eigen::Map<const Eigen::Quaterniond> _q(x + 3);

    Eigen::Map<const Eigen::Vector3d> dp(delta);

    Eigen::Quaterniond dq = Utility::deltaQ(Eigen::Map<const Eigen::Vector3d>(delta + 3));

    Eigen::Map<Eigen::Vector3d> p(x_plus_delta);
    Eigen::Map<Eigen::Quaterniond> q(x_plus_delta + 3);

#if 1
    std::cout << "执行到这里了" << std::endl;
    std::cout << "_p = " << _p.transpose() << std::endl;
    std::cout << "_q = " << _q.coeffs().transpose() << std::endl;
    std::cout << "dp = " << dp.transpose() << std::endl;
    std::cout << "dq = " << dq.coeffs().transpose() << std::endl;
    std::cout << "p = " << p.transpose() << std::endl;
    std::cout << "q = " << q.coeffs().transpose() << std::endl;
#endif

    p = _p + dp;
    q = (_q * dq).normalized();
#if 1
    std::cout << "p = " << p.transpose() << std::endl;
    std::cout << "q = " << q.coeffs().transpose() << std::endl;
#endif

    return true;
}
bool PoseLocalParameterization::ComputeJacobian(const double *x, double *jacobian) const
{
#if 1
    std::cout << "执行到这里了***" << std::endl;
#endif
    Eigen::Map<Eigen::Matrix<double, 7, 6, Eigen::RowMajor>> j(jacobian);
    j.topRows<6>().setIdentity();
    j.bottomRows<1>().setZero();

    return true;
}
