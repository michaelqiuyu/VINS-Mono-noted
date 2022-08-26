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

    // xc's todo: 理论上dq的实部应该是1，但是在测试过程中发现，有时候dq的实部并不是1，而且其虚部也并不是近似零向量
    // 理论上讲，此时得到的迭代项是不准的，这可能与优化过程中的函数值相关，也就是不是在局部最优附近优化

#if 0
    std::cout << "执行到这里了" << std::endl;
    std::cout << "_p = " << _p.transpose() << std::endl;
    std::cout << "_q = " << _q.coeffs().transpose() << std::endl;
    std::cout << "dp = " << dp.transpose() << std::endl;
    std::cout << "dq = " << dq.coeffs().transpose() << std::endl;
    std::cout << "p = " << p.transpose() << std::endl;
    std::cout << "q = " << q.coeffs().transpose() << std::endl;
#endif

    p = _p + dp;
    q = (_q * dq).normalized();  // 右乘扰动，注意jacobian的求解需要与这里的扰动方向一致
#if 0
    std::cout << "p = " << p.transpose() << std::endl;
    std::cout << "q = " << q.coeffs().transpose() << std::endl;
#endif

    return true;
}

bool PoseLocalParameterization::ComputeJacobian(const double *x, double *jacobian) const
{
#if 0
    std::cout << "执行到这里了***" << std::endl;
#endif
    /**
     * 参考：https://blog.csdn.net/hzwwpgmwy/article/details/86490556
     *
     * 真实的jacobian是evaluate中计算的global_jacobian * 这里的jacobian
     * 也就是LocalParameterization的成员函数MultiplyByJacobian所做的事情：local_jacobian = global_jacobian * jacobian
     * jacobian存在的意义是对于过参数化的优化变量进行调整，得到没有过参数化的表达，比如在这个地方的位姿是7个维度(global_size)，实际上有效的维度是6个(local_size)
     * 位姿的前三个维度是平移，后4个维度是旋转，其中旋转的前三维是四元数的虚部，最后一维是四元数的实部
     * global_jacobian是15 * 7，jabocian是7 * 6，需要注意的是，这里并没有遵循一般的global_jacobian = e对q求导，jacobian = q对正切空间维度的变量(实际旋转的维度，在这里就是角度的一半，四元数的double cover效应)
     *
     * 此处的特殊性在于，直接使用e对正切空间维度的优化变量(角度)进行求解，但是为了使得global_jacobian的维度是15 * 7，因此global_jacobian的最后一列是0，且jacobian的最后一行是0，jacobian的前6行是一个单位阵
     */
    Eigen::Map<Eigen::Matrix<double, 7, 6, Eigen::RowMajor>> j(jacobian);
    j.topRows<6>().setIdentity();
    j.bottomRows<1>().setZero();

    return true;
}
