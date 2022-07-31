#include <iostream>
#include <opencv2/core/persistence.hpp>
#include <Eigen/Core>
#include <eigen3/Eigen/Dense>

#define precision 1e-5


const double fx = 190.978477;
const double fy = 190.973307;
const double cx = 254.931706;
const double cy = 256.897442;
const double k1 = 0.003482389402;
const double k2 = 0.000715034845;
const double k3 = -0.002053236141;
const double k4 = 0.000202936736;

void UnprojectByORB3(const cv::Point2f &p2D, double &theta);

void UnprojectByVINS(const cv::Point2f &p2D, double &theta);

void Check(const cv::Point2f &p2D, double &theta);

int main() {
    cv::Point2f p2D(100, 100);
    double theta1, theta2;

//    // test: ORB-SLAM3的逻辑与VINS-MONO的逻辑的结果是否一致
//    UnprojectByORB3(p2D, theta1);
//    std::cout << "使用ORB的代码求解入射角为：" << theta1 << std::endl;
//    Check(p2D, theta1);
//    std::cout << std::endl;
//
//    UnprojectByVINS(p2D, theta2);
//    std::cout << "使用VINS的代码求解入射角为：" << theta2 << std::endl;
//    Check(p2D, theta2);


    // test: 是否实特征值总是只有一个
    for (int i = 0; i < 600; i = i + 10) {
        for (int j = 0; j < 500; j = j + 10) {
            p2D = cv::Point2f(i, j);
            UnprojectByVINS(p2D, theta2);
        }
    }

    return 0;
}




void UnprojectByORB3(const cv::Point2f &p2D, double &theta) {
    cv::Point2f pw((p2D.x - cx) / fx, (p2D.y - cy) / fy);
    float theta_d = sqrtf(pw.x * pw.x + pw.y * pw.y);
    // theta_d = fminf(fmaxf(-CV_PI / 2.f, theta_d), CV_PI / 2.f);

    if (theta_d > 1e-8) {
        theta = theta_d;

        // 开始迭代
        for (int j = 0; j < 10; j++) {
            float theta2 = theta * theta,
                    theta4 = theta2 * theta2,
                    theta6 = theta4 * theta2,
                    theta8 = theta4 * theta4;
            float k0_theta2 = k1 * theta2,
                    k1_theta4 = k2 * theta4;
            float k2_theta6 = k3 * theta6,
                    k3_theta8 = k4 * theta8;
            float theta_fix = (theta * (1 + k0_theta2 + k1_theta4 + k2_theta6 + k3_theta8) - theta_d) /
                              (1 + 3 * k0_theta2 + 5 * k1_theta4 + 7 * k2_theta6 + 9 * k3_theta8);
            theta = theta - theta_fix;  // Gauss-Newton迭代公式
            if (fabsf(theta_fix) < precision)  // 如果更新量变得很小，表示接近最终值
                break;
        }
    }

}

void UnprojectByVINS(const cv::Point2f &p2D, double &theta) {
    double tol = 1e-10;
    cv::Point2f pw((p2D.x - cx) / fx, (p2D.y - cy) / fy);
    double p_u_norm = std::sqrt(pw.x * pw.x + pw.y * pw.y);

    int npow = 9;
    Eigen::MatrixXd coeffs(npow + 1, 1);
    coeffs.setZero();
    coeffs(0) = -p_u_norm;
    coeffs(1) = 1.0;

    coeffs(3) = k1;
    coeffs(5) = k2;
    coeffs(7) = k3;
    coeffs(9) = k4;

    Eigen::MatrixXd A(npow, npow);
    A.setZero();
    A.block(1, 0, npow - 1, npow - 1).setIdentity();
    A.col(npow - 1) = - coeffs.block(0, 0, npow, 1) / coeffs(npow);

    Eigen::EigenSolver<Eigen::MatrixXd> es(A);
    Eigen::MatrixXcd eigval = es.eigenvalues();

    std::vector<double> thetas;
    for (int i = 0; i < eigval.rows(); ++i) {
        if (fabs(eigval(i).imag()) > tol) {
            continue;
        }

        double t = eigval(i).real();

        if (t < -tol) {
            std::cout << "出现负数" << std::endl;
            continue;
        }
        else if (t < 0.0) {
            t = 0.0;
        }

        thetas.push_back(t);
    }

    if (thetas.empty()) {
        theta = p_u_norm;
    }
    else {
        theta = *std::min_element(thetas.begin(), thetas.end());
    }

    if (thetas.size() > 1) {
        std::cout << "所有的theta为：";
        for (int i = 0; i < thetas.size(); i++) {
            std::cout << thetas[i] << " ";
        }
        std::cout << std::endl;
    }

}

void Check(const cv::Point2f &p2D, double &theta) {
    cv::Point2f pw((p2D.x - cx) / fx, (p2D.y - cy) / fy);
    float theta_d = sqrtf(pw.x * pw.x + pw.y * pw.y);
    std::cout << "theta_d = " << theta_d << std::endl;

    float proj_theta_d = theta * (1 + k1 * std::pow(theta, 2) + k2 * std::pow(theta, 4) + k3 * std::pow(theta, 6) + k4 * std::pow(theta, 8));
    std::cout << "error = " << theta_d - proj_theta_d << std::endl;


}