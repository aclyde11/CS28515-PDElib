//
// Created by Austin Clyde on 5/28/18.
//

#include "Problem2d.h"


Problem2d::Problem2d() {
    mesh = VariMesh(1, 1, 1, 1);
}

Eigen::VectorXd multiplyStiff(Eigen::VectorXd v, VariMesh mesh) {
    double r = mesh.dx[0] / mesh.dy[0];
    double a = 0.5 * (r + 1.0 / r);
    double b = 0.5 * r;

    Eigen::VectorXd result(v.size());

    for (int i = 0; i < 2; i++) {
        result(i * 4) = a * v(i * 4) + b * v(i * 4 + 1) + b * v(i * 4 + 3);
        result(i * 4 + 1) = -1.0 * b * v(i * 4) + a * v(i * 4 + 1) - b * v(i * 4 + 2);
        result(i * 4 + 2) = -1.0 * b * v(i * 4 + 1) + a * v(i * 4 + 2) - b * v(i * 4 + 3);
        result(i * 4 + 3) = -1.0 * b * v(i * 4) - b * v(i * 4 + 2) + a * v(i * 4 + 3);
    }
    return result;
}

Eigen::VectorXd conjugateGradient(Eigen::MatrixXd C, Eigen::VectorXd g, Eigen::VectorXd y0) {
    Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    double eps = 1e-3;

    int N_max = 2;
    Eigen::VectorXd s = C * y0 - g;
    Eigen::VectorXd r(s), e(s.size());
    double alpha, beta;

    while (r.norm() > eps) {
        alpha = -1.0 * (r.dot(s) / (C * s).dot(s));
        std::cout << "alpha = " << alpha << std::endl;

        y0 = y0 + alpha * s;
        std::cout << "y0 = \n" << y0.format(OctaveFmt) << std::endl;

        r = r + alpha * C * s;
        std::cout << "r0 = \n" << r.format(OctaveFmt) << std::endl;

        beta = -1.0 * r.dot(C * s) / (C * s).dot(s);
        std::cout << "beta = " << beta << std::endl;

        s = r + beta * s;
        std::cout << "s = " << s.format(OctaveFmt) << std::endl;
        std::cout << "rnorm = " << r.norm() << std::endl;
    }

    return y0;
}

Eigen::VectorXd conjugateGradientPp(Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd x0) {
    Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    double eps = 1e-3;
    int N_max = 2;

    Eigen::MatrixXd M = A.diagonal().asDiagonal();
    Eigen::MatrixXd Minv = A.diagonal().asDiagonal().inverse();

    Eigen::VectorXd g = Minv.cwiseSqrt() * b;
    Eigen::VectorXd y = M.cwiseSqrt() * x0;

    Eigen::VectorXd s = (Minv.cwiseSqrt() * A * Minv.cwiseSqrt()) * y - g;

    Eigen::VectorXd x = Minv.cwiseSqrt() * y;
    Eigen::VectorXd sigma = Minv.cwiseSqrt() * s;
    Eigen::VectorXd rho_0 = M.cwiseSqrt() * s;
    Eigen::VectorXd rho_k(rho_0.size());

    double alpha, beta;
    for (int i = 0; i < N_max; i++) {
        alpha = -1.0 * (Minv * rho_0).dot(rho_0) / (A * sigma).dot(sigma);
        x = x + alpha * sigma;
        rho_k = rho_0 + alpha * A * sigma;
        beta = -1.0 * (Minv * rho_k).dot(rho_k) / (Minv * rho_0).dot(rho_0);
        sigma = Minv * rho_k + beta * sigma;
        rho_0 = rho_k;
        std::cout << x.format(OctaveFmt) << std::endl;
        std::cout << (Minv.cwiseSqrt().inverse() * x).format(OctaveFmt) << std::endl;
    }

    return Minv.cwiseSqrt().inverse() * x;
}
