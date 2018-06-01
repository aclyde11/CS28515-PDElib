//
// Created by Austin Clyde on 5/28/18.
//

#include "Problem2d.h"

Problem2d::Problem2d() {
    mesh = VariMesh(1, 1, 1, 1);

}

Eigen::VectorXd multiplyStiff(Eigen::VectorXd v, VariMesh mesh) {
    Eigen::MatrixXd Y(mesh.x_nodes, mesh.y_nodes);
    for (int i = 0; i < mesh.x_nodes; i++) {
        for (int j = 0; j < mesh.y_nodes; j++)
            Y(i, j) = 0;
    }

    Eigen::MatrixXd Z = Eigen::Map<Eigen::MatrixXd>(v.data(), mesh.x_nodes, mesh.y_nodes);
    double r, a, b;
    for (int i = 0; i < mesh.x_nodes - 1; i++) {
        for (int j = 0; j < mesh.y_nodes - 1; j++) {
            r = mesh.dy[j] / mesh.dx[i];
            a = 1.0 / (2.0 * r);
            b = r / 2.0;

            Y(i, j) += a * (Z(i, j) - Z(i, j + 1)) + b * (Z(i, j) - Z(i + 1, j + 1));
            Y(i + 1, j) += a * (Z(i + 1, j) - Z(i + 1, j + 1)) + b * (Z(i + 1, j) - Z(i, j));
            Y(i, j + 1) += a * (Z(i, j + 1) - Z(i + 1, j + 1)) + b * (Z(i, j + 1) - Z(i, j));
            Y(i + 1, j + 1) += a * (Z(i + 1, j + 1) - Z(i + 1, j + 1)) + b * (Z(i + 1, j + 1) - Z(i + 1, j + 1));
        }
    }
    return Eigen::Map<Eigen::VectorXd>(Y.data(), Y.size());
}

Eigen::VectorXd generateStiffD(VariMesh mesh) {
    Eigen::MatrixXd Y(mesh.x_nodes, mesh.y_nodes);
    for (int i = 0; i < mesh.x_nodes; i++) {
        for (int j = 0; j < mesh.y_nodes; j++)
            Y(i, j) = 0;
    }

    double r, a, b;
    for (int i = 0; i < mesh.x_nodes - 1; i++) {
        for (int j = 0; j < mesh.y_nodes - 1; j++) {
            r = mesh.dy[j] / mesh.dx[i];
            a = 1.0 / (2.0 * r);
            b = r / 2.0;

            Y(i, j) += 1.0 / (2.0 * r) + r / 2.0;
            Y(i + 1, j) += 1.0 / (2.0 * r) + r / 2.0;
            Y(i, j + 1) += 1.0 / (2.0 * r) + r / 2.0;
            Y(i + 1, j + 1) += 1.0 / (2.0 * r) + r / 2.0;
        }
    }
    return Eigen::Map<Eigen::VectorXd>(Y.data(), Y.size());
}

Eigen::VectorXd conjugateGradient(Eigen::MatrixXd C, Eigen::VectorXd g, Eigen::VectorXd y0) {
    Eigen::VectorXd s = C * y0 - g;
    Eigen::VectorXd r(s), e(s.size());
    double alpha, beta;

    for (int i = 0; i < ITER_MAX && r.norm() > EPS; i++) {
        alpha = -1.0 * (r.dot(s) / (C * s).dot(s));
        y0 = y0 + alpha * s;
        r = r + alpha * C * s;
        beta = -1.0 * r.dot(C * s) / (C * s).dot(s);
        s = r + beta * s;
    }

    return y0;
}

Eigen::VectorXd conjugateGradientPreconditioning(Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd x0) {
    Eigen::MatrixXd M = A.diagonal().asDiagonal();
    Eigen::MatrixXd Minv = M.inverse();

    Eigen::VectorXd g = Minv.cwiseSqrt() * b;
    Eigen::VectorXd y = M.cwiseSqrt() * x0;

    Eigen::VectorXd s = Minv * y - g;

    Eigen::VectorXd x = Minv.cwiseSqrt() * y;
    Eigen::VectorXd sigma = Minv.cwiseSqrt() * s;
    Eigen::VectorXd rho_0 = M.cwiseSqrt() * s;
    Eigen::VectorXd rho_k(rho_0.size());

    double alpha, beta;
    int i = 0;
    for (i = 0; i < ITER_MAX && rho_0.norm() > EPS; i++) {
        alpha = -1.0 * (Minv * rho_0).dot(rho_0) / (A * sigma).dot(sigma);
        x = x + alpha * sigma;
        rho_k = rho_0 + alpha * A * sigma;
        beta = (Minv * rho_k).dot(rho_k) / (Minv * rho_0).dot(rho_0);
        sigma = Minv * rho_k + beta * sigma;
        rho_0 = rho_k;
    }
    return x;
}

Eigen::VectorXd conjugateGradientPreconditioningNM(Eigen::VectorXd b, Eigen::VectorXd x0,
                                                   std::function<Eigen::VectorXd(Eigen::VectorXd, VariMesh)> Amul,
                                                   Eigen::VectorXd D, VariMesh mesh) {
    Eigen::VectorXd M = D;
    Eigen::VectorXd Minv = D.cwiseInverse();

    Eigen::VectorXd g = Minv.cwiseSqrt().cwiseProduct(b);
    Eigen::VectorXd y = M.cwiseSqrt().cwiseProduct(x0);
    Eigen::VectorXd s = Minv.cwiseProduct(y) - g;
    Eigen::VectorXd x = Minv.cwiseSqrt().cwiseProduct(y);
    Eigen::VectorXd sigma = Minv.cwiseSqrt().cwiseProduct(s);
    Eigen::VectorXd rho_0 = M.cwiseSqrt().cwiseProduct(s);
    Eigen::VectorXd rho_k(rho_0.size());

    double alpha, beta;
    int i = 0;
    for (i = 0; i < ITER_MAX && rho_0.norm() > EPS; i++) {
        alpha = -1.0 * (Minv.cwiseProduct(rho_0)).dot(rho_0) / (Amul(sigma, mesh)).dot(sigma);
        x = x + alpha * sigma;
        rho_k = rho_0 + alpha * Amul(sigma, mesh);
        beta = (Minv.cwiseProduct(rho_k)).dot(rho_k) / (Minv.cwiseProduct(rho_0)).dot(rho_0);
        sigma = Minv.cwiseProduct(rho_k) + beta * sigma;
        rho_0 = rho_k;
    }
    return x;
}

void write_coords(Eigen::VectorXd X, VariMesh mesh, std::string file) {
    std::ofstream ofs(file, std::fstream::out);

    if (!ofs) {
        std::cout << "Error opening file for output" << std::endl;
        return;
    }

    Eigen::MatrixXd out = Eigen::Map<Eigen::MatrixXd>(X.data(), mesh.x_nodes, mesh.y_nodes);
    for (int i = 0; i < mesh.x_nodes; i++) {
        for (int j = 0; j < mesh.y_nodes; j++) {
            ofs << i * mesh.dx[0] << "," << j * mesh.dy[0] << "," << out(i, j) << "\n";
        }
    }
    ofs.close();
}

