//
// Created by Austin Clyde on 5/28/18.
//

#include "FEMLaplacian.h"

FEMLaplacian::FEMLaplacian() {
    ;
}

FEMLaplacian::FEMLaplacian(VariMesh mesh, std::function<double(int, int, VariMesh)> u_init, std::string file,
                           std::function<Eigen::VectorXd(Eigen::VectorXd, VariMesh)> Amul,
                           std::function<Eigen::VectorXd(VariMesh)> GenDiag) {
    this->mesh = mesh;
    this->file = file;
    this->Amul = Amul;

    Eigen::MatrixXd init_m(mesh.x_nodes, mesh.y_nodes);
    for (int i = 0; i < mesh.x_nodes; i++) {
        for (int j = 0; j < mesh.y_nodes; j++) {
            init_m(i, j) = u_init(i, j, mesh);
        }
    }

    b = Eigen::Map<Eigen::VectorXd>(init_m.data(), init_m.size());
    U = Eigen::VectorXd(b.size());
    D = GenDiag(mesh);
}

void FEMLaplacian::run() {
    Eigen::VectorXd x0(mesh.x_nodes * mesh.y_nodes);
    Eigen::VectorXd U = conjugateGradientPreconditioningNM(x0);
    write_coords(U, mesh, file);
}

void FEMLaplacian::write_out_mesh(std::string file) {
    write_coords(b, mesh, file);
}

Eigen::VectorXd FEMLaplacian::conjugateGradientPreconditioningNM(Eigen::VectorXd x0) {
    Eigen::VectorXd M = D;
    Eigen::VectorXd Minv = D.cwiseInverse();

    Eigen::VectorXd x = x0;
    Eigen::VectorXd rho_0 = Amul(x0, mesh) - b;
    Eigen::VectorXd sigma = Minv.cwiseProduct(rho_0);
    Eigen::VectorXd rho_k(rho_0.size());
    double init_norm = rho_0.norm();
    std::cout << "init norm: " << rho_0.norm() << std::endl;
    double alpha, beta;
    int i = 0;
    for (i = 0; i < ITER_MAX && rho_0.norm() >= init_norm / 1000; i++) {
        alpha = -1.0 * (Minv.cwiseProduct(rho_0)).dot(rho_0) / (Amul(sigma, mesh)).dot(sigma);
        x = x + alpha * sigma;
        rho_k = rho_0 + alpha * Amul(sigma, mesh);
        beta = (Minv.cwiseProduct(rho_k)).dot(rho_k) / (Minv.cwiseProduct(rho_0)).dot(rho_0);
        sigma = Minv.cwiseProduct(rho_k) + beta * sigma;
        rho_0 = rho_k;
        std::cout << "i = " << i << ", rho norm = " << rho_k.norm() << std::endl;
    }
    std::cout << "ITERS = " << i << std::endl;
    return x;
}

Eigen::VectorXd multiplyStiffExample0(Eigen::VectorXd v, VariMesh mesh) {
    Eigen::MatrixXd Y(mesh.x_nodes, mesh.y_nodes);
    Eigen::MatrixXd Z = Eigen::Map<Eigen::MatrixXd>(v.data(), mesh.x_nodes, mesh.y_nodes);
    double r, a, b;
    for (int i = 0; i < mesh.x_nodes - 1; i++) {
        for (int j = 0; j < mesh.y_nodes - 1; j++) {
            r = mesh.dy[j] / mesh.dx[i];
            a = 1.0 / (2.0 * r);
            b = r / 2.0;

            Y(i, j) += a * (Z(i, j) - Z(i, j + 1)) + b * (Z(i, j) - Z(i + 1, j));
            Y(i + 1, j) += a * (Z(i + 1, j) - Z(i + 1, j + 1)) + b * (Z(i + 1, j) - Z(i, j));
            Y(i, j + 1) += a * (Z(i, j + 1) - Z(i, j)) + b * (Z(i, j + 1) - Z(i + 1, j + 1));
            Y(i + 1, j + 1) += a * (Z(i + 1, j + 1) - Z(i + 1, j)) + b * (Z(i + 1, j + 1) - Z(i, j + 1));
        }
    }


    for (int i = 0; i < mesh.x_nodes; i++) {
        Y(i, 0) = Z(i, 0);
        Y(i, mesh.y_nodes - 1) = Z(i, mesh.y_nodes - 1);
    }
    for (int j = 0; j < mesh.y_nodes; j++) {
        Y(0, j) = Z(0, j);
        Y(mesh.x_nodes - 1, j) = Z(mesh.x_nodes - 1, j);
    }

    return Eigen::Map<Eigen::VectorXd>(Y.data(), Y.size());
}

Eigen::VectorXd multiplyStiffExample1(Eigen::VectorXd v, VariMesh mesh) {
    Eigen::MatrixXd Y(mesh.x_nodes, mesh.y_nodes);
    Eigen::MatrixXd Z = Eigen::Map<Eigen::MatrixXd>(v.data(), mesh.x_nodes, mesh.y_nodes);
    double r, a, b;
    for (int i = 0; i < mesh.x_nodes - 1; i++) {
        for (int j = 0; j < mesh.y_nodes - 1; j++) {
            r = mesh.dy[j] / mesh.dx[i];
            a = 1.0 / (2.0 * r);
            b = r / 2.0;

            Y(i, j) += a * (Z(i, j) - Z(i, j + 1)) + b * (Z(i, j) - Z(i + 1, j));
            Y(i + 1, j) += a * (Z(i + 1, j) - Z(i + 1, j + 1)) + b * (Z(i + 1, j) - Z(i, j));
            Y(i, j + 1) += a * (Z(i, j + 1) - Z(i, j)) + b * (Z(i, j + 1) - Z(i + 1, j + 1));
            Y(i + 1, j + 1) += a * (Z(i + 1, j + 1) - Z(i + 1, j)) + b * (Z(i + 1, j + 1) - Z(i, j + 1));
        }
    }


    for (int i = 0; i < mesh.x_nodes; i++) {
        Y(i, 0) = Z(i, 0);
        Y(i, mesh.y_nodes - 1) = Z(i, mesh.y_nodes - 1);
    }
    for (int j = 0; j < mesh.y_nodes; j++) {
        Y(0, j) = Z(0, j);
        Y(mesh.x_nodes - 1, j) = Z(mesh.x_nodes - 1, j);
    }

    return Eigen::Map<Eigen::VectorXd>(Y.data(), Y.size());
}

Eigen::VectorXd generateStiffD(VariMesh mesh) {
    Eigen::MatrixXd Y(mesh.x_nodes, mesh.y_nodes);

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