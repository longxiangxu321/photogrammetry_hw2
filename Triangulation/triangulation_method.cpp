#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>

using namespace easy3d;

class RT {
public:
    Matrix33 R;
    Vector3D t;

    RT(Matrix matrix, Vector3D vector) {
        R = matrix;
        t = vector;
    }
};

bool get_normalization_matrix(const unsigned long long &num_pairs,
                              const std::vector<Vector2D> &points,
                              Matrix &N) {
    // normalize the 2d points
    // centroid of points
    Vector2D c(0, 0);
    for (const Vector2D &point: points) {
        c = c + point;
    }
    c = c / num_pairs;
    // average distance from centroid
    double d = 0;
    for (const Vector2D &point: points) {
        d = d + distance(c, point);
    }
    d = d / num_pairs;
    // scaling factor
    double s = sqrt(2) / d;
    // scaling matrix
    Matrix33 S(s, 0, 0,
               0, s, 0,
               0, 0, 1);
    // translation matrix
    Matrix33 T(1, 0, -c.x(),
               0, 1, -c.y(),
               0, 0, 1);
    // normalization matrix
    N = S * T;
    return true;
}

bool get_fundamental_matrix(const unsigned long long &num_pairs,
                            const std::vector<Vector2D> &points_0,
                            const std::vector<Vector2D> &points_1,
                            Matrix33 &F) {

    Matrix33 N_0;
    Matrix33 N_1;
    get_normalization_matrix(num_pairs, points_0, N_0);
    get_normalization_matrix(num_pairs, points_1, N_1);
    // std::cout << "Normalization matrix for first camera is: " << N_0 << std::endl;
    // std::cout << "Normalization matrix for second camera is: " << N_1 << std::endl;

    // normalization execution
    std::vector<Vector2D> q_0;
    std::vector<Vector2D> q_1;
    for (int i = 0; i < num_pairs; i++) {
        q_0.emplace_back(N_0 * (points_0[i].homogeneous()));
        q_1.emplace_back(N_1 * (points_1[i].homogeneous()));
    }

    // use svd to get the fundamental matrix with full rank
    Matrix W(num_pairs, 9);
    for (int i = 0; i < num_pairs; i++) {
        double arr[] = {q_0[i].x() * q_1[i].x(),
                        q_0[i].y() * q_1[i].x(),
                        q_1[i].x(),
                        q_0[i].x() * q_1[i].y(),
                        q_0[i].y() * q_1[i].y(),
                        q_1[i].y(),
                        q_0[i].x(),
                        q_0[i].y(),
                        1
        };
        Vector row_i(9, arr);
        W.set_row(i, row_i);
    }
    Matrix U(num_pairs, num_pairs, 0.0);
    Matrix D(num_pairs, 9, 0.0);
    Matrix V(9, 9, 0.0);
    svd_decompose(W, U, D, V);
    Vector f = V.get_column(8);
    Matrix33 F_full(f[0], f[1], f[2],
                    f[3], f[4], f[5],
                    f[6], f[7], f[8]);

    // use svd to get the fundamental matrix with rank two
    Matrix U_full(3, 3, 0.0);
    Matrix D_full(3, 3, 0.0);
    Matrix V_full(3, 3, 0.0);
    svd_decompose(F_full, U_full, D_full, V_full);
    auto D_two = D_full;
    D_two.set_row(2, {0, 0, 0});
    Matrix F_two = U_full * D_two * V_full.transpose();

    // de-normalize the fundamental matrix
    F = N_1.transpose() * F_two * N_0;
    Matrix ff(9, 1, 0.0);
    return true;
}

bool get_relative_pose(const unsigned long long &num_pairs,
                       const std::vector<Vector2D> &points_0,
                       const std::vector<Vector2D> &points_1,
                       const Matrix33 &F,
                       const Matrix33 &K_0,
                       const Matrix33 &K_1,
                       std::vector<RT> &Rts) {
    Matrix33 E = K_1.transpose() * F * K_0;

    // use svd to get the relative poses
    Matrix U(3, 3, 0.0);
    Matrix D(3, 3, 0.0);
    Matrix V(3, 3, 0.0);
    svd_decompose(E, U, D, V);
    Matrix33 W(0, -1, 0, 1, 0, 0, 0, 0, 1);
    Matrix R_0 = determinant(U * W * V.transpose()) * (U * W * V.transpose());
    Matrix R_1 = determinant(U * W.transpose() * V.transpose()) * (U * W.transpose() * V.transpose());
    // std::cout << R_0 << std::endl;
    // std::cout << R_1 << std::endl;
    Vector3D t_0 = U * Vector3D(0, 0, 1);
    Vector3D t_1 = -t_0;

    // four candidate poses
    Rts.emplace_back(RT(R_0, t_0));
    Rts.emplace_back(RT(R_0, t_1));
    Rts.emplace_back(RT(R_1, t_0));
    Rts.emplace_back(RT(R_1, t_1));
    // std::cout << "R_0: " << R_0 << std::endl;
    // std::cout << "R_1: " << R_1 << std::endl;
    // std::cout << "t_0: " << t_0 << std::endl;
    // std::cout << "t_1: " << t_1 << std::endl;
    return true;
};

bool get_3d_coordinates(const unsigned long long &num_pairs,
                        const std::vector<Vector2D> &points_0,
                        const std::vector<Vector2D> &points_1,
                        const Matrix33 &K_0,
                        const Matrix33 &K_1,
                        const std::vector<RT> &Rts,
                        std::vector<Vector3D> &points_3d,
                        Matrix34 &M_0,
                        Matrix34 &M_1,
                        Matrix33 &R,
                        Vector3D &t) {
    Matrix34 Rt_0(1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1, 0);
    M_0 = K_0 * Rt_0;

    int max_count = 0;

    // try each pose
    for (int i = 0; i < 4; i++) {
        int count = 0;
        auto R_i = Rts[i].R;
        auto t_i = Rts[i].t;
        Matrix34 Rt_1;

        // the extrinsic parameters of the second camera
        Rt_1.set_column(0, R_i.get_column(0));
        Rt_1.set_column(1, R_i.get_column(1));
        Rt_1.set_column(2, R_i.get_column(2));
        Rt_1.set_column(3, t_i);
        Matrix34 M_temp = K_1 * Rt_1;
        std::vector<Vector3D> points_3d_temp;
        for (int j = 0; j < num_pairs; j++) {
            // use svd to do triangulation
            Matrix A(4, 4);
            A.set_row(0, points_0[j].x() * M_0.get_row(2) - M_0.get_row(0));
            A.set_row(1, points_0[j].y() * M_0.get_row(2) - M_0.get_row(1));
            A.set_row(2, points_1[j].x() * M_temp.get_row(2) - M_temp.get_row(0));
            A.set_row(3, points_1[j].y() * M_temp.get_row(2) - M_temp.get_row(1));
            Matrix U(4, 4, 0.0);
            Matrix D(4, 4, 0.0);
            Matrix V(4, 4, 0.0);
            svd_decompose(A, U, D, V);
            Vector4D P_0_h = V.get_column(3);
            Vector3D P_0_c = P_0_h.cartesian();
            points_3d_temp.emplace_back(P_0_c);
            Vector3D P_1_c = R_i * P_0_c + t_i;

            // accumulate the case that meet the constraint
            if (P_0_c.z() > 0 && P_1_c.z() > 0) {
                count++;
            }
        }
        std::cout << "For combination " << i << " The count of correct pairs is: " << count << std::endl;
        // replace the pose with the best one
        if (max_count < count) {
            max_count = count;
            points_3d = points_3d_temp;
            R = R_i;
            t = t_i;
            M_1 = M_temp;
        }
    }
    return true;
};

bool get_difference(const unsigned long long &num_pairs,
                    const Matrix34 &M_0,
                    const Matrix34 &M_1,
                    const std::vector<Vector2D> &points_0,
                    const std::vector<Vector2D> &points_1,
                    const std::vector<Vector3D> &points_3d
) {

    double MSE_0;
    double MSE_1;
    for (int i = 0; i < num_pairs; i++) {
        // reproject the 3d points on image planes
        Vector3D p_0_h = M_0 * points_3d[i].homogeneous();
        Vector3D p_1_h = M_1 * points_3d[i].homogeneous();
        Vector2D p_0_c = p_0_h.cartesian();
        Vector2D p_1_c = p_1_h.cartesian();
        double dist_0 = distance(p_0_c, points_0[i]);
        double dist_1 = distance(p_1_c, points_1[i]);
        // accumulate the square differences
        MSE_0 = MSE_0 + dist_0 * dist_0 / num_pairs;
        MSE_1 = MSE_1 + dist_1 * dist_1 / num_pairs;
    }
    // calculate the RMSE
    double RMSE_0 = sqrt(MSE_0);
    double RMSE_1 = sqrt(MSE_1);
    double RMSE = sqrt(pow(RMSE_0, 2) + pow(RMSE_1, 2));
    std::cout << "RMSE_0 " << RMSE_0 << std::endl;
    std::cout << "RMSE_1 " << RMSE_1 << std::endl;
    std::cout << "RMSE " << RMSE << std::endl;
    return true;
}

bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const {

    Matrix33 K_0(fx, 0, cx, 0, fy, cy, 0, 0, 1);
    Matrix33 K_1 = K_0;

    // std::cout << "Intrinsic parameters matrix is: " << K_0 << std::endl;

    for (int i = 8; i < 161; i++) {
        std::vector<Vector2D> points_00;
        std::vector<Vector2D> points_11;
        for (int j = 0; j < i; j++) {
            points_00.emplace_back(points_0[j]);
            points_11.emplace_back(points_1[j]);
        }
        unsigned long long num_pairs = points_00.size();
        std::cout << "The numer of points used is: " << num_pairs << std::endl;

        // get fundamental matrix
        Matrix33 F;
        get_fundamental_matrix(num_pairs, points_00, points_11, F);

        // std::cout << "Fundamental matrix is: " << F << std::endl;

        // get candidate relative poses
        std::vector<RT> Rts;
        get_relative_pose(num_pairs, points_00, points_11, F, K_0, K_1, Rts);

        // get 3d coordinates
        Matrix34 M_0;
        Matrix34 M_1;
        get_3d_coordinates(num_pairs, points_00, points_11, K_0, K_1, Rts, points_3d, M_0, M_1,
                           R, t);

        // std::cout << "The chosen relative pose is,R: " << R << " t: " << t << std::endl;

        // evaluate the re-projection error
        get_difference(num_pairs, M_0, M_1, points_00, points_11, points_3d);

    }

    return true;
}