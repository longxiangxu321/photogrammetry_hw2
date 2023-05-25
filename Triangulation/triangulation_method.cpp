/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;


/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{

    /// define a 2D vector/point
    Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    Vector3D a(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    Vector2D p = a.cartesian();

    /// get the Homogeneous coordinates of p
    Vector3D q = p.homogeneous();

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    Matrix33 A;

    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    Matrix34 M(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    Matrix WW(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    WW.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    int num_rows = WW.rows();

    /// get the number of columns.
    int num_cols = WW.cols();

    /// get the the element at row 1 and column 2
    double value = WW(1, 2);

    /// get the last column of a matrix
    Vector last_column = WW.get_column(WW.cols() - 1);

    /// define a 3 by 3 identity matrix
    Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.

    if (points_0.size()<8)
    {
        std::cout<<"At least 8 points are needed for recover fundamental matrix"<<std::endl;
        return false;
    }

    if (points_0.size()!=points_1.size())
    {
        std::cout<<"Point number not equal for two images"<<std::endl;
    }

    /// Normalization
    Vector3D sample_center_0(0, 0, 0);   ///c
    Vector3D sample_center_1(0, 0, 0);
    double avg_dist_0 = 0; ///d
    double avg_dist_1 = 0;
    double scale_fac_0, scale_fac_1;

    for (int i = 0; i<points_0.size(); i++) {
        sample_center_0 += points_0[i];
        sample_center_1 += points_1[i];
    }


    for (int i = 0; i< points_0.size(); i++) {
        double temp_dist_0 = distance(sample_center_0, points_0[i]);
        double temp_dist_1 = distance(sample_center_1, points_1[i]);
        avg_dist_0 += temp_dist_0;
        avg_dist_1 += temp_dist_1;
    }


    scale_fac_0 = sqrt(2) / avg_dist_0; ///s
    scale_fac_1 = sqrt(2) / avg_dist_1;
    Matrix33 normalization_T_0 = (scale_fac_0, 0, -sample_center_0.x(),
                                        0, scale_fac_0, -sample_center_0.y(),
                                        0, 0, 1);
    Matrix33 normalization_T_1 = (scale_fac_1, 0, -sample_center_1.x(),
                                        0, scale_fac_1, -sample_center_1.y(),
                                        0, 0, 1);

    std::vector<Vector2D> normalized_points_0;
    std::vector<Vector2D> normalized_points_1;

    for (int i = 0; i< points_0.size(); i++) {
        Vector2D temp_pt_0 = normalization_T_0 * (points_0[i]);
        Vector2D temp_pt_1 = normalization_T_1 * (points_1[i]);
        normalized_points_0.emplace_back(temp_pt_0);
        normalized_points_1.emplace_back(temp_pt_1);
    }

    Matrix W(points_0.size(), 9);
    for (int i = 0; i<points_0.size(); i++) {
        double arr[] = {normalized_points_0[i].x() * normalized_points_1[i].x(), //u0 u1
                          normalized_points_0[i].y() * normalized_points_1[i].x(), //v0 u1
                          normalized_points_1[i].x(), //u1
                          normalized_points_0[i].x() * normalized_points_1[i].y(), //u0 v1
                          normalized_points_0[i].y() * normalized_points_1[i].y(), //v0 v1
                          normalized_points_1[i].y(), //v1
                          normalized_points_0[i].x(), //u0
                          normalized_points_0[i].y(), //v0
                          1
                          };
        Vector row_0(9, arr);
        W.set_row(i, row_0);
    }

    Matrix U, D, V;
    svd_decompose(W, U, D, V);
    Vector f = V.get_column(8);

    Matrix33 F_hat = (f[0], f[1], f[2],
                    f[3], f[4], f[5],
                    f[6], f[7], f[8]);


    Matrix U_F, D_F, V_F;
    svd_decompose(F_hat, U_F, D_F, V_F);
    D_F.set_row(2, {0,0,0});
    Matrix F = U_F * D_F *V_F;

    return points_3d.size() > 0;
}