#pragma once

#include "../Utilities/EngineMath.h"

namespace EngineMathLib {

    /**
    * @class Matrix2x2
    * @brief Represents a 2x2 matrix for linear algebra operations in 2D space.
    *
    * This class provides support for common matrix operations such as addition,
    * subtraction, multiplication (by scalar and matrix), transposition, inversion,
    * and determinant calculation. It is typically used for
    * 2D linear transformations
    * like scaling and rotation. Internally, it uses `double` precision and a small
    * EPSILON value for floating-point comparisons.
    */

    class
        Matrix2x2 {
    public:
        double m00, m01, m10, m11;

        /**
         * @brief Default constructor.
         * Initializes the matrix as an identity matrix.
         */
        Matrix2x2()
            : m00(1.0), m01(0.0), m10(0.0), m11(1.0) {
        }

        /**
         * @brief Constructs a matrix with specific values.
         * @param a Element at row 0, column 0.
         * @param b Element at row 0, column 1.
         * @param c Element at row 1, column 0.
         * @param d Element at row 1, column 1.
         */
        Matrix2x2(double a, double b, double c, double d)
            : m00(a), m01(b), m10(c), m11(d) {
        }

        /**
         * @brief Creates a scaling matrix.
         * @param scaleX Scaling factor on the X axis.
         * @param scaleY Scaling factor on the Y axis.
         * @return Matrix2x2 A matrix that scales vectors by scaleX and scaleY.
         */
        static Matrix2x2
            Scale(double scaleX, double scaleY) {
            return Matrix2x2(scaleX, 0.0, 0.0, scaleY);
        }

        /**
         * @brief Creates a rotation matrix.
         * @param angle Angle in radians to rotate counter-clockwise.
         * @return Matrix2x2 A matrix that rotates vectors by the given angle.
         */
        static Matrix2x2
            Rotate(double angle) {
            double cosA = cos(angle);
            double sinA = sin(angle);
            return Matrix2x2(cosA, -sinA, sinA, cosA);
        }

        /**
         * @brief Returns the transpose of the matrix.
         * @return Matrix2x2 The transposed matrix.
         */
        Matrix2x2
            transpose() const {
            return Matrix2x2(m00, m10, m01, m11);
        }

        /**
         * @brief Computes the determinant of the matrix.
         * @return double The determinant value.
         */
        double
            determinant() const {
            return (m00 * m11) - (m01 * m10);
        }

        /**
         * @brief Computes the inverse of the matrix.
         * @return Matrix2x2 The inverse matrix if invertible,
         * identity matrix otherwise.
         */
        Matrix2x2
            inverse() const {
            double det = determinant();
            if (fabs(det) < EPSILON) return Matrix2x2();  // Fallback to identity
            double invDet = 1.0 / det;
            return Matrix2x2(
                m11 * invDet, -m01 * invDet,
                -m10 * invDet, m00 * invDet
            );
        }

        /**
         * @brief Matrix multiplication.
         * @param other The matrix to multiply with.
         * @return Matrix2x2 The resulting matrix after multiplication.
         */
        Matrix2x2
            operator*(const Matrix2x2& other) const {
            return Matrix2x2(
                m00 * other.m00 + m01 * other.m10,
                m00 * other.m01 + m01 * other.m11,
                m10 * other.m00 + m11 * other.m10,
                m10 * other.m01 + m11 * other.m11
            );
        }

        /**
         * @brief Matrix addition.
         * @param other The matrix to add.
         * @return Matrix2x2 The resulting matrix after addition.
         */
        Matrix2x2
            operator+(const Matrix2x2& other) const {
            return Matrix2x2(
                m00 + other.m00, m01 + other.m01,
                m10 + other.m10, m11 + other.m11
            );
        }

        /**
         * @brief Matrix subtraction.
         * @param other The matrix to subtract.
         * @return Matrix2x2 The resulting matrix after subtraction.
         */
        Matrix2x2
            operator-(const Matrix2x2& other) const {
            return Matrix2x2(
                m00 - other.m00, m01 - other.m01,
                m10 - other.m10, m11 - other.m11
            );
        }

        /**
         * @brief Scalar multiplication.
         * @param scalar The scalar to multiply each element by.
         * @return Matrix2x2 The resulting scaled matrix.
         */
        Matrix2x2
            operator*(double scalar) const {
            return Matrix2x2(
                m00 * scalar, m01 * scalar,
                m10 * scalar, m11 * scalar
            );
        }

        /**
         * @brief Scalar division.
         * @param scalar The scalar to divide each element by.
         * @return Matrix2x2 The resulting matrix after division,
         ( or identity if scalar is near zero.
         */
        Matrix2x2
            operator/(double scalar) const {
            if (fabs(scalar) < EPSILON) return Matrix2x2(); // Fallback to identity
            double inv = 1.0 / scalar;
            return (*this) * inv;
        }

        /**
         * @brief Equality comparison between two matrices.
         * @param other The matrix to compare with.
         * @return true If all elements are approximately equal.
         * @return false Otherwise.
         */
        bool
            operator==(const Matrix2x2& other) const {
            return fabs(m00 - other.m00) < EPSILON &&
                fabs(m01 - other.m01) < EPSILON &&
                fabs(m10 - other.m10) < EPSILON &&
                fabs(m11 - other.m11) < EPSILON;
        }

        /**
         * @brief Inequality comparison between two matrices.
         * @param other The matrix to compare with.
         * @return true If any element differs beyond EPSILON.
         * @return false If all elements are approximately equal.
         */
        bool
            operator!=(const Matrix2x2& other) const {
            return !(*this == other);
        }
    };

}