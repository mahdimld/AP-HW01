#include <iostream>
#include <random>
#include "hw1.h"
#include <iomanip>
#include <cmath>

Matrix algebra::zeros(size_t n, size_t m)
{

    Matrix output{};

    std::vector<double> op{};

    for (size_t i{}; i < n; i++)
    {

        for (size_t j{}; j < m; j++)
        {

            op.push_back(0);
        }

        output.push_back(op);
        op.clear();
    }

    return output;
}

Matrix algebra::ones(size_t n, size_t m)
{

    Matrix output{};

    std::vector<double> op{};

    for (size_t i{}; i < n; i++)
    {

        for (size_t j{}; j < m; j++)
        {

            op.push_back(1);
        }

        output.push_back(op);
        op.clear();
    }

    return output;
}

Matrix algebra::random(size_t n, size_t m, double min, double max)
{

    Matrix output{};
    std::vector<double> op{};

    if (min > max)
    {
        throw std::logic_error("Caution: min cannot be greater than max");
    }

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(min, max);

    for (size_t i{}; i < n; i++)
    {
        for (size_t j{}; j < m; j++)
        {

            op.push_back(dist(mt));
        }

        output.push_back(op);
        op.clear();
    }

    return output;
}

void algebra::show(const Matrix &matrix)
{

    for (size_t i{}; i < matrix.size(); i++)
    {
        for (size_t j{}; j < matrix[0].size(); j++)
        {

            std::cout << std::setw(3) << std::setprecision(3) << matrix[i][j] << " ";
        }

        std::cout << std::endl;
    }
}

Matrix algebra::multiply(const Matrix &matrix, double c)
{

    std::vector<double> op{};
    Matrix output{};
    if (matrix.size() == 0)
    {
        return output;
    }

    for (size_t i{}; i < matrix.size(); i++)
    {
        for (size_t j{}; j < matrix[0].size(); j++)
        {

            op.push_back(c * matrix[i][j]);
        }

        output.push_back(op);
        op.clear();
    }

    return output;
}

Matrix algebra::multiply(const Matrix &matrix1, const Matrix &matrix2)
{
    std::vector<double> op{};
    Matrix output{};

    if (matrix1.size() == 0 && matrix2.size() == 0)
    {
        return output;
    }
    size_t n{matrix1.size()};
    size_t m{matrix1[0].size()};
    size_t q{matrix2.size()};
    size_t p{matrix2[0].size()};

    double sum{};

    if (n == 0 || m == 0 || q == 0 || p == 0)
    {

        throw std::logic_error("Caution: multiplication of 2 empty matrix");
    }

    if (m != q)
    {

        throw std::logic_error("Caution: matrices with wrong dimensions cannot be multiplied");
    }
    else
    {

        for (size_t i{}; i < n; i++)
        {
            for (size_t j{}; j < p; j++)
            {
                for (size_t k{}; k < m; k++)
                {
                    sum += (matrix1[i][k] * matrix2[k][j]);
                }
                op.push_back(sum);
                sum = 0;
            }

            output.push_back(op);
            op.clear();
        }
    }

    return output;
}

Matrix algebra::sum(const Matrix &matrix, double c)
{

    std::vector<double> op{};
    Matrix output{};

    if (matrix.size() == 0)
    {

        return output;
    }

    for (size_t i{}; i < matrix.size(); i++)
    {
        for (size_t j{}; j < matrix[0].size(); j++)
        {

            op.push_back(c + matrix[i][j]);
        }

        output.push_back(op);
        op.clear();
    }

    return output;
}

Matrix algebra::sum(const Matrix &matrix1, const Matrix &matrix2)
{

    std::vector<double> op{};
    Matrix output{};
    if (matrix1.empty() && matrix2.empty())
    {
        return output;
    }
    size_t height_of_matrix1{matrix1.size()};

    size_t height_of_matrix2{matrix2.size()};

    if (height_of_matrix1 != height_of_matrix2)
    {

        throw std::logic_error("Caution: matrices with wrong dimensions cannot be summed");
    }

    size_t length_of_matrix1{matrix1[0].size()};
    size_t length_of_matrix2{matrix2[0].size()};
    size_t n{height_of_matrix1}, m{length_of_matrix1};

    if (length_of_matrix1 != length_of_matrix2)
    {
        throw std::logic_error("Caution: matrices with wrong dimensions cannot be summed");
    }

    for (size_t i{}; i < n; i++)
    {
        for (size_t j{}; j < m; j++)
        {

            op.push_back(matrix1[i][j] + matrix2[i][j]);
        }

        output.push_back(op);
        op.clear();
    }

    return output;
}

Matrix algebra::transpose(const Matrix &matrix)
{
    std::vector<double> op{};
    Matrix output{};

    if (matrix.empty())
    {

        return output;
    }

    for (size_t i{}; i < matrix[0].size(); i++)
    {
        for (size_t j{}; j < matrix.size(); j++)
        {

            op.push_back(matrix[j][i]);
        }

        output.push_back(op);
        op.clear();
    }

    return output;
}

Matrix algebra::minor(const Matrix &matrix, size_t n, size_t m)
{

    std::vector<double> op{};
    Matrix output{};

    for (size_t i{}; i < matrix.size(); i++)
    {

        if (i != n)
        {

            for (size_t j{}; j < matrix[0].size(); j++)
            {

                if (j != m)
                {
                    op.push_back(matrix[i][j]);
                }
            }

            output.push_back(op);
            op.clear();
        }
    }

    return output;
}

double algebra::determinant(const Matrix &matrix)
{
    if (matrix.empty())
    {
        return 1.0;
    }

    double sum{};

    if (matrix.size() != matrix[0].size())
    {

        throw std::logic_error("Caution: non-square matrices have no determinant");
    }

    if (matrix.size() == 2 && matrix[0].size() == 2)
    {

        sum = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }
    else if (matrix.size() == 1 && matrix[0].size() == 1)
    {
        sum = matrix[0][0];
    }
    else
    {
        for (size_t j{}; j < matrix[0].size(); j++)
        {

            sum += matrix[0][j] * pow(-1, 0 + j) * determinant(minor(matrix, 0, j));
        }
    }
    return sum;
}

Matrix algebra::inverse(const Matrix &matrix)
{

    std::vector<double> op{};
    Matrix adj{}, adjT{}, inv{};
    double det{determinant(matrix)};
    if (matrix.empty())
    {
        return adj;
    }

    if (matrix.size() != matrix[0].size())
    {

        throw std::logic_error("Caution: non-square matrices have no inverse");
    }

    if (det == 0)
    {

        throw std::logic_error("Caution: singular matrices have no inverse");
    }

    for (size_t i{}; i < matrix.size(); i++)
    {
        for (size_t j{}; j < matrix[0].size(); j++)
        {

            op.push_back(pow(-1, i + j) * determinant(minor(matrix, i, j)));
        }

        adj.push_back(op);
        op.clear();
    }

    adjT = transpose(adj);
    inv = multiply(adjT, 1 / det);

    return inv;
}
Matrix algebra::concatenate(const Matrix &matrix1, const Matrix &matrix2, int axis = 0)
{

    std::vector<double> op{};
    Matrix output{};

    if (axis == 0)
    {

        if (matrix1[0].size() != matrix2[0].size())
        {

            throw std::logic_error("Caution: matrices with wrong dimensions cannot be concatenated");
        }
        for (size_t i{}; i < matrix1.size() + matrix2.size(); i++)
        {
            for (size_t j{}; j < matrix1[0].size(); j++)
            {

                if (i < matrix1.size())
                {
                    op.push_back(matrix1[i][j]);
                }
                else if (i >= matrix1.size() && i < matrix1.size() + matrix2.size())
                {
                    op.push_back(matrix2[i - matrix1.size()][j]);
                }
            }

            output.push_back(op);
            op.clear();
        }
    }

    else if (axis == 1)
    {

        if (matrix1.size() != matrix2.size())
        {

            throw std::logic_error("Caution: matrices with wrong dimensions cannot be concatenated");
        }

        for (size_t i{}; i < matrix1.size(); i++)
        {
            for (size_t j{}; j < matrix1[0].size() + matrix2[0].size(); j++)
            {

                if (j < matrix1[0].size())
                {
                    op.push_back(matrix1[i][j]);
                }
                else if (j >= matrix1[0].size() && j < matrix1[0].size() + matrix2[0].size())
                {
                    op.push_back(matrix2[i][j - matrix1[0].size()]);
                }
            }

            output.push_back(op);
            op.clear();
        }
    }

    return output;
}

Matrix algebra::ero_swap(const Matrix &matrix, size_t r1, size_t r2)
{
    std::vector<double> op{};
    Matrix output{matrix};

    if (r1 >= output.size() || r2 >= output.size())
    {

        throw std::logic_error("Caution: r1 or r2 inputs are out of range");
    }

    for (size_t j{}; j < output[0].size(); j++)
    {

        op.push_back(output[r1][j]);
    }
    output[r1].clear();
    for (size_t j{}; j < output[0].size(); j++)
    {

        output[r1].push_back(output[r2][j]);
    }
    output[r2].clear();
    for (size_t j{}; j < output[0].size(); j++)
    {

        output[r2].push_back(op[j]);
    }

    return output;
}

Matrix algebra::ero_multiply(const Matrix &matrix, size_t r, double c)
{
    std::vector<double> op{};
    Matrix output{};

    if (r >= matrix.size())
    {

        throw std::logic_error("Caution: r1 or r2 inputs are out of range");
    }

    for (size_t i{}; i < matrix.size(); i++)
    {
        if (i == r)
        {
            for (size_t j{}; j < matrix[0].size(); j++)
            {

                op.push_back(matrix[i][j] * c);
            }

            output.push_back(op);
            op.clear();
        }

        else
        {

            for (size_t j{}; j < matrix[0].size(); j++)
            {

                op.push_back(matrix[i][j]);
            }

            output.push_back(op);
            op.clear();
        }
    }

    return output;
}

Matrix algebra::ero_sum(const Matrix &matrix, size_t r1, double c, size_t r2)
{
    Matrix output{matrix};
    std::vector<double> op;

    if (r1 >= output.size() || r2 >= output.size())
    {

        throw std::logic_error("Caution: r1 or r2 inputs are out of range");
    }

    for (size_t j{0}; j < output[0].size(); j++)
    {

        op.push_back(c * output[r1][j]);
    }

    for (size_t j{0}; j < output[0].size(); j++)
    {

        output[r2][j] += op[j];
    }

    return output;
}

Matrix algebra::upper_triangular(const Matrix &matrix)
{
    Matrix output{matrix};
    if (output.empty())
    {
        return output;
    }
    if (output.size() != output[0].size())
    {

        throw std::logic_error("Caution: non-square matrices have no upper triangular form");
    }

    for (size_t i{}; i < output.size(); i++)
    {

        for (size_t j{i}; j < output.size() - 1; j++)
        {

            output = ero_sum(output, i, -1 * output[j + 1][i] / output[i][i], j + 1);
        }
    }

    return output;
}