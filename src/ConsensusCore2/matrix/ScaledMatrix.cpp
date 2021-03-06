
#include <algorithm>
#include <cmath>
#include <iomanip>

#include "ScaledMatrix.h"

namespace PacBio {
namespace Consensus {

ScaledMatrix::ScaledMatrix(int rows, int cols, Direction dir)
    : SparseMatrix(rows, cols), logScalars_(cols, 0.0), dir_{dir}
{
}
ScaledMatrix::ScaledMatrix(const ScaledMatrix& other)
    : SparseMatrix(other), logScalars_(other.logScalars_)
{
}

void ScaledMatrix::Reset(size_t rows, size_t cols)
{
    std::vector<double>(cols, 0.0).swap(logScalars_);
    SparseMatrix::Reset(rows, cols);
}

ScaledMatrix::Direction ScaledMatrix::SetDirection(const Direction dir)
{
    const Direction res = dir_;
    dir_ = dir;
    std::fill(logScalars_.begin(), logScalars_.end(), 0.0);
    return res;
}

std::ostream& operator<<(std::ostream& os, const ScaledMatrix& mat)
{
    os << "MATRIX (" << mat.Rows() << ", " << mat.Columns() << ") BEGIN" << std::endl;
    os << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < mat.Rows(); ++i) {
        os << ' ';
        for (size_t j = 0; j < mat.Columns(); ++j) {
            os << ' ' << std::setw(9) << std::log(mat.Get(i, j)) + mat.GetLogScale(j);
        }
        os << std::endl;
    }
    os << "END" << std::endl;
    return os;
}

}  // namespace Consensus
}  // namespace PacBio
