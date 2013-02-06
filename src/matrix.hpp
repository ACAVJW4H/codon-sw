#ifndef CODONALIGN_MATRIX_H
#define CODONALIGN_MATRIX_H

#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace codonalign
{

template<typename T>
class Matrix
{
public:
    typedef typename std::vector<T>::const_iterator const_iterator;
    Matrix();
    Matrix(const size_t nrows, const size_t ncols);
    inline const T& operator()(const size_t row, const size_t col) const;
    inline T& operator()(const size_t row, const size_t col);
    inline T& operator()(const std::pair<size_t,size_t>& p) { return operator()(p.first, p.second); };
    inline const T& operator()(const std::pair<size_t,size_t>& p) const { return operator()(p.first, p.second); };
    inline T get_or_default(const int row, const int col, T def);
    void fill(const T& value);
    inline size_t nrows() const { return nrows_; };
    inline size_t ncols() const { return ncols_; };
    inline size_t size() const { return size_; };
    inline std::pair<size_t, size_t> rev_index(const size_t s) const;
    inline size_t index(const size_t row, const size_t col) const;
    inline const_iterator begin() const { return arr.begin(); };
    inline const_iterator end() const { return arr.end(); };
private:
    size_t nrows_, ncols_, size_;
    std::vector<T> arr;
};

template<typename T>
Matrix<T>::Matrix() :
    nrows_(0),
    ncols_(0),
    size_(0),
    arr(size_)
{
}

template<typename T>
Matrix<T>::Matrix(const size_t nrows, const size_t ncols) :
    nrows_(nrows),
    ncols_(ncols),
    size_(nrows*ncols),
    arr(size_)
{
}

template<typename T>
inline const T& Matrix<T>::operator()(const size_t row, const size_t col) const
{
    return arr[index(row, col)];
}

template<typename T>
inline T& Matrix<T>::operator()(const size_t row, const size_t col)
{
    return arr[index(row, col)];
}

template<typename T>
inline size_t Matrix<T>::index(const size_t row, const size_t col) const
{
    if(row > nrows())
        throw std::runtime_error("row index out of bounds: " + std::to_string(row));
    else if(col > ncols())
        throw std::runtime_error("column index out of bounds: " + std::to_string(col));
    const size_t idx = ncols() * row + col;
    assert(idx < size());
    return idx;
}

template<typename T>
inline std::pair<size_t, size_t>  Matrix<T>::rev_index(const size_t s) const
{
    if(s > size())
        throw std::runtime_error("index out of bounds: " + std::to_string(s) + " > " + std::to_string(size()));
    return { s / ncols(), s % ncols() };
}

template<typename T>
void Matrix<T>::fill(const T& value)
{
    std::fill(std::begin(arr), std::end(arr), value);
}

/// Get the value at <c>[row,col]</c> or \c def if either index is out of bounds.
template<typename T>
inline T Matrix<T>::get_or_default(const int row, const int col, T def)
{
    if(row < 0 || col < 0 || col > ncols() || row > nrows()) {
        return def;
    }
    return arr[index(row, col)];
}

} // namespace codonalign

/// Write matrix to stream
template<typename T>
::std::ostream& operator<<(std::ostream& output, const codonalign::Matrix<T>& m)
{
    for(size_t i = 0; i < m.nrows(); i++) {
        for(size_t j = 0; j < m.ncols(); j++) {
            output << m(i, j) << '\t';
        }
        output << std::endl;
    }
    return output;  // for multiple << operators.
}

#endif
