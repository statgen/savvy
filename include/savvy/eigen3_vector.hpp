
#ifndef LIBSAVVY_EIGEN3_VECTOR_HPP
#define LIBSAVVY_EIGEN3_VECTOR_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace savvy
{
  template <typename T>
  class eigen3_sparse_vector : public Eigen::SparseVector<T>
  {
  public:
    typedef T value_type;
    using Eigen::SparseVector<T>::SparseVector;
    T& operator[](std::size_t idx)
    {
      return Eigen::SparseVector<T>::coeffRef(idx);
    }

    T operator[](std::size_t idx) const
    {
      return Eigen::SparseVector<T>::coeff(idx);
    }
  };


  template <typename T>
  class eigen3_dense_vector : public Eigen::Matrix<T, 1, Eigen::Dynamic>
  {
  public:
    typedef T value_type;
    eigen3_dense_vector()
    {
    }

    eigen3_dense_vector(std::size_t size)
    {
      eigen3_dense_vector::resize(size);
    }

    T& operator[](std::size_t idx)
    {
      return Eigen::Matrix<T, 1, Eigen::Dynamic>::coeffRef(idx);
    }

    const T& operator[](std::size_t idx) const
    {
      return Eigen::Matrix<T, 1, Eigen::Dynamic>::coeffRef(idx);
    }

    void resize(std::size_t size)
    {
      Eigen::Matrix<T, 1, Eigen::Dynamic>::setConstant(size, 0);
    }
  };

  template <typename T>
  using sparse_eigen3_allele_vector = allele_vector<eigen3_sparse_vector<T>>;

  template <typename T>
  using dense_eigen3_allele_vector = allele_vector<dense_sparse_vector<T>>;
}
#endif //LIBSAVVY_EIGEN3_VECTOR_HPP
