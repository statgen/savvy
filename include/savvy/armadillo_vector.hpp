#ifndef LIBSAVVY_ARMADILLO_VECTOR_HPP
#define LIBSAVVY_ARMADILLO_VECTOR_HPP

#include "allele_vector.hpp"

#include <armadillo>
#include <cstddef>

namespace savvy
{
  namespace armadillo
  {
    template <typename T>
    class sparse_vector : public arma::SpCol<T>
    {
    public:
      typedef T value_type;
      using arma::SpCol<T>::SpCol;

      void resize(std::size_t sz)
      {
        arma::SpCol<T>::resize(sz, 1);
      }
    };

    template <typename T>
    class dense_vector : public arma::Col<T>
    {
    public:
      typedef T value_type;
      using arma::Col<T>::Col;

      void resize(std::size_t sz)
      {
        std::size_t before_size = arma::Col<T>::size();
        arma::Col<T>::resize(sz);
        if (arma::Col<T>::size() > before_size)
        {
          for (std::size_t i = before_size; i < arma::Col<T>::size(); ++i)
            (*this)[i] = value_type();
        }
      }
    };

    template <typename T>
    using sparse_allele_vector = allele_vector<sparse_vector<T>>;

    template <typename T>
    using dense_allele_vector = allele_vector<dense_vector<T>>;
  }
}

#endif //LIBSAVVY_ARMADILLO_VECTOR_HPP