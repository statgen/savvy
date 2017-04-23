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
    class dense_vector : public arma::Mat<T>
    {
    public:
      typedef T value_type;
      using arma::Mat<T>::Mat;

      void resize(std::size_t sz)
      {
        std::size_t before_size = arma::Mat<T>::size();
        arma::Mat<T>::resize(sz);
        if (arma::Mat<T>::size() > before_size)
        {
          for (std::size_t i = before_size; i < arma::Mat<T>::size(); ++i)
            (*this)[i] = value_type();
        }
      }
    };

    //template <typename T>
    //using sparse_allele_vector = allele_vector<arma::sp_mat<T>>;

    template <typename T>
    using dense_allele_vector = allele_vector<dense_vector<T>>;
  }
}

#endif //LIBSAVVY_ARMADILLO_VECTOR_HPP