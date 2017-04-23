#ifndef LIBSAVVY_UBLAS_VECTOR_HPP
#define LIBSAVVY_UBLAS_VECTOR_HPP

#include "allele_vector.hpp"

#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace savvy
{
  template <typename T>
  using sparse_ublas_allele_vector = allele_vector<boost::numeric::ublas::compressed_vector<T>>;

  template <typename T>
  using dense_ublas_allele_vector = allele_vector<boost::numeric::ublas::vector<T>>;
}

#endif //LIBSAVVY_UBLAS_VECTOR_HPP