#ifndef LIBSAVVY_ALLELE_STATUS_HPP
#define LIBSAVVY_ALLELE_STATUS_HPP

#include <cstdint>

namespace savvy
{
  enum class allele_status : std::int8_t { is_missing = -1, has_ref, has_alt };
}

#endif //LIBSAVVY_ALLELE_STATUS_HPP
