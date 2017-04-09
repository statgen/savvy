#ifndef LIBVC_ALLELE_STATUS_HPP
#define LIBVC_ALLELE_STATUS_HPP

#include <cstdint>

namespace savvy
{
  enum class allele_status : std::int8_t { is_missing = -1, has_ref, has_alt };
}

#endif //LIBVC_ALLELE_STATUS_HPP
