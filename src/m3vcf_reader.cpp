#include "m3vcf_reader.hpp"

#include <cmath>

std::uint64_t ceil_divide(std::uint64_t dividend, std::uint64_t divisor)
{
  return (std::uint64_t)(1) + ((dividend - (std::uint64_t)(1)) / divisor);
}

namespace vc
{

  m3vcf_reader::m3vcf_reader(const std::string& file_path)
    :
    file_path_(file_path)
  {

  }
}