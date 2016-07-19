#include <iostream>
#include <algorithm>
#include "cvcf_reader.hpp"

int main(int argc, char** argv)
{

  vc::cvcf::marker m;



  std::uint64_t ploidy_level = 2;
  std::uint64_t sample_size = 1000;
  std::vector<int> zero_one_two_vec(sample_size, 0);

  std::for_each(m.non_ref_begin(), m.non_ref_end(), [&zero_one_two_vec, ploidy_level](const vc::cvcf::marker::sparse_allele& a)
  {
    if (a.status == vc::allele_status::has_alt)
      ++(zero_one_two_vec[a.offset / ploidy_level]);
  });

  std::size_t i = 0;
  std::for_each(m.begin(), m.end(), [&zero_one_two_vec, &i, ploidy_level](const vc::allele_status& s)
  {
    if (s == vc::allele_status::has_alt)
      ++(zero_one_two_vec[i / ploidy_level]);
  });

  i = 0;
  for (auto const &s : m)
  {
    if (s == vc::allele_status::has_alt)
      ++(zero_one_two_vec[i / ploidy_level]);
  }

  return 0;
}


