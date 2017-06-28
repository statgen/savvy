
#include "savvy/vcf_reader.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/savvy.hpp"

#include <fstream>

int main(int argc, char** argv)
{
  savvy::vcf::reader input(argv[1]);
  savvy::vcf::dense_allele_variant_iterator<float> cur(input);
  savvy::vcf::dense_allele_variant_iterator<float> eof{};

  std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
  std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
  savvy::sav::writer compact_output(argv[2], sample_ids.begin(), sample_ids.end());

  while (cur != eof)
  {
    compact_output << *cur;

    ++cur;
  }

  return 0;
}