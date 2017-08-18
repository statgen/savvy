#include "savvy/vcf_reader.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/savvy.hpp"

#include <fstream>

int main(int argc, char** argv)
{
  savvy::sav::reader input(argv[1]);
  savvy::basic_variant_iterator<savvy::sav::reader, savvy::dense_allele_vector<float>> cur(input);
  savvy::basic_variant_iterator<savvy::sav::reader, savvy::dense_allele_vector<float>> eof{};

  std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
  std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());

  savvy::vcf::writer vcf_output(argc > 2 ? std::string(argv[2]) : std::string("/dev/stdout"), sample_ids.begin(), sample_ids.end());

  while (cur != eof)
  {
    vcf_output << *cur;

    ++cur;
  }

  return 0;
}