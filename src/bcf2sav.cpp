
#include "savvy/vcf_reader.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/savvy.hpp"

#include <fstream>

int main(int argc, char** argv)
{
  savvy::vcf::reader input(argv[1]);
  savvy::basic_variant_iterator<savvy::vcf::reader, savvy::dense_allele_vector<float>> cur(input);
  savvy::basic_variant_iterator<savvy::vcf::reader, savvy::dense_allele_vector<float>> eof{};

  std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
  std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
  auto variant_metadata = input.prop_fields();
  auto file_metadata = input.metadata();
  savvy::sav::writer compact_output(argv[2], sample_ids.begin(), sample_ids.end(), file_metadata.begin(), file_metadata.end(), variant_metadata.begin(), variant_metadata.end());

  while (cur != eof)
  {
    compact_output << *cur;

    ++cur;
  }

  return 0;
}