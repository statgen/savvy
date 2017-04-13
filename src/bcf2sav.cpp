
#include "savvy/vcf_reader.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/savvy.hpp"

#include <fstream>

int main(int argc, char** argv)
{
  savvy::vcf::reader input(argv[1]);
  savvy::vcf::dense_variant_iterator<float> cur(input);
  savvy::vcf::dense_variant_iterator<float> eof{};
  const std::string chrom = cur != eof ? cur->chromosome() : "";
  const int ploidy =  cur != eof ? savvy::get_ploidy(input, *cur) : 0;

  std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
  std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
  savvy::sav::writer compact_output(argv[2], chrom, ploidy, sample_ids.begin(), sample_ids.end());

  while (cur != eof)
  {
    if (cur->chromosome() != chrom)
    {
      std::cerr << "Multiple chromosomes encountered. CMF files can only contain one chromosome." << std::endl;
      return -1;
    }
    compact_output << *cur;

    ++cur;
  }

  return 0;
}