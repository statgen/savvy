
#include "vcf_reader.hpp"
#include "cmf_reader.hpp"

#include <fstream>

int main(int argc, char** argv)
{


  vc::vcf::block buff;
  vc::vcf::reader input(argv[1]);
  vc::vcf::reader::input_iterator eof;
  vc::vcf::reader::input_iterator cur(input, buff);
  const std::string chrom = cur != eof ? vc::vcf::reader::get_chromosome(input, *cur) : "";
  const int ploidy =  cur != eof ? cur->ploidy() : 0;

  std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
  std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
  vc::cmf::writer compact_output(argv[2], vc::vcf::reader::get_chromosome(input, *cur), ploidy, sample_ids.begin(), sample_ids.end());

  while (cur != eof)
  {
    if (vc::vcf::reader::get_chromosome(input, *cur) != chrom)
    {
      std::cerr << "Multiple chromosomes encountered. CMF files can only contain one chromosome." << std::endl;
      return -1;
    }
    compact_output << vc::cmf::marker(cur->pos(), cur->ref(), cur->alt(), cur->begin(), cur->end());

    ++cur;
  }

  return 0;
}