
#include "vcf_reader.hpp"
#include "cvcf_reader.hpp"

#include <fstream>

int main(int argc, char** argv)
{

  std::ofstream ofs(argv[2], std::ios::binary);

  vc::vcf::block buff;
  vc::vcf::reader input(argv[1]);
  vc::vcf::reader::input_iterator eof;
  vc::vcf::reader::input_iterator cur(input, buff);

  std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
  std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
  vc::cvcf::writer compact_output(ofs, "20", 2, sample_ids.begin(), sample_ids.end());

  while (cur != eof)
  {
    compact_output << vc::cvcf::marker(cur->pos(), cur->ref(), cur->alt(), cur->begin(), cur->end());

    ++cur;
  }

  return 0;
}