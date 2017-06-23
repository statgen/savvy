
#include <savvy/reader.hpp>
#include "savvy/vcf_reader.hpp"
#include "savvy/m3vcf_reader.hpp"

int main(int argc, char** argv)
{

  std::ofstream ofs(argv[2], std::ios::binary);


  savvy::vcf::reader input(argv[1]);
  savvy::vcf::dense_allele_variant_iterator<float> eof;
  savvy::vcf::dense_allele_variant_iterator<float> cur(input);

  if (cur != eof)
  {
    std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
    std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
    savvy::m3vcf::writer compact_output(ofs, cur->chromosome(), savvy::get_ploidy(input, *cur), sample_ids.begin(), sample_ids.end());

    savvy::m3vcf::block output_block(sample_ids.size(), 2);
    while (cur != eof)
    {
      if (!output_block.add_marker(cur->locus(), cur->ref(), cur->alt(), cur->begin(), cur->end()))
      {
        compact_output << output_block;
        output_block = savvy::m3vcf::block(sample_ids.size(), 2);
        output_block.add_marker(cur->locus(), cur->ref(), cur->alt(), cur->begin(), cur->end());
      }

      ++cur;
    }
    if (output_block.marker_count())
      compact_output << output_block;
  }

  return 0;
}