/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <savvy/reader.hpp>
#include "savvy/vcf_reader.hpp"
#include "savvy/m3vcf_reader.hpp"

int main(int argc, char** argv)
{

  std::ofstream ofs(argv[2], std::ios::binary);


  savvy::vcf::reader<1> input(argv[1], savvy::fmt::allele);
  savvy::site_info variant;
  std::vector<float> genotypes;

  if (input.read(variant, genotypes))
  {
    std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
    std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
    savvy::m3vcf::writer compact_output(ofs, variant.chromosome(), savvy::get_ploidy(input, genotypes), sample_ids.begin(), sample_ids.end());

    savvy::m3vcf::block output_block(sample_ids.size(), 2);
    while (input.good())
    {
      if (!output_block.add_marker(variant.position(), variant.ref(), variant.alt(), genotypes.begin(), genotypes.end()))
      {
        compact_output << output_block;
        output_block = savvy::m3vcf::block(sample_ids.size(), 2);
        output_block.add_marker(variant.position(), variant.ref(), variant.alt(), genotypes.begin(), genotypes.end());
      }

      input.read(variant, genotypes);
    }
    if (output_block.marker_count())
      compact_output << output_block;
  }

  return 0;
}