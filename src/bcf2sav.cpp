
#include "savvy/vcf_reader.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/savvy.hpp"

#include <fstream>

int main(int argc, char** argv)
{
  savvy::vcf::reader<1> input(argv[1], savvy::fmt::allele);
  savvy::site_info variant;
  std::vector<float> genotypes;

  std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
  std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
  auto headers = input.headers();
  headers.reserve(headers.size() + 3);
  headers.insert(headers.begin(), {"INFO","<ID=FILTER,Description=\"Variant filter\">"});
  headers.insert(headers.begin(), {"INFO","<ID=QUAL,Description=\"Variant quality\">"});
  headers.insert(headers.begin(), {"INFO","<ID=ID,Description=\"Variant ID\">"});
  savvy::sav::writer compact_output(argv[2], sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end());

  while (input.read(variant, genotypes))
    compact_output.write(variant, genotypes);

  return 0;
}