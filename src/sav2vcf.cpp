#include <cmath>
#include "savvy/vcf_reader.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/savvy.hpp"
#include "savvy/utility.hpp"

#include <fstream>

int main(int argc, char** argv)
{
  savvy::sav::reader<1> input(argv[1], savvy::fmt::allele);
  savvy::site_info variant;
  std::vector<float> genotypes;

  std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
  std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
  auto variant_metadata = input.prop_fields();

  auto headers = input.headers();

  std::string output_path = argc > 2 ? argv[2] : "/dev/stdout";
  if (output_path == "-")
    output_path = "/dev/stdout";

  for (auto it = headers.begin(); it != headers.end(); )
  {
    std::string header_id = savvy::parse_header_id(it->second);
    if (it->first == "INFO" && (header_id == "ID" || header_id == "QUAL" || header_id == "FILTER"))
      it = headers.erase(it);
    else
    {
      if (it->first == "FORMAT")
      {
        if (header_id == "GT")
          it->second = "<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
        else if (header_id == "GP")
          it->second = "<ID=GP,Number=3,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1\">";
      }
      else if (it->first == "fileDate")
      {
        // TODO: update date
      }

      ++it;
    }
  }

  savvy::vcf::writer vcf_output(output_path, sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end());

  while (input.read(variant, genotypes))
    vcf_output.write(variant, genotypes);

  return 0;
}