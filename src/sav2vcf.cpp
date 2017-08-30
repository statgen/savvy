#include "savvy/vcf_reader.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/savvy.hpp"
#include "savvy/utility.hpp"

#include <fstream>

int main(int argc, char** argv)
{
  savvy::sav::reader input(argv[1]);
  savvy::basic_variant_iterator<savvy::sav::reader, savvy::dense_allele_vector<float>> cur(input);
  savvy::basic_variant_iterator<savvy::sav::reader, savvy::dense_allele_vector<float>> eof{};

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
      ++it;
  }

  savvy::vcf::writer vcf_output(output_path, sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end());

  while (cur != eof)
  {
    vcf_output << *cur;

    ++cur;
  }

  return 0;
}