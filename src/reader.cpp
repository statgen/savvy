
#include "reader.hpp"

namespace vc
{
  reader::reader(const std::string& file_path)
  {
    if (vc::detail::has_extension(file_path, ".cmf"))
      cmf_reader_ = std::make_unique<cmf::reader>(file_path);
    else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, ".vcf.gz") || detail::has_extension(file_path, ".bcf"))
      vcf_reader_ = std::make_unique<vcf::reader>(file_path);
  }

  region_reader::region_reader(const std::string& file_path, region reg)
  {
    if (vc::detail::has_extension(file_path, ".cmf"))
      cmf_reader_ = std::make_unique<cmf::region_reader>(file_path, reg);
    else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, ".vcf.gz") || detail::has_extension(file_path, ".bcf"))
      vcf_reader_ = std::make_unique<vcf::region_reader>(file_path, reg);
  }
}