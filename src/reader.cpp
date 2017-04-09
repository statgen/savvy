
#include "reader.hpp"

namespace savvy
{
  reader::reader(const std::string& file_path)
  {
    if (savvy::detail::has_extension(file_path, ".sav"))
      sav_reader_ = std::make_unique<sav::reader>(file_path);
    else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, ".vcf.gz") || detail::has_extension(file_path, ".bcf"))
      vcf_reader_ = std::make_unique<vcf::reader>(file_path);
  }

  indexed_reader::indexed_reader(const std::string& file_path, const region& reg)
  {
    if (savvy::detail::has_extension(file_path, ".sav"))
      sav_reader_ = std::make_unique<sav::indexed_reader>(file_path, reg);
    else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, ".vcf.gz") || detail::has_extension(file_path, ".bcf"))
      vcf_reader_ = std::make_unique<vcf::indexed_reader>(file_path, reg);
  }

  void indexed_reader::reset_region(const region& reg)
  {
    if (sav_reader_)
      dynamic_cast<sav::indexed_reader*>(sav_reader_.get())->reset_region(reg);
    else if (vcf_reader_)
      dynamic_cast<vcf::indexed_reader*>(vcf_reader_.get())->reset_region(reg);
  }
}