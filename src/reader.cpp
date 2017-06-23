
#include "savvy/reader.hpp"

namespace savvy
{
  std::vector<std::string> reader_base::prop_fields() const
  {
    if (sav_reader_)
      return sav_reader_->prop_fields();
    else if (vcf_reader_)
      return vcf_reader_->prop_fields();
    return {};
  }

  reader_base::sample_iterator reader_base::samples_begin() const
  {
    reader_base::sample_iterator ret;
    if (sav_reader_)
      ret = reader_base::sample_iterator(sav_reader_->samples_begin());
    else if (vcf_reader_)
      ret = reader_base::sample_iterator(vcf_reader_->samples_begin());
    return ret;
  }

  reader_base::sample_iterator reader_base::samples_end() const
  {
    reader_base::sample_iterator ret;
    if (sav_reader_)
      ret = reader_base::sample_iterator(sav_reader_->samples_end());
    else if (vcf_reader_)
      ret = reader_base::sample_iterator(vcf_reader_->samples_end());
    return ret;
  }

  std::size_t reader_base::sample_size() const
  {
    std::size_t ret{};
    if (sav_reader_)
      ret = static_cast<std::size_t>(sav_reader_->samples_end() - sav_reader_->samples_begin());
    else if (vcf_reader_)
      ret = static_cast<std::size_t>(vcf_reader_->samples_end() - vcf_reader_->samples_begin());
    return ret;
  }

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