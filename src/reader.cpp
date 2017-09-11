
#include "savvy/reader.hpp"

namespace savvy
{
  std::vector<std::string> reader_base::prop_fields() const
  {
    if (sav_impl())
      return sav_impl()->prop_fields();
    else if (vcf_impl())
      return vcf_impl()->prop_fields();
    return {};
  }

  reader_base::sample_iterator reader_base::samples_begin() const
  {
    reader_base::sample_iterator ret;
    if (sav_impl())
      ret = reader_base::sample_iterator(sav_impl()->samples_begin());
    else if (vcf_impl())
      ret = reader_base::sample_iterator(vcf_impl()->samples_begin());
    return ret;
  }

  reader_base::sample_iterator reader_base::samples_end() const
  {
    reader_base::sample_iterator ret;
    if (sav_impl())
      ret = reader_base::sample_iterator(sav_impl()->samples_end());
    else if (vcf_impl())
      ret = reader_base::sample_iterator(vcf_impl()->samples_end());
    return ret;
  }

  std::size_t reader_base::sample_size() const
  {
    std::size_t ret{};
    if (sav_impl())
      ret = static_cast<std::size_t>(sav_impl()->samples_end() - sav_impl()->samples_begin());
    else if (vcf_impl())
      ret = static_cast<std::size_t>(vcf_impl()->samples_end() - vcf_impl()->samples_begin());
    return ret;
  }

  reader::reader(const std::string& file_path, fmt data_format)
  {
    if (savvy::detail::has_extension(file_path, ".sav"))
      sav_reader_ = detail::make_unique<sav::reader>(file_path, data_format);
    else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, ".vcf.gz") || detail::has_extension(file_path, ".bcf"))
      vcf_reader_ = detail::make_unique<vcf::reader>(file_path, data_format);
  }

  std::vector<std::string> indexed_reader::chromosomes() const
  {
    if (sav_reader_)
      return sav_reader_->chromosomes();
    else if (vcf_reader_)
      return vcf_reader_->chromosomes();
    return {};
  }

  indexed_reader::indexed_reader(const std::string& file_path, const region& reg, fmt data_format)
  {
    if (savvy::detail::has_extension(file_path, ".sav"))
      sav_reader_ = detail::make_unique<sav::indexed_reader>(file_path, reg, data_format);
    else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, ".vcf.gz") || detail::has_extension(file_path, ".bcf"))
      vcf_reader_ = detail::make_unique<vcf::indexed_reader>(file_path, reg, data_format);
  }

  void indexed_reader::reset_region(const region& reg)
  {
    if (sav_reader_)
      sav_reader_->reset_region(reg);
    else if (vcf_reader_)
      vcf_reader_->reset_region(reg);
  }
}