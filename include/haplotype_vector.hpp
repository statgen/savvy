
#ifndef LIBVC_HAPLOTYPE_VECTOR_HPP
#define LIBVC_HAPLOTYPE_VECTOR_HPP

#include "compressed_vector.hpp"

#include <string>
#include <vector>
#include <unordered_map>

namespace vc
{
  template<typename T>
  class haplotype_vector : public T
  {
  public:
    haplotype_vector()
    {
    }

    haplotype_vector(
      std::string&& chromosome,
      std::uint64_t locus,
      std::string&& ref,
      std::string&& alt,
      std::unordered_map<std::string, std::string>&& properties,
      T&& mixin_vector)
      :
      T(std::move(mixin_vector)),
      chromosome_(std::move(chromosome)),
      locus_(locus),
      ref_(std::move(ref)),
      alt_(std::move(alt)),
      properties_(std::move(properties))
    {

    }

//    template <typename RandAccessIterType>
//    haplotype_vector(
//      const std::string& chromosome,
//      std::uint64_t locus,
//      const std::string& ref,
//      const std::string& alt,
//      std::uint64_t sample_count,
//      std::uint8_t ploidy,
//      RandAccessIterType gt_beg,
//      RandAccessIterType gt_end)
//      :
//      T(std::move(mixin_vector)),
//      chromosome_(chromosome),
//      locus_(locus),
//      ref_(ref),
//      alt_(alt),
//      sample_count_(sample_count),
//      ploidy_(ploidy)
//    {
//      T::resize(sample_count_ * ploidy_, 0.0);
//      for (it = gt_beg; it != gt_end; ++it)
//      {
//        if (*it != 0.0)
//        {
//          (*this)[std::distance(gt_beg, it)] = *it;
//        }
//      }
//    }

    const std::string& chromosome() const { return chromosome_; }
    const std::string& ref() const { return ref_; }
    const std::string& alt() const { return alt_; }
    std::uint64_t locus() const { return locus_; }
    const std::string& prop(const std::string& key) const
    {
      auto it = properties_.find(key);
      if (it == properties_.end())
        return empty_string;
      return it->second;
    }
  private:
    std::string chromosome_;
    std::string ref_;
    std::string alt_;
    std::unordered_map<std::string, std::string> properties_;
    std::uint64_t locus_;
    static const std::string empty_string;
  };

  template<typename T>
  const std::string haplotype_vector<T>::empty_string = {};

  template <typename T>
  using dense_haplotype_vector = haplotype_vector<std::vector<T>>;
  template <typename T>
  using sparse_haplotype_vector = haplotype_vector<compressed_vector<T>>;
}
#endif //LIBVC_HAPLOTYPE_VECTOR_HPP
