
#ifndef LIBSAVVY_ALLELE_VECTOR_HPP
#define LIBSAVVY_ALLELE_VECTOR_HPP

#include "compressed_vector.hpp"

#include <string>
#include <vector>
#include <unordered_map>

namespace savvy
{
  struct variant_details
  {
    variant_details(std::string&& chromosome, std::uint64_t locus, std::string&& ref, std::string&& alt, std::unordered_map<std::string, std::string>&& properties) :
      chromosome(std::move(chromosome)),
      locus(locus),
      ref(std::move(ref)),
      alt(std::move(alt)),
      properties(std::move(properties))
    {

    }

    std::string chromosome;
    std::uint64_t locus;
    std::string ref;
    std::string alt;
    std::unordered_map<std::string, std::string> properties;
  };

  template<typename T>
  class variant_vector : public T
  {
  public:
    variant_vector()
    {
    }

    variant_vector(
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

    virtual ~variant_vector() {}

    void operator<<(variant_details&& details)
    {
      chromosome_ = std::move(details.chromosome);
      locus_ = details.locus;
      ref_ = std::move(details.ref);
      alt_ = std::move(details.alt);
      properties_ = std::move(details.properties);
    }

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

  template <typename T>
  class allele_vector : public variant_vector<T>
  {
  public:
    using variant_vector<T>::variant_vector;
  };

  template <typename T>
  class genotype_vector : public variant_vector<T>
  {
  public:
    using variant_vector<T>::variant_vector;
  };

  template<typename T>
  const std::string variant_vector<T>::empty_string = {};

  template <typename T>
  using dense_allele_vector = allele_vector<std::vector<T>>;
  template <typename T>
  using sparse_allele_vector = allele_vector<compressed_vector<T>>;

  template <typename T>
  using dense_genotype_vector = genotype_vector<std::vector<T>>;
  template <typename T>
  using sparse_genotype_vector = genotype_vector<compressed_vector<T>>;
}
#endif //LIBSAVVY_ALLELE_VECTOR_HPP
