#ifndef SAVVY_VARIANT_GROUP_ITERATOR_HPP
#define SAVVY_VARIANT_GROUP_ITERATOR_HPP

#include "reader.hpp"

#include <regex>
#include <tuple>
#include <iterator>
#include <cstddef>


namespace savvy
{
  template<typename RdrT, typename VecT>
  class basic_variant_group_iterator
  {
  public:
    typedef basic_variant_group_iterator<RdrT, VecT> self_type;
    typedef std::ptrdiff_t difference_type;
    typedef savvy::variant<VecT> value_type;
    typedef const value_type& reference;
    typedef const value_type *pointer;
    typedef std::input_iterator_tag iterator_category;

    basic_variant_group_iterator(savvy::indexed_reader& rdr, std::string marker_group_file_line) :
      rdr_(&rdr)
    {
      std::tie(group_id_, sites_) = parse_marker_group_line(marker_group_file_line);
      init();
    }

    basic_variant_group_iterator(savvy::indexed_reader& rdr, std::string group_id, std::list<site_info> sites) :
      rdr_(&rdr),
      group_id_(std::move(group_id)),
      sites_(std::move(sites))
    {
      init();
    }

    basic_variant_group_iterator() :
      rdr_(nullptr)
    {
    }

    self_type& operator++()
    {
      increment();
      return *this;
    }

    void operator++(int)
    { increment(); }

    reference operator*()
    { return variant_; }

    pointer operator->()
    { return &variant_; }

    bool operator==(const self_type& rhs)
    { return (rdr_ == rhs.rdr_); }

    bool operator!=(const self_type& rhs)
    { return (rdr_ != rhs.rdr_); }

    const std::string& group_id() const
    { return group_id_; }

    const std::list<site_info>& sites() const
    { return sites_; }

  private:
    void init()
    {
      std::vector<savvy::region> un_merged_regions(sites_.size(), {""});
      auto out_it = un_merged_regions.begin();
      for (auto in_it = sites_.begin(); in_it != sites_.end(); ++in_it, ++out_it)
        *out_it = site_info_to_region(*in_it);

      merged_regions_ = merge_regions(un_merged_regions);
      site_it_ = sites_.begin();
      region_it_ = merged_regions_.begin();

      rdr_->reset_region(*region_it_);

      increment();
    }

    void increment()
    {
      while (region_it_ != merged_regions_.end())
      {
        while (rdr_->read(variant_, variant_.data()) && site_it_ != sites_.end())
        {
          while (site_it_ != sites_.end())
          {
            if (site_it_->position() >= variant_.position() || site_it_->chromosome() != variant_.chromosome())
              break;
            ++site_it_;
          }

          if (site_it_->chromosome() != variant_.chromosome())
            break;

          for (auto pos_it = site_it_; pos_it != sites_.end() && site_it_->position() == pos_it->position(); ++pos_it)
          {
            std::string target_id = pos_it->chromosome() + ":" + std::to_string(pos_it->position()) + "_" + pos_it->ref() + "/" + pos_it->alt();
            std::string current_id = variant_.chromosome() + ":" + std::to_string(variant_.position()) + "_" + variant_.ref() + "/" + variant_.alt();
            if (
              pos_it->chromosome() == variant_.chromosome() &&
              pos_it->position() == variant_.position() &&
              pos_it->ref() == variant_.ref() &&
              pos_it->alt() == variant_.alt())
            {
              return;
            }
          }
        }

        ++region_it_;
        if (region_it_ != merged_regions_.end())
          rdr_->reset_region(*region_it_);
      }

      rdr_ = nullptr;
    }

    std::vector<region> merge_regions(const std::vector<region>& un_merged_regions)
    {
      std::vector<region> ret;

      for (auto it = un_merged_regions.begin(); it != un_merged_regions.end(); ++it)
      {
        if (ret.empty() || ret.back().chromosome() != it->chromosome())
        {
          ret.emplace_back(*it);
        }
        else
        {
          std::uint64_t from = std::min(ret.back().from(), it->from());
          std::uint64_t to = std::max(ret.back().to(), it->to());
          ret.back() = region(ret.back().chromosome(), from, to);
        }
      }

      return ret;
    }

    // [CHROM]:[POS]_[REF]/[ALT]
    static savvy::site_info marker_id_to_site_info(std::string::const_iterator beg, std::string::const_iterator end)
    {
      auto colon_it = std::find(beg, end, ':');
      std::string chrom(beg, colon_it);
      if (colon_it != end)
      {
        auto underscore_it = std::find(++colon_it, end, '_');
        std::uint64_t pos = static_cast<std::uint64_t>(std::atoll(std::string(colon_it, underscore_it).c_str()));
        if (underscore_it != end)
        {
          auto slash_it = std::find(++underscore_it, end, '/');
          std::string ref(underscore_it, slash_it);
          if (slash_it != end)
          {
            std::string alt(++slash_it, end);
            return savvy::site_info{std::move(chrom), pos, std::move(ref), std::move(alt), {}};
          }
        }
      }

      return savvy::site_info{};
    }

    // [CHROM]:[POS]_[REF]/[ALT]
    static savvy::region site_info_to_region(const savvy::site_info& site)
    {
      std::size_t length = std::max(site.ref().size(), site.alt().size());
      if (length > 0)
      {
        return savvy::region{site.chromosome(), site.position(), site.position() + length - 1};
      }

      return savvy::region{""};
    }

    static std::tuple<std::string, std::list<savvy::site_info>> parse_marker_group_line(const std::string& input)
    {
      std::tuple<std::string, std::list<savvy::site_info>> ret;
      auto delim_it = std::find(input.begin(), input.end(), '\t');
      if (delim_it != input.end())
      {
        std::get<0>(ret) = std::string(input.begin(), delim_it);
        ++delim_it;

        std::string::const_iterator next_delim_it;
        while ((next_delim_it = std::find(delim_it, input.end(), '\t')) != input.end())
        {
          std::get<1>(ret).emplace_back(marker_id_to_site_info(delim_it, next_delim_it));
          delim_it = next_delim_it + 1;
        }

        std::get<1>(ret).emplace_back(marker_id_to_site_info(delim_it, input.end()));
      }

      return ret;
    }

  private:
    RdrT* rdr_;
    std::string group_id_;
    std::list<site_info> sites_;
    std::list<site_info>::iterator site_it_;
    std::vector<region> merged_regions_;
    std::vector<region>::iterator region_it_;
    value_type variant_;
  };

  template <typename VecT>
  using variant_group_iterator = basic_variant_group_iterator<indexed_reader, VecT>;

  namespace sav
  {
    template <typename VecT>
    using variant_group_iterator = basic_variant_group_iterator<indexed_reader, VecT>;
  }

  namespace vcf
  {
    template <typename VecT>
    using variant_group_iterator = basic_variant_group_iterator<indexed_reader<1>, VecT>;
  }
}

#endif //SAVVY_VARIANT_GROUP_ITERATOR_HPP