# Marker Group Iterator

```c++
#include <regex>
#include <tuple>
#include <iterator>
#include <cstddef>


template <typename VecT>
class marker_group_iterator
{
public:
  typedef marker_group_iterator<VecT> self_type;
  typedef std::ptrdiff_t difference_type;
  typedef savvy::variant<VecT> value_type;
  typedef const value_type& reference;
  typedef const value_type* pointer;
  typedef std::input_iterator_tag iterator_category;

  marker_group_iterator(savvy::indexed_reader& rdr, std::string marker_group_file_line) :
    rdr_(&rdr)
  {
    std::tie(group_id_, marker_ids_) = parse_marker_group_line(marker_group_file_line);
    std::vector<savvy::region> un_merged_regions(marker_ids_.size(), {""});
    auto out_it = un_merged_regions.begin();
    for (auto in_it = marker_ids_.begin(); in_it != marker_ids_.end(); ++in_it,++out_it)
      *out_it = marker_id_to_region(*in_it);

    merged_regions_ = savvy::region::merge(un_merged_regions.begin(), un_merged_regions.end());
    marker_id_it_ = marker_ids_.begin();
    region_it_ = merged_regions_.begin();

    rdr_->reset_region(*region_it_);

    increment();
  }

  marker_group_iterator() :
    rdr_(nullptr)
  {
  }

  self_type& operator++(){ increment(); return *this; }
  void operator++(int) { increment(); }
  reference operator*() { return variant_; }
  pointer operator->() { return &variant_; }
  bool operator==(const self_type& rhs) { return (rdr_ == rhs.rdr_); }
  bool operator!=(const self_type& rhs) { return (rdr_ != rhs.rdr_); }

  const std::string& group_id() const { return group_id_; }
  const std::list<std::string>& marker_ids() const { return marker_ids_; }
private:
  void increment()
  {
    while (region_it_ != merged_regions_.end())
    {
      if (!region_it_->chromosome().empty())
      {
        while (rdr_->read(variant_, variant_.data()) && marker_id_it_ != marker_ids_.end())
        {
          while (marker_id_it_ != marker_ids_.end())
          {
            savvy::region reg = marker_id_to_region(*marker_id_it_);
            if (reg.from() >= variant_.position() || reg.chromosome() != variant_.chromosome())
              break;
            ++marker_id_it_;
          }

          std::string target_id = *marker_id_it_;
          std::string current_id = variant_.chromosome() + ":" + std::to_string(variant_.position()) + "_" + variant_.ref() + "/" + variant_.alt();
          if ((*marker_id_it_) == variant_.chromosome() + ":" + std::to_string(variant_.position()) + "_" + variant_.ref() + "/" + variant_.alt())
          {
            ++marker_id_it_;
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

  // [CHROM]:[POS]_[REF]/[ALT]
  static savvy::region marker_id_to_region(const std::string& marker_id)
  {
    auto colon_it = std::find(marker_id.begin(), marker_id.end(), ':');
    std::string chrom(marker_id.begin(), colon_it);
    if (colon_it != marker_id.end())
    {
      auto underscore_it = std::find(++colon_it, marker_id.end(), '_');
      std::string str_pos(colon_it, underscore_it);
      if (underscore_it != marker_id.end())
      {
        auto slash_it = std::find(++underscore_it, marker_id.end(), '/');
        std::string ref(underscore_it, slash_it);
        if (slash_it != marker_id.end())
        {
          std::string alt(++slash_it, marker_id.end());
          std::size_t length = std::max(ref.size(), alt.size());
          if (length > 0)
          {
            std::uint64_t pos = std::atoll(str_pos.c_str());
            return savvy::region{chrom, pos, pos + length - 1};
          }
        }
      }
    }

    return savvy::region{""};
  }

  static std::tuple<std::string, std::list<std::string>> parse_marker_group_line(const std::string& input)
  {
    std::tuple<std::string, std::list<std::string>> ret;
    auto delim_it = std::find(input.begin(), input.end(), '\t');
    if (delim_it != input.end())
    {
      std::get<0>(ret) = std::string(input.begin(), delim_it);
      ++delim_it;

      std::string::const_iterator next_delim_it;
      while ((next_delim_it = std::find(delim_it, input.end(), '\t')) != input.end())
      {
        std::get<1>(ret).emplace_back(delim_it, next_delim_it);
        delim_it = next_delim_it + 1;
      }

      std::get<1>(ret).emplace_back(delim_it, input.end());
    }

    return ret;
  }
private:
  savvy::indexed_reader* rdr_;
  std::string group_id_;
  std::list<std::string> marker_ids_;
  std::list<std::string>::iterator marker_id_it_;
  std::vector<savvy::region> merged_regions_;
  std::vector<savvy::region>::iterator region_it_;
  value_type variant_;
};

std::ifstream marker_group_file("/Users/lefaivej/Developer/projects/savvy/test_file_groups.grp", std::ios::binary);
savvy::indexed_reader marker_file("marker_file.sav", {""}, savvy::fmt::dosage);
std::string marker_group_line;
std::vector<float> group_dose;
while (std::getline(marker_group_file, marker_group_line))
{
  marker_group_iterator<savvy::compressed_vector<float>> it(marker_file, marker_group_line);
  marker_group_iterator<savvy::compressed_vector<float>> end{};

  group_dose.resize(0);
  if (it != end)
    group_dose.resize(it->data().size(), 0.f);

  for ( ; it != end; ++it)
  {
    it.group_id();
    it.marker_ids();

    std::cerr << it->chromosome() << ":" << it->position() << "_" << it->ref() << "/" << it->alt() << std::endl;
    it->data(); // dosages

    assert(group_dose.size() == it->data().size());
    for (auto dose_it = it->data().begin(); dose_it != it->data().end(); ++dose_it)
      group_dose[dose_it.offset()] += *dose_it;
  }
}
```