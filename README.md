# libvc
Interface to various variant calling formats.


#Example Usage
Some of the functions and classes in these examples are not yet implemented but are here to show direction of where library is going.

## Iterate Markers in Single File 
```c++
vc::open_marker_file("chr1.bcf", [](auto&& file_reader)
{
  for (const auto& marker: make_iterable_marker_stream(file_reader))
  {
    for (const vc::allele_status& haplotype: marker)
    {
        
    }
  }
});
```

## Open Multiple Files
```c++
class has_low_af_filter
{
public:
  has_low_af_filter(double threshold) : threshold_(threshold) {}
  template <typename T>
  bool operator()(const T& m) const
  {
    return (vc::calculate_allele_frequency(m) < threshold_);
  }
private:
  const double threshold_;
};

vc::open_marker_files(std::make_tuple("chr1.bcf", "chr1.m3vcf", "chr1.cmf"), [](auto&& file_reader1, auto&& file_reader2, auto&& file_reader3)
{
  {
    auto markers = make_iterable_marker_stream(file_reader1);
    std::size_t marker_count = std::count(markers.begin(), markers.end());
  }
  
  {
    auto markers = make_iterable_marker_stream(file_reader2);
    std::size_t markers_with_low_af_count = std::count_if(markers.begin(), markers.end(), has_low_af_filter(0.01));
  }
  
  {
    auto markers = make_iterable_marker_stream(file_reader3);
    auto it = std::find_if(markers.begin(), markers.end(), has_low_af_filter(0.0001));
  }
});
```

## Alternate Variadic Syntax for Multiple Files
```c++
struct file_handler
{
  template <typename T, typename T2>
  void operator()(T&& input_file_reader, T2&& input_file_reader2)
  {

  }
};

vc::open_marker_files(triple_file_handler(), "chr1.bcf", "chr1.m3vcf");
```

## Converting Files
```c++
auto bcf_markers = vc::make_iterable_marker_stream(vc::vcf::reader("chr1.bcf"));

vc::cmf::writer cmf_file;
vc::cmf::output_iterator out_it(cmf_file);

if (subset_file)
  std::copy_if(bcf_markers.begin(), bcf_markers.end(), out_it, [](const vc::vcf::marker& m) { return (vc::calculate_allele_frequency(m) < 0.1); });
else
  std::copy(bcf_markers.begin(), bcf_markers.end(), out_it);
```

## Specializing File Formats
```c++
template <typename T>
std::vector<float> init_haplotypes(const T& marker,
  const float missing_value = std::numeric_limits<float>::quiet_NaN(),
  const float alt_value = 1.0,
  const float ref_value = 0.0)
{
  std::vector<float> ret(marker.haplotype_count(), ref_value);
  
  std::size_t i = 0;
  for (auto it = marker.begin(); it != marker.end(); ++it,++i)
  {
    if (*it == vc::allele_status::has_alt)
      ret[i] = alt_value;
    else if (*it == vc::allele_status::is_missing)
      ret[i] = missing_value;
  }
  
  return ret;
}

template <>
std::vector<float> init_haplotypes<vc::cmf::marker>(const vc::cmf::marker& marker, const float missing_value, const float alt_value, const float ref_value)
{
  std::vector<float> ret(marker.haplotype_count(), ref_value);
  
  for (auto it = marker.non_ref_begin(); it != marker.non_ref_end(); ++it)
    ret[it->offset] = (it->status == vc::allele_status::has_alt ? alt_value : missing_value);
  
  return ret;
}


vc::open_marker_file(argv[1], [](auto&& file_reader)
{
  std::vector<float> phenotypes = init_phenotypes(file_reader.samples_begin(), file_reader.samples.end());
  for (const auto& marker: make_iterable_marker_stream(file_reader))
  {
    std::vector<float> haplotypes = init_haplotypes(marker, std::numeric_limits<float>::epsilon());
    float avg = std::accumulate(haplotypes.begin(), haplotypes.end(), 0.0);
    float dot_product = std::inner_product(haplotypes.begin(), haplotypes.end(), phenotypes.begin(), 0.0);
  }
});
```

## Indexed Files
```c++
vc::open_indexed_marker_file(argv[1], [](auto&& indexed_file_reader)
{
  std::size_t start_pos = 100000;
  std::size_t end_pos = 200000;

  {
    // Iterate chromosome position range with range_query class.
    auto marker_range = vc::make_range_query(indexed_file_reader, start_pos, end_pos);
    std::for_each(marker_range.begin(), marker_range.end(), [](const auto& marker)
    {
      std::for_each(marker.begin(), marker.end(), [](const vc::allele_status& haplotype)
      {

      });
    });
  }

  {
    // Seek to lower bound of chromosome position and iterate from there.
    indexed_file_reader.lower_bound_seek(start_pos);
    for (const auto& marker: make_iterable_marker_stream(indexed_file_reader))
    {
      for (const vc::allele_status& hap: marker)
      {

      }

      if (marker.pos() >= end_pos)
        break;
    }
  }
});
```