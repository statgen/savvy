# libvc
Interface to various variant calling formats.


#Examples
Some of the functions/classes in these examples are not yet implemented but are here to show direction of where library is going.

## Iterate Markers in Single File 
```c++
vc::open_marker_file("chr1.bcf", [](auto&& file_reader)
{
  for (const auto& marker: make_iterable_marker_stream(file_reader))
  {
    for (const auto& hap: marker)
    {
        
    }
  }
});
```

## Open Multiple Files and Use Standard Algorithms
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