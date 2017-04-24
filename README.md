# libsavvy
Interface to various variant calling formats.


## Read Variants from File 
```c++
savvy::reader f("chr1.cmf");
savvy::sparse_allele_vector<float> variant;
while (f >> variant)
{
  variant.locus();
  variant.chromosome();
  variant.ref();
  variant.alt();
  for (const float& haplotype : variant)
  {
    ...
  }
}
```

## Indexed Files
```c++
savvy::indexed_reader f("chr1.cmf", {"X", 100000, 199999});
savvy::dense_allele_vector<float> variant;
while (f >> variant)
{
  ...
}

f.reset_region({"X", 200000, 299999});
while (f >> variant)
{
  ...
}
```

## Input Iterators 
```c++
savvy::reader f("chr1.cmf");
savvy::sparse_variant_iterator<float> it(f);
while (it != savvy::sparse_variant_iterator<float>{})
{
  it->locus();
  it->chromosome();
  it->ref();
  it->alt();
  for (const float& haplotype : *it)
  {
    ...
  }
  ++it;
}
```

## Custom Vectors
The allele_vector class utilizes the "mixin" pattern to efficiently support 3rd-party linear algebra libraries. 
```c++
savvy::allele_vector<std::vector<float>> std_vector;
savvy::allele_vector<savvy::compressed_vector<double>> savvy_sparse_vector;
savvy::allele_vector<boost::numeric::ublas::compressed_vector<float>> ublas_sparse_vector;
```

## Read Predicates
Reading genotypes can be bypassed when using a read predicate.
```c++
savvy::indexed_reader f("chr1.cmf");
savvy::sparse_allele_vector<float> variant;

while (f.read_if(variant, [](const auto& v) { return std::stof(v.prop("AF")) < 0.1; }))
{
  ...
}
```
```c++
savvy::indexed_reader f("chr1.cmf");
savvy::dense_allele_vector<float> buf;

{
  struct 
  {
    template <typename T>
    bool operator()(const T& v)
    {
      ++variant_count;
      if (v.ref().size() == v.alt().size())
        ++snp_count;
      genotype_count += std::stoi(v.prop("NS"));
      return false;
    }
    
    std::size_t variant_count = 0;
    std::size_t genotype_count = 0;
    std::size_t snp_count = 0;
  } file_statistics;
  
  while (f.read_if(buf, file_statistics)) { }
  
  // Consume file statisicts ...
}
```

## Simple Linear Regression Example
```c++
auto lin_reg = [](const std::vector<float>& x, const std::vector<float>& y)
{
  const std::size_t n = x.size();
  const float s_x     = std::accumulate(x.begin(), x.end(), 0.0f);
  const float s_y     = std::accumulate(y.begin(), y.end(), 0.0f);
  const float s_xx    = std::inner_product(x.begin(), x.end(), x.begin(), 0.0f);
  const float s_xy    = std::inner_product(x.begin(), x.end(), y.begin(), 0.0f);
  const float m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const float b       = (s_y - m * s_x) / n;
  auto fx             = [m,b](float x) { return m * x + b; };
  float se_line       = 0.0f; for (std::size_t i = 0; i < n; ++i) se_line += y[i] - fx(x[i]);
  const float y_mean  = s_y / n;
  float se_y_mean     = 0.0f; for (std::size_t i = 0; i < n; ++i) se_y_mean += y[i] - y_mean;
  const float r2      = 1 - se_line / se_y_mean;

  return std::make_tuple(m, b, r2); // slope, y-intercept, r-squared
};

savvy::dense_allele_vector<float> variant;
std::vector<float> pheno((r.samples_end() - r.samples_begin()) * 2);

savvy::reader f("chr1.cmf");
while (f >> variant)
{
  auto [ m, b, r2 ] = lin_reg(variant, pheno);
  // ...
}
```

## Custom Missing Allele Value
The default value for missing genotypes is `std::numeric_limits<T::value_type>::quiet_NaN()`. This can be overwritten when using the read method.
```c++
savvy::reader f("chr1.cmf");
savvy::dense_allele_vector<float> variant;

const float missing_value = std::numeric_values<float>::epsilon();
while (f.read(variant, missing_value))
{
  ...
}
```

## Converting Files
```c++
savvy::vcf::reader bcf_file("file.bcf");
savvy::vcf::dense_variant_iterator it(bcf_file);
savvy::vcf::dense_variant_iterator eof;

if (it != eof)
{
  savvy::cmf::writer output("file.cmf", it->chromosome(), savvy::get_ploidy(bcf_file, *it), bcf_file.samples_begin(), bcf_file.samples_end());
  savvy::cmf::output_iterator out_it(cmf_file);
  
  if (subset_file)
    std::copy_if(std::move(it), eof, out_it, 
      [](const auto& v) { return (savvy::calculate_allele_frequency(v) < 0.1); });
  else
    std::copy(std::move(it), eof, out_it);
}
```