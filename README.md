# libvc
Interface to various variant calling formats.


## Read Variants from File 
```c++
vc::reader f("chr1.cmf");
vc::sparse_allele_vector<float> variant;
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
vc::indexed_reader f("chr1.cmf", {"X", 100000, 199999});
vc::dense_allele_vector<float> variant;
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
vc::reader f("chr1.cmf");
vc::sparse_variant_iterator<float> it(f);
while (it != vc::sparse_variant_iterator<float>{})
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
vc::allele_vector<std::vector<float>> std_vector;
vc::allele_vector<vc::compressed_vector<double>> vc_sparse_vector;
vc::allele_vector<boost::numeric::ublas::compressed_vector<float>> ublas_sparse_vector;
```

## Read Predicates
Reading genotypes can be bypassed when using a read predicate.
```c++
vc::indexed_reader f("chr1.cmf");
vc::sparse_allele_vector<float> variant;

while (f.read_if(variant, [](const auto& v) { return std::stof(v.prop("AF")) < 0.1; }))
{
  ...
}
```
```c++
vc::indexed_reader f("chr1.cmf");
vc::dense_allele_vector<float> buf;

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

## Custom Allele Status Values
Values used to load vectors can be customized.
```c++
vc::reader f("chr1.cmf");
vc::dense_allele_vector<float> variant;

const float is_missing = std::numeric_values<float>::epsilon();
const float has_alt = 1;
const float has_ref = 0;

while (f.read(variant, is_missing, has_alt, has_ref))
{
  ...
}
```

## Converting Files
```c++
vc::vcf::reader bcf_file("file.bcf");
vc::vcf::dense_variant_iterator it(bcf_file);
vc::vcf::dense_variant_iterator eof;

if (it != eof)
{
  vc::cmf::writer output("file.cmf", it->chromosome(), vc::get_ploidy(bcf_file, *it), bcf_file.samples_begin(), bcf_file.samples_end());
  vc::cmf::output_iterator out_it(cmf_file);
  
  if (subset_file)
    std::copy_if(std::move(it), eof, out_it, 
      [](const auto& v) { return (vc::calculate_allele_frequency(v) < 0.1); });
  else
    std::copy(std::move(it), eof, out_it);
}
```