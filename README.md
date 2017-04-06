# libvc
Interface to various variant calling formats.


## Read Variants from File 
```c++
vc::reader f("chr1.cmf");
vc::sparse_haplotype_vector<float> variant;
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
vc::dense_haplotype_vector<float> variant;
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

## Custom Haplotype Vectors
The hapotype_vector class utilizes the "mixin" pattern to efficiently support 3rd-party linear algebra libraries. 
```c++
vc::hapotype_vector<std::vector<float>> std_vector;
vc::hapotype_vector<vc::compressed_vector<double>> vc_sparse_vector;
vc::hapotype_vector<boost::numeric::ublas::compressed_vector<float>> ublas_sparse_vector;
```

## Read Predicates
Reading genotypes can be bypassed when using a read predicate.
```c++
vc::indexed_reader f("chr1.cmf");
vc::sparse_haplotype_vector<float> variant;

while (f.read_if(variant, [](const auto& v) { return std::stof(v.prop("AF")) < 0.1; }))
{
  ...
}
```
```c++
vc::indexed_reader f("chr1.cmf");
vc::dense_haplotype_vector<float> buf;

{
  struct file_statistics_functor
  {
    template <typename T>
    bool operator()(const T& v)
    {
      ++variant_count;
      genotype_count += std::stoi(v.prop("NS"));
      return false;
    }
    
    std::size_t variant_count = 0;
    std::size_t genotype_count = 0;
  } file_statistics;
  
  while (f.read_if(buf, file_statistics)) { }
  
  // Consume file statisicts ...
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