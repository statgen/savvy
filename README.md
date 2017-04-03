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
vc::indexed_reader f("chr1.cmf", "X", 100000, 199999);
vc::dense_haplotype_vector<float> variant;
while (f >> variant)
{
  ...
}

f.reset_region("X", 200000, 299999);
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

## Converting Files
```c++
vc::vcf::reader bcf_file("file.bcf");

vc::cmf::writer cmf_file("file.cmf");
vc::cmf::output_iterator out_it(cmf_file);

if (subset_file)
  std::copy_if(vc::vcf::dense_variant_iterator{bcf_file}, vc::vcf::dense_variant_iterator{}, out_it, 
    [](const auto& v) { return (vc::calculate_allele_frequency(v) < 0.1); });
else
  std::copy(vc::vcf::dense_variant_iterator{bcf_file}, vc::vcf::dense_variant_iterator{}, out_it);
```