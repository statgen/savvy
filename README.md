# Savvy Library
Interface to various variant calling formats.

## Installing
The easiest way to install savvy and its dependencies is to use [cget](http://cget.readthedocs.io/en/latest/src/intro.html#installing-cget).
```bash
cget install --prefix <install_prefix> statgen/savvy
```
CMakeLists.txt:
```cmake
add_executable(prog main.cpp)
target_link_libraries(prog savvy hts z zstd)
```

## Read Variants from File 
```c++
savvy::reader f("chr1.sav", savvy::fmt::gt);
savvy::variant<std::vector<float>> var;

while (f >> var)
{
  var.locus();
  var.chromosome();
  var.ref();
  var.alt();
  for (const float& allele : var.data())
  {
    ...
  }
}
```

## Indexed Files
```c++
savvy::indexed_reader f("chr1.sav", {"X", 100000, 199999}, savvy::fmt::gt);
savvy::variant<std::vector<float>> var;

while (f >> v)
{
  ...
}

f.reset_region({"X", 200000, 299999});
while (f >> v)
{
  ...
}
```

## Read annotations and genotypes separately.
```c++
savvy::reader f("chr1.sav", savvy::fmt::gt);
savvy::site_info anno;
std::vector<float> alleles;

while (f.read(anno, alleles))
{
  anno.locus();
  anno.chromosome();
  anno.ref();
  anno.alt();
  for (const float& allele : alleles)
  {
    ...
  }
}
```

## Subsetting Samples
```c++
savvy::reader f("chr1.sav", savvy::fmt::gt);
std::vector<std::string> requested = {"ID001","ID002","ID003"};
std::vector<std::string> intersect = f.subset_samples({requested.begin(), requested.end()});

savvy::site_info anno;
std::vector<float> alleles;
while (f.read(anno, alleles))
{
  ...
}
```

## 3rd-party Vectors
The reader classes utilize generic programming to efficiently support 3rd-party linear algebra libraries. 
```c++
savvy::site_info anno;
std::vector<float> std_vector;
savvy::compressed_vector<double> savvy_sparse_vector;
boost::numeric::ublas::compressed_vector<float> ublas_sparse_vector;

savvy::reader f("chr1.sav", savvy::fmt::ds);
f.read(anno, std_vector);
f.read(anno, savvy_sparse_vector);
f.read(anno, ublas_sparse_vector);
```

## Read Predicates
Reading genotypes can be bypassed when using a read predicate.
```c++
savvy::indexed_reader< f("chr1.sav", savvy::fmt::gt);
savvy::site_info anno;
savvy::compressed_vector<float> gt;

// Read only if allele frequency is less than 0.1.
while (f.read_if([](const site_info& v) { return std::stof(v.prop("AF")) < 0.1; }, anno, gt))
{
  ...
}
```
```c++
savvy::indexed_reader f("chr1.sav", savvy::fmt::gt);
savvy::anno;
std::vector<float> buf;

{
  struct 
  {
    bool operator()(const site_info& v)
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
  
  while (f.read_if(file_statistics, anno, buf)) { }
  
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
  float se_line       = 0.0f; for (std::size_t i = 0; i < n; ++i) se_line += std::pow(y[i] - fx(x[i]), 2);
  const float y_mean  = s_y / n;
  float se_y_mean     = 0.0f; for (std::size_t i = 0; i < n; ++i) se_y_mean += std::pow(y[i] - y_mean, 2);
  const float r2      = 1 - se_line / se_y_mean;

  return std::make_tuple(m, b, r2); // slope, y-intercept, r-squared
};

savvy::reader f("chr1.sav", savvy::fmt::ac);
savvy::site_info anno;
std::vector<float> geno;
std::vector<float> pheno(f.sample_size());

while (f.read(anno, geno))
{
  auto [ m, b, r2 ] = lin_reg(geno, pheno);
  // ...
}
```

## Armadillo Example
```c++
auto arma_lin_reg = [](const arma::Col<float>& x, const arma::Col<float>& y)
{
  float m = arma::as_scalar(arma::cov(x, y) / arma::cov(x, x));
  float b = arma::as_scalar(arma::mean(y) - m * arma::mean(x));
  float r2 = arma::as_scalar(arma::square(arma::cor(x, y)));

  return std::make_tuple(m, b, r2); // slope, y-intercept, r-squared
};

savvy::reader f("chr1.sav", savvy::fmt::gt);
savvy::site_info anno;
savvy::armadillo::dense_vector<float> geno;
arma::Col<float> pheno(f.sample_size() * 2);

while (f.read(anno, geno))
{
  auto [ m, b, r2 ] = arma_lin_reg(geno, pheno);
  // ...
}
```

## Multiple Data Vectors
```c++
savvy::vcf::reader<2> f("chr1.bcf", savvy::fmt::ac, savvy::fmt::ds);
savvy::site_info anno;
std::vector<float> genotypes;
std::vector<float> dosages;
while (f.read(anno, genotypes, dosages))
{
  anno.locus();
  anno.chromosome();
  anno.ref();
  anno.alt();
  
  for (const float& gt : genotypes)
  {
    ...
  }
  
  for (const float& ds : dosages)
  {
    ...
  }
}
```

## C++ 17 Class Template Argument Deduction
C++ 17 supports class template argument deduction, which means that template arguments can be deduced by constructor arguments. Compilers that do not support this must specify the number of data vectors to read as a template argument to the reader.
```c++
// With C++ 17
savvy::vcf::reader f("chr1.bcf", savvy::fmt::gt);
savvy::vcf::reader f("chr1.bcf", savvy::fmt::gt, savvy::fmt::gl);

// Without C++ 17
savvy::vcf::reader<1> f("chr1.bcf", savvy::fmt::gt);
savvy::vcf::reader<2> f("chr1.bcf", savvy::fmt::gt, savvy::fmt::gl);
``` 

# SAV Command Line Interface
File manipulation for SAV format.

## Import
```shell
sav import --sort --index file.bcf file.sav
```

## Merge
```shell
sav merge file1.sav file2.sav > merged.sav
```

## Export
```shell
sav export --regions chr1,chr2:10000-20000 --sample-ids ID1,ID2,ID3 file.sav > file.vcf
```

#Packaging
```shell
docker build -t savvy-packaging - < packaging-dockerfile-ubuntu16
mkdir -p packages
docker run -v $(pwd):/savvy-src -v $(pwd)/packages:/out savvy-packaging /savvy-src/package-linux.sh /savvy-src /out
```