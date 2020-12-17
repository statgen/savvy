This branch contains version 2.0 progress. The latest stable release can be found in the [releases](https://github.com/statgen/savvy/releases) section.

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
savvy::reader f("chr1.sav");
savvy::variant var;
std::vector<int> geno;
  
while (f >> var)
{
  var.position();
  var.chromosome();
  var.ref();
  var.alt();
  
  var.get_format("GT", geno);
  for (const int& allele : geno)
  {
    ...
  }
}
```

## Random Access
```c++
savvy::indexed_reader f("chr1.sav");
savvy::variant var;

f.reset_bounds({"X", 60001, 2699520});
while (f >> var)
{
  ...
}

f.reset_bounds({"X", 154931044, 155260560});
while (f >> var)
{
  ...
}
```

## Subsetting Samples
```c++
savvy::reader f("chr1.sav");
std::vector<std::string> requested = {"ID001","ID002","ID003"};
std::vector<std::string> intersect = f.subset_samples({requested.begin(), requested.end()});

savvy::variant var;
while (f.read(var))
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
savvy::variant var;
std::vector<float> geno;
std::vector<float> pheno(f.sample_size());

while (f.read(var))
{
  var.get_format("GT", geno);
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
savvy::reader f("chr1.sav");
savvy::variant var;
std::vector<int> genotypes;
std::vector<float> dosages;
while (f.read(var))
{
  var.position();
  var.chromosome();
  var.ref();
  var.alt();
  
  var.get_format("GT", genotypes);
  for (const int& gt : genotypes)
  {
    ...
  }
  
  var.get_format("DS", dosages);
  for (const float& ds : dosages)
  {
    ...
  }
}
```

# SAV Command Line Interface
File manipulation for SAV format.

## Import
```shell
sav import --sort --index file.bcf file.sav
```

## Concatenate
```shell
sav concat file1.sav file2.sav > concat.sav
```

## Export
```shell
sav export --regions chr1,chr2:10000-20000 --sample-ids ID1,ID2,ID3 file.sav > file.vcf
```

## Slice Queries
In addition to querying genomic regions, S1R indexes can be used to quickly subset records by their offset within a file.
```shell
# export first 1,000 records
sav export --slice 0:1000 file.sav > file.vcf

# export second 1,000 records (1,000-1,999)
sav export --slice 1000:2000 file.sav > file.vcf
```

# Packaging
```shell
docker build -t savvy-packaging - < packaging-dockerfile-ubuntu16
mkdir -p packages
docker run -v $(pwd):/savvy-src -v $(pwd)/packages:/out savvy-packaging /savvy-src/package-linux.sh /savvy-src /out
```