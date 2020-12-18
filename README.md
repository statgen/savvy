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
  var.alts();

  int ac;
  var.get_info("AC", ac);
  
  var.get_format("GT", geno);
  for (int allele : geno)
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

## Multiple Data Vectors
```c++
savvy::reader f("chr1.sav");
savvy::variant var;
std::vector<int> genotypes;
std::vector<float> dosages;
while (f.read(var))
{  
  var.get_format("GT", genotypes);
  for (int gt : genotypes)
  {
    ...
  }
  
  var.get_format("DS", dosages);
  for (float ds : dosages)
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