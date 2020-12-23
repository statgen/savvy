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

## Copying files
```c++
#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

savvy::reader in("in.sav");
savvy::writer out("out.bcf", savvy::fmt::bcf, in.headers(), in.samples());
savvy::variant var;
while (in >> var)
  out << var;
```

## Creating New Files
```c++
#include <savvy/writer.hpp>

std::vector<std::string> sample_ids = {"ID1", "ID2", "ID3"};
std::vector<std::pair<std::string, std::string>> headers = {};

savvy::writer out("out.sav", savvy::fmt::sav2, headers, sample_ids);

std::vector<int8_t> geno = {0,0,1,0,0,1};

savvy::variant var("chr1", 10000000, "A", {"AC"}); // chrom, pos, ref, alts
var.set_info("AC", std::accumulate(geno.begin(), geno.end(), 0));
var.set_info("AN", geno.size());
var.set_format("GT", geno);

out.write(var);
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

## Parameter Trade-offs
| Action | Pro | Con |
|:-------|:----|:----|
|Increasing block size|Smaller file size (especially with pbwt)|Reduces precision of random access|
|Increasing compression level|Smaller file size|Slower compression speed (decompression not affected)|
|Enabling PBWT|Smaller file size when used with some fields|Slower compression and decompression|

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