[![Anaconda-Server Badge](https://anaconda.org/bioconda/savvy/badges/installer/conda.svg)](https://anaconda.org/bioconda/savvy)


All branches in this repository are development branches. The latest stable release can be found in the [releases](https://github.com/statgen/savvy/releases) section.

# Savvy Library
Savvy is the official C++ interface for the [SAV file format](sav_spec_v2.md) and offers seamless support for BCF and VCF files.

Since the release of version 2.0, Savvy no longer supports writing of SAV 1.x files but will continue to support reading of existing 1.x files.  

## Installing
The easiest way to install Savvy and its dependencies from source is to use [cget](http://cget.readthedocs.io/en/latest/src/intro.html#installing-cget).
```shell
cget install --prefix <install_prefix> statgen/savvy # default <install_prefix> is ./cget/
```

Installing binaries of Savvy and its dependencies can be done with [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html).
```shell
conda install -c conda-forge -c bioconda savvy
```

## Including in Projects
CMakeLists.txt:
```cmake
# Configure with cmake option: -DCMAKE_TOOLCHAIN_FILE=<install_prefix>/cget/cget.cmake
find_package(savvy REQUIRED)
add_executable(prog main.cpp)
target_link_libraries(prog savvy)
```

## Reading Variants from File 
```c++
#include <savvy/reader.hpp>

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
In addition to the genomic region queries that CSI indices enable for VCF/BCF files, S1R indices also enable SAV files to be queried by record offset.  

### Genomic Queries
```c++
#include <savvy/reader.hpp>

savvy::reader f("chrX.sav");
savvy::variant var;

f.reset_bounds(savvy::genomic_region("X", 60001, 2699520));
while (f >> var)
{
  ...
}

// Shorthand
f.reset_bounds({"X", 154931044, 155260560});
while (f >> var)
{
  ...
}
```

### Slice Queries
```c++
#include <savvy/reader.hpp>
savvy::reader f("chr1.sav");

// Get the 10,000th record through the 19,999th record (0-based amd non-inclusive)
f.reset_bounds(savvy::slice_bounds(10000, 20000));

// Shorthand
f.reset_bounds({20000, 30000});
```

## Subsetting Samples
```c++
#include <savvy/reader.hpp>
savvy::reader f("chr1.sav");
std::vector<std::string> requested = {"ID001","ID002","ID003"};
std::vector<std::string> intersect = f.subset_samples({requested.begin(), requested.end()});

savvy::variant var;
while (f.read(var))
{
  ...
}
```

## Copying Files
```c++
#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

savvy::reader in("in.sav");
savvy::writer out("out.bcf", savvy::file::format::bcf, in.headers(), in.samples());
savvy::variant var;
while (in >> var)
  out << var;
```

## Creating New Files
```c++
#include <savvy/writer.hpp>

std::vector<std::string> sample_ids = {"ID1", "ID2", "ID3"};

std::vector<std::pair<std::string, std::string>> headers = {
  {"fileformat", "VCFv4.2"},
  {"FILTER", "<ID=PASS,Description=\"All filters passed\">"},
  {"contig", "<ID=chr1,length=248956422>"},
  {"INFO",   "<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Counts\">"},
  {"INFO",   "<ID=AN,Number=1,Type=Integer,Description=\"Total Number Allele Counts\">"},
  {"FORMAT", "<ID=GT,Type=Integer,Description=\"Genotype\">"} 
};

savvy::writer out("out.sav", savvy::file::format::sav2, headers, sample_ids);

std::vector<std::int8_t> geno = {0,0,1,0,0,1};

savvy::variant var("chr1", 10000000, "A", {"AC"}); // chrom, pos, ref, alts
var.set_info("AC", std::count(geno.begin(), geno.end(), 1));
var.set_info("AN", geno.size());
var.set_format("GT", geno);

out.write(var);
```

# SAV Command Line Interface
File manipulation for SAV format.

## Import
The `import` sub-command generates a SAV file from a BCF or VCF file. An S1R index is automatically generated and appended to the end of the resulting SAV file.
```shell
sav import file.bcf file.sav
```

## Export
The `export` sub-command can be used to manipulate SAV files and/or convert between file formats.
```shell
sav export --regions chr1,chr2:10000-20000 --sample-ids ID1,ID2,ID3 file.sav > file.vcf
```

## Concatenate
Fast concatenation of SAV files (similar to `bcftools concat --naive`) can be achieved with the `concat` sub-command. This command avoids deserialization of variant data by performing a byte-for-byte copy of compressed variant blocks. The S1R index is also quickly concatenated without having to parse records in the SAV file. 
```shell
sav concat file1.sav file2.sav > concat.sav
```

## Slice Queries
In addition to querying genomic regions, S1R indices can be used to quickly subset records by their offset within a file.
```shell
# export first 1,000 records
sav export --slice 0:1000 file.sav > file.vcf

# export second 1,000 records (1,000-1,999)
sav export --slice 1000:2000 file.sav > file.vcf
```

## Statistics
There are two sub-commands for gathering statistics on sav files. The `stat` command parses the entire file to calculate statistics. The `stat-index` sub-command only parsed the S1R index, making it a faster alternative for some statistics (e.g., number of variant records, chromosomes, etc.).
```shell
sav stat file.sav
sav stat-index file.sav
```

## Sort
The `sort` sub-command sorts variant records by chromosome and position.  It can also be used to sort in descending order, which is supported by S1R indices.
```shell
sav sort unsorted.sav > sorted.sav
sav sort --direction desc unsorted.sav > reversed.sav
```

## Header
The `head` and `rehead` sub-commands are used for retrieving and manipulating header information.
```shell
sav head file.sav > header.txt
sav head --sample-ids file.sav > sample_ids.txt
sav rehead --sample-ids new_ids_file.txt old.sav new.sav
```

## Parameter Trade-offs
| Action | Pro | Con |
|:-------|:----|:----|
|Increasing block size|Smaller file size (especially with PBWT)|Reduces precision of random access|
|Increasing compression level|Smaller file size|Slower compression speed (decompression not affected)|
|Enabling PBWT|Smaller file size when used with some fields|Slower compression and decompression|

# Packaging
```shell
docker build -t savvy-packaging - < packaging-dockerfile-ubuntu16
mkdir -p packages
docker run -v $(pwd):/savvy-src -v $(pwd)/packages:/out savvy-packaging /savvy-src/package-linux.sh /savvy-src /out
```

# Optional Build Targets
* `-DBUILD_TESTS=ON` allows running of tests with `make test`
* `-DBUILD_EVAL=ON` enables building of sav-eval executable used to evaluate deserialization performance
* `-DBUILD_SPARSE_REGRESSION=ON` enables building of sav-at executable