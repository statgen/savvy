#include <cmath>
#include "sav/merge.hpp"
#include "savvy/reader.hpp"
#include "savvy/utility.hpp"

#include <cstdlib>
#include <getopt.h>

#include <fstream>
#include <vector>
#include <set>

class merge_prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;

  std::vector<option> long_options_;
  std::vector<std::string> input_paths_;
  std::string output_path_;
  int compression_level_ = -1;
  std::uint16_t block_size_ = default_block_size;
  savvy::fmt format_ = savvy::fmt::allele;
  bool help_ = false;
public:
  merge_prog_args() :
    long_options_(
      {
        {"block-size", required_argument, 0, 'b'},
        {"format", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {"output", no_argument, 0, 'o'},
        {0, 0, 0, 0}
      }),
    output_path_("/dev/stdout")
  {
  }

  const std::vector<std::string>& input_paths() { return input_paths_; }
  const std::string& output_path() { return output_path_; }
  std::uint8_t compression_level() const { return std::uint8_t(compression_level_); }
  std::uint16_t block_size() const { return block_size_; }
  savvy::fmt format() const { return format_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "----------------------------------------------\n";
    os << "Usage: sav merge [opts] input.{sav,vcf,vcf.gz,bcf} input2.{sav,vcf,vcf.gz,bcf} [additional_input.{sav,vcf,vcf.gz,bcf} ...] \n";
    os << "\n";
    os << " -#               : # compression level (1-19, default: " << default_compression_level << ")\n";
    os << " -b, --block-size : Number of markers in compression block (0-65535, default: " << default_block_size << ")\n";
    os << " -f, --format     : Format field to copy (GT or HDS, default: GT)\n";
    os << " -h, --help       : Print usage\n";
    os << " -o, --output     : Output file (default: stdout)\n";
    os << "----------------------------------------------\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "0123456789b:f:ho:", long_options_.data(), &long_index )) != -1)
    {
      std::string str_opt_arg(optarg ? optarg : "");
      char copt = char(opt & 0xFF);
      switch (copt)
      {
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
          if (compression_level_ < 0)
            compression_level_ = 0;
          compression_level_ *= 10;
          compression_level_ += copt - '0';
          break;
        case 'b':
          block_size_ = std::uint16_t(atoi(optarg) > 0xFFFF ? 0xFFFF : atoi(optarg));
          break;
        case 'f':
          if (str_opt_arg == "HDS")
          {
            format_ = savvy::fmt::haplotype_dosage;
          }
          else if (str_opt_arg != "GT")
          {
            std::cerr << "Invalid format field value (" << str_opt_arg << ")\n";
            return false;
          }
          break;
        case 'h':
          help_ = true;
          break;
        case 'o':
          output_path_ = std::string(optarg);
        default:
          return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count < 2)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else
    {
      input_paths_.reserve((unsigned)remaining_arg_count);
      for (int i = optind; i < argc; ++i)
      {
        input_paths_.emplace_back(argv[i]);
      }
    }


    if (compression_level_ < 0)
      compression_level_ = default_compression_level;
    else if (compression_level_ > 19)
      compression_level_ = 19;

    return true;
  }
};

int merge_main(int argc, char** argv)
{
  merge_prog_args args;
  if (!args.parse(argc, argv))
  {
    args.print_usage(std::cerr);
    return EXIT_FAILURE;
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }

  std::vector<savvy::vcf::reader<1>> input_files;
  input_files.reserve(args.input_paths().size());
  for (auto it = args.input_paths().begin(); it != args.input_paths().end(); ++it)
    input_files.emplace_back(*it, args.format());
  std::vector<savvy::site_info> sites(args.input_paths().size());
  std::vector<std::vector<float>> genos(args.input_paths().size()); //std::vector<savvy::compressed_vector<float>> genos(args.input_paths().size());

  {
    std::size_t num_samples = 0;
    for (auto& f : input_files)
      num_samples += f.sample_size();
    std::vector<std::string> sample_ids;
    sample_ids.reserve(num_samples);

    std::vector<std::pair<std::string,std::string>> merged_headers;

    std::set<std::string> info_fields;

    for (auto& f : input_files)
    {
      for (auto it = f.samples_begin(); it != f.samples_end(); ++it)
        sample_ids.emplace_back(*it);

      merged_headers.reserve(merged_headers.size() + (f.headers().end() - f.headers().begin()));
      for (auto it = f.headers().begin(); it != f.headers().end(); ++it)
      {
        if (it->first != "INFO" || info_fields.insert(savvy::parse_header_id(it->second)).second)
        {
          merged_headers.push_back(*it);
        }
      }
    }


    if (info_fields.insert("FILTER").second)
      merged_headers.insert(merged_headers.begin(), {"INFO","<ID=FILTER,Description=\"Variant filter\">"});
    if (info_fields.insert("QUAL").second)
      merged_headers.insert(merged_headers.begin(), {"INFO","<ID=QUAL,Description=\"Variant quality\">"});
    if (info_fields.insert("ID").second)
      merged_headers.insert(merged_headers.begin(), {"INFO","<ID=ID,Description=\"Variant ID\">"});

    savvy::sav::writer::options opts;
    opts.compression_level = args.compression_level();
    opts.block_size = args.block_size();
    opts.data_format = args.format();

    savvy::sav::writer output(args.output_path(), sample_ids.begin(), sample_ids.end(), merged_headers.begin(), merged_headers.end(), opts);

    // TODO: This will only work for single chromosome.

    std::vector<std::size_t> matching_indices;
    matching_indices.reserve(input_files.size());

    std::size_t min_pos_index = input_files.size();
    std::uint64_t min_pos = std::numeric_limits<std::uint64_t>::max();
    for (std::size_t i = 0; i < input_files.size(); ++i)
    {
      if (!input_files[i].read(sites[i], genos[i]))
        genos[i].resize(0);

      if (genos[i].size() && sites[i].position() < min_pos)
      {
        min_pos = sites[i].position();
        min_pos_index = i;
      }
    }

    while (min_pos_index < input_files.size())
    {
      matching_indices.resize(0);
      for (std::size_t i = 0; i < input_files.size(); ++i)
      {
        if (i == min_pos_index
          || (genos[i].size()
          && sites[i].position() == sites[min_pos_index].position()
          && sites[i].ref() == sites[min_pos_index].ref()
          && sites[i].alt() == sites[min_pos_index].alt()))
        {
          matching_indices.emplace_back(i);
        }
      }

      //output.write(sites[min_pos_index]);
      for (auto it = matching_indices.begin(); it != matching_indices.end(); ++it)
      {
        //output.write(genos[*it]);
        if (!input_files[*it].read(sites[*it], genos[*it]))
          genos[*it].resize(0);
      }

      min_pos_index = input_files.size();
      min_pos = std::numeric_limits<std::uint64_t>::max();

      for (std::size_t i = 0; i < input_files.size(); ++i)
      {
        if (genos[i].size() && sites[i].position() < min_pos)
        {
          min_pos = sites[i].position();
          min_pos_index = i;
        }
      }
    }
  }

  return EXIT_FAILURE;
}