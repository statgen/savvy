#include "sav/sort.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/utility.hpp"

#include <cstdlib>
#include <getopt.h>

#include <fstream>
#include <deque>
#include <vector>
#include <set>
#include <random>
#include <algorithm>

class sort_prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;

  std::vector<option> long_options_;
  std::string input_path_;
  std::string output_path_;
  int compression_level_ = -1;
  std::uint16_t block_size_ = default_block_size;
  bool help_ = false;
  std::uint32_t ploidy_ = 2;
public:
  sort_prog_args() :
    long_options_(
      {
        {"block-size", required_argument, 0, 'b'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
      }),
    output_path_("/dev/stdout")
  {
  }

  const std::string& input_path() const { return input_path_; }
  const std::string& output_path() const { return output_path_; }
  std::uint8_t compression_level() const { return std::uint8_t(compression_level_); }
  std::uint16_t block_size() const { return block_size_; }
  std::uint32_t ploidy() const { return ploidy_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "----------------------------------------------\n";
    os << "Usage: sav sort [opts ...] [in.sav] [out.sav]\n";
    os << "\n";
    os << " -#               : # compression level (1-19, default: " << default_compression_level << ")\n";
    os << " -b, --block-size : Number of markers in compression block (0-65535, default: " << default_block_size << ")\n";
    os << " -h, --help       : Print usage\n";
    os << "----------------------------------------------\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "0123456789b:h", long_options_.data(), &long_index )) != -1)
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
        case 'h':
          help_ = true;
          break;
        default:
          return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 0)
    {
      input_path_ = "/dev/stdin";
      output_path_ = "/dev/stdout";
    }
    else if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
      output_path_ = "/dev/stdout";
    }
    else if (remaining_arg_count == 2)
    {
      input_path_ = argv[optind];
      output_path_ = argv[optind + 1];
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    if (compression_level_ < 0)
      compression_level_ = default_compression_level;
    else if (compression_level_ > 19)
      compression_level_ = 19;

    return true;
  }
};

typedef savvy::variant<savvy::compressed_vector<float>> var_type;

struct left_comparator
{
  bool operator()(const savvy::site_info& a, const savvy::site_info& b)
  {
    if (a.chromosome() == b.chromosome())
      return a.position() < b.position();
    return a.chromosome() < b.chromosome();
  }
};

struct right_comparator
{
  bool operator()(const savvy::site_info& a, const savvy::site_info& b)
  {
    if (a.chromosome() == b.chromosome())
      return (a.position() + std::max(a.ref().size(), a.alt().size())) < (b.position() + std::max(b.ref().size(), b.alt().size()));
    return a.chromosome() < b.chromosome();
  }
};

struct mid_comparator
{
  bool operator()(const savvy::site_info& a, const savvy::site_info& b)
  {
    if (a.chromosome() == b.chromosome())
    {
      double mid_a = static_cast<double>(a.position()) + (static_cast<double>(std::max(a.ref().size(), a.alt().size())) / 2.0);
      double mid_b = static_cast<double>(b.position()) + (static_cast<double>(std::max(b.ref().size(), b.alt().size())) / 2.0);

      return mid_a < mid_b;
    }

    return a.chromosome() < b.chromosome();
  }
};

class random_string_generator
{
public:
  random_string_generator() :
    rg_(std::random_device{}()),
    dist_(0, char_array_.size() - 1)
  {
  }
  std::string operator()(std::size_t length)
  {
    std::string ret;
    ret.reserve(length);

    while (length--)
      ret += char_array_[dist_(rg_)];

    return ret;
  }
private:
  std::vector<char> char_array_ = {'0','1','2','3','4','5','6','7','8','9',
                                   'A','B','C','D','E','F','G','H','I','J',
                                   'K','L','M','N','O','P','Q','R','S','T',
                                   'U','V','W','X','Y','Z'};
  std::mt19937 rg_;
  std::uniform_int_distribution<std::size_t> dist_;
};

class headless_writer : public savvy::sav::writer
{
public:
  template<typename RandAccessStringIterator, typename RandAccessKVPIterator>
  headless_writer(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, options opts = options())
    :
    savvy::sav::writer("", samples_beg, samples_end, headers_beg, headers_end, opts)
  {
    output_buf_ = (opts.compression_level > 0 ? std::unique_ptr<std::streambuf>(new shrinkwrap::zstd::obuf(file_path, opts.compression_level)) : std::unique_ptr<std::streambuf>(savvy::sav::writer::create_std_filebuf(file_path, std::ios::binary | std::ios::out)));
    output_stream_.rdbuf(output_buf_.get());
    file_path_ = file_path;
  }

  const std::string& file_path() const { return this->file_path_; }
};

class headless_reader : public savvy::sav::reader
{
public:
  template<typename RandAccessStringIterator, typename RandAccessKVPIterator>
  headless_reader(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, savvy::fmt format) :
    savvy::sav::reader("", format)
  {
    file_data_format_ = format;
    file_path_ = file_path;
    input_stream_ = shrinkwrap::zstd::istream(file_path);
    for (auto it = headers_beg; it != headers_end; ++it)
    {
      if (it->first == "INFO")
      {
        std::string info_field = savvy::parse_header_id(it->second);
        metadata_fields_.push_back(std::move(info_field));
      }
      else if (it->first == "FORMAT")
      {
        std::string format_field = savvy::parse_header_id(it->second);
        if (format_field == "GT")
          file_data_format_ = savvy::fmt::allele;
//                    else if (format_field == "GP")
//                      file_data_format_ = fmt::genotype_probability;
        else if (format_field == "HDS")
          file_data_format_ = savvy::fmt::haplotype_dosage;
      }
      headers_.emplace_back(it->first, it->second);
    }

    this->sample_ids_.assign(samples_beg, samples_end);

  }
};


template <typename Cmp>
int run_sort(const sort_prog_args& args)
{
  Cmp cmp;

  random_string_generator str_gen;
  std::size_t temp_file_size = 7;
  savvy::sav::reader r(args.input_path());
  std::deque<headless_writer> temp_writers;
  std::deque<headless_reader> temp_readers;

  savvy::sav::writer::options write_opts;
  write_opts.data_format = r.data_format();
  savvy::sav::writer output(args.output_path(), r.samples_begin(), r.samples_end(), r.headers().begin(), r.headers().end(), write_opts);

  {
    std::vector<var_type> variants(temp_file_size);

    std::size_t read_counter = temp_file_size;
    while (read_counter == temp_file_size)
    {
      read_counter = 0;
      for (; read_counter < temp_file_size; ++read_counter)
      {
        r >> variants[read_counter];
        if (!r.good())
          break;
      }

      if (read_counter)
      {
        std::vector<std::reference_wrapper<var_type>> variant_refs(variants.begin(), variants.end());
        std::sort(variant_refs.begin(), variant_refs.begin() + read_counter, cmp);
        std::string temp_path = "/tmp/tmp-" + str_gen(8) + ".sav";
        temp_writers.emplace_back(temp_path, r.samples_begin(), r.samples_end(), r.headers().begin(), r.headers().end(), write_opts);
        temp_readers.emplace_back(temp_path, r.samples_begin(), r.samples_end(), r.headers().begin(), r.headers().end(), write_opts.data_format);
        std::remove(temp_path.c_str());
        for (std::size_t i = 0; i < read_counter; ++i)
        {
          temp_writers.back() << variant_refs[i].get();
        }
      }
    }
  }



  {
    temp_writers.clear();
    std::vector<savvy::variant<savvy::compressed_vector<float>>> variants(temp_readers.size());
    std::size_t i = 0;
    for (auto it = temp_readers.begin(); it != temp_readers.end(); ++it)
    {
      it->read(variants[i], variants[i].data());
      ++i;
    }


    std::size_t min_index;

    do
    {
      min_index = variants.size();
      for (std::size_t rdr_index = 0; rdr_index < variants.size(); ++rdr_index)
      {
        if (temp_readers[rdr_index].good())
        {
          if (min_index < variants.size())
          {
            if (cmp(variants[rdr_index], variants[min_index]))
              min_index = rdr_index;
          }
          else
          {
            min_index = rdr_index;
          }
        }
      }

      if (min_index < variants.size())
      {
        output << variants[min_index];
        temp_readers[min_index] >> variants[min_index];
      }

    } while (min_index < variants.size());

  }

  return EXIT_FAILURE;
}


int sort_main(int argc, char** argv)
{
  sort_prog_args args;
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



  return run_sort<left_comparator>(args);
}