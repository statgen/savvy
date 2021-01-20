/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/sort.hpp"
#include "savvy/reader.hpp"
#include "savvy/writer.hpp"

#include <getopt.h>
#include <chrono>

//================================================================//
less_than_comparator::less_than_comparator(savvy::s1r::sort_point type, std::unordered_map<std::string, std::size_t> contig_order_map) :
  contig_order_map_(std::move(contig_order_map)),
  sort_type_(type)
{
}

bool less_than_comparator::operator()(const savvy::site_info& a, const savvy::site_info& b) const
{
  switch (sort_type_)
  {
  case savvy::s1r::sort_point::mid: return mid(a, b);
  case savvy::s1r::sort_point::beg: return left(a, b);
  default: return right(a, b);
  }
}

bool less_than_comparator::contig_compare(const std::string& a, const std::string& b) const
{
  auto a_res = contig_order_map_.find(a);
  auto b_res = contig_order_map_.find(b);

  if (a_res != contig_order_map_.end() && b_res != contig_order_map_.end())
  {
    return (a_res->second < b_res->second);
  }
  else if (a_res == contig_order_map_.end() && b_res == contig_order_map_.end())
  {
    return a < b;
  }
  else
  {
    return a_res != contig_order_map_.end(); // put a before b if b is not in headers, or put b before a if a is not in headers.
  }
}

bool less_than_comparator::left(const savvy::site_info& a, const savvy::site_info& b) const
{
  if (a.chromosome() == b.chromosome())
    return a.position() < b.position();
  return contig_compare(a.chromosome(), b.chromosome());
}

bool less_than_comparator::right(const savvy::site_info& a, const savvy::site_info& b) const
{
  if (a.chromosome() == b.chromosome())
  {
    std::size_t a_max = a.ref().size();
    for (auto it = a.alts().begin(); it != a.alts().end(); ++it)
      a_max = std::max(a_max, it->size());

    std::size_t b_max = b.ref().size();
    for (auto it = b.alts().begin(); it != b.alts().end(); ++it)
      b_max = std::max(b_max, it->size());

    return (a.position() + a_max) < (b.position() + b_max);
  }
  return contig_compare(a.chromosome(), b.chromosome());
}

bool less_than_comparator::mid(const savvy::site_info& a, const savvy::site_info& b) const
{
  if (a.chromosome() == b.chromosome())
  {
    std::size_t a_max = a.ref().size();
    for (auto it = a.alts().begin(); it != a.alts().end(); ++it)
      a_max = std::max(a_max, it->size());

    std::size_t b_max = b.ref().size();
    for (auto it = b.alts().begin(); it != b.alts().end(); ++it)
      b_max = std::max(b_max, it->size());

    double mid_a = static_cast<double>(a.position()) + (static_cast<double>(a_max) / 2.0);
    double mid_b = static_cast<double>(b.position()) + (static_cast<double>(b_max) / 2.0);

    return mid_a < mid_b;
  }

  return contig_compare(a.chromosome(), b.chromosome());
}
//================================================================//

//================================================================//
greater_than_comparator::greater_than_comparator(savvy::s1r::sort_point type, std::unordered_map<std::string, std::size_t> contig_order_map) :
  contig_order_map_(std::move(contig_order_map)),
  sort_type_(type)
{
}

bool greater_than_comparator::operator()(const savvy::site_info& a, const savvy::site_info& b) const
{
  switch (sort_type_)
  {
  case savvy::s1r::sort_point::mid: return mid(a, b);
  case savvy::s1r::sort_point::beg: return left(a, b);
  default: return right(a, b);
  }
}

bool greater_than_comparator::contig_compare(const std::string& a, const std::string& b) const
{
  auto a_res = contig_order_map_.find(a);
  auto b_res = contig_order_map_.find(b);

  if (a_res != contig_order_map_.end() && b_res != contig_order_map_.end())
  {
    return (a_res->second > b_res->second);
  }
  else if (a_res == contig_order_map_.end() && b_res == contig_order_map_.end())
  {
    return a > b;
  }
  else
  {
    return b_res != contig_order_map_.end(); // put b before a if a is not in headers, or put a before b if b is not in headers.
  }
}

bool greater_than_comparator::left(const savvy::site_info& a, const savvy::site_info& b) const
{
  if (a.chromosome() == b.chromosome())
    return a.position() > b.position();
  return contig_compare(a.chromosome(), b.chromosome());
}

bool greater_than_comparator::right(const savvy::site_info& a, const savvy::site_info& b) const
{
  if (a.chromosome() == b.chromosome())
  {
    std::size_t a_max = a.ref().size();
    for (auto it = a.alts().begin(); it != a.alts().end(); ++it)
      a_max = std::max(a_max, it->size());

    std::size_t b_max = b.ref().size();
    for (auto it = b.alts().begin(); it != b.alts().end(); ++it)
      b_max = std::max(b_max, it->size());

    return (a.position() + a_max) > (b.position() + b_max);
  }
  return contig_compare(a.chromosome(), b.chromosome());
}

bool greater_than_comparator::mid(const savvy::site_info& a, const savvy::site_info& b) const
{
  if (a.chromosome() == b.chromosome())
  {
    std::size_t a_max = a.ref().size();
    for (auto it = a.alts().begin(); it != a.alts().end(); ++it)
      a_max = std::max(a_max, it->size());

    std::size_t b_max = b.ref().size();
    for (auto it = b.alts().begin(); it != b.alts().end(); ++it)
      b_max = std::max(b_max, it->size());

    double mid_a = static_cast<double>(a.position()) + (static_cast<double>(a_max) / 2.0);
    double mid_b = static_cast<double>(b.position()) + (static_cast<double>(b_max) / 2.0);

    return mid_a > mid_b;
  }

  return contig_compare(a.chromosome(), b.chromosome());
}
//================================================================//

//================================================================//
random_string_generator::random_string_generator() :
  rg_(std::random_device{}()^std::chrono::high_resolution_clock().now().time_since_epoch().count()),
  dist_(0, char_array_.size() - 1)
{
}

std::string random_string_generator::operator()(std::size_t length)
{
  std::string ret;
  ret.reserve(length);

  while (length--)
    ret += char_array_[dist_(rg_)];

  return ret;
}
//================================================================//


class sort_prog_args
{
private:
  std::vector<option> long_options_;
  std::string input_path_;
  std::string output_path_ = "/dev/stdout";
  std::string direction_ = "asc";
  savvy::s1r::sort_point point_ = savvy::s1r::sort_point::beg;
  bool help_ = false;
public:
  sort_prog_args() :
    long_options_(
      {
        {"direction", required_argument, 0, 'd'},
        {"help", no_argument, 0, 'h'},
        {"output", required_argument, 0, 'o'},
        {"point", required_argument, 0, 'p'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& direction() const { return direction_; }
  const std::string& input_path() const { return input_path_; }
  const std::string& output_path() const { return output_path_; }
  savvy::s1r::sort_point point() const { return point_; }

  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav sort [opts ...] <in.sav> \n";
    os << "\n";
    os << " -d, --direction   Specifies whether to sort in ascending or descending order (asc or desc; default: asc)\n";
    os << " -h, --help        Print usage\n";
    os << " -o, --output      Path to output SAV file (default: /dev/stdout).\n";
    os << " -p, --point       Specifies which allele position to sort by (beg, mid or end; default: beg)\n";

    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "d:ho:p:", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'd':
        direction_ = std::string(optarg ? optarg : "");
        if (direction_ != "asc" && direction_ != "desc")
        {
          std::cerr << "Invalid --direction argument (" << direction_ << ")." << std::endl;
          return false;
        }
        break;
      case 'h':
        help_ = true;
        return true;
      case 'o':
        output_path_ = std::string(optarg ? optarg : "");
        break;
      case 'p':
      {
        std::string sort_str = std::string(optarg ? optarg : "");
        if (sort_str.front()=='b')
        {
          point_ = savvy::s1r::sort_point::beg;
        }
        else if (sort_str.front()=='e')
        {
          point_ = savvy::s1r::sort_point::end;
        }
        else if (sort_str.front()=='m')
        {
          point_ = savvy::s1r::sort_point::mid;
        }
        else
        {
          std::cerr << "Invalid --sort-point argument (" << sort_str << ")." << std::endl;
          return false;
        }
        break;
      }
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count < 1)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else if (remaining_arg_count > 1)
    {
      std::cerr << "Too many arguments\n";
      return false;
    }
    else
    {
      input_path_ = argv[optind];
    }

    return true;
  }
};

template <typename SiteCompare>
int run(savvy::reader& in, savvy::writer& out, const SiteCompare& compare_fn)
{
  random_string_generator str_gen;
  const std::size_t temp_file_size = 8192;
  std::deque<std::string> temp_file_paths;
  std::deque<savvy::writer> temp_writers;
  std::deque<savvy::reader> temp_readers;


  std::vector<savvy::variant> in_mem_variants(temp_file_size);
  std::vector<std::reference_wrapper<savvy::variant>> in_mem_variant_refs(in_mem_variants.begin(), in_mem_variants.end());

  std::size_t read_counter = temp_file_size;
  while (read_counter == temp_file_size)
  {
    read_counter = 0;
    for (; read_counter < temp_file_size; )
    {
      in >> in_mem_variant_refs[read_counter].get();
      if (!in.good())
      {
        break;
      }
      else
      {
        ++read_counter;
      }
    }

    if (read_counter)
    {
      std::sort(in_mem_variant_refs.begin(), in_mem_variant_refs.begin() + read_counter, compare_fn);

      temp_file_paths.emplace_back("/tmp/tmp-" + str_gen(32) + ".sav");
      temp_writers.emplace_back(temp_file_paths.back(), savvy::file::format::sav2, in.headers(), in.samples());
      // TODO: open FILE*  and unlink file
      // temp_readers.emplace_back(temp_path, /*in.samples().begin(), in.samples().end(), in.headers().begin(), in.headers().end(),*/ out_format);
      // std::remove(temp_path.c_str());
      for (std::size_t i = 0; i<read_counter; ++i)
      {
        temp_writers.back() << in_mem_variant_refs[i].get();
      }
    }
  }



  temp_writers.clear();

  for (auto it = temp_file_paths.begin(); it != temp_file_paths.end(); ++it)
  {
    temp_readers.emplace_back(*it);
    std::remove(it->c_str());
  }

  std::vector<savvy::variant> write_variants(temp_readers.size());
  std::size_t i = 0;
  for (auto it = temp_readers.begin(); it != temp_readers.end(); ++it)
  {
    it->read(write_variants[i]);
    ++i;
  }


  std::size_t min_index;

  do
  {
    min_index = write_variants.size();
    for (std::size_t rdr_index = 0; rdr_index < write_variants.size(); ++rdr_index)
    {
      if (temp_readers[rdr_index].good())
      {
        if (min_index < write_variants.size())
        {
          if (compare_fn(write_variants[rdr_index], write_variants[min_index]))
            min_index = rdr_index;
        }
        else
        {
          min_index = rdr_index;
        }
      }
      else if (temp_readers[rdr_index].bad())
      {
        std::cerr << "Error: read failure with temp reader " << rdr_index << std::endl;
        return false;
      }
    }

    if (min_index < write_variants.size())
    {
      out << write_variants[min_index];
      temp_readers[min_index] >> write_variants[min_index];
    }

  } while (min_index < write_variants.size());

  return true;
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

  savvy::reader rdr(args.input_path());
  if (!rdr)
  {
    std::cerr << "Error: failed to open input file" << std::endl;
    return EXIT_FAILURE;
  }

  savvy::writer wtr(args.output_path(), savvy::file::format::sav2, rdr.headers(), rdr.samples());

  std::unordered_map<std::string, std::size_t> contig_order_map;
  contig_order_map.reserve(rdr.headers().size()); // std::count_if(in.headers().begin(), in.headers().end(), [](const std::pair<std::string,std::string>& e) { return e.first == "contig"; }));

  std::size_t contig_counter = 0;
  for (auto it = rdr.headers().begin(); it != rdr.headers().end(); ++it)
  {
    if (it->first == "contig")
    {
      std::string contig = savvy::parse_header_sub_field(it->second, "ID");
      contig_order_map[contig] = contig_counter++;
    }
  }


  if (args.direction() == "desc")
  {
    greater_than_comparator greater_than(args.point(), std::move(contig_order_map));
    return run(rdr, wtr, greater_than) ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  else
  {
    less_than_comparator less_than(args.point(), std::move(contig_order_map));
    return run(rdr, wtr, less_than) ? EXIT_SUCCESS : EXIT_FAILURE;
  }
}

