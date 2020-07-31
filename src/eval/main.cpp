/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/reader.hpp"
#include "savvy/bcf_writer.hpp"

#include <htslib/vcf.h>

#include <chrono>
#include <getopt.h>
#include <sys/stat.h>

class sav_eval_prog_args
{
private:
  std::vector<option> long_options_;
  std::list<std::string> input_paths_;
  std::string fmt_field_ = "GT";
  bool generate_ = false;
  bool help_ = false;
public:
  sav_eval_prog_args() :
    long_options_(
      {
        {"format", no_argument, 0, 'f'},
        {"generate", no_argument, 0, 'g'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::list<std::string>& input_paths() const { return input_paths_; }
  std::string input_path() const { return input_paths_.size() ? input_paths_.front() : ""; }
  const std::string& fmt_field() const { return fmt_field_; }
  bool run_generate() const { return generate_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav-eval [opts ...] <in.{bcf,sav,vcf.gz}> [ <in2.{bcf,sav,vcf.gz}> ...] \n";
    os << "Or: sav-eval --generate [opts ...] <in.{bcf,sav,sav2,vcf.gz}> \n";
    os << "\n";
    os << " -f, --format    Format field to use with --generate (default: GT)\n";
    os << " -g, --generate  Generate test files instead of running evaluation\n";
    os << " -h, --help      Print usage\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "f:gh", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'f':
        fmt_field_ = optarg ? optarg : "";
        break;
      case 'g':
        generate_ = true;
        break;
      case 'h':
        help_ = true;
        return true;
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
    else if (generate_ && remaining_arg_count > 1)
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    for (int i = 0; i < remaining_arg_count; ++i)
    {
      input_paths_.emplace_back(argv[optind + i]);
    }

    return true;
  }
};

template<typename T>
std::uint32_t adler32(T beg_it, T end_it, std::uint32_t cs = 1)
{
  uint32_t r0 = cs, r1 = 0;

  for (T it = beg_it; it != end_it; ++it)
  {
    r0 = (r0 + *it) % 65521;
    r1 = (r1 + r0) % 65521;
  }

  return (r1 << 16) | r0;
}

std::string get_prefix(const std::string& file_path)
{
  std::string prefix = file_path;
  if (savvy::detail::has_extension(prefix, ".sav") || savvy::detail::has_extension(prefix, ".bcf") || savvy::detail::has_extension(prefix, ".vcf"))
    prefix.erase(prefix.end() - 4, prefix.end());
  else if (savvy::detail::has_extension(prefix, ".vcf.gz"))
    prefix.erase(prefix.end() - 7);
  else
    return "";
  return prefix;
}

template <typename Itr>
int eval_gt(Itr paths_beg, Itr paths_end)
{
  using namespace std::chrono;

  std::size_t n_iterations = 3;

  for (auto it = paths_beg; it != paths_end; ++it)
  {
    std::string input_path = *it;


    if (savvy::detail::has_extension(input_path, ".bcf") || savvy::detail::has_extension(input_path, ".vcf.gz") || savvy::detail::has_extension(input_path, ".vcf"))
    {
      struct stat st;
      std::cout << (stat(input_path.c_str(), &st) == 0 ? st.st_size : -1);

      std::vector<std::int64_t> read_times;
      read_times.reserve(n_iterations);
      for (std::size_t i = 0; i <= n_iterations; ++i)
      {
        auto start = steady_clock::now();
        htsFile* fp = bcf_open(input_path.c_str(), "r");
        bcf_hdr_t *h = bcf_hdr_read(fp);
        bcf1_t *b = bcf_init1();
        int* geno{};
        int geno_sz{};

        std::size_t cnt = 0;
        int r;

        if (i == 0)
        {
          std::uint32_t cs = 1;
          while ((r=bcf_read(fp, h, b)) >= 0)
          {
            bcf_get_genotypes(h, b, &geno, &geno_sz);
            for (std::size_t j = 0; j < geno_sz; ++j)
            {
              if (geno[j] == bcf_int32_vector_end ) continue;

              //int is_phased = bcf_gt_is_phased(ptr[j]);

              geno[j] = (geno[j] >> 1) - 1; //bcf_gt_allele(geno[j]);
              if (geno[j] == -1)
                geno[j] = bcf_int32_missing;
            }
            cs = adler32(geno, geno + geno_sz, cs);
            ++cnt;
          }
          std::cout << "\t" << cs;
        }
        else
        {
          while ((r=bcf_read(fp, h, b)) >= 0)
          {
            bcf_get_genotypes(h, b, &geno, &geno_sz);
            for (std::size_t j = 0; j < geno_sz; ++j)
            {
              if (geno[j] == bcf_int32_vector_end ) continue;

              //int is_phased = bcf_gt_is_phased(ptr[j]);

              geno[j] = (geno[j] >> 1) - 1; //bcf_gt_allele(geno[j]);
              if (geno[j] == -1)
                geno[j] = bcf_int32_missing;
            }
            ++cnt;
          }

          read_times.emplace_back(duration_cast<milliseconds>(steady_clock::now() - start).count());
          std::cout << "\t" << read_times.back();
        }

        if (r < -1)
        {
          std::cerr << "Error reading hts file" << std::endl;
          exit(-1);
        }

        bcf_destroy1(b);
        bcf_hdr_destroy(h);
        hts_close(fp);
        free(geno);
      }

      float avg = std::accumulate(read_times.begin(), read_times.end(), 0.) / read_times.size();
      std::cout << "\t" << avg << "\t" << input_path << "\thts" << std::endl;
    }
//    else if (savvy::detail::has_extension(input_path, ".sav"))
//    {
//    }
//    else
    {
      struct stat st;
      std::cout << (stat(input_path.c_str(), &st) == 0 ? st.st_size : -1);

      std::vector<std::int64_t> read_times;
      read_times.reserve(n_iterations);
      for (std::size_t i = 0; i <= n_iterations; ++i)
      {
        auto start = steady_clock::now();
        savvy::v2::variant var;
        savvy::v2::reader rdr(input_path);
        std::size_t cnt = 0;
        std::vector<int> geno;

        if (i == 0)
        {
          std::uint32_t cs = 1;
          while (rdr.read(var))
          {
            var.get_format("GT", geno);
            cs = adler32(geno.begin(), geno.end(), cs);
            ++cnt;
          }
          std::cout << "\t" << cs;
        }
        else
        {

          while (rdr.read(var))
          {
            var.get_format("GT", geno);
            ++cnt;
          }

          read_times.emplace_back(duration_cast<milliseconds>(steady_clock::now() - start).count());
          std::cout << "\t" << read_times.back();
        }
      }

      float avg = std::accumulate(read_times.begin(), read_times.end(), 0.) / read_times.size();
      std::cout << "\t" << avg << "\t" << input_path << "\tsavvy" << std::endl;
    }
  }

  return EXIT_SUCCESS;
}

//void init_headers(std::vector<std::pair<std::string,std::string>>& destination, bcf_hdr_t* hdr)
//{
//  if (hdr)
//  {
//    destination.reserve(std::size_t(hdr->nhrec - 1));
//    for (int i = 0; i < hdr->nhrec; ++i)
//    {
//      std::string key, val;
//      if (hdr->hrec[i]->key && hdr->hrec[i]->value)
//      {
//        key = hdr->hrec[i]->key;
//        val = hdr->hrec[i]->value;
//      }
//      else if (hdr->hrec[i]->key && hdr->hrec[i]->nkeys) // (hdr->hrec[i]->type == BCF_HL_INFO || hdr->hrec[i]->type == BCF_HL_FLT || hdr->hrec[i]->type == BCF_HL_STR))
//      {
//        bcf_hrec_t* r = hdr->hrec[i];
//        key = r->key;
//        std::stringstream ss_val;
//
//        ss_val << "<";
//        for (int j = 0; j < r->nkeys - 1; ++j) // minus 1 to remove IDX;
//        {
//          if (j > 0)
//            ss_val << ",";
//          if (r->keys[j])
//            ss_val << r->keys[j];
//          ss_val << "=";
//          if (r->vals[j])
//            ss_val << r->vals[j];
//        }
//        ss_val << ">";
//        val = ss_val.str();
//      }
//
//      if (key.size())
//        destination.emplace_back(std::move(key), std::move(val));
//      //ret.insert(std::upper_bound(ret.begin(), ret.end(), std::make_pair(key, std::string()), [](const auto& a, const auto& b) { return a.first < b.first; }), {std::move(key), std::move(val)});
//    }
//  }
//}

savvy::v2::site_info get_site_info(std::size_t allele_index, bcf_hdr_t* hdr, bcf1_t* rec)
{
  bcf_unpack(rec, BCF_UN_ALL);
  std::size_t n_info = rec->n_info;
  std::size_t n_flt = rec->d.n_flt;
  bcf_info_t* info = rec->d.info;
  std::unordered_map<std::string, std::string> props;
  props.reserve(n_info + 2);

  if (std::isnan(rec->qual))
  {
    props["QUAL"] = ".";
  }
  else
  {
    std::string qual(std::to_string(rec->qual));
    qual.erase(qual.find_last_not_of('0') + 1); // rtrim zeros.
    qual.erase(qual.find_last_not_of('.') + 1);
    props["QUAL"] = std::move(qual);
  }

  std::stringstream ss;
  for (std::size_t i = 0; i < n_flt; ++i)
  {
    if (i > 0)
      ss << ";";
    ss << bcf_hdr_int2id(hdr, BCF_DT_ID, rec->d.flt[i]);
  }
  std::string fltr(ss.str());
  if (fltr == ".")
    fltr.clear();
  props["FILTER"] = std::move(fltr);
  props["ID"] = rec->d.id;


  for (std::size_t i = 0; i < n_info; ++i)
  {
    // bcf_hdr_t::id[BCF_DT_ID][$key].key
    const char* key = hdr->id[BCF_DT_ID][info[i].key].key;
    if (key)
    {
      switch (info[i].type)
      {
      case BCF_BT_NULL:
        props[key] = "1"; // Flag present so should be true.
        break;
      case BCF_BT_INT8:
      case BCF_BT_INT16:
      case BCF_BT_INT32:
        props[key] = std::to_string(info[i].v1.i);
        break;
      case BCF_BT_FLOAT:
        props[key] = std::to_string(info[i].v1.f);
        props[key].erase(props[key].find_last_not_of('0') + 1); // rtrim zeros.
        props[key].erase(props[key].find_last_not_of('.') + 1);
        break;
      case BCF_BT_CHAR:
        props[key] = std::string((char*)info[i].vptr, info[i].vptr_len);
        break;
      }
    }
  }

  return savvy::v2::site_info(
    std::string(bcf_hdr_id2name(hdr, rec->rid)),
    static_cast<std::uint64_t>(rec->pos + 1),
    std::string(rec->d.allele[0]),
    {std::string(rec->n_allele > 1 ? rec->d.allele[allele_index] : "")},
    "",
    rec->qual);
}
//
//int generate_sav2_dense(const std:: string& input_path, const std::string& fmt_field)
//{
//  std::string prefix_path = get_prefix(input_path);
//
//  std::size_t failed_checksum_count = 0;
//
//  std::size_t cl = 19;
//
//  for (std::size_t i = 0; i < 2; ++i)
//  {
//    std::string sav2_file_path;
//    if (i)
//      sav2_file_path = prefix_path + ".pbwt" + (cl == 19 ? ".c19.sav2" : ".sav2");
//    else
//      sav2_file_path = prefix_path + (cl == 19 ? ".c19.sav2" : ".sav2");
//
//    std::uint32_t orig_checksum = 1, sav2_checksum = 1;
//    std::size_t orig_cnt = 0, sav2_cnt = 0;
//
//    if (true)
//    {
//      htsFile* fp = bcf_open(input_path.c_str(), "r");
//      bcf_hdr_t *h = bcf_hdr_read(fp);
//      bcf1_t *b = bcf_init1();
//
//      int* row_buf{};
//      int row_buf_n{};
//
//      std::vector<std::pair<std::string,std::string>> hdrs;
//      init_headers(hdrs, h);
//      if (i)
//      {
//        hdrs.emplace_back("INFO", "<ID=_PBWT_SORT_" + fmt_field + ", Type=Flag, Format=\"" + fmt_field + "\">");
//        hdrs.emplace_back("INFO", "<ID=_PBWT_RESET, Type=Flag");
//      }
//      std::vector<std::string> samples(h->samples, h->samples + bcf_hdr_nsamples(h));
//
//      savvy::sav2::writer wrt(sav2_file_path, hdrs, samples, cl);
//      savvy::sav2::variant var;
//      std::vector<int> vec;
//
//      int r;
//      while ((r=bcf_read(fp, h, b)) >= 0)
//      {
//        ++orig_cnt;
//        if (r < -1)
//        {
//          std::cerr << "Error reading hts file" << std::endl;
//          return EXIT_FAILURE;
//        }
//
//        bcf_get_format_int32(h, b, fmt_field.c_str(), &row_buf, &row_buf_n);
//        vec.assign(row_buf, row_buf + row_buf_n);
//
//        orig_checksum = adler32(vec.data(), vec.data() + vec.size(), orig_checksum);
//
//        savvy::sav2::site_info* p = &var;
//        *p = get_site_info(1, h, b);
//        var.set_format(fmt_field, vec);
//
//        wrt.write_record(var);
//      }
//
//      bcf_destroy1(b);
//      bcf_hdr_destroy(h);
//      hts_close(fp);
//      free(row_buf);
//    }
//
//    {
//      std::vector<int> vec;
//      //std::vector<std::int16_t> vec;
//      savvy::sav2::variant var;
//      savvy::sav2::reader rdr(sav2_file_path);
//      //rdr.reset_bounds({"chr20", 500000, 500030});
//      while (rdr.read_record(var))
//      {
//        ++sav2_cnt;
//        var.get_format(fmt_field, vec);
//
//        sav2_checksum = adler32(vec.data(), vec.data() + vec.size(), sav2_checksum);
//      }
//
//    }
//
//    std::cout << orig_checksum << "\t" << sav2_checksum << "\t" << orig_cnt << "\t" << sav2_cnt << "\t" << sav2_file_path << std::endl;
//    if (orig_checksum != sav2_checksum)
//      ++failed_checksum_count;
//  }
//
//  return failed_checksum_count;
//}
//
//int generate_sav2_sparse(const std:: string& input_path)
//{
//  std::string prefix_path = get_prefix(input_path);
//  if (prefix_path.empty())
//    return EXIT_FAILURE;
//
//  std::size_t failed_checksum_count = 0;
//
//  for (float t : {-1.f, 0.0f, 0.0001f, 0.001f, 0.01f, 0.1f, 1.f})
//  {
//    std::uint32_t orig_checksum = 1, sav2_checksum = 1;
//    std::size_t orig_cnt{}, sav2_cnt{};
//    std::string sav2_file_path;
//    if (t < 1.f && t > 0.f)
//    {
//      auto str_t = std::to_string(t);
//      sav2_file_path = prefix_path + ".pbwt.t" + savvy::detail::rtrim(str_t, "0") + ".c19.sav2";
//    }
//    else if (t == 1.f)
//    {
//      sav2_file_path = prefix_path  + ".c19.sav2";
//    }
//    else if (t == 0.f)
//    {
//      sav2_file_path = prefix_path  + ".pbwt.c19.sav2";
//    }
//    else
//    {
//      sav2_file_path = prefix_path  + ".dense.c19.sav2";
//    }
//
//    if (true)
//    {
//      savvy::reader rdr(input_path, savvy::fmt::gt);
//
//      auto hdrs = rdr.headers();
//      hdrs.emplace_back("contig", "<ID=chr20>");
//      if (t < 1.f && t > -1.f)
//      {
//        hdrs.emplace_back("INFO", "<ID=_PBWT_SORT_GT, Type=Flag, Format=\"GT\">");
//        hdrs.emplace_back("INFO", "<ID=_PBWT_RESET, Type=Flag");
//      }
//      savvy::sav2::writer wrt(sav2_file_path, hdrs, rdr.samples(), 19);
//      savvy::compressed_vector<std::int16_t> vec;
//      std::vector<std::int16_t> dense_vec;
//      savvy::site_info site;
//      savvy::sav2::variant var;
//
//      while (rdr.read(site, vec))
//      {
//        ++orig_cnt;
//        savvy::sav2::site_info* p = &var;
//        *p = savvy::sav2::site_info(site.chromosome(), site.position(), site.ref(), {site.alt()}, site.prop("ID"), std::atof(site.prop("QUAL").c_str()));
//        for (auto it = vec.begin(); it != vec.end(); ++it)
//        {
//          if (!(*it))
//          {
//            *it = savvy::sav2::typed_value::missing_value<decltype(vec)::value_type>();
//          }
//        }
//
////        std::cout << var.pos();
////        for (auto it = vec.begin(); it != vec.end(); ++it)
////        {
////          std::cout << "\t" << it.offset() << ":"<< (*it);
////        }
////        std::cout << std::endl;
//
//        dense_vec.resize(0);
//        dense_vec.resize(vec.size());
//        for (auto it = vec.begin(); it != vec.end(); ++it)
//          dense_vec[it.offset()] = *it;
//
//        orig_checksum = adler32(vec.value_data(), vec.value_data() + vec.non_zero_size(), orig_checksum);
//        orig_checksum = adler32(vec.index_data(), vec.index_data() + vec.non_zero_size(), orig_checksum);
//
//        if ((float(vec.non_zero_size()) / float(vec.size())) >= t)
//          var.set_format("GT", dense_vec);
//        else
//          var.set_format("GT", vec);
//
//        wrt.write_record(var);
//      }
//
//      //std::cerr << cnt << std::endl;
//
//      if (rdr.bad() || !wrt)
//        return EXIT_FAILURE;
//    }
//
//    if (true)
//    {
//      savvy::compressed_vector<std::int16_t> vec;
//      //std::vector<std::int16_t> vec;
//      savvy::sav2::variant var;
//      savvy::sav2::reader rdr(sav2_file_path);
//      //rdr.reset_bounds({"chr20", 500000, 500030});
//      while (rdr.read_record(var))
//      {
//        ++sav2_cnt;
//        var.get_format("GT", vec);
//
////        std::cout << var.pos();
////        for (auto it = vec.begin(); it != vec.end(); ++it)
////        {
////          std::cout << "\t" << it.offset() << ":" << (*it);
////        }
////        std::cout << std::endl;
//
//        sav2_checksum = adler32(vec.value_data(), vec.value_data() + vec.non_zero_size(), sav2_checksum);
//        sav2_checksum = adler32(vec.index_data(), vec.index_data() + vec.non_zero_size(), sav2_checksum);
//      }
//    }
//
//    std::cout << orig_checksum << "\t" << sav2_checksum << "\t" << orig_cnt << "\t" << sav2_cnt << std::endl;
//    if (orig_checksum != sav2_checksum)
//      ++failed_checksum_count;
//  }
//
//  return failed_checksum_count;
//}

int main(int argc, char** argv)
{
  sav_eval_prog_args args;
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

//  if (args.run_generate())
//  {
//    if (args.fmt_field() == "GT")
//      return generate_sav2_sparse(args.input_path());
//    else
//      return generate_sav2_dense(args.input_path(), args.fmt_field());
//  }
  return eval_gt(args.input_paths().begin(), args.input_paths().end());
}

// sav stat-index ./test_file_hard.sav.s1r | grep "^marker count" | cut -f 2- | xargs echo | awk 't=0; {for(i=1;i<=NF;i++) t+=$i; print t}'
