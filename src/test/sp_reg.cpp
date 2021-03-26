/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "tcdf.hpp"
#include "savvy/reader.hpp"

//#ifndef __cpp_lib_as_const
//#define __cpp_lib_as_const
//#endif
#include <xtensor/xarray.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <random>
#include <getopt.h>

std::vector<std::string> split_string_to_vector(const char* in, char delim)
{
  std::vector<std::string> ret;
  const char* d = nullptr;
  std::string token;
  const char* s = in;
  const char*const e = in + strlen(in);
  while ((d = std::find(s, e,  delim)) != e)
  {
    ret.emplace_back(std::string(s, d));
    s = d ? d + 1 : d;
  }
  ret.emplace_back(std::string(s,d));
  return ret;
}

savvy::genomic_region string_to_region(const std::string& s)
{
  const std::size_t colon_pos = s.find(':');
  if (colon_pos == std::string::npos)
  {
    return savvy::genomic_region(s);
  }
  else
  {
    std::string chr = s.substr(0, colon_pos);
    const std::size_t hyphen_pos = s.find('-', colon_pos + 1);
    if (hyphen_pos == std::string::npos)
    {
      std::string slocus = s.substr(colon_pos + 1);
      std::uint64_t ilocus = std::uint64_t(std::atoll(slocus.c_str()));
      return savvy::genomic_region(chr, ilocus, ilocus);
    }
    else
    {
      std::string sbeg = s.substr(colon_pos + 1, hyphen_pos - chr.size() - 1);
      std::string send = s.substr(hyphen_pos + 1);
      if (send.empty())
      {
        return savvy::genomic_region(chr, std::uint64_t(std::atoll(sbeg.c_str())));
      }
      else
      {
        return savvy::genomic_region(chr, std::uint64_t(std::atoll(sbeg.c_str())), std::uint64_t(std::atoll(send.c_str())));
      }
    }
  }

}

class prog_args
{
private:
  std::vector<option> long_options_;
  std::vector<std::string> covariate_fields_;
  std::string id_field_;
  std::string phenotype_field_;
  std::string geno_path_;
  std::string pheno_path_;
  std::string output_path_ = "/dev/stdout";
  std::string fmt_field_ = "";
  std::unique_ptr<savvy::genomic_region> region_;
  double min_mac_ = 1.0;
  bool no_sparse_ = false;
  bool always_sparse_ = false;
  bool logit_ = false;
  bool help_ = false;
public:
  prog_args() :
    long_options_(
      {
        {"cov", required_argument, 0, 'c'},
        {"fmt-field", required_argument, 0, '\x02'},
        {"help", no_argument, 0, 'h'},
        {"id", required_argument, 0, 'i'},
        {"logit", no_argument, 0, 'b'},
        {"min-mac", required_argument, 0, '\x02'},
        {"no-sparse", no_argument, 0, '\x01'},
        {"always-sparse", no_argument, 0, '\x01'},
        {"output", required_argument, 0, 'o'},
        {"pheno", required_argument, 0, 'p'},
        {"region", required_argument, 0, 'r'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::vector<std::string>& cov_columns() const { return covariate_fields_; }
  const std::string& id_column() const { return id_field_; }
  const std::string& pheno_column() const { return phenotype_field_; }
  const std::string& geno_path() const { return geno_path_; }
  const std::string& pheno_path() const { return pheno_path_; }
  const std::string& output_path() const { return output_path_; }
  const std::string& fmt_field() const { return fmt_field_; }
  const std::unique_ptr<savvy::genomic_region>& region() const { return region_; }
  double min_mac() const { return min_mac_; }
  bool sparse_disabled() { return no_sparse_; }
  bool force_sparse() { return always_sparse_; }
  bool logit_enabled() const { return logit_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sp-reg [opts ...] <geno_file> <pheno_file> \n";
    os << "\n";
    os << " -c, --cov            Comma separated list of covariate columns\n";
    os << " -h, --help           Print usage\n";
    os << " -i, --id             Sample ID column (defaults to first column)\n";
    os << " -b, --logit          Enable logistic model\n";
    os << " -o, --output         Output path (default: /dev/stdout)\n";
    os << " -p, --pheno          Phenotype column\n";
    os << " -r, --region         Genomic region to test (chrom:beg-end)\n";
    os << "     --min-mac        Minimum minor allele count (default: 1)\n";
    os << "     --no-sparse      Disables sparse optimizations\n";
    os << "     --always-sparse  Forces sparse optimizations even for dense file records\n";
    os << "     --fmt-field      Format field to use (DS, HDS, or GT)\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "\x01\x02:bhc:o:p:r:", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case '\x01':
        if (std::string("no-sparse") == long_options_[long_index].name)
        {
          no_sparse_ = true;
        }
        else if (std::string("always-sparse") == long_options_[long_index].name)
        {
          always_sparse_ = true;
        }
        else
        {
          return std::cerr << "Error: invalid option " << long_options_[long_index].name << std::endl, false;
        }
        break;
      case '\x02':
        if (std::string("min-mac") == long_options_[long_index].name)
        {
          min_mac_ = std::atof(optarg ? optarg : "");
        }
        else if (std::string("fmt-field") == long_options_[long_index].name)
        {
          fmt_field_ = optarg ? optarg : "";
          if (fmt_field_ != "DS" && fmt_field_ != "HDS" && fmt_field_ != "GT")
            return std::cerr << "Error: --fmt-field must be DS, HDS, or GT\n", false;
        }
        else
        {
          return std::cerr << "Error: invalid option " << long_options_[long_index].name << std::endl, false;
        }
        break;
      case 'b':
        logit_ = true;
        break;
      case 'h':
        help_ = true;
        return true;
      case 'c':
        covariate_fields_ = split_string_to_vector(optarg ? optarg : "", ',');
        break;
      case 'o':
        output_path_ = optarg ? optarg : "";
        break;
      case 'p':
        phenotype_field_ = optarg ? optarg : "";
        break;
      case 'r':
        region_.reset(new savvy::genomic_region(string_to_region(optarg ? optarg : "")));
        break;
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 2)
    {
      geno_path_ = argv[optind];
      pheno_path_ = argv[optind + 1];
      //phenotype_field_ = argv[optind + 2];
    }
    else if (remaining_arg_count < 2)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    return true;
  }
};

int test_xtensor()
{
  using namespace xt;
  using namespace xt::linalg;

  xt::xarray<double> a = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  xt::xtensor<double, 1> y = {0., 1., 1., 0., 0.};
  xt::xtensor<double, 2> x = {
    {0.05, 0.12},
    {0.18, 0.22},
    {0.31, 0.35},
    {0.42, 0.38},
    {0.5, 0.49}};

  std::cerr << y.shape()[0] << std::endl;
  std::cerr << xt::transpose(y).shape()[0] << std::endl;
  std::cerr << x << std::endl;
  x = xt::concatenate(xt::xtuple(xt::ones<double>({y.size(),std::size_t(1)}), x), 1);
  std::cerr << x << std::endl;
  std::cerr << xt::adapt(y.shape()) << std::endl;
  std::cerr << y << std::endl;
  std::cerr << xt::adapt(xt::col(x, 0).shape()) << std::endl;
  std::cerr << xt::col(x, 0) << std::endl;


  std::cerr << xt::linalg::dot(xt::transpose(x), x) << std::endl;


  //B' = (X!X)−1X!y
  //b = inv(X.T.dot(X)).dot(X.T).dot(y)
  // b = (x'x)-1 X
  xt::xarray<double> betas = dot(dot(inv(dot(transpose(x), x)), transpose(x)), y);
  xt::xarray<double> pbetas = dot(dot(pinv(dot(transpose(x), x)), transpose(x)), y);
  std::cerr << betas << std::endl;
  std::cerr << pbetas << std::endl;
  xt::xarray<double> residuals = y - dot(x, betas);
  std::cerr << residuals << std::endl;
  double sos_res = 0.;
  for (auto e : residuals)
    sos_res += e*e;

  double psos_res = 0.;
  for (auto e : (y - dot(x, pbetas)))
    psos_res += e*e;

  {
    //auto [solution, sos_residuals, rank, s] = lstsq(x, y); // xt::transpose(y));
    xarray<double> solution, sos_residuals, rank, s;
    std::tie(solution, sos_residuals, rank, s) = lstsq(x, y);
    double sos_res2 = sos_residuals[0];

    std::cerr << solution << std::endl;
    std::cerr << sos_residuals << std::endl;
    std::cerr << rank << std::endl;
    std::cerr << s << std::endl;
  }
  return 0;
}

typedef double scalar_type;

struct phenotype_file_data
{
  std::vector<std::string> ids;
  std::vector<scalar_type> resp_data;
  std::vector<std::vector<scalar_type>> cov_data;
};

auto parse_pheno_file(const prog_args& args, phenotype_file_data& dest)
{
  dest = phenotype_file_data();
//  xt::xtensor<scalar_type, 1> resp_data;
//  xt::xtensor<scalar_type, 2> cov_data;

  //std::size_t id_idx = args.id_column().empty() ? 0 : std::size_t(-1);
  //std::vector<std::size_t> covariate_idxs(args.cov_columns().size(), std::size_t(-1));

  std::ifstream pheno_file(args.pheno_path(), std::ios::binary);

  const std::size_t id_code = 1;
  const std::size_t resp_code = 2;
  const std::size_t cov_code = 3;

  std::string line;
  if (std::getline(pheno_file, line))
  {
    auto header_names = split_string_to_vector(line.c_str(), '\t');
    if (header_names.empty())
      return std::cerr << "Error: empty header\n", false;

    if (header_names[0].size() && header_names[0][0] == '#')
      header_names[0].erase(header_names[0].begin());

    std::vector<std::size_t> mask(header_names.size());
    if (args.id_column().empty())
    {
      std::size_t default_id_idx = args.pheno_path().rfind(".ped") == (args.pheno_path().size() - 4) ? 1 : 0;
      std::cerr << "Notice: using column " << (default_id_idx + 1) << " for sample ID column since --id not specified\n";
      mask[default_id_idx] = id_code;
    }

    for (std::size_t i = 0; i < header_names.size(); ++i)
    {
      if (header_names[i] == args.id_column())
      {
        mask[i] = id_code;
      }
      else if (header_names[i] == args.pheno_column())
      {
        mask[i] = resp_code;
      }
      else
      {
        for (std::size_t j = 0; j < args.cov_columns().size(); ++j)
        {
          if (header_names[i] == args.cov_columns()[j])
          {
            mask[i] = cov_code;
            break;
          }
        }
      }
    }

    if (std::count(mask.begin(), mask.end(), id_code) == 0)
      return std::cerr << "Error: missing identifier column\n", false; // TODO: better error message
    if (std::count(mask.begin(), mask.end(), resp_code) == 0)
      return std::cerr << "Error: missing response column\n", false; // TODO: better error message
    if (std::count(mask.begin(), mask.end(), cov_code) != args.cov_columns().size())
      return std::cerr << "Error: could not find all covariate columns\n", false; // TODO: better error message

    char* p = nullptr;
    while (std::getline(pheno_file, line))
    {
      auto str_fields = split_string_to_vector(line.c_str(), '\t');
      dest.cov_data.emplace_back(args.cov_columns().size(), std::numeric_limits<scalar_type>::quiet_NaN());
      dest.resp_data.emplace_back(std::numeric_limits<scalar_type>::quiet_NaN());
      std::size_t j = 0;
      for (std::size_t i = 0; i < str_fields.size(); ++i)
      {
        if (mask[i] == id_code)
        {
          dest.ids.emplace_back(std::move(str_fields[i]));
        }
        else if (mask[i] == resp_code)
        {
          scalar_type v = std::strtod(str_fields[i].c_str(), &p);
          if (p == str_fields[i].c_str() && !str_fields[i].empty() && str_fields[i][0] != '.' && std::tolower(str_fields[i][0]) != 'n')
            return std::cerr << "Error: encountered non-numeric phenotype\n", false;
          else
            dest.resp_data.back() = v;
        }
        else if (mask[i] == cov_code)
        {
          scalar_type v = std::strtod(str_fields[i].c_str(), &p);
          if (p == str_fields[i].c_str() && !str_fields[i].empty() && str_fields[i][0] != '.' && std::tolower(str_fields[i][0]) != 'n')
            return std::cerr << "Error: encountered non-numeric covariate\n", false;
          else
            dest.cov_data.back()[j] = v;
          ++j;
        }
      }
    }
  }

  return true;
}

//typedef xt::xtensor<scalar_type, 1> residuals_type;
typedef xt::xarray<scalar_type> residuals_type;

template <typename T, typename T2>
residuals_type compute_residuals(const T& y, const T2& x_orig)
{
  using namespace xt;
  using namespace xt::linalg;
  T2 x = concatenate(xtuple(xt::ones<scalar_type>({y.size(), std::size_t(1)}), x_orig), 1);
//  auto a = dot(transpose(x), x);
//  std::cerr << a << std::endl;
//  auto b = pinv(a);
//  std::cerr << "END A ------------------------------" << std::endl;
//  std::cerr << b << std::endl;
//  std::cerr << "END B ------------------------------" << std::endl;
//  auto c = dot(b, transpose(x));
//  std::cerr << c << std::endl;
  auto pbetas = dot(dot(pinv(dot(transpose(x), x)), transpose(x)), y);
  std::cerr << pbetas << std::endl;
  residuals_type residuals = y - dot(x, pbetas);
  std::cerr << "sum(y): " << sum(y) << std::endl;
  return residuals;
}

template <typename T, typename T2>
residuals_type compute_residuals_logit(const T& y, const T2& x_orig)
{
  using namespace xt;
  using namespace xt::linalg;
  const scalar_type epsilon = 0.00001;

  T2 x = concatenate(xtuple(xt::ones<scalar_type>({y.size(), std::size_t(1)}), x_orig), 1);

  T y2 = xt::maximum(epsilon, xt::minimum(y, 1. - epsilon));
  T div_y = xt::operator/(1., y2);
  T logit_y = -xt::log(1. / y2 - 1.);

//  auto a = dot(transpose(x), x);
//  std::cerr << a << std::endl;
//  auto b = pinv(a);
//  std::cerr << "END A ------------------------------" << std::endl;
//  std::cerr << b << std::endl;
//  std::cerr << "END B ------------------------------" << std::endl;
//  auto c = dot(b, transpose(x));
//  std::cerr << c << std::endl;
  auto pbetas = dot(dot(pinv(dot(transpose(x), x)), transpose(x)), logit_y);
  std::cerr << pbetas << std::endl;
  auto xw = dot(x, pbetas);
  T y_hat = 1. / (1. + xt::exp(-xw));
  residuals_type residuals = y - y_hat;
  scalar_type se = xt::mean(residuals * residuals)();

  std::cerr << "sum(y): " << sum(y) << std::endl;
  return residuals;
}

template <typename T>
T square(T v)
{
  return v * v;
}
#if 0
auto linreg_ttest_old(const std::vector<float>& y, const std::vector<float>& x)
{
  const std::size_t n = x.size();
  const float s_x     = std::accumulate(x.begin(), x.end(), 0.0f);
  const float s_y     = std::accumulate(y.begin(), y.end(), 0.0f);
  const float s_xx    = std::inner_product(x.begin(), x.end(), x.begin(), 0.0f);
  const float s_xy    = std::inner_product(x.begin(), x.end(), y.begin(), 0.0f);
  const float m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const float b       = (s_y - m * s_x) / n;
  auto fx             = [m,b](float x) { return m * x + b; };
  float se_line       = 0.0f; for (std::size_t i = 0; i < n; ++i) se_line += square(y[i] - fx(x[i]));
  const float x_mean  = s_x / n;
  float se_x_mean     = 0.0f; for (std::size_t i = 0; i < n; ++i) se_x_mean += square(x[i] - x_mean);
  const float dof     = n - 2;
  const float std_err = std::sqrt(se_line / dof) / std::sqrt(se_x_mean);
  float t = m / std_err;
  //boost::math::students_t_distribution<float> dist(dof);
  float pval =  tcdf(t, n - 1); //cdf(complement(dist, std::fabs(std::isnan(t) ? 0 : t))) * 2;
  return std::make_tuple(m, std_err, t, pval); // slope, std error, t statistic, p value
}
#endif
auto linreg_ttest(const std::vector<scalar_type>& y, const std::vector<scalar_type>& x, const scalar_type s_y)
{
  assert(y.size() == x.size());
  const std::size_t n = x.size();
  scalar_type s_x{}; //     = std::accumulate(x.begin(), x.end(), 0.0f);
  //scalar_type s_y{}; //     = std::accumulate(y.begin(), y.end(), 0.0f);
  scalar_type s_xx{}; //    = std::inner_product(x.begin(), x.end(), x.begin(), 0.0f);
  scalar_type s_xy{}; //    = std::inner_product(x.begin(), x.end(), y.begin(), 0.0f);

  for (std::size_t i = 0; i < n; ++i)
  {
    s_x += x[i];
    //s_y += y[i];
    s_xx += x[i] * x[i];
    s_xy += x[i] * y[i];
  }

  const scalar_type m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const scalar_type b       = (s_y - m * s_x) / n;
  auto fx             = [m,b](scalar_type x) { return m * x + b; };
  const scalar_type x_mean  = s_x / n;

  double se_line{};
  double se_x_mean{};
  for (std::size_t i = 0; i < n; ++i)
  {
    se_line += square(y[i] - fx(x[i]));
    se_x_mean += square(x[i] - x_mean);
  }

  const scalar_type dof     = n - 2;
  const scalar_type std_err = std::sqrt(se_line / dof) / std::sqrt(se_x_mean);
  scalar_type t = m / std_err;
  //boost::math::students_t_distribution<scalar-type> dist(dof);
  scalar_type pval =  tcdf(t, dof); //cdf(complement(dist, std::fabs(std::isnan(t) ? 0 : t))) * 2;

  return std::make_tuple(m, std_err, t, pval); // slope, std error, t statistic, p value
}
#if 0
auto sp_lin_reg_old(const std::vector<float>& y, const savvy::compressed_vector<float>& x)
{
  const std::size_t n = x.size();
  const float s_x     = std::accumulate(x.begin(), x.end(), 0.0f);
  const float s_y     = std::accumulate(y.begin(), y.end(), 0.0f);
  const float s_xx    = std::inner_product(x.begin(), x.end(), x.begin(), 0.0f);
  float s_xy    = 0.0f; for (auto it = x.begin(); it != x.end(); ++it) s_xy += (*it * y[x.index_data()[it - x.begin()]]);
  const float m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const float b       = (s_y - m * s_x) / n;
  auto fx             = [m,b](float x) { return m * x + b; };
  float se_line       = 0.0f; for (std::size_t i = 0; i < n; ++i) se_line += square(y[i] - fx(x[i]));
  const float x_mean  = s_x / n;
  float se_x_mean     = 0.0f; for (std::size_t i = 0; i < n; ++i) se_x_mean += square(x[i] - x_mean);
  const float dof     = n - 2;
  const float std_err = std::sqrt(se_line / dof) / std::sqrt(se_x_mean);
  float t = m / std_err;
  //std::students_t_distribution<float> dist(dof);
  float pval = tcdf(t, n - 1); //cdf(complement(dist, std::fabs(t))) * 2;
  /*
  beta = ((c+1)*sxy-sx*sy)/((c+1)*sxx-sx*sx);
  varE = 1/(c+1.)/(c-1.)*((c+1)*syy-sy*sy-beta*beta*((c+1)*sxx-sx*sx));
  sebeta = sqrt((c+1)*varE/((c+1)*sxx-sx*sx));
  r = ((c+1)*sxy-sx*sy)/sqrt(((c+1)*sxx-sx*sx)*((c+1)*syy-sy*sy));
  t = r * sqrt((c-1)/(1-r*r+pEmmaxHelper::ZEPS));
  pval = pEmmaxHelper::tcdf(t, c-1);
  */
  return std::make_tuple(m, std_err, t, pval); // slope, std error, t statistic, p value
}
#endif

auto linreg_ttest(const std::vector<scalar_type>& y, const savvy::compressed_vector<scalar_type>& x, const scalar_type& s_y, const scalar_type& s_yy)
{
  assert(y.size() == x.size());
  const std::size_t n = x.size();
  scalar_type s_x{}; //     = std::accumulate(x.begin(), x.end(), 0.0f);
  scalar_type s_xx{}; //    = std::inner_product(x.begin(), x.end(), x.begin(), 0.0f);
  scalar_type s_xy{}; //    = std::inner_product(x.begin(), x.end(), y.begin(), 0.0f);

  const auto x_beg = x.begin();
  const auto x_end = x.end();
//  const scalar_type* x_values = x.value_data();
//  const std::size_t* x_indices = x.index_data();
  for (auto it = x_beg; it != x_end; ++it)
  {
    s_x += *it;
    s_xx += (*it) * (*it);
    s_xy += (*it) * y[it.offset()];
  }

  //const float s_y     = std::accumulate(y.begin(), y.end(), 0.0f);
  const scalar_type m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const scalar_type x_mean  = s_x / n;

  float se_x_mean{};

  if (false)
  {
    const scalar_type b       = (s_y - m * s_x) / n;
    auto fx             = [m,b](scalar_type x) { return m * x + b; };
    const scalar_type f_of_zero = fx(0.0f);

    float se_line{};
    std::size_t i = 0;
    for (auto it = x.begin(); it != x.end(); ++i)
    {
      if (i == it.offset())
      {
        se_line += square(y[i] - fx(*it));
        se_x_mean += square(*it - x_mean);
        ++it;
      }
      else
      {
        se_line += square(y[i] - f_of_zero);
      }
    }

    for ( ; i < n; ++i)
    {
      se_line += square(y[i] - f_of_zero);
    }

    se_x_mean += (square(0.0f - x_mean) * scalar_type(n - x.non_zero_size()));

    const scalar_type dof = n - 2;
    const scalar_type std_err = std::sqrt(se_line / dof) / std::sqrt(se_x_mean);
    scalar_type t = m / std_err;
    //std::students_t_distribution<float> dist(dof);
    scalar_type pval = tcdf(t, dof); //cdf(complement(dist, std::fabs(std::isnan(t) ? 0 : t))) * 2;
    /*
    beta = ((c+1)*sxy-sx*sy)/((c+1)*sxx-sx*sx);
    varE = 1/(c+1.)/(c-1.)*((c+1)*syy-sy*sy-beta*beta*((c+1)*sxx-sx*sx));
    sebeta = sqrt((c+1)*varE/((c+1)*sxx-sx*sx));
    r = ((c+1)*sxy-sx*sy)/sqrt(((c+1)*sxx-sx*sx)*((c+1)*syy-sy*sy));
    t = r * sqrt((c-1)/(1-r*r+pEmmaxHelper::ZEPS));
    pval = pEmmaxHelper::tcdf(t, c-1);
    */

    return std::make_tuple(m, std_err, t, pval); // slope, std error, t statistic, p value
  }
  else
  {
    for (auto it = x.begin(); it != x.end(); ++it)
    {
      se_x_mean += square(*it - x_mean);
    }

    se_x_mean += (square(0.0f - x_mean) * scalar_type(n - x.non_zero_size()));
    double se2 = 1./(n*(n-2)) * (n*s_yy - s_y*s_y - square(m)*(n*s_xx - square(s_x)));

    const scalar_type dof = n - 2;
    const scalar_type std_err = std::sqrt(se2) / std::sqrt(se_x_mean);
    scalar_type t = m / std_err;
    //std::students_t_distribution<float> dist(dof);
    scalar_type pval = tcdf(t, dof); //cdf(complement(dist, std::fabs(std::isnan(t) ? 0 : t))) * 2;

    return std::make_tuple(m, std_err, t, pval); // slope, std error, t statistic, p value
  }
}

//void lin_reg(const residuals_type& y, const std::vector<scalar_type>& x)
//{
//}
//
//void lin_reg(const residuals_type& y, const savvy::compressed_vector<scalar_type>& x)
//{
//}

bool load_phenotypes(const prog_args& args, savvy::reader& geno_file, xt::xtensor<scalar_type, 1>& pheno_vec, xt::xtensor<scalar_type, 2>& cov_mat)
{
  phenotype_file_data full_pheno;
  if (!parse_pheno_file(args, full_pheno))
    return false;



  std::unordered_set<std::string> samples_with_phenotypes;
  std::unordered_map<std::string, std::size_t> id_map;
  samples_with_phenotypes.reserve(full_pheno.ids.size());
  id_map.reserve(full_pheno.ids.size());
  for (std::size_t i = 0; i < full_pheno.resp_data.size(); ++i)
  {
    if (std::isnan(full_pheno.resp_data[i]) /*|| std::find_if(covariate_data[i].begin(); covariate_data[i].end(), std::isnan) != covariate_data[i].end())*/)
    {
      // missing

    }
    else
    {
      id_map[full_pheno.ids[i]] = i;
      samples_with_phenotypes.emplace(full_pheno.ids[i]);
    }
  }

  auto sample_intersection = geno_file.subset_samples(samples_with_phenotypes);

  pheno_vec = xt::xtensor<scalar_type, 1>::from_shape({sample_intersection.size()});
  cov_mat = xt::xtensor<scalar_type, 2>::from_shape({sample_intersection.size(), args.cov_columns().size()});
  for (std::size_t i = 0; i < sample_intersection.size(); ++i)
  {
    std::size_t src_idx = id_map[sample_intersection[i]];
    pheno_vec(i) = full_pheno.resp_data[src_idx];
    for (std::size_t j = 0; j < args.cov_columns().size() /*TODO: change if bias added*/; ++j)
      cov_mat(i, j) = full_pheno.cov_data[src_idx][j];
  }

  return true;
}

void slope_test()
{
  using namespace xt;
  using namespace xt::linalg;

  xarray<double> x = {1., 2., 3.};
  x.reshape({3,1});
  xtensor<double, 1> y = {3.1, 2.9, 3.2};
  //y.reshape({3,1});

  auto pbetas = dot(dot(pinv(dot(transpose(x), x)), transpose(x)), y);
  std::cerr << pbetas << std::endl;

  xtensor<double, 2> dmat = xt::concatenate(xtuple(xt::ones<double>({3, 1}), x), 1);
  std::cerr << dmat << std::endl;
  xarray<double> i = pinv(dot(transpose(dmat), dmat));
  auto pbetas2 = dot(dot(i, transpose(dmat)), y);
  std::cerr << pbetas2 << std::endl;

  //-----------------------//
  const std::size_t n = x.size();
  const double s_x     = std::accumulate(x.begin(), x.end(), 0.0);
  const double s_y     = std::accumulate(y.begin(), y.end(), 0.0);
  const double s_xx    = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
  const double s_xy    = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
  const double m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);

  const double b       = (s_y - m * s_x) / n;
  std::cerr << b << " , " << m << std::endl;
  auto fx              = [m,b](double x) { return m * x + b; };
  double se_line       = 0.0f; for (std::size_t i = 0; i < n; ++i) se_line += ::square(y[i] - fx(x[i]));
  const double x_mean  = s_x / n;
  double se_x_mean     = 0.0f; for (std::size_t i = 0; i < n; ++i) se_x_mean += ::square(x[i] - x_mean);
  const double dof     = n - 2;
  const double std_err = std::sqrt(se_line / dof) / std::sqrt(se_x_mean);
  float t = m / std_err;

}

void test()
{
  // Example from https://en.wikipedia.org/wiki/Simple_linear_regression#Numerical_example
  std::vector<double> x = {1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63, 1.65, 1.68, 1.70, 1.73, 1.75, 1.78, 1.80, 1.83};
//  std::vector<double> y = {52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47, 66.28, 68.10, 69.92, 72.19, 74.46};
  std::vector<double> y = {52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 29.93, 61.29, 163.11, 164.47, 66.28, 68.10, 69.92, 72.19, 74.46};

  std::size_t n = x.size();
  double s_x = std::accumulate(x.begin(), x.end(), 0.);
  double s_y = std::accumulate(y.begin(), y.end(), 0.);
  double s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.);
  double s_yy = std::inner_product(y.begin(), y.end(), y.begin(), 0.);
  double s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.);
  std::cerr << s_x << " " << s_y << std::endl;
  std::cerr << s_xx << " " << s_yy << std::endl;
  std::cerr << s_xy << std::endl;

  double beta = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  double alpha = (1./n) * s_y - beta * (1./n) * s_x;

  //se2 = se_line / dof
  double se2 = 1./(n*(n-2)) * (n*s_yy - s_y*s_y - beta*beta*(n*s_xx - s_x*s_x));
  double sbeta2 = n*se2 / (n*s_xx - s_x*s_x);
  double salpha2 = sbeta2*(1./n)*s_xx;

  double r = (n*s_xy - s_x*s_y) / std::sqrt((n*s_xx - s_x * s_x) * (n*s_yy - s_y * s_y));

  double m, std_err, t, pval;
  std::tie(m, std_err, t, pval) = linreg_ttest(y, x, s_y);

  return;
}

void challenger_test()
{
  using namespace xt;
  using namespace xt::linalg;

  xarray<double> x = {53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81};
  xtensor<double, 1> y = {1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0};
  assert(y.size() == x.size());
  x.reshape({y.size(), 1});
  x = xt::concatenate(xt::xtuple(xt::ones<double>({y.size(), std::size_t(1)}), x), 1);
  xarray<double> x_transpose_copy; // We use this later for the result of (X^t)W. W is too large (n x n), so we populate the result without generating W.

  xtensor<double, 1> beta = {2.90476190, -0.03738095}; // TODO: use linear model to produce initial betas

  std::size_t n_iter = 8;
  for (std::size_t i = 0; i < n_iter; ++i)
  {
    //std::cerr << xt::transpose(beta) << std::endl;
    //std::cerr << x << std::endl;

    xarray<double> z = dot(x, beta);
    //std::cerr << z << std::endl;
    xarray<double> p = 1. / (1. + xt::exp(-z));
    //std::cerr << p << std::endl;

    xarray<double> F = dot(transpose(x), y - p);
    //std::cerr << F << std::endl;
    x_transpose_copy = transpose(x);
    for (std::size_t i = 0; i < p.size(); ++i)
    {
      xt::col(x_transpose_copy, i) *= p(i) * (1. - p(i));
    }
    //std::cerr << transpose(x) << std::endl;
    //std::cerr << x_transpose_copy << std::endl;

    //xtensor<double, 2> I = -dot(x_transpose_copy, x);
    xtensor<double, 2> I = dot(x_transpose_copy, x);
    //std::cerr << I << std::endl;
    //std::cerr << pinv(I) << std::endl;

    beta = beta + dot(pinv(I), F);
    std::cerr << beta << std::endl;
  }

  auto a = 0;
}

int main(int argc, char** argv)
{
  //challenger_test();
  //test();
  //slope_test();
  //return test_xtensor();
  prog_args args;
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

//  std::string pheno_file_path = "../test-data/pheno_2000.tsv";
//  std::string geno_file_path = "../test-data/rand2000.freeze9.merged.chr20.filtered.anno.gtonly.minDP10.passonly.w280000.ligated.10000001-20000000.b8192.c19.sav";

  savvy::reader geno_file(args.geno_path());
  if (!geno_file)
    return std::cerr << "Could not open geno file\n", EXIT_FAILURE;

  if (args.region() && !geno_file.reset_bounds(*args.region()))
    return std::cerr << "Could not open genomic region\n", EXIT_FAILURE;

  std::string format_field = args.fmt_field();
  std::unordered_set<std::string> fmt_avail;

  for (const auto& h : geno_file.format_headers())
    fmt_avail.insert(h.id);

  if (format_field.empty())
  {
    if (fmt_avail.find("DS") != fmt_avail.end()) format_field = "DS";
    else if (fmt_avail.find("HDS") != fmt_avail.end()) format_field = "HDS";
    else if (fmt_avail.find("GT") != fmt_avail.end()) format_field = "GT";
    else return std::cerr << "Error: file must contain DS, HDS, or GT format fields\n", EXIT_FAILURE;
    std::cerr << "Notice: --fmt-field not specified so auto selecting " << format_field << std::endl;
  }
  else
  {
    if (fmt_avail.find(format_field) == fmt_avail.end())
      return std::cerr << "Error: requested format field (" << format_field << ") not found in file headers\n", EXIT_FAILURE;
  }


  xt::xtensor<scalar_type, 1> xresp;
  xt::xtensor<scalar_type, 2> xcov;
  if (!load_phenotypes(args, geno_file, xresp, xcov))
    return std::cerr << "Could not load phenotypes\n", EXIT_FAILURE;



//  auto xresp_data = xt::adapt(response_data, {response_data.size()});
//  auto xcov_data = xt::adapt(covariate_data, {std::size_t(covariate_data.size() / 2), std::size_t(2)});
  //residuals_type res = compute_residuals(xt::adapt(response_data, {response_data.size()}), xt::adapt(covariate_data, {std::size_t(covariate_data.size() / 2), std::size_t(cov_fields.size())}));
  residuals_type res = args.logit_enabled() ? compute_residuals_logit(xresp, xcov) : compute_residuals(xresp, xcov);
  std::cerr << res << std::endl;
  std::vector<scalar_type> res_std(res.begin(), res.end());
  scalar_type res_sum = std::accumulate(res_std.begin(), res_std.end(), scalar_type());
  scalar_type rss = std::inner_product(res_std.begin(), res_std.end(), res_std.begin(), scalar_type());

  //std::size_t n_samples = geno_file.samples().size();

  std::ofstream output_file(args.output_path(), std::ios::binary);

  output_file <<  "chrom\tpos\tmaf\tmac\tbeta\tse\tt\tpval\n";

  savvy::variant var;
  savvy::compressed_vector<scalar_type> sparse_geno;
  std::vector<scalar_type> dense_geno;
  while (geno_file >> var)
  {
    std::size_t ploidy = 0;
    bool is_sparse = false;
    bool found = false;
    for (const auto& f : var.format_fields())
    {
      if (f.first == format_field)
      {
        found = true;
        is_sparse = args.force_sparse() || (!args.sparse_disabled() && f.second.is_sparse());
        is_sparse ? f.second.get(sparse_geno) : f.second.get(dense_geno);
        ploidy = is_sparse ? sparse_geno.size() / res_std.size() : dense_geno.size() / res_std.size();
        is_sparse ? savvy::stride_reduce(sparse_geno, ploidy) : savvy::stride_reduce(dense_geno, ploidy);
        break;
      }
    }

    if (!found)
    {
      std::cerr << "Warning: skipping variant with not GT field\n";
      continue;
    }

    assert(ploidy != 0);

    float ac = 0.f, af = 0.f;
    std::int64_t an = 0;
    // For now, we are pulling from INFO fields but will likely always compute AC (along with case/ctrl AC) in the future.
    if (var.get_info("AC", ac) && var.get_info("AN", an) && an > 0)
    {
      af = float(ac) / an;
    }
    else if (!var.get_info("AF", af))
    {
      // For computing ac and af we use AN of sample subset.
      an = (res_std.size() * ploidy);
      ac = is_sparse ? std::accumulate(sparse_geno.begin(), sparse_geno.end(), 0.f) : std::accumulate(dense_geno.begin(), dense_geno.end(), 0.f);
      af = ac / an;
    }
    else
    {
      an = (geno_file.samples().size() * ploidy);
      ac = af * an;
    }

    float mac = (ac > (an/2) ? an - ac : ac);
    float maf = (af > 0.5 ? 1.f - af : af);

    if (mac < args.min_mac()) continue;

#if 0
    if (is_sparse)
    {
      float beta, se, t, pval;
      std::tie(beta, se, t, pval) = linreg_ttest(res_std, sparse_geno, res_sum);
      //std::tie(beta, se, t, pval) = linreg_ttest_fast(res_std, sparse_geno, res_sum);

      dense_geno.clear();
      dense_geno.resize(sparse_geno.size());
      for (auto it = sparse_geno.begin(); it != sparse_geno.end(); ++it)
        dense_geno[it.offset()] = *it;
      std::tie(beta, se, t, pval) = linreg_ttest(res_std, dense_geno, res_sum);
      auto a = 0;
    }
#endif

    float beta, se, t, pval;
    std::tie(beta, se, t, pval) = is_sparse ? linreg_ttest(res_std, sparse_geno, res_sum, rss) : linreg_ttest(res_std, dense_geno, res_sum);
    output_file << var.chromosome()
        << "\t" << var.position()
        << "\t" << maf
        << "\t" << mac
        << "\t" << beta
        << "\t" << se
        << "\t" << t
        << "\t" << pval << "\n";
  }

  // gzip -cd sp-reg-results-chr19-ldl.tsv | tail -n+2 | awk -F'\t' '$4>5  {print $2"\t"$8}' | gnuplot --persist -e "set logscale y; set yrange [0.99:5e-32] reverse; set xrange [1:65000000]; plot '-' using 1:2 w points"

  return 0;
}