/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/site_info.hpp"
//#include "savvy/armadillo_vector.hpp"
//#include "savvy/ublas_vector.hpp"
#include "savvy/reader.hpp"

#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>
#include <assert.h>
#include <tuple>
#include <type_traits>

#include <boost/math/distributions.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <OpenCL/opencl.h>
#include <boost/compute.hpp>
#include <savvy/ublas_vector.hpp>

namespace ublas = boost::numeric::ublas;

void tab_delimited_write_floats(std::ostream& os, std::tuple<float, float, float, float>&& fields)
{
  std::cout << std::get<0>(fields) << "\t" << std::get<1>(fields) << "\t" << std::get<2>(fields) << "\t" << std::get<3>(fields) << std::endl;
}

void print(std::ostream& os, const ublas::matrix<float>& m)
{
  for (std::size_t i = 0; i < m.size1(); ++i)
  {
    for (std::size_t j = 0; j < m.size2(); ++i)
    {
      os << m(i, j);
    }
    os << std::endl;
  }
}

template <typename T>
T square(const T& v)
{
  return v * v;
}

auto linreg_ttest_gpu(const std::vector<float>& x, const std::vector<float>& y, const boost::compute::vector<float>& device_x, const boost::compute::vector<float>& device_y, const boost::compute::vector<float>& device_coefs, boost::compute::command_queue& queue)
{
  const std::size_t n = x.size();
  auto start = std::chrono::high_resolution_clock().now();
  boost::compute::copy(x.begin(), x.end(), device_x.begin(), queue);
  auto op1 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();

  start = std::chrono::high_resolution_clock().now();
  float sum_x = std::accumulate(x.begin(), x.end(), 0.0f);
  float sum_y = std::accumulate(y.begin(), y.end(), 0.0f);
  auto op2 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();

  start = std::chrono::high_resolution_clock().now();
  float sum_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0f);
  float sum_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0f);
  auto op3 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();

  start = std::chrono::high_resolution_clock().now();
  float s_x    = {}; boost::compute::reduce(device_x.begin(), device_x.end(), device_coefs.begin(), queue);
  float s_y    = {}; boost::compute::reduce(device_y.begin(), device_y.end(), device_coefs.begin() + 1, queue);
  auto op4 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();

  start = std::chrono::high_resolution_clock().now();
  float s_xx    = boost::compute::inner_product(device_x.begin(), device_x.end(), device_x.begin(), 0.0f, queue);
  float s_xy    = boost::compute::inner_product(device_x.begin(), device_x.end(), device_y.begin(), 0.0f, queue);
  auto op5 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();
  auto h = 0;

  start = std::chrono::high_resolution_clock().now();
  boost::compute::transform_reduce(device_x.begin(), device_x.end(), device_x.begin(), device_coefs.begin() + 2, boost::compute::multiplies<float>(), boost::compute::plus<float>(), queue);
  boost::compute::transform_reduce(device_x.begin(), device_x.end(), device_y.begin(), device_coefs.begin() + 3, boost::compute::multiplies<float>(), boost::compute::plus<float>(), queue);
  auto op6 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();

  h = 0;

  const float m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const float b       = (s_y - m * s_x) / n;
  //auto fx             = [m,b](float x) { return m * x + b; };
  BOOST_COMPUTE_CLOSURE(float, lin_func, (float x, float y), (m, b),
  {
    float ret = y - (m * x + b);
    return ret * ret;
  });

  start = std::chrono::high_resolution_clock().now();
  float se_line{}; for (std::size_t i = 0; i < n; ++i) se_line += square(y[i] - (m * x[i] + b));
  auto op7 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();
  se_line = 0.0f;

  start = std::chrono::high_resolution_clock().now();
  boost::compute::transform_reduce(device_x.begin(), device_x.end(), device_y.begin(), &se_line, lin_func, boost::compute::plus<float>(), queue);
  auto op8 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();

  h = 0;
  std::cout << "    copy: " << op1 << std::endl;
  std::cout << " std sum: " << op2 << std::endl;
  std::cout << "std prod: " << op3 << std::endl;
  std::cout << " bst sum: " << op4 << std::endl;
  std::cout << "bst prod: " << op5 << std::endl;
  std::cout << "bst t_rd: " << op6 << std::endl;
  std::cout << "  std fx: " << op7 << std::endl;
  std::cout << "bst clsr: " << op8 << std::endl;
  std::cout << sum_x << sum_y << sum_xx << sum_xy << std::endl;
  std::cout << std::endl;

//
//  const float x_mean  = s_x / n;
//  float se_x_mean[1]     = {}; //for (std::size_t i = 0; i < n; ++i) se_x_mean += square(x[i] - x_mean);
//  boost::compute::transform_reduce(x.begin(), x.end(), std::begin(se_x_mean), boost::compute::placeholders::_1 - x_mean, boost::compute::plus<float>());
//  const float dof     = n - 2;
//  const float std_err = std::sqrt(se_line[0] / dof) / std::sqrt(se_x_mean[0]);
//  float t = m / std_err;
//  boost::math::students_t_distribution<float> dist(dof);
//  float pval = cdf(complement(dist, std::fabs(std::isnan(t) ? 0 : t))) * 2;
  return std::make_tuple(1.0, 1.0, 1.0, 1.0); //return std::make_tuple(m, std_err, t, pval); // slope, std error, t statistic, p value
}

auto linreg_ttest_old(const std::vector<float>& x, const std::vector<float>& y)
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
  boost::math::students_t_distribution<float> dist(dof);
  float pval = cdf(complement(dist, std::fabs(std::isnan(t) ? 0 : t))) * 2;
  return std::make_tuple(m, std_err, t, pval); // slope, std error, t statistic, p value
}

auto linreg_ttest(const std::vector<float>& x, const std::vector<float>& y, const float s_y)
{
  const std::size_t n = x.size();
  float s_x{}; //     = std::accumulate(x.begin(), x.end(), 0.0f);
  //float s_y{}; //     = std::accumulate(y.begin(), y.end(), 0.0f);
  float s_xx{}; //    = std::inner_product(x.begin(), x.end(), x.begin(), 0.0f);
  float s_xy{}; //    = std::inner_product(x.begin(), x.end(), y.begin(), 0.0f);

  for (std::size_t i = 0; i < n; ++i)
  {
    s_x += x[i];
    //s_y += y[i];
    s_xx += x[i] * x[i];
    s_xy += x[i] * y[i];
  }

  const float m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const float b       = (s_y - m * s_x) / n;
  auto fx             = [m,b](float x) { return m * x + b; };
  const float x_mean  = s_x / n;

  float se_line{};
  float se_x_mean{};
  for (std::size_t i = 0; i < n; ++i)
  {
    se_line += square(y[i] - fx(x[i]));
    se_x_mean += square(x[i] - x_mean);
  }

  const float dof     = n - 2;
  const float std_err = std::sqrt(se_line / dof) / std::sqrt(se_x_mean);
  float t = m / std_err;
  boost::math::students_t_distribution<float> dist(dof);
  float pval = cdf(complement(dist, std::fabs(std::isnan(t) ? 0 : t))) * 2;

  return std::make_tuple(m, std_err, t, pval); // slope, std error, t statistic, p value
}

auto sp_lin_reg_old(const savvy::compressed_vector<float>& x, const std::vector<float>& y)
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
  boost::math::students_t_distribution<float> dist(dof);
  float pval = cdf(complement(dist, std::fabs(t))) * 2;
  return std::make_tuple(m, std_err, t, pval); // slope, std error, t statistic, p value
}

auto linreg_ttest(const savvy::compressed_vector<float>& x, const std::vector<float>& y, const float s_y)
{
  const std::size_t n = x.size();
  float s_x{}; //     = std::accumulate(x.begin(), x.end(), 0.0f);
  float s_xx{}; //    = std::inner_product(x.begin(), x.end(), x.begin(), 0.0f);
  float s_xy{}; //    = std::inner_product(x.begin(), x.end(), y.begin(), 0.0f);

  const auto x_beg = x.begin();
  const auto x_end = x.end();
  const float* x_values = x.value_data();
  const std::size_t* x_indices = x.index_data();
  for (auto it = x_beg; it != x_end; ++it)
  {
    s_x += *it;
    s_xx += (*it) * (*it);
    s_xy += (*it) * y[x_indices[it - x_beg]];
  }

  //const float s_y     = std::accumulate(y.begin(), y.end(), 0.0f);
  const float m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const float b       = (s_y - m * s_x) / n;
  auto fx             = [m,b](float x) { return m * x + b; };
  const float x_mean  = s_x / n;

  float se_line{};
  float se_x_mean{};

  const float f_of_zero = fx(0.0f);
  const auto x_idx_beg = x.index_data();
  const auto x_idx_end = x.index_data() + x.non_zero_size();
  std::size_t i = 0;
  for (auto it = x_idx_beg; it != x_idx_end; ++i)
  {
    if (i == (*it))
    {
      float x_val = x_values[it - x_idx_beg];
      se_line += square(y[i] - fx(x_val));
      se_x_mean += square(x_val - x_mean);
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

  se_x_mean += (square(0.0f - x_mean) * float(n - x.non_zero_size()));

  const float dof     = n - 2;
  const float std_err = std::sqrt(se_line / dof) / std::sqrt(se_x_mean);
  float t = m / std_err;
  boost::math::students_t_distribution<float> dist(dof);
  float pval = cdf(complement(dist, std::fabs(std::isnan(t) ? 0 : t))) * 2;

  return std::make_tuple(m, std_err, t, pval); // slope, std error, t statistic, p value
}

auto multi_lin_reg_ttest(const ublas::vector<float>& geno, const ublas::vector<float>& pheno, const ublas::matrix<double>& covariates)
{
  const ublas::matrix<float> beta;
  return true;
}

//template <typename T1, typename T2>
//auto ublas_lin_reg(T1& x, const T2& y)
//{
//  typedef typename T1::value_type flt_type;
//  const std::size_t n    = x.size();
//  const flt_type s_x     = ubl::sum(x);
//  const flt_type s_y     = ubl::sum(y);
//  const flt_type s_xx    = ubl::inner_prod(x, x);
//  const flt_type s_xy    = ubl::inner_prod(x, y);
//  const flt_type m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
//  const flt_type b       = (s_y - m * s_x) / n;
//  auto fx                = [m,b](float x) { return m * x + b; };
//  flt_type se_line       = 0.0f; for (std::size_t i = 0; i < n; ++i) se_line += square(y[i] - fx(x[i]));
//  const flt_type y_mean  = s_y / n;
//  flt_type se_y_mean     = 0.0f; for (std::size_t i = 0; i < n; ++i) se_y_mean += square(y[i] - y_mean);
//  const flt_type r2      = 1 - se_line / se_y_mean;
//
//  return std::make_tuple(m, b, r2); // slope, y-intercept, r-squared
//}

void run_simple(const std::string& file_path)
{
  auto start = std::chrono::high_resolution_clock().now();
  std::int64_t compute_time = 0;

//  boost::compute::device device = boost::compute::system::default_device();
//  boost::compute::context context(device);
//  boost::compute::command_queue q(context, device);

  savvy::site_info anno;
  savvy::compressed_vector<float> x;
  savvy::reader<1> r(file_path, savvy::fmt::allele);

  std::random_device rnd_device;
  std::mt19937 mersenne_engine(rnd_device());
  std::uniform_int_distribution<int> dist(0, 100);
  auto gen = std::bind(dist, mersenne_engine);

  std::vector<float> y(r.sample_size() * 2);
  generate(y.begin(), y.end(), gen);
  const float y_sum = std::accumulate(y.begin(), y.end(), 0.0f);

//  std::vector<float> fake_x(150000, 1.0f);
//  std::vector<float> fake_y(150000, 1.0f);
//  std::generate(fake_x.begin(), fake_x.end(), gen);
//  std::generate(fake_y.begin(), fake_y.end(), gen);
//
//  boost::compute::vector<float> device_y(fake_y.begin(), fake_y.end(), q);
//  boost::compute::vector<float> device_x(fake_y.size(), 0.0f, q);
//  boost::compute::vector<float> device_coefs(4, 0.0f, q);

  std::cout << "pos\tref\talt\tslope\tse\ttstat\tpval\n";
  std::cout << std::fixed << std::setprecision( 4 );
  int i = 0;
  while (r.read(anno, x))
  {
    float m, se, tstat, pval;
    //auto compute_start = std::chrono::high_resolution_clock().now();
    std::tie(m, se, tstat, pval) = linreg_ttest(x, y, y_sum); //linreg_ttest_gpu(fake_x, fake_y, device_x, device_y, device_coefs, q);
    //compute_time += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - compute_start).count();

//    ++i;
//    if (i == 5)
//      break;

//    std::cout << x.locus() << "\t" << x.ref() << "\t" << x.alt() << "\t";
//    std::cout << m << "\t" << se << "\t" << tstat << "\t" << pval << std::endl;
  }
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock().now() - start).count();
  std::cout << "wall time: " << elapsed << "s" << std::endl;
  //std::cout << "compute time: : " <<  std::chrono::duration_cast<std::chrono::seconds>(std::chrono::microseconds(compute_time)).count() << "s" << std::endl;
}

void print_matrix(const ublas::matrix<float>& m)
{
  for (auto i = 0; i < m.size1(); i++)
  {
    for (auto j = 0; j < m.size2(); j++)
      std::cout << m(i, j) << " ";
    std::cout << std::endl;
  }
}

void print_matrix(const ublas::vector<float>& m)
{
  for (auto i = 0; i < m.size(); i++)
    std::cout << m(i) << " ";
  std::cout << std::endl;
}

template <typename T>
auto invert_matrix(const ublas::matrix<T>& input)
{
  // create a working copy of the input
  ublas::matrix<T> a(input);

  // create a permutation matrix for the LU-factorization
  ublas::permutation_matrix<std::size_t> pm(a.size1());

  // perform LU-factorization
  unsigned long res = lu_factorize(a, pm);
  if (res != 0)
    throw std::invalid_argument("Not invertable");

  ublas::matrix<T> ret(a.size1(), a.size2());

  // create identity matrix of "inverse"
  ret.assign(ublas::identity_matrix<T> (a.size1()));

  // backsubstitute to get the inverse
  lu_substitute(a, pm, ret);

  return ret;
}

static_assert(sizeof(std::array<char, 1>) == sizeof(char), "Compiler not supported. std::array has extra data. This is allowed by standard but not expected.");

template <std::size_t NumCovariates>
class gpu_multi_reg
{
public:
  static const std::size_t num_cols = NumCovariates + 2;
public:
  gpu_multi_reg(boost::compute::command_queue& queue, const std::vector<float>& observed_responses, const std::vector<std::array<float, NumCovariates>>& covariates) :
    queue_(queue),
    observed_responses_(observed_responses.begin(), observed_responses.end(), queue_),
    covariates_(covariates.size())
  {
//    for (auto & row : covariates_)
//      row[0] = 1.0;
//    for (std::size_t i = 0; i < covariates.size(); ++i)
//      std::copy(covariates[i].begin(), covariates[i].end(), covariates_[i + 1].begin());
  }

  auto operator()(const std::vector<float>& genotypes)
  {

//    const std::size_t num_rows = observed_responses_.size();
//    for (std::size_t i = 0; i < num_rows; ++i)
//      covariates_[i][num_cols - 1] = genotypes[i];
//


//    ublas::matrix<float> beta_variances;
//    try
//    {
//      beta_variances = invert_matrix<float>(ublas::prod(ublas::trans(covariates_), covariates_));
//    }
//    catch (std::exception& e)
//    {
//      std::cerr << e.what() << std::endl;
//      return std::make_tuple(ublas::vector<float>(covariates_.size2(), 0.0), ublas::vector<float>(covariates_.size2(), 0.0));
//    }
//
//    ublas::vector<float> beta = ublas::prod(ublas::prod(beta_variances, ublas::trans(covariates_)), observed_responses_);
//    ublas::vector<float> residuals = observed_responses_ - ublas::prod(covariates_, beta);
//
//    float square_error{};
//    for (const float& r : residuals)
//      square_error += square(r);
//
//    const float se_line_mean = square_error / num_rows;
//
//    const float dof = num_rows - num_cols; //n - (k + 1)
//    const float variance = square_error / dof;
//
//    ublas::matrix<float> var_cov_mat = variance * beta_variances;
//
//    ublas::vector<float> standard_errors(num_cols);
//
//    for (std::size_t i = 0; i < num_cols; ++i)
//      standard_errors[i] = std::sqrt(var_cov_mat(i,i));
//
//    return std::make_tuple(std::move(beta), std::move(standard_errors));
  }
private:
  template <typename T, std::size_t Sz>
  struct required_float_type
  {
    typedef typename std::conditional<std::is_same<T, float>::value && Sz == 1, float,
      typename std::conditional<std::is_same<T, float>::value && Sz == 2, boost::compute::float2_,
        typename std::conditional<std::is_same<T, float>::value && (Sz == 3 || Sz == 4), boost::compute::float4_,
          typename std::conditional<std::is_same<T, float>::value && (Sz > 4 && Sz <= 8), boost::compute::float8_,
            typename std::conditional<std::is_same<T, float>::value && (Sz > 8 && Sz <= 16), boost::compute::float16_,
              typename std::conditional<std::is_same<T, double>::value && Sz == 1, double,
                typename std::conditional<std::is_same<T, double>::value && Sz == 2, boost::compute::double2_,
                  typename std::conditional<std::is_same<T, double>::value && (Sz == 3 || Sz == 4), boost::compute::double4_,
                    typename std::conditional<std::is_same<T, double>::value && (Sz > 4 && Sz <= 8), boost::compute::double8_,
                      typename std::conditional<std::is_same<T, double>::value && (Sz > 8 && Sz <= 16), boost::compute::double16_, void>::type>::type>::type>::type>::type>::type>::type>::type>::type>::type type;

//    typedef std::conditional<
//      std::is_same<T, float>::value && Sz == 1,
//      float,
//      std::conditional<
//        std::is_same<T, float>::value && Sz == 2,
//        boost::compute::float2_,
//        boost::compute::float4_
//      >::type
//    >::type type;
  };
private:
  boost::compute::command_queue & queue_;
  boost::compute::vector<float> observed_responses_;
  typename required_float_type<float, 1>::type foo_;
  boost::compute::vector<typename required_float_type<float, 16>::type> covariates_;

};

class slow_multi_reg
{
public:
  slow_multi_reg(const ublas::vector<float>& observed_responses, const ublas::matrix<float, ublas::column_major>& covariates) :
    observed_responses_(observed_responses),
    covariates_(covariates.size1(), covariates.size2() + 2)
  {
    std::fill((covariates_.begin2() + 0).begin(), (covariates_.begin2() + 0).end(), 1.0);
    for (std::size_t i = 0; i < covariates.size2(); ++i)
      std::copy((covariates.begin2() + i).begin(), (covariates.begin2() + i).end(), (covariates_.begin2() + (i + 1)).begin());
  }

  auto operator()(const ublas::vector<float>& genotypes)
  {
    std::copy(genotypes.begin(), genotypes.end(), (covariates_.begin2() + covariates_.size2() - 1).begin());

    const std::size_t num_rows = observed_responses_.size();
    const std::size_t num_cols = covariates_.size2();

    ublas::matrix<float> beta_variances;
    try
    {
      beta_variances = invert_matrix<float>(ublas::prod(ublas::trans(covariates_), covariates_));
    }
    catch (std::exception& e)
    {
      //std::cerr << e.what() << std::endl;
      return std::make_tuple(ublas::vector<float>(covariates_.size2(), 0.0), ublas::vector<float>(covariates_.size2(), 0.0));
    }

    ublas::vector<float> beta = ublas::prod(ublas::prod(beta_variances, ublas::trans(covariates_)), observed_responses_);
    ublas::vector<float> residuals = observed_responses_ - ublas::prod(covariates_, beta);

    float square_error{};
    for (const float& r : residuals)
      square_error += square(r);

    const float se_line_mean = square_error / num_rows;

    const float dof = num_rows - num_cols; //n - (k + 1)
    const float variance = square_error / dof;

    ublas::matrix<float> var_cov_mat = variance * beta_variances;

    ublas::vector<float> standard_errors(num_cols);

    for (std::size_t i = 0; i < num_cols; ++i)
      standard_errors[i] = std::sqrt(var_cov_mat(i,i));

    return std::make_tuple(std::move(beta), std::move(standard_errors));
  }
private:
  const ublas::vector<float>& observed_responses_;
  ublas::matrix<float> covariates_;

};

class fast_multi_reg
{
public:
  fast_multi_reg(const ublas::vector<float>& observed_responses, const ublas::matrix<float, ublas::column_major>& covariates) :
    observed_responses_(observed_responses),
    covariates_(covariates.size1(), covariates.size2() + 1),
    row_maj_covariates_(covariates.size1(), covariates.size2() + 1),
    transpose_product_matrix_(covariates_.size2() + 1, covariates_.size2() + 1)
  {
    std::fill((covariates_.begin2() + 0).begin(), (covariates_.begin2() + 0).end(), 1.0);
    for (std::size_t i = 0; i < covariates.size2(); ++i)
      std::copy((covariates.begin2() + i).begin(), (covariates.begin2() + i).end(), (covariates_.begin2() + (i + 1)).begin());
    ublas::matrix<float> tmp = ublas::prod(ublas::trans(covariates_), covariates_);

    std::size_t i = 0;
    for (auto row = tmp.begin1(); row != tmp.end1(); ++row, ++i)
      std::copy(row.begin(), row.end(), (transpose_product_matrix_.begin1() + i).begin());

    row_maj_covariates_ = covariates_;
  }

  auto operator()(const ublas::vector<float>& genotypes)
  {
    const std::size_t num_rows = observed_responses_.size();
    const std::size_t num_cols = covariates_.size2() + 1;
    const std::size_t num_covs = covariates_.size2();
    const std::size_t last_index_of_result_matrix = num_covs;

    // zeros bottom row
    for (std::size_t i = 0; i < num_cols; ++i)
    {
      transpose_product_matrix_(last_index_of_result_matrix, i) = 0.0;
      //transpose_product_matrix_(i, last_index_of_result_matrix) = 0.0;
    }

    // Sets bottom row of result matrix.
    for (std::size_t i = 0; i < num_rows; ++i) // TODO: test performance of flipping loop order.
    {
      for (std::size_t j = 0; j < num_covs; ++j)
      {
        transpose_product_matrix_(last_index_of_result_matrix, j) += (genotypes[i] * covariates_(i, j));
      }
      transpose_product_matrix_(last_index_of_result_matrix, last_index_of_result_matrix) += (genotypes[i] * genotypes[i]);
    }

    // Fills reflected column.
    for (std::size_t j = 0; j < num_covs; ++j)
      transpose_product_matrix_(j, last_index_of_result_matrix) = transpose_product_matrix_(last_index_of_result_matrix, j);

    ublas::matrix<float> beta_variances;
    try
    {
      beta_variances = invert_matrix<float>(transpose_product_matrix_);
    }
    catch (std::exception& e)
    {
      return std::make_tuple(ublas::vector<float>(num_cols, 0.0), ublas::vector<float>(num_cols, 0.0));
    }

    ublas::matrix<float> beta_var_and_x_trans_prod(beta_variances.size1(), num_rows, 0.0);

    for (std::size_t i = 0; i < num_cols; ++i)
    {
      for (std::size_t j = 0; j < num_rows; ++j)
      {
        for (std::size_t k = 0; k < num_covs; ++k)
        {
          beta_var_and_x_trans_prod(i, j) += (beta_variances(i, k) * row_maj_covariates_(j, k));
        }
        beta_var_and_x_trans_prod(i, j) += (beta_variances(i, last_index_of_result_matrix) * genotypes[j]);
      }
    }

    ublas::vector<float> beta = ublas::prod(beta_var_and_x_trans_prod, observed_responses_);

    ublas::vector<float> y_hat(num_rows, 0.0);
    for (std::size_t i = 0; i < num_rows; ++i)
    {
      for (std::size_t j = 0; j < num_covs; ++j)
      {
        y_hat[i] += (covariates_(i, j) * beta[j]);
      }
      y_hat[i] += (genotypes[i] * beta[last_index_of_result_matrix]);
    }


    ublas::vector<float> residuals = observed_responses_ - y_hat;

    float square_error{};
    for (const float& r : residuals)
      square_error += square(r);

    const float se_line_mean = square_error / num_rows;

    const float dof = num_rows - num_cols; //n - (k + 1)
    const float variance = square_error / dof;

    ublas::matrix<float> var_cov_mat = variance * beta_variances;

    ublas::vector<float> standard_errors(num_cols);

    for (std::size_t i = 0; i < num_cols; ++i)
      standard_errors[i] = std::sqrt(var_cov_mat(i,i));

    return std::make_tuple(std::move(beta), std::move(standard_errors));
  }
private:
  const ublas::vector<float>& observed_responses_;
  ublas::matrix<float, ublas::column_major> covariates_;
  ublas::matrix<float, ublas::row_major> row_maj_covariates_;
  ublas::matrix<float> transpose_product_matrix_;
};

const std::vector<float> predictor_1 = {41.9, 43.4, 43.9, 44.5, 47.3, 47.5, 47.9, 50.2, 52.8, 53.2, 56.7, 57.0, 63.5, 65.3, 71.1, 77.0, 77.8};
const std::vector<float> predictor_2 = {29.1, 29.3, 29.5, 29.7, 29.9, 30.3, 30.5, 30.7, 30.8, 30.9, 31.5, 31.7, 31.9, 32.0, 32.1, 32.5, 32.9};
const std::vector<float> response = {251.3, 251.3, 248.3, 267.5, 273.0, 276.5, 270.3, 274.9, 285.0, 290.0, 297.0, 302.5, 304.5, 309.3, 321.7, 330.7, 349.0};

void run_multi(const std::string& file_path)
{
//  ublas::vector<float> ty(response.size());
//  std::copy(response.begin(), response.end(), ty.begin());
//
//  ublas::matrix<float, ublas::column_major> tcov(predictor_2.size(), 1);
//  std::copy(predictor_2.begin(), predictor_2.end(), (tcov.begin2()).begin());
//
//  ublas::vector<float> tgeno(predictor_1.size());
//  std::copy(predictor_1.begin(), predictor_1.end(), tgeno.begin());

  std::random_device rnd_device;
  std::mt19937 mersenne_engine(rnd_device());
  std::uniform_real_distribution<float> dist(0.0, 1.0);
  auto gen = std::bind(dist, mersenne_engine);

  savvy::reader<1> r(file_path, savvy::fmt::allele);

  ublas::vector<float> y(r.sample_size() * 2);
  generate(y.begin(), y.end(), gen);

  ublas::matrix<float> cov(r.sample_size() * 2, 8);
  for (auto row_it = cov.begin1(); row_it != cov.end1(); ++row_it)
    generate(row_it.begin(), row_it.end(), gen);

  auto start = std::chrono::high_resolution_clock().now();
  //slow_multi_reg slow_reg(y, cov);
  fast_multi_reg fast_reg(y, cov);

  savvy::site_info anno;
  boost::numeric::ublas::vector<float> variant;
  while (r.read(anno, variant))
  {
    ublas::vector<float> coefs, se;
    //std::tie(coefs, se) = slow_reg(variant);
    std::tie(coefs, se) = fast_reg(variant);


    boost::math::students_t_distribution<float> dist(y.size() - (2 + 1));
    float t = coefs[2] / se[2];
    float pval = cdf(complement(dist, std::fabs(std::isnan(t) ? 0 : t))) * 2;
  }
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock().now() - start).count();
  std::cout << "wall time: " << elapsed << "s" << std::endl;
}

//void run_multi(const std::string& file_path)
//{
//  auto start = std::chrono::high_resolution_clock().now();
//
//  const std::size_t num_rows = response.size();
//
//  ublas::vector<float> y(num_rows);
//  std::copy(response.begin(), response.end(), y.begin());
//
//
//  ublas::matrix<float, ublas::column_major> x(num_rows, 3);
//  const std::size_t num_cols = x.size2();
//  std::fill((x.begin2() + 0).begin(), (x.begin2() + 0).end(), 1.0);
//  std::copy(predictor_1.begin(), predictor_1.end(), (x.begin2() + 1).begin());
//  std::copy(predictor_2.begin(), predictor_2.end(), (x.begin2() + 2).begin());
//
//
//  ublas::matrix<float> beta_variances = invert_matrix<float>(ublas::prod(ublas::trans(x), x));
//
//  ublas::vector<float> beta = ublas::prod(ublas::prod(beta_variances, ublas::trans(x)), y); // ublas::trans(beta) ;
//  auto fx             = [&beta](const std::vector<float>& x) { return beta[0] + beta[1] * x[0] + beta[2] * x[1]; };
//  std::cout << fx({47, 31}) << std::endl;
//
//  ublas::vector<float> residuals = y - ublas::prod(x, beta);
//  float square_error{};
//  for (const float& r : residuals)
//  {
//    std::cout << r << " ";
//    square_error += square(r);
//  }
//  std::cout << std::endl;
//  const float mean_square_error = square_error / residuals.size();
//  std::cout << mean_square_error << std::endl;
//
//  float se_line{};
//  float se_x_mean{};
//  for (std::size_t i = 0; i < num_rows; ++i)
//  {
//    std::cout << (y[i] - fx({x(i, 1), x(i, 2)})) << " ";
//    se_line += square(y[i] - fx({x(i, 1), x(i, 2)}));
//  }
//  std::cout << std::endl;
//  const float se_line_mean = se_line / num_rows;
//
//  const float dof = num_rows - num_cols; //n - (k + 1)
//  const float variance = se_line / dof;
//
//  ublas::matrix<float> var_cov_mat = variance * beta_variances;
//
//
//
//  const float std_err_beta_1 = std::sqrt(var_cov_mat(1, 1));
//  const float std_err_beta_2 = std::sqrt(var_cov_mat(2, 2));
//
//  const float b1_t = beta[1] / std_err_beta_1;
//  const float b2_t = beta[2] / std_err_beta_2;
//
//  boost::math::students_t_distribution<float> dist(dof);
//  float pval1 = cdf(complement(dist, std::fabs(b1_t))) * 2;
//  float pval2 = cdf(complement(dist, std::fabs(b2_t))) * 2;
//
//  std::cout << beta[0] << " + " << beta[1] << "x + " << beta[2] << "x + " << ublas::sum(residuals) << std::endl;
//
//
//  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock().now() - start).count();
//  std::cout << "elapsed: " << elapsed << "s" << std::endl;
//}


//int main(int argc, char** argv)
//{
//  boost::compute::device device = boost::compute::system::default_device();
//  boost::compute::context context(device);
//  boost::compute::command_queue queue(context, device);
//
//  std::vector<float> a;
//  std::vector<std::array<float, 2>> b;
//  gpu_multi_reg<2> test(queue, a, b);
//
//  boost::compute::float4_ foo;
//  std::vector<std::array<float, 4>> bar;
//  std::vector<std::tuple<float, float, float, float>> man;
//  //foo = man;
//
//  boost::compute::vector<float> res(4, 0.0, queue);
//  boost::compute::vector<boost::compute::float4_> res_4(1, boost::compute::float4_(0, 0, 0, 0), queue);
//
//
//  std::vector<float> v(1000);
//  boost::compute::vector<float> device_vec(1000 * 4, queue.get_context());
//  boost::compute::vector<boost::compute::float4_> device_vec_4(1000, boost::compute::float4_(1, 2, 3, 0), queue);
//
//
//
//
//  {
//    auto start = std::chrono::high_resolution_clock().now();
//    boost::compute::copy(v.begin(), v.end(), boost::compute::make_transform_iterator(device_vec_4.begin(), boost::compute::get<0>()), queue);
//    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();
//    std::cout << elapsed << std::endl;
//  }
//
//  {
//    auto start = std::chrono::high_resolution_clock().now();
//    boost::compute::copy(v.begin(), v.end(), boost::compute::make_transform_iterator(device_vec_4.begin(), boost::compute::get<0>()), queue);
//    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();
//    std::cout << elapsed << std::endl;
//  }
//
//  {
//    auto start = std::chrono::high_resolution_clock().now();
//    boost::compute::copy(v.begin(), v.end(), boost::compute::make_transform_iterator(device_vec_4.begin(), boost::compute::get<0>()), queue);
//    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();
//    std::cout << elapsed << std::endl;
//  }
//
//  {
//    auto start = std::chrono::high_resolution_clock().now();
//    boost::compute::copy(v.begin(), v.end(), boost::compute::make_transform_iterator(device_vec_4.begin(), boost::compute::get<0>()), queue);
//    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();
//    std::cout << elapsed << std::endl;
//  }
//
//  std::cout << std::endl;
//
//
//  {
//    auto start = std::chrono::high_resolution_clock().now();
//    boost::compute::copy(v.begin(), v.end(), device_vec.begin(), queue);
//    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();
//    std::cout << elapsed << std::endl;
//  }
//
//  {
//    auto start = std::chrono::high_resolution_clock().now();
//    boost::compute::copy(v.begin(), v.end(), device_vec.begin(), queue);
//    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();
//    std::cout << elapsed << std::endl;
//  }
//
//  {
//    auto start = std::chrono::high_resolution_clock().now();
//    boost::compute::copy(v.begin(), v.end(), device_vec.begin(), queue);
//    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();
//    std::cout << elapsed << std::endl;
//  }
//
//  {
//    auto start = std::chrono::high_resolution_clock().now();
//    boost::compute::copy(v.begin(), v.end(), device_vec.begin(), queue);
//    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();
//    std::cout << elapsed << std::endl;
//  }
//
//  return 0;
//
//
//  run_multi(argv[1]);
//
//
//
//
//  return 0;
//}

int main()
{
//  for (const auto& device : boost::compute::system::devices())
//  {
//    std::cout << "NAME: " << device.name() << std::endl;
//    std::cout << "compute_units: " << device.compute_units() << std::endl;
//    std::cout << "max_work_group_size: " << device.max_work_group_size() << std::endl;
//    std::cout << "max_work_item_dimensions: " << device.max_work_item_dimensions() << std::endl;
//    std::cout << "preferred_vector_width: " << device.preferred_vector_width<float>() << std::endl;
//    std::cout << "profile: " << device.profile() << std::endl;
//    std::cout << "max_memory_alloc_size: " << device.max_memory_alloc_size() << std::endl;
//    std::cout << "local_memory_size: " << device.local_memory_size() << std::endl;
//    std::cout << "global_memory_size: " << device.global_memory_size() / 1024 / 1024 << std::endl;
//    std::cout << "clock_frequency: " << device.clock_frequency() << std::endl;
//    std::cout << std::endl;
//  }


  auto device = boost::compute::system::default_device();
  boost::compute::context context(device);
  boost::compute::command_queue queue(context, device);

  std::vector<float> x(3000000, 1.0f);
  std::cout << x.size() << std::endl;

  float cpu_res{};
  std::vector<float> gpu_res(100);

  boost::compute::vector<float> device_x(x.begin(), x.end(), queue);
  boost::compute::vector<float> device_result(gpu_res.size(), 0.0f, queue);

  auto start = std::chrono::high_resolution_clock().now();;
  for (std::size_t i = 0; i < gpu_res.size(); ++i)
    cpu_res = std::accumulate(x.begin(), x.end(), 0.0f);
  auto cpu_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();


  start = std::chrono::high_resolution_clock().now();
  boost::compute::copy_async(x.begin(), x.end(), device_x.begin(), queue);
  for (std::size_t i = 0; i < gpu_res.size(); ++i)
    boost::compute::reduce(device_x.begin(), device_x.end(), device_result.begin() + i, queue);
  boost::compute::copy(device_result.begin(), device_result.end(), gpu_res.begin(), queue);
  auto gpu_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock().now() - start).count();

  std::cout << "cpu res: " << cpu_res << std::endl;
  std::cout << "gpu res: " << gpu_res.back() << std::endl;

  std::cout << "cpu time: " << cpu_time << "us" << std::endl;
  std::cout << "gpu time: " << gpu_time << "us" << std::endl;



  return 0;
}