#include "savvy/allele_vector.hpp"
//#include "savvy/armadillo_vector.hpp"
//#include "savvy/ublas_vector.hpp"
#include "savvy/reader.hpp"

#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <iostream>

#include <boost/math/distributions.hpp>


//namespace ubl = boost::numeric::ublas;

template <typename T>
T square(const T& v) { return v * v; }

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
  float pval = cdf(complement(dist, std::fabs(t))) * 2;
  return std::make_tuple(m, std_err, t, pval); // slope, std error, t statistic, p value
}

auto linreg_ttest(const std::vector<float>& x, const std::vector<float>& y)
{
  const std::size_t n = x.size();
  float s_x{}; //     = std::accumulate(x.begin(), x.end(), 0.0f);
  float s_y{}; //     = std::accumulate(y.begin(), y.end(), 0.0f);
  float s_xx{}; //    = std::inner_product(x.begin(), x.end(), x.begin(), 0.0f);
  float s_xy{}; //    = std::inner_product(x.begin(), x.end(), y.begin(), 0.0f);
  for (std::size_t i = 0; i < n; ++i)
  {
    s_x += x[i];
    s_y += y[i];
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
  float pval = cdf(complement(dist, std::fabs(t))) * 2;

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

auto linreg_ttest(const savvy::compressed_vector<float>& x, const std::vector<float>& y)
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
    s_xy += (*it * y[x_indices[it - x_beg]]);
  }


  const float s_y     = std::accumulate(y.begin(), y.end(), 0.0f);
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
  float pval = cdf(complement(dist, std::fabs(t))) * 2;

  return std::make_tuple(m, std_err, t, pval); // slope, std error, t statistic, p value
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


int main(int argc, char** argv)
{

//  {
//    savvy::ublas::dense_allele_vector <float> x;
//    x.resize(4);
//    x[0] = -2;
//    //x[1] = 0;
//    //x[2] = 0;
//    x[3] = 4;
//    ubl::vector<float> y;
//    y.resize(4);
//    y[0] = -3;
//    y[1] = -1;
//    y[2] = 2;
//    y[3] = 3;
//
//    float m, b, r2;
//    std::tie(m, b, r2) = ublas_lin_reg(x, y);
//    std::cout << m << "," << b << "," << r2 << std::endl;
//  }

  {
    auto start = std::chrono::high_resolution_clock().now();
    savvy::dense_allele_vector<float> x;
    savvy::reader r(argv[1]);

    std::random_device rnd_device;
    std::mt19937 mersenne_engine(rnd_device());
    std::uniform_int_distribution<int> dist(0, 100);
    auto gen = std::bind(dist, mersenne_engine);

    std::vector<float> y(r.sample_size() * 2);
    generate(y.begin(), y.end(), gen);

    std::cout << "pos\tref\talt\tslope\tse\ttstat\tpval\n";
    std::cout << std::fixed << std::setprecision( 4 );
    while (r.read(x, std::numeric_limits<float>::epsilon()))
    {
      float m, se, tstat, pval;
      std::tie(m, se, tstat, pval) = linreg_ttest(x, y);

      std::cout << x.locus() << "\t" << x.ref() << "\t" << x.alt() << "\t";
      std::cout << m << "\t" << se << "\t" << tstat << "\t" << pval << std::endl;
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock().now() - start).count();
    std::cout << "elapsed: " << elapsed << "s" << std::endl;
  }

  return 0;
}