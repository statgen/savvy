/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/site_info.hpp"
#include "savvy/armadillo_vector.hpp"
#include "savvy/ublas_vector.hpp"
#include "savvy/reader.hpp"

#include <algorithm>
#include <numeric>
#include <iostream>


namespace ubl = boost::numeric::ublas;

template <typename T>
T square(const T& v) { return v * v; }

auto lin_reg(const std::vector<float>& x, const std::vector<float>& y)
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
  const float y_mean  = s_y / n;
  float se_y_mean     = 0.0f; for (std::size_t i = 0; i < n; ++i) se_y_mean += square(y[i] - y_mean);
  const float r2      = 1 - se_line / se_y_mean;

  return std::make_tuple(m, b, r2); // slope, y-intercept, r-squared
}

auto sp_lin_reg(const savvy::compressed_vector<float>& x, const std::vector<float>& y)
{
  const std::size_t n = x.size();
  const float s_x     = std::accumulate(x.begin(), x.end(), 0.0f);
  const float s_y     = std::accumulate(y.begin(), y.end(), 0.0f);
  const float s_xx    = std::inner_product(x.begin(), x.end(), x.begin(), 0.0f);
  float s_xy          = 0.0f; for (auto it = x.begin(); it != x.end(); ++it) s_xy += (*it * y[x.index_data()[it - x.begin()]]);
  const float m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const float b       = (s_y - m * s_x) / n;
  auto fx             = [m,b](float x) { return m * x + b; };
  float se_line       = 0.0f; for (std::size_t i = 0; i < n; ++i) se_line += square(y[i] - fx(x[i]));
  const float y_mean  = s_y / n;
  float se_y_mean     = 0.0f; for (std::size_t i = 0; i < n; ++i) se_y_mean += square(y[i] - y_mean);
  const float r2      = 1 - se_line / se_y_mean;

  return std::make_tuple(m, b, r2); // slope, y-intercept, r-squared
}

template <typename T1, typename T2>
auto ublas_lin_reg(T1& x, const T2& y)
{
  typedef typename T1::value_type flt_type;
  const std::size_t n    = x.size();
  const flt_type s_x     = ubl::sum(x);
  const flt_type s_y     = ubl::sum(y);
  const flt_type s_xx    = ubl::inner_prod(x, x);
  const flt_type s_xy    = ubl::inner_prod(x, y);
  const flt_type m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const flt_type b       = (s_y - m * s_x) / n;
  auto fx                = [m,b](float x) { return m * x + b; };
  flt_type se_line       = 0.0f; for (std::size_t i = 0; i < n; ++i) se_line += square(y[i] - fx(x[i]));
  const flt_type y_mean  = s_y / n;
  flt_type se_y_mean     = 0.0f; for (std::size_t i = 0; i < n; ++i) se_y_mean += square(y[i] - y_mean);
  const flt_type r2      = 1 - se_line / se_y_mean;

  return std::make_tuple(m, b, r2); // slope, y-intercept, r-squared
}

template <typename T>
auto arma_lin_reg(const arma::Col<T>& x, const arma::Col<T>& y)
{
  //typedef typename T1::value_type flt_type;
  typedef T flt_type;
  const std::size_t n    = x.size();
  const flt_type s_x     = arma::sum(x);
  const flt_type s_y     = arma::sum(y);
  const flt_type s_xx    = arma::dot(x, x);
  const flt_type s_xy    = arma::dot(x, y);
  const flt_type m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const flt_type b       = (s_y - m * s_x) / n;
  auto fx                = [m,b](float x) { return m * x + b; };
  flt_type se_line       = 0.0f; for (std::size_t i = 0; i < n; ++i) se_line += square(y[i] - fx(x[i]));
  const flt_type y_mean  = s_y / n;
  flt_type se_y_mean     = 0.0f; for (std::size_t i = 0; i < n; ++i) se_y_mean += square(y[i] - y_mean);
  const flt_type r2      = 1 - se_line / se_y_mean;

  return std::make_tuple(m, b, r2); // slope, y-intercept, r-squared
}

template <typename T>
auto arma_sp_lin_reg(const arma::SpCol<T>& x, const arma::SpCol<T>& y)
{
  //typedef typename T1::value_type flt_type;
  typedef T flt_type;
  const std::size_t n    = x.size();
  const flt_type s_x     = arma::sum(x);
  const flt_type s_y     = arma::sum(y);
  const flt_type s_xx    = arma::dot(x, x);
  const flt_type s_xy    = arma::dot(x, y);
  const flt_type m       = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const flt_type b       = (s_y - m * s_x) / n;
  auto fx                = [m,b](float x) { return m * x + b; };
  flt_type se_line       = 0.0f; for (std::size_t i = 0; i < n; ++i) se_line += square(y[i] - fx(x[i]));
  const flt_type y_mean  = s_y / n;
  flt_type se_y_mean     = 0.0f; for (std::size_t i = 0; i < n; ++i) se_y_mean += square(y[i] - y_mean);
  const flt_type r2      = 1 - se_line / se_y_mean;

  return std::make_tuple(m, b, r2); // slope, y-intercept, r-squared
}

template <typename T1, typename T2>
auto arma_r_squared(T1& x, const T2& y)
{
  typedef typename T1::value_type flt_type;
  flt_type ret = std::pow(arma::as_scalar(arma::cor(x, y)), 2);
  return ret;
}


int main()
{
  {
    savvy::armadillo::dense_allele_vector<float> x;
    x.resize(4);
    x[0] = -2;
    x[1] = 0;
    x[2] = 0;
    x[3] = 4;
    arma::Col<float> y = {-3, -1, 2, 3};
    float r2 = arma_r_squared(x, y);
    r2 = arma::as_scalar(arma::square(arma::cor(x, y)));
    std::cout << r2 << std::endl;

  }

  {
    savvy::armadillo::dense_allele_vector<float> x;
    x.resize(4);
    x[0] = -2;
    x[1] = 0;
    x[2] = 0;
    x[3] = 4;
    arma::Col<float> y = {-3, -1, 2, 3};
    float m, b, r2;
    std::tie(m, b, r2) = arma_lin_reg(x, y);
    std::cout << m << "," << b << "," << r2 << std::endl;

  }

  {
    savvy::armadillo::sparse_allele_vector<float> x;
    x.resize(4);
    x[0] = -2;
    //x[1] = 0;
    //x[2] = 0;
    x[3] = 4;
    arma::SpCol<float> y;
    y.resize(4, 1);
    y[0] = -3;
    y[1] = -1;
    y[2] = 2;
    y[3] = 3;

    float m, b, r2;
    std::tie(m, b, r2) = arma_sp_lin_reg(x, y);
    std::cout << m << "," << b << "," << r2 << std::endl;

  }

  {
    savvy::ublas::dense_allele_vector <float> x;
    x.resize(4);
    x[0] = -2;
    //x[1] = 0;
    //x[2] = 0;
    x[3] = 4;
    ubl::vector<float> y;
    y.resize(4);
    y[0] = -3;
    y[1] = -1;
    y[2] = 2;
    y[3] = 3;

    float m, b, r2;
    std::tie(m, b, r2) = ublas_lin_reg(x, y);
    std::cout << m << "," << b << "," << r2 << std::endl;
  }

  {
    savvy::dense_allele_vector<float> x;
    x.assign({-2, 0, 0, 4});
    std::vector<float> y = {-3, -1, 2, 3};

    float m, b, r2;
    std::tie(m, b, r2) = lin_reg(x, y);
    std::cout << m << "," << b << "," << r2 << std::endl;
  }

  {
    savvy::sparse_allele_vector<float> x;
    x.resize(4);
    x[0] = -2;
    //x[1] = 0;
    //x[2] = 0;
    x[3] = 4;
    std::vector<float> y = {-3, -1, 2, 3};

    float m, b, r2;
    std::tie(m, b, r2) = sp_lin_reg(x, y);
    std::cout << m << "," << b << "," << r2 << std::endl;
  }
  return 0;
}