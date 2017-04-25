#include "savvy/allele_vector.hpp"

#include <algorithm>
#include <numeric>
#include <iostream>
#include <cmath>

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
  float se_line       = 0.0f; for (std::size_t i = 0; i < n; ++i) se_line += std::pow(y[i] - fx(x[i]), 2);
  const float y_mean  = s_y / n;
  float se_y_mean     = 0.0f; for (std::size_t i = 0; i < n; ++i) se_y_mean += std::pow(y[i] - y_mean, 2);
  const float r2      = 1 - se_line / se_y_mean;

  return std::make_tuple(m, b, r2); // slope, y-intercept, r-squared
}

int main()
{
  savvy::dense_allele_vector<float> x;
  x.assign({-2, -1, 1, 4});
  std::vector<float> y = {-3, -1, 2, 3};

  float m, b, r2;
  std::tie(m, b, r2) = lin_reg(x, y);
  std::cout << m << "," << b << "," << r2 << std::endl;

  return 0;
}