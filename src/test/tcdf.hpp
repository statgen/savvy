/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TCDF_HPP
#define TCDF_HPP

#include <cmath>
#include <system_error>
#include <sstream>

static double betacf(double a, double b, double x) {
  int m,m2;
  double aa,c,d,del,h,qab,qam,qap;
  double FPMIN = 1e-30;
  int MAXIT = 1000;
  double EPS = 1e-10;

  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap; if (std::fabs(d) < FPMIN) d=FPMIN; d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {
    m2=2*m; aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if (std::fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (std::fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;
    if (std::fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (std::fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (std::fabs(del-1.0) < EPS) break;
  }
  if (m > MAXIT) {
    std::stringstream ss;
    ss << "a or b too big, or MAXIT too small in betacf" << a << b << x;
    throw std::runtime_error(ss.str());
  }
  return h;
}

inline static double gammln(double xx) {
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*std::log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++)
    ser += cof[j]/++y;
  return -tmp+std::log(2.5066282746310005*ser/x);
}

inline static double betai(double a, double b, double x) {
  double bt;
  if (x < 0.0 || x > 1.0) {
    throw std::runtime_error("Bad x in routine betai");
  }
  if (x == 0.0 || x == 1.0) bt=0.0;
  else bt=exp((gammln(a+b)-gammln(a)-gammln(b))+(a*std::log(x))+(b*std::log(1.0-x)));
  if (x < (a+1.0)/(a+b+2.0))
    return bt*betacf(a,b,x)/a;
  else
    return 1.0-bt*betacf(b,a,1.0-x)/b;
}

inline static double tcdf(double t, double nu) {
  if ( std::isnan(t) ) return 1.;
  else return betai(nu/2.,0.5,nu/(nu+t*t));
}

#endif