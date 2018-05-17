#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

void rk4(double ta, double tb, double h, std::vector<double> & y);
double f(double t, const std::vector<double> & y, int id);

const double O=0.25;
const double A=1;
const double C=1;
const double a=1.25;
const double b=1;
const double c=0.75;
int main (void)
{
  const double N  = 6;
  const double TA = 0;
  const double TB = 10000;
  const double H  = 0.1;
  std::vector<double> y = {2.5, 0, 0, 0, 2.355, 0}; // {x0, y0, z0, vx0, vy0, vz0}  ---1.707---

   rk4(TA,TB, H, y);

   return 0;
 }
 void rk4(double ta, double tb, double h, std::vector<double> & y)
 {
   std::vector<double> k1, k2, k3, k4, aux;
   k1.resize(y.size());
   k2.resize(y.size());
   k3.resize(y.size());
   k4.resize(y.size());
   aux.resize(y.size());

   const int N = (tb-ta)/h;

   for (int nt = 0; nt < N; ++nt) {
     double t = ta + h*nt;
     // k1
     for(int ii = 0; ii < y.size(); ++ii) {
       k1[ii] = h*f(t, y, ii);
     }
     // k2 aux
     for(int ii = 0; ii < y.size(); ++ii) {
       aux[ii] = y[ii] + k1[ii]/2;
     }
     //k2
     for(int ii = 0; ii < y.size(); ++ii) {
       k2[ii] = h*f(t + h/2, aux, ii);
     }
     // k3 aux
     for(int ii = 0; ii < y.size(); ++ii) {
       aux[ii] = y[ii] + k2[ii]/2;
     }
     //k3
     for(int ii = 0; ii < y.size(); ++ii) {
       k3[ii] = h*f(t + h/2, aux, ii);
     }
     // k4 aux
     for(int ii = 0; ii < y.size(); ++ii) {
       aux[ii] = y[ii] + k3[ii];
     }
     //k4
     for(int ii = 0; ii < y.size(); ++ii) {
       k4[ii] = h*f(t + h, aux, ii);
     }
     // write new y
     for(int ii = 0; ii < y.size(); ++ii) {
       y[ii] = y[ii] + (k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii])/6.0;
     }
     std::cout << "\t" << y[0] << "\t" << y[1] << "\t" << y[2]<<std::endl;
   }
 }
 double f(double t, const std::vector<double> & y, int id)
 {
   if (0 == id) {
     return y[3]+O*y[1];
   }
   else if (1 == id) {
     return y[4] -O*y[0];
   }
   else if (2 == id){
     return y[5];
   }
   else if (3 == id){
     return (O*y[4])-A*(2*y[0]/(a*a))*(1/(C+((y[0]*y[0])/(a*a))+((y[1]*y[1])/(b*b))+((y[2]*y[2])/(c*c))))*std::log10(std::exp(1));
   }
   else if (4 == id){
     return -(O*y[3])-A*(2*y[1]/(b*b))*(1/(C+((y[0]*y[0])/(a*a))+((y[1]*y[1])/(b*b))+((y[2]*y[2])/(c*c))))*std::log10(std::exp(1));
   }
   else if (5 == id){
     return -A*(2*y[2]/(c*c))*(1/(C+((y[0]*y[0])/(a*a))+((y[1]*y[1])/(b*b))+((y[2]*y[2])/(c*c))))*std::log10(std::exp(1));
   }
 }
