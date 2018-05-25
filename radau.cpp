#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

void rk4(double ta, double tb, double h, std::vector<double> & y);
double f(double t, const std::vector<double> & y, int id);
//void euler(double ta, double tb, double h, std::vector<double> & y);
void radau55(double ta, double tb, double h, std::vector<double> & y);
void initialize(std::vector<double> & a);

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
  const double TB = 10;
  const double H  = 0.1;
  std::vector<double> y = {2.5, 0, 0, 0, 2.355, 0}; // {x0, y0, z0, vx0, vy0, vz0}  ---1.707---
    rk4(TA,TB, H, y);
    std::cout << "YAAAAAAAAAAAAAAAAAAAAAAAAAAA" <<std::endl;
   // euler(TA,TB, H, y);
   radau55(TA,TB,H,y);
   return 0;
 }
void initialize(std::vector<double> & a)
{
  for(int ii = 0; ii < a.size(); ++ii){
    a[ii] = 0.0;
  }
}

 void radau55(double ta, double tb, double h, std::vector<double> & y)
 {
   const double b2[3] ={(16.0-std::sqrt(6.0))/36.0,(16.0+std::sqrt(6.0))/36.0, 1.0/9.0};
   const double c2[3] ={(4.0-std::sqrt(6.0))/10.0, (4.0+std::sqrt(6.0))/10.0, 1.0};
   const double a2[9] ={/*11 */(88.0-(7.0*std::sqrt(6.0)))/360.0,/* 12  */(296.0-(169.0*std::sqrt(6.0)))/1800.0,/* 13 */(-2.0+(3.0*std::sqrt(6.0)))/225.0,/*  21  */(296+(196.0*std::sqrt(6.0)))/1800,/*  22 */(88.0+(7.0*std::sqrt(6.0)))/360.0,/*  23  */(-2.0-(3.0*std::sqrt(6.0)))/225.0,/* 31   */(16.0-(std::sqrt(6.0)))/36.0,/* 32*/(16.0+(std::sqrt(6.0)))/36.0, /*  33  */1.0/9.0};
   
   const int N = (tb-ta)/h;
   for (int nt = 0; nt < N; ++nt) {
     double t = ta + h*nt;
     double sumatoriakj =0;
     double sumatoriaki =0;
     std::vector<double>  valorki(3);
     std::vector<double> aux(3);
     for(int i=0; i<3; i++)
       {
        aux[i]= y[i]+h*sumatoriakj;
       }
     for(int i=0; i<3; i++){
        valorki[i]=f(t+(c2[i]*h),aux,i);
     }
     for(int i=0; i<3; i++){
       sumatoriaki+=b2[i]*valorki[i];
     }
     for(int i=0; i<3;i++){
       for(int j=0; j<3;j++){
	 sumatoriakj+=a2[j+3*i]*valorki[j];
       }
     }
     
     for(int ii = 0; ii < 3; ++ii) {
       y[ii] =  y[ii]+ h*sumatoriaki;
     }
     
     std::cout << "\t" << y[0] << "\t" << y[1] << "\t" << y[2]<<std::endl;
   }
 }



/*void euler(double ta, double tb, double h, std::vector<double> & y)
  {
  const int N = (tb-ta)/h;
  for (int nt = 0; nt < N; ++nt) {
  double t = ta + h*nt;
  for(int ii = 0; ii < y.size(); ++ii) {
  y[ii] += h*f(t,y,ii);
  }
  std::cout << "\t" << y[0] << "\t" << y[1] << "\t" << y[2]<<std::endl;
  }
  }*/
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
    std::cout << "\t" << y[0] << "\t" << y[1] << "\t" << y[2]<< "\t" << y[3] << "\t" << y[4] << "\t" << y[5]<<std::endl;
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
