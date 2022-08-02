// Filename: ml_v2.cpp
// Version:  Standalone C++ Version 2.0
// Date:     August 2, 2022
//
// Author:   Martin Moehle
//           Mathematical Institute
//           University of Tuebingen
//           Germany
//           martin.moehle@uni-tuebingen.de
//
// Random number generator for the Mittag-Leffler distribution
//
// Please compile this file using the command "g++ ml_v2.cpp -o ml.exe"
// or "g++ ml_v2.cpp -std=c++11 -o ml.exe".
//
// This C++ program demonstrates random number generation for the (second
// type) two-parameter Mittag-Leffler distribution ML(a,b). For the
// mathematical background we refer the reader to (Section 8 of)
//
// Moehle, M.: A restaurant process with cocktail bar and relations
// to the three-parameter Mittag-Leffler distribution, J. Appl. Probab. 58,
// 978-1006 (2021).
//
// In contrast to the standalone ANSI C Version "mitlef_v2.c" this C++
// program "ml_v2.cpp" uses the C++11 random number generation component
// <random> in order to generate pseudo random numbers.
//
// Version Date             Comments, changes etc.
// ----------------------------------------------------------------------------
// 1.0     May 10, 2022     Version which works for parameters 0 < a < 1 and
//                          b >= 0.
//
// 2.0     August 2, 2022   Program now works for parameters 0 < a < 1 and
//                          b > -a. The new functions "T", "S", "v_rand_3"
//                          and "v_rand_4" cover the case of negative b (>-a).
//                          Moreover, function "loggamma" and functions for
//                          the mean number of steps of the four algorithms
//                          and sample steps added.

#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

// global variables of the program
unsigned int steps; // counter for the number of steps of the used rejection method

double PI = 3.141592653589793238462643;

// coefficients for the gamma function approximation

double cc[] = {1.000000000190015, 76.18009172947146, -86.50532032941677,
              24.01409824083091, -1.231739572450155, 0.001208650973866179,
	          -0.000005395239384953};

double gamma(double x) {
// gamma function appproximation via C. Lanczos, A precision approximation of
// the gamma function, SIAM Journal on Numerical Analysis, Series B, Vol. 1,
// 1964, pp. 86-96.
   double h;
   int i;
   if (x<0.5) return PI/(sin(PI*x)*gamma(1.0-x)); // reflection formula
   h=cc[0]; for (i=5;i>=0;i--) h+=cc[i+1]/(x+i);
   return sqrt(2.0*PI)*exp((x-0.5)*log(x+4.5)-(x+4.5))*h;
} // end gamma

double loggamma(double x) {
// logarithm of the gamma function for x > 0
   double h;
   int i;
   if (x<0.5) return log(PI)-log(sin(PI*x))-loggamma(1.0-x); // reflection formula
   h=cc[0]; for (i=5;i>=0;i--) h+=cc[i+1]/(x+i);
   return 0.5*log(2.0*PI)+log(h)+(x-0.5)*log(x+4.5)-(x+4.5);
} // end loggamma

// gobal_urng() shares a single URNG with all other functions in the toolkit

std::default_random_engine & global_urng() {
   static std::default_random_engine u{};
   return u;
}

void randomize(unsigned timeseed) {
// sets the shared URNG to an unpredictable state.
//
//   old version:
//     static std::random_device rd{};
//     global_urng().seed(rd());
//
   global_urng().seed(timeseed);
} // end randomize

double uniform_rand(void) {
// returns a double variate uniformly distributed in [0,1)
   static std::uniform_real_distribution<> d{};
   using parm_t = decltype(d)::param_type;
   return d( global_urng(), parm_t{0.0,1.0});
} // end uniform_rand

double exponential_rand(double a) {
// returns a double variate exponentially distributed with parameter a
   static std::exponential_distribution<> d{};
   using parm_t = decltype(d)::param_type;
   return d( global_urng(), parm_t{a});
} // end exponential_rand

double gamma_rand(double a) {
// returns a double variate gamma distributed with parameter a and 1
// corresponding to the density x^(a-1) exp(-x) / gamma(a)
   static std::gamma_distribution<> d{};
   using parm_t = decltype(d)::param_type;
   return d( global_urng(), parm_t{a,1.0});
} // end gamma_rand

double beta_rand(double p, double q) {
// returns a double variate beta distributed with parameter p and q
   double x,y;
   x=gamma_rand(p);
   y=gamma_rand(q);
   return x/(x+y);
} // end beta_rand

double h(double a, double v) {
// This auxiliary function is used by the following procedures v_rand_1,
// v_rand_2 and mittag_leffler_type2_two_parameter_rand. For more details
// see Eq. (36) in Moehle (2021).
   return sin((1.0-a)*v)/sin(a*v)*exp(log(sin(a*v)/sin(v))/(1.0-a));
} // end h

double v_rand_1(double a, double b) {
// Pseudo random number generator for the (0,PI)-valued random variable
// V = V(a,b), 0 < a < 1, b >= 0, having density v |-> c(a,b)/(h(a,v))^(b*(1-a)/a),
// with normalizing constant
//
// c(a,b) := gamma(b+1)*gamma(b*(1-a)/a+1)/PI/gamma(b/a+1).
//
// Algorithm 1 of Section 8 in Moehle (2021) is used. The function v_rand_1
// can be alternatively used (instead of v_rand_2) by the function
// mittag_leffler_type2_two_parameter_rand.

   double u,  // local variable for uniform random number
          v,  // variable for random number of V
          h0, // variable for h(a,0+) = (1-a)*a^(a/(1-a))
          p;  // variable for b*(1-a)/a
   p=b*(1.0-a)/a;
   h0=(1.0-a)*exp(a/(1.0-a)*log(a));
   steps=0; // initialize the global variable steps
   // rejection Algorithm 1 of Moehle (2021)
   do {
      u=uniform_rand();    // sample u uniformly on (0,1)
      v=PI*uniform_rand(); // sample v uniformly on (0,PI)
      steps+=1; // count this step
   } while (u>exp(p*log(h0/h(a,v))));
   return v;
} // end v_rand_1

double v_rand_2(double a, double b) {
// same as v_rand_1, but Algorithm 2 in Section 2 of Moehle (2021) is used.
// v_rand_2 is efficient (uniformly for all parameters 0 < a < 1 and b >= 0)
// and (for most parameter values) more efficient than v_rand_1, in
// particular for large values of b.
   double u,  // local variable for uniform random number
          v,  // variable for random number of V
          h0, // variable for h(a,0+) = (1-a)*a^(a/(1-a))
          p,  // variable for b*(1-a)/a
          t;  // variable for (cos(v/2))^(2*b) * (h(a,v)/h(a,0+))^p
   if (a>0.5) { b=b*(1.0-a)/a; a=1.0-a; } // transformation to a <= 0.5
   p=b*(1.0-a)/a;
   h0=(1.0-a)*exp(a/(1.0-a)*log(a));
   steps=0; // initialize the global variable steps
   // rejection Algorithm 2 of Moehle (2021)
   do {
      u=uniform_rand(); // sample u uniformly on (0,1)
      // sample v from V(1/2,b) via Eq. (45) of Moehle (2021)
      v=2.0*acos(sqrt(beta_rand(b+0.5,0.5)));
      t=exp(2.0*b*log(cos(v/2.0))+p*log(h(a,v)/h0));
      steps+=1; // count this step
   } while (u*t>1.0);
   return v;
} // end v_rand_2

double T(double a, double b, double v) {
// The transformation T(v) := (1-v/PI)^(b/a+1), 0 < v < PI.
   return exp((b/a+1.0)*log(1.0-v/PI));
} // end T

double S(double a, double b, double w) {
// The inverse transformation S := T^{-1} of T is given by
// S(w) = PI(1-w^(a/(a+b))), 0 < w < 1. This function is used by
// the following function "v_rand_3".
   return PI*(1.0-exp(a/(a+b)*log(w)));   
} // end S
double v_rand_3(double a, double b) {
// Algorithm to sample from V if the parameter b (>-a) is negative. Note that
// the corresponding rejection method is not documented in Moehle (2021). For
// negative parameter b, the density of V is unbounded. Thus, the transformed
// random variable W := T(V) with transformation T(v) := (1-v/PI)^(b/a+1) is
// considered. The transformation T is chosen such that the density of W is
// bounded. Thus, a standard rejection algorithm is applicable to sample from W.
// The back transformed value v:=S(w), where S denotes the inverse of T, is
// then a sample from V. This algorithm theoretically works for any 0 < a < 1
// and -a < b < 0, however it is not very efficient (and also numerically unstable)
// for (a,b) close to (0,0) or (1,-1) due to the fact that in this region, V
// is very concentrated around PI and hence the density of V is highly peaked.
   const double eps=1e-8; // a small number greater than the maschine exactness
   double u,   // local variable for uniform random number
          v,   // variable for random number of V := S(a,b,w)
          w,   // variable for random number of W
          h0,  // variable for h(a,0+) = (1-a)*a^(a/(1-a))
          t,   // variable for (h0/h(a,v))^(b/a(1-a)) w^(-b/(a+b))
          C;   // variable for (sin(a*PI))^(1/(1-a))*sin((1-a)*PI)/sin(a*PI)
   C=exp(log(sin(a*PI))/(1.0-a))*sin((1.0-a)*PI)/sin(a*PI);
   h0=(1.0-a)*exp(a/(1.0-a)*log(a));
   steps=0; // initialize the global variable steps
   do {
      u=uniform_rand(); // sample u and w
      w=uniform_rand(); // uniformly on (0,1)
      v=S(a,b,w);       // compute the back transformed v := S(a,b,w) 
      // compute t = t(w)
      if (v<PI-eps)
         // use the usual formula for t = t(w) if v is not too close to PI
         t=exp(b/a*(1.0-a)*log(h0/h(a,v))-b/(a+b)*log(w));
      else // v is close to PI
         // use t(0+) instead of t(w) to avoid numerical problems
         t=exp(b/a*log(PI)+b*(1.0-a)/a*log(h0/C));
      steps+=1; // count this step
   } while (u>t);
   return v; // return the back transformed value v = S(a,b,w)
} // end v_rand_3

double v_rand_4(double a, double b) {
// This algorithm again covers the case of negative b (>-a), but "v_rand_4"
// should be prefered compared to "v_rand_3", since it is more efficient. This
// algorithm uses the same transformation T as v_rand_3, but in addition,
// similar as in "v_rand_2" compared to "v_rand_1", a head function with
// a = 1/2 is used. The method theoretically works for all parameter values
// 0 < a < 1  and -a < b < 0. The efficiency is not very good and the
// algorithm becomes numerically unstable if a is close to 0 or 1 and
// b is close to -a.
   const double eps=1e-8; // a small number greater than the maschine exactness
   double u,  // local variable for uniform random number
          v,  // variable for random number of V := S(a,b,w)
          w,  // variable for random number of W
          h0, // variable for h(a,0+) = (1-a)*a^(a/(1-a))
          p,  // variable for b*(1-a)/a
          t,  // variable for (cos(S(0.5,b,w)/2))^(2*b) * (h(a,v)/h0)^p * w^(b/(a+b)-b/(0.5+b))
          C;  // variable for (sin(a*PI))^(1/(1-a))*sin((1-a)*PI)/sin(a*PI)
   if (a>0.5) { b*=(1.0-a)/a; a=1.0-a; } // transformation to a <= 0.5
   C=exp(log(sin(a*PI))/(1.0-a))*sin((1.0-a)*PI)/sin(a*PI);
   p=b*(1.0-a)/a;
   h0=(1.0-a)*exp(a/(1.0-a)*log(a));
   steps=0; // initialize the global variable steps
   do {
      u=uniform_rand();
      // sample w from W = T(1/2,b,V(1/2,b)) via (45) of Moehle (2021)
      w=T(0.5,b,2.0*acos(sqrt(beta_rand(b+0.5,0.5))));
      v=S(a,b,w); // compute the back transformed v := S(a,b,w)
         if (v<PI-eps)
         // use the usual formula for t = t(w) if v is not too close to PI
         t=exp(
           2.0*b*log(cos(S(0.5,b,w)/2.0))
           +p*log(h(a,v)/h0)
           +(b/(a+b)-b/(b+0.5))*log(w)
         );
      else // v is close to PI
         // use t(0+) instead of t(w) to avoid numerical problems
         t=exp(2.0*b*log(PI/2.0)+p*log(C/h0)-b/a*log(PI));
      steps+=1; // count this step
   } while (u*t>1.0);
   return v; // return the back transformed value v = S(a,b,w)
} // end v_rand_4

double mittag_leffler_type2_two_parameter_rand(double a, double b) {
// Let L be second type two-parameter Mittag-Leffler distributed
// with parameters 0 < a < 1 and b > 0. This function is based on the fact
// that L has the same distribution as (G/h(V))^(1-a), where G is gamma
// distributed with parameters b*(1-a)/a+1 and 1 having density
// x |-> x^(b*(1-a)/a)*exp(-x)/gamma(b*(1-a)/a+1), the function h(a,.) is
// defined via
//
// h(a,u) := (sin(a*u)/sin(u))^(1/(1-a)) * sin((1-a)*u)/sin(a*u)
//
// and the (0,PI)-valued random variable V is independent of G and has
// density v |-> c(a,b)/(h(a,v))^(b*(1-a)/a), with normalizing constant
//
// c(a,b) := gamma(b+1)*gamma(b*(1-a)/a+1)/PI/gamma(b/a+1).
//
// Pseudo random numbers of V are generated via the rejection
// method (Algorithm 2) (or alternatively Algorithm 1) described
// in Section 8 of Moehle (2021).

   double g,  // variable for random numbers of G (gamma distributed)
          v;  // variable for random number of V
   // sample g gamma distributed with parameters b*(1-a)/a+1 and 1
   g=gamma_rand(b*(1.0-a)/a+1.0);
   // sample v from the random variable V(a,b)
   if (b<0) v=v_rand_4(a,b); // particular algorithm for negative b; you may
   // replace this by "v_rand_3(a,b)", but note that "v_rand_3" is never better
   // than "v_rand_4".
   else // b >=0
   v=v_rand_2(a,b); // you may replace this by "v_rand_1(a,b)", but note
   // that "v_rand_1" is not very efficient for large values of beta and
   // never better than "v_rand_2".
   // return the value (g/h(a,v))^(1-a)
   return exp((1.0-a)*log(g/h(a,v)));
} // end mittag_leffler_type2_two_parameter_rand

double mean_1(double a, double b) {
// Expected number E(N) of steps of the rejection function v_rand_1;
//   see the formula in Remark 24 of Moehle (2021).
   double h; // local variable
   h=loggamma(b+1.0)+loggamma(b*(1.0-a)/a+1.0)-loggamma(b/a+1.0)
     -b*log(a)-b*(1.0-a)/a*log(1.0-a);
   return exp(h);
} // end mean_1

double mean_2(double a, double b) {
// Expected number E(N) of steps of the rejection function v_rand_2;
// see the formula for d(a,b) at the bottom of page 1004 of Moehle (2021)
   double h; // local variable
   h=loggamma(b+0.5)+loggamma(b/a*(1.0-a)+1.0)-0.5*log(PI)-loggamma(b/a+1.0)   
     -b*log(a)-b*(1.0-a)/a*log(1.0-a);
   return exp(h);
} // end mean_2

double mean_3(double a, double b) {
// Expected number E(N) of steps of the rejection function v_rand_3.
// This function is not documented in Moehle (2021). The formula is
//
//          gamma(b+1) * gamma(b(1-a)/a+1)
// E(N) = -------------------------------------.
//        gamma(b/a+2) * a^b * (1-a)^(b(1-a)/a)
//
   double h;
   h=loggamma(b+1)+loggamma(b*(1.0-a)/a+1.0)-loggamma(b/a+2.0)
     -b*log(a)-b*(1.0-a)/a*log(1.0-a);
   return exp(h);
} // end mean_3

double mean_4(double a, double b) {
// Expected number E(N) of steps of the rejection function v_rand_4.
// This function is not documented in Moehle (2021). The formula is
//
//        2 * gamma(b+3/2) * gamma(b(1-a)/a+1)
// E(N) = ------------------------------------------------.
//        sqrt(PI) * gamma(b/a+2) * a^b * (1-a)^(b(1-a)/a)
   double h;   
   h=log(2.0)+loggamma(b+1.5)+loggamma(b*(1.0-a)/a+1.0)-0.5*log(PI)
   -loggamma(b/a+2.0)-b*log(a)-b*(1.0-a)/a*log(1.0-a);
   return exp(h);   
} // end mean_4

int main(int argc, char *argv[]) {
   int n,              // number of pseudo random numbers to be generated
       i;              // local index
   double a,b,         // parameters of the ML(a,b) distribution
          x,           // variable for pseudo random numbers
          m,m2,v,      // mean, second moment and variance of ML(a,b)
          sm,sm2,svar, // sample mean ...
          ssteps;      // sample number of steps

   // time-based seed
   unsigned timeseed
   = std::chrono::system_clock::now().time_since_epoch().count();

//   std::cout << "timeseed = " << timeseed << std::endl;

   if (argc != 4) {
      std::cout << "This is ml, Version 2, August 2, 2022."
                << std::endl << std::endl;
      std::cout << "USAGE: ml a b n" << std::endl << std::endl;
      std::cout << "where 0 < a < 1, b > -a and n is a positive integer."
                << std::endl << std::endl;
      std::cout << "Example: ml 0.5 1 50";
      exit(0);
   }
   a=atof(argv[1]);
   b=atof(argv[2]);
   n=atoi(argv[3]);
   if ((a<=0)||(a>=1)) {
      std::cout << "Parameter a should satisfy 0 < a < 1." << std::endl;
      std::cout << "Program aborted." << std::endl;
      exit(1); }
   if (b<=-a) {
      std::cout << "Parameter b should satisfy b > -a." << std::endl;
      std::cout << "Program aborted." << std::endl;
      exit(1); }
   if (n<=0) {
      std::cout << "Parameter n should be a positive integer." << std::endl;
      std::cout << "Program aborted." << std::endl;
      exit(1); }
   std::cout << std::setprecision(3) << std::fixed;
   std::cout << "Second type two-parameter Mittag-Leffler distribution ML(a,b)"
   << std::endl;
   std::cout << "Parameters: a = " << a << " b = " << b << std::endl;
   m=gamma(b+1.0)/a/gamma(a+b);              // mean of ML(a,b)
   m2=(a+b)*gamma(b+1.0)/a/a/gamma(2.0*a+b); // second moment of ML(a,b)
   v=m2-m*m;                                 // variance of ML(a,b)
   std::cout << "Mean:          " << m << std::endl;
   std::cout << "Second moment: " << m2 << std::endl;
   std::cout << "Variance:      " << v << std::endl;
   // randomized initialization
   randomize(timeseed);
   std::cout << "Generating n = " << n << " samples from ML(a,b) via ";
   if (b>=0) { 
      std::cout << "Algorithm 2." << std::endl; 
      std::cout << "Expected number of steps of rejection Algorithm 2: E(N) = ";
      if (a<=0.5) std::cout << mean_2(a,b); else std::cout << mean_2(1.0-a,b*(1.0-a)/a);
      std::cout << std::endl;
      std::cout << "(Expected number of steps of rejection Algorithm 1: E(N) = "
      << mean_1(a,b) << ")" << std::endl;
   }
   else { // b < 0
      std::cout << "Algorithm 4, since b < 0." << std::endl;
      std::cout << "Expected number of steps of rejection Algorithm 4: E(N) = ";
      if (a<=0.5) std::cout << mean_4(a,b); else std::cout << mean_4(1.0-a,b*(1.0-a)/a);
      std::cout << std::endl;
      std::cout << "(Expected number of steps of rejection Algorithm 3: E(N) = "
      << mean_3(a,b) << ")" << std::endl;
   }
   sm=0.0;sm2=0.0;ssteps=0;
   for (i=1;i<=n;++i) {
      x=mittag_leffler_type2_two_parameter_rand(a,b);
      sm+=x; sm2+=x*x;ssteps+=steps;
      std::cout << std::setprecision(3) << std::fixed << x << ' ';
   }
   std::cout << std::endl;
   sm=sm/n; sm2=sm2/n;
   svar=sm2-sm*sm;
   ssteps=ssteps/n;
   std::cout << "Sample mean:            " << sm << std::endl;
   std::cout << "Second sample moment:   " << sm2 << std::endl;
   std::cout << "Sample variance:        " << svar << std::endl;
   std::cout << "Sample number of steps: " << ssteps << " (compare with E(N))" << std::endl;
   return 0;
} // end main
