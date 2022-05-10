// Filename: ml.cpp
// Version:  Standalone C++ Version
// Date:     May 10, 2022
//
// Author:   Martin Moehle
//           Mathematical Institute
//           University of Tuebingen
//           Germany
//           martin.moehle@uni-tuebingen.de
//
// Random number generator for the Mittag-Leffler distribution
//
// Please compile this file using the command "g++ ml.cpp -o ml.exe"
// or "g++ ml.cpp -std=c++11 -o ml.exe".
//
// This C++ program demonstrates random number generation for the (second
// type) two-parameter Mittag-Leffler distribution ML(a,b). For the
// mathematical background we refer the reader to (Section 8 of)
//
// Moehle, M.: A restaurant process with cocktail bar and relations
// to the three-parameter Mittag-Leffler distribution, J. Appl. Probab. 58,
// 978-1006 (2021).
//
// In contrast to the standalone ANSI C Version "mitlef.c" this C++
// program "ml.cpp" uses the C++11 random number generation component
// <random> in order to generate pseudo random numbers.

#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

double PI = 3.14159265358979323846;

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
// V = V(a,b), 0<a<1, b>0, having density v |-> c(a,b)/(h(a,v))^(b*(1-a)/a),
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
   // rejection Algorithm 1 of Moehle (2021)
   do {
      u=uniform_rand();    // sample u uniformly on (0,1)
      v=PI*uniform_rand(); // sample v uniformly on (0,PI)
   } while (u>exp(p*log(h0/h(a,v))));
   return v;
} // end v_rand_1

double v_rand_2(double a, double b) {
// same as v_rand_1, but Algorithm 2 in Section 2 of Moehle (2021) is used.
// v_rand_2 is more efficient than v_rand_1 in particular for large values
// of b.
   double u,  // local variable for uniform random number
          v,  // variable for random number of V
          h0, // variable for h(a,0+) = (1-a)*a^(a/(1-a))
          p,  // variable for b*(1-a)/a
          t;  // variable for (cos(v/2))^(2*b) * (h(a,v)/h(a,0+))^p
   if (a>0.5) { b=b*(1.0-a)/a; a=1.0-a; } // transformation to a <= 0.5
   p=b*(1.0-a)/a;
   h0=(1.0-a)*exp(a/(1.0-a)*log(a));
   // rejection Algorithm 2 of Moehle (2021)
   do {
      u=uniform_rand(); // sample u uniformly on (0,1)
      // sample v from V(1/2,b) via Eq. (45) of Moehle (2021)
      v=2.0*acos(sqrt(beta_rand(b+0.5,0.5)));
      t=exp(2.0*b*log(cos(v/2.0))+p*log(h(a,v)/h0));
   } while (u*t>1.0);
   return v;
} // end v_rand_2

double mittag_leffler_type2_two_parameter_rand(double a, double b) {
// Let L be second type two-parameter Mittag-Leffler distributed
// with parameters 0<a<1 and b>0. This function is based on the fact
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
   v=v_rand_2(a,b); // you may replace this by v_rand_1(a,b), but note
   // that v_rand_1 is slower than v_rand_2 for large values of b.
   // return the value (g/h(a,v))^(1-a)
   return exp((1.0-a)*log(g/h(a,v)));
} // end mittag_leffler_type2_two_parameter_rand

int main(int argc, char *argv[]) {
   int n,              // number of pseudo random numbers to be generated
       i;              // local index
   double a,b,         // parameters of the ML(a,b) distribution
          x,           // variable for pseudo random numbers
          m,m2,v,      // mean, second moment and variance of ML(a,b)
          sm,sm2,svar; // sample mean ...

   // time-based seed
   unsigned timeseed
   = std::chrono::system_clock::now().time_since_epoch().count();

//   std::cout << "timeseed = " << timeseed << std::endl;

   if (argc != 4) {
      printf("USAGE: ml a b n\n\n");
      printf("where 0 < a < 1, b >= 0 and n is a positive integer.\n\n");
      printf("Example: ml 0.5 1 50");
      exit(0);
   }
   a=atof(argv[1]);
   b=atof(argv[2]);
   n=atoi(argv[3]);
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
   std::cout << "Generating n = " << n << " pseudo random numbers for ML(a,b):"
   << std::endl;
   sm=0.0;sm2=0.0;
   for (i=1;i<=n;++i) {
      x=mittag_leffler_type2_two_parameter_rand(a,b);
      sm+=x; sm2+=x*x;
      std::cout << std::setprecision(3) << std::fixed << x << ' ';
   }
   std::cout << std::endl;
   sm=sm/n; sm2=sm2/n;
   svar=sm2-sm*sm;
   std::cout << "Sample mean:          " << sm << std::endl;
   std::cout << "Second sample moment: " << sm2 << std::endl;
   std::cout << "Sample variance:      " << svar << std::endl;
   return 0;
} // end main
