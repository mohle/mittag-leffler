/* Filename: mitlef.c
   Version:  Standalone ANSI C Version
   Date:     May 22, 2022

   Author:   Martin Moehle
             Mathematical Institute
             University of Tuebingen
             Germany
             martin.moehle@uni-tuebingen.de

   Please compile this file using the command "gcc mitlef.c -o mitlef.exe".

   This C program demonstrates random number generation for the (second
   type) two-parameter Mittag-Leffler distribution ML(a,b). For the
   mathematical background we refer the reader to (Section 8 of)

   Moehle, M.: A restaurant process with cocktail bar and relations
   to the three-parameter Mittag-Leffler distribution, J. Appl. Probab. 58,
   978-1006 (2021).

   For an analog program written in C++11 see "ml.cpp". */

#include <stdlib.h> /* RAND_MAX, rand, srand, ... */
#include <stdio.h>  /* printf, ... */
#include <time.h>   /* time, ... */
#include <math.h>   /* exp, log, ... */

#define PI 3.14159265358979323846

/* coefficients for the gamma function approximation */

double cc[] = {1.000000000190015, 76.18009172947146, -86.50532032941677,
              24.01409824083091, -1.231739572450155, 0.001208650973866179,
	      -0.000005395239384953};

double gamma(double x) {
/* gamma function approximation via C. Lanczos, A precision approximation of
the gamma function, SIAM Journal on Numerical Analysis, Series B, Vol. 1,
1964, pp. 86-96. */
   double h;
   int i;
   if (x<0.5) return PI/(sin(PI*x)*gamma(1-0-x)); /* reflection formula */
   h=cc[0]; for (i=5;i>=0;i--) h+=cc[i+1]/(x+i);
   return sqrt(2.0*PI)*exp((x-0.5)*log(x+4.5)-(x+4.5))*h;
} /* end gamma */

/* Remark: For simplicity this program uses the (historical and nowadays
outdated) C standard functions to generate pseudo random numbers uniformly
distributed on the open interval (0,1). The program code has the advantage
that it compiles with any C standard compiler. It should work standalone,
so no extra include files are needed. In serious applications the outdated
functions my_randomize() and uniform_rand() should be implemented using a
suitable modern random number generator. For academic demonstration however
the following function uniform_rand() should be sufficient. */

void my_randomize(void) {
/* initialize the random number generator */
   srand((unsigned) time(NULL));
} /* my_randomize */

double uniform_rand(void) {
/* simple pseudo uniform random number in the open unit interval (0,1) */
   return (rand()+1.0) / ((double)RAND_MAX+2.0);
} /* end uniform_rand */

/* The following four functions are based on some standard methods to
generate exponential, Erlang, gamma and beta-distributed pseudo random
numbers. There exist possibly better ways to generate such random numbers,
but the code of the provided functions is short and useful for academic
demonstration. */

double exponential_rand(double b) {
/* inverse transform method */
   return -b*log(uniform_rand());
} /* end exponential_rand */

double erlang_rand(double b, int c) {
/* Choose U_1,...,U_c independent and all uniformly distributed on (0,1).
Then X := -b*log(U_1...U_c) is Erlang distributed with parameter c */
   double h;
   int i;
   h=1.0; for(i=1;i<=c;i++) h*=uniform_rand();
   return -b*log(h);
} /* end erlang_rand */

double gamma_rand(double b, double c) {
/* If c<1 then choose V and W independent with V exponentially
distributed with parameter b and W beta distributed with parameters
c and 1-c. Then VW is gamma distributed with parameters b and c.

Assume now that c >= 1. Then decompose c=m+delta in its integer
part m and the fractional part delta. Then generate a gamma(b,delta)
random variable X and take an independent random variable Y being
Erlang(b,m) distributed. Then X+Y is gamma(b,c)-distributed. */
   double v,w,y1,y2,
          delta,  /* variable for the fractional part of c */
          m;      /* variable for the integer part of c */
   /* compute the fractional part and the integer part of c */
   delta=modf(c,&m); /* delta := fractional part of c
                        m := integer part of c */
   do {
      if (delta>0.0) { y1=exp(log(uniform_rand())/delta); } else { y1=0.0; }
      y2=exp(log(uniform_rand())/(1.0-delta));
   }
   while (y1+y2>1.0);
   w=y1/(y1+y2); /* w := beta_rand(delta,1-delta) */
   v=exponential_rand(b);
   return v*w+erlang_rand(b,m);
} /* end gamma_rand */

double beta_rand(double p, double q) {
/* Choose X and Y independent with X gamma distributed with parameters
1 and p and Y gamma distributed with parameters 1 and q. Then X/(X+Y)
is beta distributed with parameters p and q. */
   double x,y;
   x=gamma_rand(1.0,p);
   y=gamma_rand(1.0,q);
   return x/(x+y);
} /* end beta_rand */

/* In the following the specific functions for the generation of (second type)
two-parameter Mittag-Leffler distributed pseudo random numbers are provided. */

double h(double a, double v) {
/* This auxiliary function is used by the following procedures v_rand_1,
   v_rand_2 and mittag_leffler_type2_two_parameter_rand. For more details
   see Eq. (36) in Moehle (2021). */
   return sin((1.0-a)*v)/sin(a*v)*exp(log(sin(a*v)/sin(v))/(1.0-a));
} /* end h */

double v_rand_1(double a, double b) {
/* Pseudo random number generator for the (0,PI)-valued random variable
V = V(a,b), 0<a<1, b>0, having density v |-> c(a,b)/(h(a,v))^(b*(1-a)/a),
with normalizing constant

c(a,b) := gamma(b+1)*gamma(b*(1-a)/a+1)/PI/gamma(b/a+1).

Algorithm 1 of Section 8 in Moehle (2021) is used. The function v_rand_1
can be alternatively used (instead of v_rand_2) by the function
mittag_leffler_type2_two_parameter_rand. */

   double u,  /* local variable for uniform random number */
          v,  /* variable for random number of V */
          h0, /* variable for h(a,0+) = (1-a)*a^(a/(1-a)) */
          p;  /* variable for b*(1-a)/a */
   p=b*(1.0-a)/a;
   h0=(1.0-a)*exp(a/(1.0-a)*log(a));
   /* rejection Algorithm 1 of Moehle (2021) */
   do {
      u=uniform_rand();    /* sample u uniformly on (0,1) */
      v=PI*uniform_rand(); /* sample v uniformly on (0,PI) */
   } while (u>exp(p*log(h0/h(a,v))));
   return v;
} /* end v_rand_1 */

double v_rand_2(double a, double b) {
/* same as v_rand_1, but Algorithm 2 in Section 8 of Moehle (2021) is used.
v_rand_2 is efficient (uniformly for all parameters 0<a<1 and b > 0 and
(for most parameter values) more efficient than v_rand_1 in particular
for large values of b */
   double u,  /* local variable for uniform random number */
          v,  /* variable for random number of V */
          h0, /* variable for h(a,0+) = (1-a)*a^(a/(1-a)) */
          p,  /* variable for b*(1-a)/a */
          t;  /* variable for (cos(v/2))^(2*b) * (h(a,v)/h(a,0+))^p */
   if (a>0.5) { b*=(1.0-a)/a; a=1.0-a; } /* transformation to a <= 0.5 */
   p=b*(1.0-a)/a;
   h0=(1.0-a)*exp(a/(1.0-a)*log(a));
   /* rejection Algorithm 2 of Moehle (2021) */
   do {
      u=uniform_rand(); /* sample u uniformly on (0,1) */
      /* sample v from V(1/2,b) via (45) of Moehle (2021) */
      v=2.0*acos(sqrt(beta_rand(b+0.5,0.5)));
      t=exp(2.0*b*log(cos(v/2.0))+p*log(h(a,v)/h0));
   } while (u*t>1.0);
   return v;
} /* end v_rand_2 */

double mittag_leffler_type2_two_parameter_rand(double a, double b) {
/* Let L be second type two-parameter Mittag-Leffler distributed
with parameters 0<a<1 and b>0. This function is based on the fact
that L has the same distribution as (G/h(V))^(1-a), where G is gamma
distributed with parameters 1 and b*(1-a)/a+1 having density
x |-> x^(b*(1-a)/a)*exp(-x)/gamma(b*(1-a)/a+1), the function h(a,.) is
defined via

h(a,u) := (sin(a*u)/sin(u))^(1/(1-a)) * sin((1-a)*u)/sin(a*u)

and the (0,PI)-valued random variable V is independent of G and has
density v |-> c(a,b)/(h(a,v))^(b*(1-a)/a), with normalizing constant

c(a,b) := gamma(b+1)*gamma(b*(1-a)/a+1)/PI/gamma(b/a+1).

Pseudo random numbers of V are generated via the rejection
method (Algorithm 2) (or alternatively Algorithm 1) described
in Section 8 of Moehle (2021). */

   double g,  /* variable for random number of G (gamma distributed) */
          v;  /* variable for random number of V */
   /* sample g gamma distributed with parameters 1 and b*(1-a)/a+1 */
   g=gamma_rand(1.0,b*(1.0-a)/a+1.0);
   /* sample v from the random variable V(a,b) */
   v=v_rand_2(a,b); /* you may replace this by v_rand_1(a,b), but note
   that v_rand_1 is slower than v_rand_2 for large values of b */
   /* return the value (g/h(a,v))^(1-a) */
   return exp((1.0-a)*log(g/h(a,v)));
} /* end mittag_leffler_type2_two_parameter_rand */

int main(int argc, char *argv[]) {

   int n,              /* number of pseudo random numbers to be generated */
       i;              /* local index */
   double a,b,         /* parameters of the ML(a,b) distribution */
          x,           /* variable for pseudo random numbers */
          m,m2,v,      /* mean, second moment and variance of ML(a,b) */
          sm,sm2,svar; /* sample mean ... */

   if (argc != 4) {
      printf("USAGE: mitlef a b n\n\n");
      printf("where 0 < a < 1, b >= 0 and n is a positive integer.\n\n");
      printf("Example: mitlef 0.5 1 50");
      exit(0);
   }
   a=atof(argv[1]);
   b=atof(argv[2]);
   n=atoi(argv[3]);
   printf("Second type two-parameter Mittag-Leffler distribution ML(a,b)\n");
   printf("Parameters: a = %5.3f b = %5.3f\n",a,b);
   m=gamma(b+1.0)/a/gamma(a+b);              /* mean of ML(a,b) */
   m2=(a+b)*gamma(b+1.0)/a/a/gamma(2.0*a+b); /* second moment of ML(a,b) */
   v=m2-m*m;                                 /* variance of ML(a,b) */
   printf("Mean:          %9.3f\n",m);
   printf("Second moment: %9.3f\n",m2);
   printf("Variance:      %9.3f\n",v);
   my_randomize();
   printf("Generating n = %d pseudo random numbers for ML(a,b):\n",n);
   sm=0.0;sm2=0.0;
   for (i=1;i<=n;i++) {
      x=mittag_leffler_type2_two_parameter_rand(a,b);
      sm+=x; sm2+=x*x;
      printf("%5.3f ",x);
   }
   printf("\n");
   sm=sm/n; sm2=sm2/n;
   svar=sm2-sm*sm;
   printf("Sample mean:          %9.3f\n",sm);
   printf("Second sample momemt: %9.3f\n",sm2);
   printf("Sample variance:      %9.3f",svar);
   return 0;
} /* end main */
