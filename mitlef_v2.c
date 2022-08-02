/* Filename: mitlef_v2.c
   Version:  Standalone ANSI C Version 2.0
   Date:     August 2, 2022

   Author:   Martin Moehle
             Mathematical Institute
             University of Tuebingen
             Germany
             martin.moehle@uni-tuebingen.de

   Please compile this file using the command "gcc mitlef_v2.c -o mitlef.exe".

   This C program demonstrates random number generation for the (second
   type) two-parameter Mittag-Leffler distribution ML(a,b). For the
   mathematical background we refer the reader to (Section 8 of)

   Moehle, M.: A restaurant process with cocktail bar and relations
   to the three-parameter Mittag-Leffler distribution, J. Appl. Probab. 58,
   978-1006 (2021).

   For an analog program written in C++11 see "ml_v2.cpp". 

   Version Date             Comments, changes etc.
   ----------------------------------------------------------------------------
   1.0     May 22, 2022     Version which works for parameters 0 < a < 1 and
                            b >= 0.

   2.0     August 2, 2022   Program now works for parameters 0 < a < 1 and
                            b > -a. The new functions "T", "S", "v_rand_3"
                            and "v_rand_4" cover the case of negative b (>-a).
                            Function "gamma_rand" fully new developed.
                            Moreover, function "loggamma" and functions for
                            the mean number of steps of the four algorithms
                            and sample steps added. */

#include <stdlib.h> /* RAND_MAX, rand, srand, ... */
#include <stdio.h>  /* printf, ... */
#include <time.h>   /* time, ... */
#include <math.h>   /* exp, log, ... */
#include <stdbool.h> /* boolean variables */

/* global variables of the program */
unsigned int steps; /* counter for the number of steps of the used rejection method */ 

#define PI 3.141592653589793238462643

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
   if (x<0.5) return PI/(sin(PI*x)*gamma(1.0-x)); /* reflection formula */
   h=cc[0]; for (i=5;i>=0;i--) h+=cc[i+1]/(x+i);
   return sqrt(2.0*PI)*exp((x-0.5)*log(x+4.5)-(x+4.5))*h;
} /* end gamma */

double loggamma(double x) {
/* logarithm of the gamma function for x > 0 */
   double h;
   int i;
   if (x<0.5) return log(PI)-log(sin(PI*x))-loggamma(1.0-x); /* reflection formula */
   h=cc[0]; for (i=5;i>=0;i--) h+=cc[i+1]/(x+i);
   return 0.5*log(2.0*PI)+log(h)+(x-0.5)*log(x+4.5)-(x+4.5);
} /* end loggamma */

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

double best(double a) {
/* Best's (1978) rejection algorithm XG to sample from the gamma distribution
   with density x -> x^(a-1) exp(-x) / gamma(a), x > 0. This algorithm works
   for shape parameter a >= 1. For more details we refer the reader to page
   410 of the book of Devroye (1986), Non-Uniform Random Variate Generation. */
   double b,c,u,v,w,x,y,z;
   bool accept;
   b=a-1.0; c=3.0*a-0.75;
   do {
      u=uniform_rand();v=uniform_rand();
      w=u*(1.0-u);y=sqrt(c/w)*(u-0.5);x=b+y;
      if (x>=0) {
         z=64*w*w*w*v*v;
         accept=(z<=1.0-2.0*y*y/x);
         if (!accept) accept=(log(z)<=2.0*(b*log(x/b)-y));    
      }
   }
   while (!accept);
   return x;
} /* end best */

double johnk(double a) {
/* Johnk's (1971) rejection algorithm to sample from gamma distribution
   with density x -> x^(a-1) exp(-x) / gamma(a), x > 0. This algorithm
   works for shape parameter 0 < a < 1. Choose E and W independent with
   E exponentially distributed with parameter 1 and W beta distributed
   with parameters a and 1-a. Then EW is gamma distributed (with scale
   parameter b=1) and shape parameter 0 < a < 1. For more details we refer
   the reader to page 418 of the book of Devroye (1986), Non-Uniform
   Random Variate Generation. */
   double x,y;
   do {
      x=exp(log(uniform_rand())/a);
      y=exp(log(uniform_rand())/(1.0-a));
   }
   while (x+y>1); /* w := x/(x+y) is now a sample from beta(a,1-a) */
   return -log(uniform_rand())*x/(x+y);
} /* end johnk */

double gamma_rand(double b, double c) {
/* Returns a random sample from the gamma distribution with scale
   parameter b > 0 and shape parameter c > 0 having density
   x -> (1/b)^c x^(c-1) exp(-x/b) / gamma(c), x > 0. For shape
   parameters c >= 1 Best's (1978) algorithm is used. For 0 < c < 1
   Johnk's (1971) algorithm is used. */
   if (c>=1) return b*best(c); else return b*johnk(c);
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
   v_rand_2, v_rand_3 and mittag_leffler_type2_two_parameter_rand. For more
   details see Eq. (36) in Moehle (2021). */
   return sin((1.0-a)*v)/sin(a*v)*exp(log(sin(a*v)/sin(v))/(1.0-a));
} /* end h */

double log_h(double a, double v) {
/* log of h */
   return (log(sin(a*v))-log(sin(v)))/(1.0-a)+log(sin((1.0-a)*v))-log(sin(a*v));
} /* end log_h */

double v_rand_1(double a, double b) {
/* Pseudo random number generator for the (0,PI)-valued random variable
V = V(a,b), 0 < a < 1, b >= 0, having density v |-> c(a,b)/(h(a,v))^(b*(1-a)/a),
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
   steps=0; /* initialize the global variable steps */
   /* rejection Algorithm 1 of Moehle (2021) */
   do {
      u=uniform_rand();    /* sample u uniformly on (0,1) */
      v=PI*uniform_rand(); /* sample v uniformly on (0,PI) */
      steps+=1; /* count this step */
   } while (u>exp(p*log(h0/h(a,v))));
   return v;
} /* end v_rand_1 */

double v_rand_2(double a, double b) {
/* same as v_rand_1, but Algorithm 2 in Section 8 of Moehle (2021) is used.
v_rand_2 is efficient (uniformly for all parameters 0 < a < 1 and b >= 0)
and (for most parameter values) more efficient than v_rand_1, in
particular for large values of b */
   double u,  /* local variable for uniform random number */
          v,  /* variable for random number of V */
          h0, /* variable for h(a,0+) = (1-a)*a^(a/(1-a)) */
          p,  /* variable for b*(1-a)/a */
          t;  /* variable for (cos(v/2))^(2*b) * (h(a,v)/h(a,0+))^p */
   if (a>0.5) { b*=(1.0-a)/a; a=1.0-a; } /* transformation to a <= 0.5 */
   p=b*(1.0-a)/a;
   h0=(1.0-a)*exp(a/(1.0-a)*log(a));
   steps=0; /* initialize the global variable steps */
   /* rejection Algorithm 2 of Moehle (2021) */
   do {
      u=uniform_rand(); /* sample u uniformly on (0,1) */
      /* sample v from V(1/2,b) via (45) of Moehle (2021) */
      v=2.0*acos(sqrt(beta_rand(b+0.5,0.5)));
      t=exp(2.0*b*log(cos(v/2.0))+p*log(h(a,v)/h0));
      steps+=1; /* count this step */
   } while (u*t>1.0);
   return v;
} /* end v_rand_2 */

double T(double a, double b, double v) {
/* The transformation T(v) := (1-v/PI)^(b/a+1), 0 < v < PI. */
   return exp((b/a+1.0)*log(1.0-v/PI));
} /* end T */

double S(double a, double b, double w) {
/* The inverse transformation S := T^{-1} of T is given by
S(w) = PI(1-w^(a/(a+b))), 0 < w < 1. This function is used by
the following function "v_rand_3". */
   return PI*(1.0-exp(a/(a+b)*log(w)));   
} /* end S */

double v_rand_3(double a, double b) {
/* Algorithm to sample from V if the parameter b (>-a) is negative. Note that
the corresponding rejection method is not documented in Moehle (2021). For
negative parameter b, the density of V is unbounded. Thus, the transformed
random variable W := T(V) with transformation T(v) := (1-v/PI)^(b/a+1) is
considered. The transformation T is chosen such that the density of W is
bounded. Thus, a standard rejection algorithm is applicable to sample from W.
The back transformed value v:=S(w), where S denotes the inverse of T, is
then a sample from V. This algorithm theoretically works for any 0 < a < 1
and -a < b < 0, however it is not very efficient (and also numerically unstable)
for (a,b) close to (0,0) or (1,-1) due to the fact that in this region, V
is very concentrated around PI and hence the density of V is highly peaked. */
   const double eps=1e-8; /* a small number greater than the maschine exactness */
   double u,  /* local variable for uniform random number */
          v,  /* variable for random number of V := S(a,b,w) */
          w,  /* variable for random number of W */
          h0, /* variable for h(a,0+) = (1-a)*a^(a/(1-a)) */
          t,  /* variable for (h0/h(a,v))^(b/a(1-a)) w^(-b/(a+b)) */
          C;  /* variable for (sin(a*PI))^(1/(1-a))*sin((1-a)*PI)/sin(a*PI) */
   C=exp(log(sin(a*PI))/(1.0-a))*sin((1.0-a)*PI)/sin(a*PI);
   h0=(1.0-a)*exp(a/(1.0-a)*log(a));
   steps=0; /* initialize the global variable steps */
   do {
      u=uniform_rand(); /* sample u and w */
      w=uniform_rand(); /* uniformly on (0,1) */
      v=S(a,b,w);       /* compute the back transformed v := S(a,b,w) */ 
      /* compute t = t(w) */
      if (v<PI-eps)
         /* use the usual formula for t = t(w) if v is not too close to PI */
         t=exp(b/a*(1.0-a)*log(h0/h(a,v))-b/(a+b)*log(w));
      else /* v is close to PI */
         /* use t(0+) instead of t(w) to avoid numerical problems */
         t=exp(b/a*log(PI)+b*(1.0-a)/a*log(h0/C));
      steps+=1; /* count this step */
   } while (u>t);
   return v; /* return the back transformed value v = S(a,b,w) */
} /* end v_rand_3 */

double v_rand_4(double a, double b) {
/* This algorithm again covers the case of negative b (>-a), but "v_rand_4"
should be prefered compared to "v_rand_3", since it is more efficient. This
algorithm uses the same transformation T as v_rand_3, but in addition,
similar as in "v_rand_2" compared to "v_rand_1", a head function with
a = 1/2 is used. The method theoretically works for all parameter values
0 < a < 1  and -a < b < 0. The efficiency is not very good and the algorithm
becomes numerically unstable if (a,b) is close to (0,0) or (1,-1). */
   const double eps=1e-8; /* a small number greater than the maschine exactness */
   double u,  /* local variable for uniform random number */
          v,  /* variable for random number of V := S(a,b,w) */
          w,  /* variable for random number of W */
          h0, /* variable for h(a,0+) = (1-a)*a^(a/(1-a)) */
          p,  /* variable for b*(1-a)/a */
          t,  /* variable for (cos(S(0.5,b,w)/2))^(2*b) * (h(a,v)/h0)^p * w^(b/(a+b)-b/(0.5+b))*/
          C;  /* variable for (sin(a*PI))^(1/(1-a))*sin((1-a)*PI)/sin(a*PI) */
   if (a>0.5) { b*=(1.0-a)/a; a=1.0-a; } /* transformation to a <= 0.5 */
   C=exp(log(sin(a*PI))/(1.0-a))*sin((1.0-a)*PI)/sin(a*PI);
   p=b*(1.0-a)/a;
   h0=(1.0-a)*exp(a/(1.0-a)*log(a));
   steps=0; /* initialize the global variable steps */
   do {
      u=uniform_rand();
      /* sample w from W = T(1/2,b,V(1/2,b)) via (45) of Moehle (2021) */
      w=T(0.5,b,2.0*acos(sqrt(beta_rand(b+0.5,0.5))));
      v=S(a,b,w); /* compute the back transformed v := S(a,b,w) */
      if (v<PI-eps)
         /* use the usual formula for t = t(w) if v is not too close to PI */
         t=exp(
           2.0*b*log(cos(S(0.5,b,w)/2.0))
           +p*log(h(a,v)/h0)
           +(b/(a+b)-b/(b+0.5))*log(w)
         );
      else /* v is close to PI */
         /* use t(0+) instead of t(w) to avoid numerical problems */
         t=exp(2.0*b*log(PI/2.0)+p*log(C/h0)-b/a*log(PI));
      steps+=1; /* count this step */
   } while (u*t>1.0);
   return v; /* return the back transformed value v = S(a,b,w) */
} /* end v_rand_4 */

double mittag_leffler_type2_two_parameter_rand(double a, double b) {
/* Let L be second type two-parameter Mittag-Leffler distributed
with parameters 0<a<1 and b>-a. This function is based on the fact
that L has the same distribution as (G/h(V))^(1-a), where G is gamma
distributed with parameters 1 and b*(1-a)/a+1 having density
x |-> x^(b*(1-a)/a)*exp(-x)/gamma(b*(1-a)/a+1), the function h(a,.) is
defined via

h(a,u) := (sin(a*u)/sin(u))^(1/(1-a)) * sin((1-a)*u)/sin(a*u)

and the (0,PI)-valued random variable V is independent of G and has
density v |-> c(a,b)/(h(a,v))^(b*(1-a)/a), with normalizing constant

c(a,b) := gamma(b+1)*gamma(b*(1-a)/a+1)/PI/gamma(b/a+1).

For b>=0 pseudo random numbers of V are generated via the rejection
method (Algorithm 2) (or alternatively Algorithm 1) described
in Section 8 of Moehle (2021). For negative b (>-a) a particular
algorithm is used, which is not documented in Moehle (2021); see
"v_rand_3" and "v_rand_4" for more details. */

   double g,  /* variable for random number of G (gamma distributed) */
          v;  /* variable for random number of V */
   /* sample g gamma distributed with parameters 1 and b*(1-a)/a+1 */
   g=gamma_rand(1.0,b*(1.0-a)/a+1.0);
/*   printf("\n\nrand: p=%f g=%f\n",b*(1.0-a)/a+1.0,g); */
   /* sample v from the random variable V(a,b) */
   if (b<0) v=v_rand_4(a,b); /* particular algorithm for negative b; you may
   replace this by "v_rand_3(a,b)", but note that "v_rand_3" is never better
   than "v_rand_4". */
   else /* b >=0 */
   v=v_rand_2(a,b); /* you may replace this by "v_rand_1(a,b)", but note
   that "v_rand_1" is not very efficient for large values of beta and
   never better than "v_rand_2". */
/* printf("rand: v=%f h(a,v)=%f (1-a)v/sin(v)*exp(-v*cot(v))=%f g/h(a,v)=%f\n", 
   v,h(a,v),(1.0-a)*v/sin(v)*exp(-v*cos(v)/sin(v)),g/h(a,v)); */
   /* return the value (g/h(a,v))^(1-a) */
   return exp((1.0-a)*log(g/h(a,v)));
} /* end mittag_leffler_type2_two_parameter_rand */

double mean_1(double a, double b) {
/* Expected number E(N) of steps of the rejection function v_rand_1;
   see the formula in Remark 24 of Moehle (2021). */
   double h; /* local variable */
   h=loggamma(b+1.0)+loggamma(b*(1.0-a)/a+1.0)-loggamma(b/a+1.0)
     -b*log(a)-b*(1.0-a)/a*log(1.0-a);
   return exp(h);
} /* end mean_1 */

double mean_2(double a, double b) {
/* Expected number E(N) of steps of the rejection function v_rand_2;
   see the formula for d(a,b) at the bottom of page 1004 of Moehle (2021). */
   double h; /* local variable */ 
   h=loggamma(b+0.5)+loggamma(b/a*(1.0-a)+1.0)-0.5*log(PI)-loggamma(b/a+1.0)   
     -b*log(a)-b*(1.0-a)/a*log(1.0-a);
   return exp(h);
} /* end mean_2 */

double mean_3(double a, double b) {
/* Expected number E(N) of steps of the rejection function v_rand_3.
This function is not documented in Moehle (2021). The formula is

          gamma(b+1) * gamma(b(1-a)/a+1)
E(N) = -------------------------------------.
       gamma(b/a+2) * a^b * (1-a)^(b(1-a)/a)
*/
   double h;
   h=loggamma(b+1)+loggamma(b*(1.0-a)/a+1.0)-loggamma(b/a+2.0)
     -b*log(a)-b*(1.0-a)/a*log(1.0-a);
   return exp(h);
} /* end mean_3 */

double mean_4(double a, double b) {
/* Expected number E(N) of steps of the rejection function v_rand_4.
This function is not documented in Moehle (2021). The formula is

       2 * gamma(b+3/2) * gamma(b(1-a)/a+1)
E(N) = ------------------------------------------------.
       sqrt(PI) * gamma(b/a+2) * a^b * (1-a)^(b(1-a)/a)
*/
   double h;   
   h=log(2.0)+loggamma(b+1.5)+loggamma(b*(1.0-a)/a+1.0)-0.5*log(PI)
   -loggamma(b/a+2.0)-b*log(a)-b*(1.0-a)/a*log(1.0-a);
   return exp(h);
} /* end mean_4 */

int main(int argc, char *argv[]) {

   int n,              /* number of pseudo random numbers to be generated */
       i;              /* local index */
   double a,b,         /* parameters of the ML(a,b) distribution */
          x,           /* variable for pseudo random numbers */
          m,m2,v,      /* mean, second moment and variance of ML(a,b) */
          sm,sm2,svar, /* sample mean, sample second moment and sample variance */
          ssteps;      /* sample number of steps */  
   if (argc!=4) {
      printf("This is mitlef, Version 2, August 2, 2022.\n\nUSAGE: mitlef a b n\n\n");
      printf("where 0 < a < 1, b > -a and n is a positive integer.\n\n");
      printf("Example: mitlef 0.5 1 50");
      exit(0);
   }
   a=atof(argv[1]);
   b=atof(argv[2]);
   n=atoi(argv[3]);
   if ((a<=0)||(a>=1)) {
      printf("Parameter a should satisfy 0 < a < 1.\n");
      printf("Program aborted.");
      exit(1); }
   if (b<=-a) {
      printf("Parameter b should satisfy b > -a.\n");
      printf("Program aborted.");
      exit(1);  }
   if (n<=0) {
      printf("Parameter n should be a positive integer.\n");
      printf("Program aborted.");
      exit(1); }
   printf("Second type two-parameter Mittag-Leffler distribution ML(a,b)\n");
   printf("Parameters: a = %5.3f b = %5.3f\n",a,b);
   m=gamma(b+1.0)/a/gamma(a+b);              /* mean of ML(a,b) */
   m2=(a+b)*gamma(b+1.0)/a/a/gamma(2.0*a+b); /* second moment of ML(a,b) */
   v=m2-m*m;                                 /* variance of ML(a,b) */
   printf("Mean:          %9.3f\n",m);
   printf("Second moment: %9.3f\n",m2);
   printf("Variance:      %9.3f\n",v);
   my_randomize();
   printf("Generating n = %d samples from ML(a,b) via ",n);
 
   if (b>=0) { 
      printf("Algorithm 2.\n");
      printf("Expected number of steps of rejection Algorithm 2: E(N) = ");
      if (a<=0.5) printf("%f\n",mean_2(a,b)); else printf("%f\n",mean_2(1.0-a,b*(1.0-a)/a));
      printf("(Expected number of steps of rejection Algorithm 1: E(N) = %f)\n",mean_1(a,b));
   }
   else { /* b < 0 */
      printf("Algorithm 4, since b < 0.\n");
      printf("Expected number of steps of rejection Algorithm 4: E(N) = ");
      if (a<=0.5) printf("%f\n",mean_4(a,b)); else printf("%f\n",mean_4(1.0-a,b*(1.0-a)/a));
      printf("(Expected number of steps of rejection Algorithm 3: E(N) = %f)\n",mean_3(a,b));
   }
   sm=0.0;sm2=0.0;ssteps=0;
   for (i=1;i<=n;i++) {
      x=mittag_leffler_type2_two_parameter_rand(a,b);
      sm+=x; sm2+=x*x;ssteps+=steps;
      printf("%5.3f ",x);
   }
   printf("\n");
   sm=sm/n;sm2=sm2/n;
   svar=sm2-sm*sm;
   ssteps=ssteps/n;
   printf("Sample mean:            %9.3f\n",sm);
   printf("Second sample momemt:   %9.3f\n",sm2);
   printf("Sample variance:        %9.3f\n",svar);
   printf("Sample number of steps: %9.3f (compare with E(N))",ssteps);
   return 0;
} /* end main */
