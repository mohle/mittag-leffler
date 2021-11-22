# mittag-leffler
Random number generator for the two-parameter Mittag-Leffler distribution.
This repository contains two files mitlef.c (C99 version) and ml.cpp (c++11 version).
Both programs (should) do exactly the same. They generate pseudo random numbers distributed according to the two-parameter Mittag-Leffler distribution ML(a,b) of type two (having moments of all orders). The algorithm is based on a rejection method, which is highly efficient and accurate for all parameter values 0<a<1 and b>0. More details are provided in Section 8 of Möhle, M. (2021) A restaurant process with cocktail bar and relations to the three-parameter Mittag-Leffler distribution, J. Appl. Probab. 58, 978-1006. (DOI: https://doi.org/10.1017/jpr.2021.10)
