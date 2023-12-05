/*
  Author:  Pate Williams (c) 1997 & James Wanless (c) 2000-23

  Algorithm Wanless Extended (orthogonal)
*/

#include <getopt.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>
#include <time.h>

#include <iostream>
#include <gmpxx.h>
#include <gmp.h>

using namespace std;

#define HAVE_GETTIMEOFDAY   1

gmp_randstate_t rstate;

long seed = 0;

mpz_class rand2()
/* returns GMP pseudo-random number of 100-bits */
{
	mpz_t temp;
	mpz_init(temp);
	mpz_urandomb(temp, rstate, 100L);
	mpz_class temp_class(temp);
	mpz_clear(temp);
	return temp_class;
}

mpz_class gcd(mpz_class a, mpz_class b)
/* returns greatest common divisor of a and b */
{
  mpz_t temp;
  mpz_init(temp);
  mpz_gcd(temp, a.get_mpz_t(), b.get_mpz_t());
  mpz_class temp_class(temp);
  mpz_clear(temp);
  return temp_class;
}

mpz_class modpos(mpz_class a, mpz_class b)
/* returns a modulo b, strictly non-negative */
{
	mpz_class temp_class;
	temp_class = a % b;
	if (temp_class < 0)
		temp_class += b;
	return temp_class;
}

mpz_class exp_mod(mpz_class x, mpz_class b, mpz_class n)
/* returns x ^ b mod n */
{
        mpz_t temp;
        mpz_init (temp);
        mpz_powm(temp, x.get_mpz_t(), b.get_mpz_t(), n.get_mpz_t());
        mpz_class temp_class(temp);
        mpz_clear (temp);
        return temp_class;
}


mpz_class WanlessEO(mpz_class N)
{
	mpz_class B = 1000l;
	mpz_class za, zd, zg = 1;
	mpz_class zq, zr;
	mpz_class base = 0;

	if (N % 2 == 0)
		return 2;
	if (N % 3 == 0)
		return 3;

	do {
		
		do {
			za = modpos(rand2(), N);

			for (mpz_class newza = 1; newza < B; newza++) {
				za = za * modpos(rand2(), N);
				za = modpos(za, N);
			}

			base++;
			cout << "B=" << B << ", base#" << base << ", za=" << za << "          \r";
			fflush(stdout);

			zq = exp_mod(za, N, N);
			zr = za;
			
			for (mpz_class newq = 1; newq < B; newq++) {
				zq = zq * zq;
				zq = modpos(zq, N);
				zr = zr * zr;
				zr = modpos(zr, N);
			}
			
			zd = modpos(zq - zr, N);
			zg = gcd(zd, N);
		
		} while (base < 50 && (zg == 1 || zg == N));
  
		B = B * 10;
		base = 0;

	} while (zg == 1);
  
	cout << "\n";
  
	return zg;
}

void factorize(mpz_class n) {
	mpz_class factor;

	factor = WanlessEO(n);

	if (factor == n)
		cout << n << "\n";

	if (factor > 1 && factor < n) {
		factorize(factor);
		factorize(n/factor);
	}
}

/* The name of this program. */
const char* program_name;

void print_usage (FILE* stream, int exit_code)
{
	fprintf(stream, "Usage: %s options [ < inputfile ]\n", program_name);
	fprintf(stream,
			"	-h	--help		Display this usage information.\n"
			"	-s	--seed		Set (pseudo-)random seed.\n");
	exit(exit_code);
}


int main(int argc, char* argv[])
{
	char answer[256];
	mpz_class N = 0;

	int next_option;
	
	const char* const short_options = "hs:";
	const struct option long_options[] = {
		{"help", 0, NULL, 'h'},
		{"seed", 1, NULL, 's'},
		{NULL,0,NULL,0}
	};
	
	program_name = argv[0];
	
	do {
		next_option = getopt_long(argc, argv, short_options, long_options, NULL);
		switch (next_option)
		{
			case 'h':
				print_usage(stdout, 0);
				
			case 's':
				seed = atol(optarg);
				break;
				
			case '?':
				print_usage(stderr, 1);
				
			case -1:
				break;
				
			default:
				abort();
		}
	}
	while (next_option != -1);
	
	gmp_randinit_default (rstate);
	
    {
#if HAVE_GETTIMEOFDAY
        struct timeval tv;
        gettimeofday (&tv, NULL);
		if (!seed)
			seed = tv.tv_sec + tv.tv_usec;
		
#else
        time_t t;
        time (&t);
		if (!seed)
			seed = t;
#endif
    }
	
	gmp_randseed_ui (rstate, seed);
	
	cout << "random seed = " << seed << "\n";
	fflush(stdout);
	
	printf("enter the number to be factored:\n");
	cin >> N;
	  
	factorize(N);

	return 0;
}
