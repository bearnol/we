/* Factoring general integers with WE method using random base(s), with residue.
 Author:  James Wanless (c) 2000-10
*/

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sys/time.h>
#include <time.h>

#include <iostream>
#include <gmpxx.h>
#include <gmp.h>

using namespace std;

#define HAVE_GETTIMEOFDAY   1

gmp_randstate_t rstate;

mpz_class base = 0;
mpz_class trials = 0;
int quiet = 0;

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

mpz_class input2(char * inputstr)
/* returns GMP number from denary string */
{
	mpz_t temp;
	mpz_init(temp);
	mpz_set_str(temp, inputstr, 10);
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

int Rabin_Miller(mpz_class n)
/* given an integer n >= 3 returns 0 composite, 1 probably prime */
{
	return mpz_probab_prime_p(n.get_mpz_t(), 10);
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

bool Wanless (mpz_class n)
{
	mpz_class basestested = 0;
	mpz_class S = 1;
	mpz_class T = 1;
	mpz_class b = 1;
	mpz_class A = 1;
	mpz_class prodA = 1;
	mpz_class wanless = 1;
	mpz_class mersenne = 1;
	mpz_class residue = 1;
	long nbits = 1;
	
	for (int i = 0; i < 127; i++)
		mersenne *= 2;
	mersenne -= 1;
	
	if (!quiet) {
		cout << "mersenne = " << mersenne << "\n";
		fflush(stdout);
	}
	
	if (Rabin_Miller(n)) {
		cout << n << "\n";
		fflush(stdout);
		return true;
	}	

	if (n < 2)
		return false;
	
	if (n % 2 == 0) {
		cout << "2\n";
		fflush(stdout);
		return false;
	}
	
	
	while (wanless < n) {
		wanless *= 2;
		nbits ++;
	}
	
	A = base;
	
	while ((T == 1 || T == n) && (basestested < trials)) {
		basestested++;
		if (!quiet) {
			cout << "base#" << basestested << "\r";
			fflush(stdout);
		}
		
		for (long i = 0; i < nbits; i++) {
			prodA *= A;
			prodA = modpos(prodA, n);
			A++;
		}
				
		b = exp_mod(prodA, wanless, n);
		S = exp_mod(b, n, n) - b;
		T = gcd(n, S);
		residue *= S;
		residue = modpos(residue, mersenne);
	}
	
	cout << "residue = " << residue << "\n";
	fflush(stdout);
	
	if (T > 1 && T < n) {
		cout << "factor = " << T << "\t[prodA = " << prodA << "]\n";
		fflush(stdout);
	}
		
	
	return (false);

}

/* The name of this program. */
const char* program_name;

void print_usage (FILE* stream, int exit_code)
{
	fprintf(stream, "Usage: %s options [ < inputfile ]\n", program_name);
	fprintf(stream,
			"	-h	--help		Display this usage information.\n"
			"	-q	--quiet		Print reduced messages.\n"
			"	-b	--base		Set initial base (default random 100-bit integer).\n"
			"	-t	--trials	Set #trials to run (default 1).\n");
	exit(exit_code);
}


int main(int argc, char* argv[])
{
	int next_option;
	
	const char* const short_options = "hqb:t:";
	const struct option long_options[] = {
		{"help", 0, NULL, 'h'},
		{"quiet", 0, NULL, 'q'},
		{"base", 1, NULL, 'b'},
		{"trials", 1, NULL, 't'},
		{NULL,0,NULL,0}
	};
	
	program_name = argv[0];
	
	do {
		next_option = getopt_long(argc, argv, short_options, long_options, NULL);
		switch (next_option)
		{
			case 'h':
				print_usage(stdout, 0);
				
			case 'q':
				quiet = 1;
				break;
				
			case 'b':
				base = input2(optarg);
				break;
				
			case 't':
				trials = input2(optarg);
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
		if (base == 0) {
			gmp_randseed_ui (rstate, tv.tv_sec + tv.tv_usec);
			if (!quiet) {
				cout << "random seed = " << tv.tv_sec + tv.tv_usec << "\n";
				fflush(stdout);
			}
		}
		
#else
        time_t t;
        time (&t);
		if (base == 0) {
			gmp_randseed_ui (rstate, t);
			if (!quiet) {
				cout << "random seed = " << t << "\n";
				fflush(stdout);
			}
		}
#endif
    }
	
	if (base == 0)
		base = rand2();
	cout << "base = " << base << "\n";
	fflush(stdout);
	
	if (trials == 0)
		trials = 1;

	cout << "trials = " << trials << "\n";
	fflush(stdout);
	
	
	
	if (!quiet) {
		cout << "number to be tested:\n";
		fflush(stdout);
	}
	
	mpz_class N;
	
	cin >> N;
	Wanless(N);
		
	return 0;
}
