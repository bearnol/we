
// Added JGW 2005-05-05
// Wanless Factorize, Wanless Extended Factorize
// copyright 2000 James Wanless. All rights reserved.

// Added JGW 2006-04-25
// Fermat probable primality test


#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>

#include <sys/time.h>
#include <sys/sysinfo.h>

#include <iostream>
#include <gmpxx.h>
#include <gmp.h>

using namespace std;

mpz_class N;

int quiet = 0;
int ecm = 0;
mpz_class base = 0;
int num_threads = 0;

pthread_t *threads;
int *thread_args;
mpz_class thread_base;
int i;
int result_code;

#define HAVE_GETTIMEOFDAY   1

gmp_randstate_t rstate;

long seed = 0;

struct point {mpz_class x, y;};

mpz_class rand2()
/* returns GMP pseudo-random number of 50-bits */
{
	mpz_t temp;
	mpz_init(temp);
	mpz_urandomb(temp, rstate, 50L);
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

mpz_class nextp(mpz_class n)
/* returns next prime after n */
{
	mpz_t temp;
	mpz_init(temp);
	mpz_nextprime(temp, n.get_mpz_t());
	mpz_class temp_class(temp);
	mpz_clear(temp);
	return temp_class;
}

mpz_class inverse(mpz_class a, mpz_class b)
/* returns inverse of a modulo b or 0 if it does not exist */
{
	mpz_t temp;
	mpz_init(temp);
	if (!mpz_invert(temp, a.get_mpz_t(), b.get_mpz_t()))
		mpz_set_si(temp, 0);
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

mpz_class Wanless (mpz_class n)
{
	mpz_class basestested = 0;
	mpz_class S = 1;
	mpz_class T = 1;
	mpz_class b = 1;
	mpz_class a = 2;
	mpz_class prodA = 1;
	mpz_class wanless = 2;
	mpz_class nbits = 1;
	mpz_class i = 1;
	
	if (n < 2)
		return 1;
	
	if (n % 2 == 0)
		return 2;	
	
	while (wanless < n) {
		wanless *= 2;
		nbits++;
	}
	
	a = thread_base;
	
	while (T == 1 || T == n) {
		basestested++;
		if (!quiet) {
			cout << "\rbase#" << basestested;
			fflush(stdout);
		}
		
		for (i = 1; i <= nbits; i++) {
			prodA *= a;
			prodA = modpos(prodA, n);
			a++;
		}
		
		b = exp_mod(prodA, wanless, n);
		S = exp_mod(b, n, n) - b;
		T = gcd(n, S);
	}
	
	if (!quiet) {
		cout << "\t[prodA = " << prodA << "]\n";
		fflush(stdout);
	}
	
	return T;
	
}

mpz_class wanless (mpz_class n)
{
	mpz_class S = 1;
	mpz_class T = 1;
	mpz_class a = 2;
	
	if (n < 2)
		return 1;
	
	if (n % 2 == 0)
		return 2;	
		
	a = 2;
	
	while (a*a <= n && (T == 1 || T == n) && a <= 2) {
		S = exp_mod(a, n, n) - a;
		T = gcd(n, S);
		a++;
	}
		
	return T;
	
}

mpz_class trialdivide(mpz_class n)
{
	mpz_class a = 2;
	
	while (a < 10000 && a*a <= n) {
		if (n % a == 0)
			return a;
		else
			a++;
	}
	
	return n;
}

int addition_1(mpz_class n, struct point P1, struct point P2, struct point *P3)
/* affine coords */
/* returns 1 if P1 = -P2 therefore P1 + P2 = O, 0 otherwise */
/* P1 != P2 */
{
	mpz_class delta_x;
	mpz_class delta_y;
	mpz_class m;
	
	delta_x = modpos(P2.x - P1.x, n);
	delta_y = modpos(P2.y - P1.y, n);
	
	if (P1.x == P2.x && ((P1.y + P2.y) == 0 || (P1.y + P2.y) == n || (P1.y + P2.y) == 2 * n)) {
		P3->x = 0, P3->y = 1;
		return 1;
	}
	
	/* calculate m = (y2 - y1)(x2 - x1) ^ -1 mod n */
	m = modpos(delta_y * inverse(delta_x, n), n);
	
	/* calculate x3 = m ^ 2 - (x1 + x2) mod n */
	P3->x = modpos(m * m - (P1.x + P2.x), n);
	
	/* calculate y3 = m(x1 - x3) - y1 mod n */
	P3->y = modpos(m * (P1.x - P3->x) - P1.y, n);
	
	return 0;
}

void addition_2(mpz_class a, mpz_class n, struct point P1, struct point *P3)
/* affine coords */
/* P1 == P2 */
{
	mpz_class m;
	
	/* calculate m = (3x1 ^ 2 + a)(2y1) ^ -1 mod n */
	m = modpos((3 * P1.x * P1.x + a) * inverse(2 * P1.y, n), n);
	
	/* calculate x3 = m ^ 2 - 2x1 mod n */
	P3->x = modpos(m * m - 2 * P1.x, n);
	
	/* calculate y3 = m(x1 - x3) - y1 mod n */
	P3->y = modpos(m * (P1.x - P3->x) - P1.y, n);
}

int multiply(mpz_class a, mpz_class k, mpz_class n, struct point P, struct point *R, mpz_class *d)
/* binary ladder */
/* returns -1 if O encountered, 0 if divisor not found, 1 otherwise */
{
	int value = 1;
	struct point A, B, C;
	
	/*  A = P */
	A = P;
	/* B = O = (0, 1) the point at infinity */
	B.x = 0, B.y = 1;
	
	while (value && k > 0) {
		if (k % 2 != 0) {
			*d = gcd(modpos(B.x - A.x, n), n);
			
			k--;
			value = (*d == 1 || *d == n);
			
			if (A.x == 0 && A.y == 1);
			else if (B.x == 0 && B.y == 1) B = A;
			else if (value) {
				addition_1(n, A, B, &C);
				B = C;
			}
		}
		else {
			*d = gcd(modpos(2 * A.y, n), n);
			
			k >>= 1;
			value = (*d == 1 || *d == n);
			
			if (value) {
				addition_2(a, n, A, &C);
				A = C;
			}
		}
	}
	
	*R = B;
	R->x = modpos(R->x, n);
	R->y = modpos(R->y, n);
	
	if (R->x == 0 && R->y == 1) return -1;
	
	return !value;
}


mpz_class LenstrasECM(mpz_class n)
{
	int found = 0;
	mpz_class B = 1000;
	mpz_class l, q, q1, newq;
	struct point x, y;
	mpz_class a, g, d;
	long curve = 0;
	
	if (n % 2 == 0)
		return 2;
	if (n % 3 == 0)
		return 3;
	
	do {
		
		do {
			a = rand2();
			x.x = 0;
			x.y = 1;
			
			curve++;
			if (!quiet) {
				cout << "\rB=" << B << ", curve#" << curve << ", a=" << a << "                    ";
				fflush(stdout);
			}
			
			newq = 2;
			for (; newq < B && found != 1;) {
				q = newq;
				q1 = q;
				l = B / q;
				while (q1 <= l)
					q1 *= q;
				
				found = multiply(a, q1, n, x, &y, &d);
				
				x.x = y.x;
				x.y = y.y;
				//  cout << "X=" << x.x << "," << x.y << "\n";
				
				newq = nextp(q);
			}
			
			g = gcd(d, n);
			
		} while (curve < 50 && (g == 1 || g == n));
		
		B = B * 10;
		curve = 0;
		
	} while (B <= 1000000000 && (g == 1 || g == n));
	
	if (!quiet) {
		cout << "\n";
		fflush(stdout);
	}
	
	return g;
}


bool fermat(mpz_class n)
{
	bool prime = false;
	mpz_class S = 1;
	mpz_class T = 1;
	mpz_class a = 2;
	
	if (n < 2)
		return false;
	
	if (n == 2)
		return true;
	
	if (n % 2 == 0)
		return false;
	
	T = n;
	
	while (a < 10 && T == n) {
		S = exp_mod(a, n, n) - a;
		T = gcd(n, S);
		a++;
	}
	
	if (T < n)
		prime = false;
	else
		prime = true;
	
	return prime;
}

void factorize(mpz_class n)
{
	if (Rabin_Miller(n) && fermat(n)) {
		cout << n << "\n";
		fflush(stdout);
		return;
	}
	
	mpz_class factor;
	
	factor = wanless(n);
	if (factor > 1 && factor < n) {
		factorize(factor);
		factorize(n/factor);
	}
	
	else {

	factor = trialdivide(n);
	if (factor > 1 && factor < n) {
		factorize(factor);
		factorize(n/factor);
	}
		
	else {
		if (ecm) {
	factor = LenstrasECM(n);
	if (factor > 1 && factor < n) {
			factorize(factor);
			factorize(n/factor);
	}
		}
		
	else {
		
	factor = Wanless(n);
	if (factor > 1 && factor < n) {
		factorize(factor);
		factorize(n/factor);
	}
		
	else {
			cout << "composite\t" << n << "\n";
			fflush(stdout);
	}
	}
	}
	}
	
	return;
}

void *perform_work(void *arguments){
  int index = *((int *)arguments);
  
  printf("THREAD %d: Started.\n", index);
  	if (base == 0)
		thread_base = rand2();
  	else
  		thread_base = base;
	if (!quiet) {
		cout << "base = " << thread_base << "\n";
		fflush(stdout);
	}
  factorize(N);
  printf("THREAD %d: Ended.\n", index);
  
  return NULL;
}

/* The name of this program. */
const char* program_name;

void print_usage (FILE* stream, int exit_code)
{
	fprintf(stream, "Usage: %s options [ < inputfile ]\n", program_name);
	fprintf(stream,
			"	-h	--help		Display this usage information.\n"
			"	-b	--base		Set initial base (default random 50-bit integer).\n"
			"	-e	--ecm		Include Lenstra's ECM factoring.\n"
			"	-q	--quiet		Print reduced messages.\n"
			"	-s	--seed		Set (pseudo-)random seed.\n"
			"	-t	--threads	Set number of threads.\n");
	exit(exit_code);
}


int main(int argc, char* argv[])
{
	int next_option;
	
	const char* const short_options = "hb:eqs:t:";
	const struct option long_options[] = {
		{"help", 0, NULL, 'h'},
		{"base", 1, NULL, 'b'},
		{"ecm", 0, NULL, 'e'},
		{"quiet", 0, NULL, 'q'},
		{"seed", 1, NULL, 's'},
		{"threads", 1, NULL, 't'},
		{NULL,0,NULL,0}
	};
	
	program_name = argv[0];
	
	do {
		next_option = getopt_long(argc, argv, short_options, long_options, NULL);
		switch (next_option)
		{
			case 'h':
				print_usage(stdout, 0);
				
			case 'b':
				base = input2(optarg);
				break;
				
			case 'e':
				ecm = 1;
				break;
				
			case 'q':
				quiet = 1;
				break;

			case 's':
				seed = atol(optarg);
				break;

			case 't':
				num_threads = atol(optarg);
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
	if (!quiet) {
		cout << "random seed = " << seed << "\n";
		fflush(stdout);
	}



	if (num_threads == 0)
		num_threads = get_nprocs();
	if (!quiet) {
		cout << "threads = " << num_threads << "\n";
		fflush(stdout);
	}

	if (!quiet) {
		cout << "number to be tested:\n";
		fflush(stdout);
	}
	cin >> N;

	threads =(pthread_t *)malloc((num_threads + 1) * sizeof(pthread_t));
	thread_args =(int *)malloc((num_threads + 1) * sizeof(int));
		
  //create all threads one by one
  for (i = 0; i < num_threads; i++) {
    printf("IN MAIN: Creating thread %d.\n", i);
    thread_args[i] = i;
    result_code = pthread_create(&threads[i], NULL, perform_work, &thread_args[i]);
    assert(!result_code);
    sleep(5);
  }
  
  printf("IN MAIN: All threads are created.\n");

  //wait for each thread to complete
  for (i = 0; i < num_threads; i++) {
    result_code = pthread_join(threads[i], NULL);
    assert(!result_code);
    printf("IN MAIN: Thread %d has ended.\n", i);
  }

  printf("MAIN program has ended.\n");

	free(thread_args);
	free(threads);
		
	return 0;
}
