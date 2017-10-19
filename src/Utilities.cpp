/*
 *  Utilities.cpp
 *  VBL
 *
 *  Created by Edoardo Milotti on 8/22/08.
 *  Copyright 2008 I.N.F.N.-Sezione di Trieste. All rights reserved.
 *
 */
#include "Utilities.h"

//this was in main
bool eventON;       // Boolean variable indicating that the special event case was selected
double EventTime;   // Time associated with a single special event, added by hand (now only one type of special event)
double pAlt;        // Probability of the special event "transition to the alternative type"

/** see http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
 * to create random numbers with c++11
 */
#if __GNUC_PREREQ(4,5)
  //from http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0, 1);
#else
  
#endif

double ran2(int &idum)
{
#if __GNUC_PREREQ(4,5)
  return dis(gen);
#else
  std::srand(std::time(0));
  return (double)(std::rand()/RAND_MAX);
#endif
}

double gammln(const double xx)
{
	int j;
	double x,y,tmp,ser;
	static const double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
		-0.5395239384953e-5};

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<6;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double factln(const int n)
{
	static double a[101];

	if (n < 0) 
		std::cout << "Negative factorial in routine factln";
	if (n <= 1) 
		return 0.0;
	if (n <= 100)
		return (a[n] != 0.0 ? a[n] : (a[n]=gammln(n+1.0)));
	else 
		return gammln(n+1.0);
}

double bico(const int n, const int k)
{
	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}


double bnldev(const double pp, const int n, int &idum)
{
	const double PI=3.141592653589793238;
	int j;
	static int nold=(-1);
	double am,em,g,angle,p,bnl,sq,t,y;
	static double pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=0;j<n;j++)
			if (ran2(idum) < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran2(idum);
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran2(idum);
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
				-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran2(idum) > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}

double gasdev(int &idum)
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (idum < 0) iset=0;
	if (iset == 0) {
		do {
			v1=2.0*ran2(idum)-1.0;
			v2=2.0*ran2(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


