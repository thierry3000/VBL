/*
 *  Utilities.cpp
 *  VBL
 *
 *  Created by Edoardo Milotti on 8/22/08.
 *  Copyright 2008 I.N.F.N.-Sezione di Trieste. All rights reserved.
 *
 */
 
#ifndef UTILITY_H
#define UTILITY_H  // header guard
// numeri casuali 

// versione standalone di ran2 (versione NR modificata)

# include <boost/random.hpp>
# include <boost/random/random_device.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>      // std::setprecision  std::scientific

#include <cmath> // log exp sin tan floor
// namespace vbl{


double ran2();
double gasdev(int &idum);
double gammln(const double xx);
double factln(const int n);
double bico(const int n, const int k);
double bnldev(const double pp, const int n, int &idum);

// }//namespace vbl{
#endif //headerguard
