//============================================================================
// Name        : CE3202-Project2-PathFinder.cpp
// Author      : Kevin, Dennis, David
// Version     : 1.0
// Copyright   : Costa Rica Institute of Technology
// Description : Project for the course Numerical Analysis of the Costa Rica Institute of Technology
//============================================================================

#include "config.h"
#include <iostream>


using namespace std;

#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE)
int main() {
	cout << "!!!Hello World SIMD is active and available!!!" << endl; // prints !!!Hello World!!!
	return 0;
}
#else
int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	return 0;
}
#endif
