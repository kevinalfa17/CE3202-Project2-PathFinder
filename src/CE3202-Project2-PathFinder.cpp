//============================================================================
// Name        : CE3202-Project2-PathFinder.cpp
// Author      : Kevin, Dennis, David
// Version     : 1.0
// Copyright   : Costa Rica Institute of Technology
// Description : Project for the course Numerical Analysis of the Costa Rica Institute of Technology
//============================================================================

#include "config.h"
#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE)
#include "Matrix/Matrix_simd.hpp"
#else
#include "Matrix/Matrix_simd.hpp"
#endif
#include <iostream>
#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE)
#include <xmmintrin.h>
#endif

using namespace std;
using namespace anpi;




/**
 * @brief Print a matrix
 * @param m: Matrix to print
 */
template <typename T>
void printMatrix(anpi::Matrix<T> &m){
	for(int i = 0; i < m.rows(); i++){
		cout << "|\t";
		for(int j = 0; j < m.cols(); j++)
			cout << "[" << m[i][j] << "]\t";
		cout << "|" <<endl;
	}
}


#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE)
int main() {
	cout << "!!!Hello World SIMD is active and available!!!!" << endl; // prints !!!Hello World!!!

	//SSE Example: calculate sqrt(x)/x for values from 1 to 400
	int length = 400;
	float* pResult = (float*) _mm_malloc (length * sizeof(float), 16);//Align 400 float spaces to 16-byte for SSE

	__m128 x; //16 bytes vector (128bits)
	__m128 xDelta = _mm_set1_ps(4.0f);// [4,4,4,4]
	__m128 *pResultSSE = (__m128*)pResult; //Destination vector for each calc

	x = _mm_set_ps(4.0f,3.0f,2.0f,1.0f); //Inicial vector [4,3,2,1]

	//End condition, for will repeat only 100 times because SSE do 4 operations at the same time
	const int SSELength = length/4;
	//100 iterations
	for(int i = 0; i<SSELength; i++){

		__m128 xSqrt = _mm_sqrt_ps(x); // xSqrt = [sqrt(x_i),sqrt(x_i+1),sqrt(x_i+2),sqrt(x_i+3)]

		pResultSSE[i] = _mm_div_ps(xSqrt,x); // pResultSSE[i] = [xSqrt_i/x,xSqrt_i+1/x,xSqrt_i+2/x,xSqrt_i/x]

		// x = [x_i + 4, x_i+1 + 4, x_i+2 + 4,x_i+3 + 4], example: [1+4,2+4,3+4,4+4] = [5,6,7,8] in the first iteration
		x = _mm_add_ps(x,xDelta);
	}

	//Print first 20 result to test
	for(int i = 0;i<20;i++){
		cout <<"Result["<<i<<"] = "<<pResult[i]<<endl;
	}

	//Matrix example
	Matrix<float> M = Matrix<float>(3, 3, float(4), Matrix<float>::Padded);
	cout <<"dcols"<<M.dcols()<<endl;
	cout <<"ecols"<<M.ecols()<<endl;

	printMatrix(M);


	return 0;
}
#else
int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	return 0;
}
#endif
