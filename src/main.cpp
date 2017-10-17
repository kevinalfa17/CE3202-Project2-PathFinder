//============================================================================
// Name        : CE3202-Project2-PathFinder.cpp
// Author      : Kevin, Dennis, David
// Version     : 1.0
// Copyright   : Costa Rica Institute of Technology
// Description : Project for the course Numerical Analysis of the Costa Rica Institute of Technology
//============================================================================

#include "config.h"
#include "Matrix/Matrix.hpp"
#include "PathFinder/NodePair.h"
#include "PathFinder/IndexMap.h"
#include "PathFinder/PathFinder.h"
#include <iostream>
#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <opencv2/opencv.hpp>
#include "UploadImage.h"


#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE)
#include <xmmintrin.h>
#endif

using namespace cv;
using namespace std;
using namespace anpi;
using namespace boost::bimaps;




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

void fillHoles(Mat &mask)
{
    /*
     This hole filling algorithm is decribed in this post
     https://www.learnopencv.com/filling-holes-in-an-image-using-opencv-python-c/
     */

    Mat maskFloodfill = mask.clone();
    floodFill(maskFloodfill, cv::Point(0,0), Scalar(255));
    Mat mask2;
    bitwise_not(maskFloodfill, mask2);
    mask = (mask2 | mask);

}

/*template< class MapType >
void print_map(const MapType & map)
{
	typedef typename MapType::const_iterator const_iterator;

	for( const_iterator i = map.begin(), iend = map.end(); i != iend; ++i )
	{
		NodePair * pair = i->second;
		cout << i->first << "-->";
		pair->printPair();
	}
}*/


#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE)
int main() {
	cout << "!!!Hello World SIMD is active and available!!!!" << endl; // prints !!!Hello World!!!
	UploadImage * uploadImage = new UploadImage();



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
	Matrix<float> M = Matrix<float>(5, 5, float(4), Matrix<float>::Padded);
	cout <<"dcols "<<M.dcols()<<endl;
	cout <<"ecols "<<M.ecols()<<endl;
	Matrix<float> M2 = Matrix<float>(5, 5, float(5), Matrix<float>::Padded);
	Matrix<float> M3 = M2 - M;
	M3(4,4) = 7;
	printMatrix(M3);

	//Index mapping example
	IndexMap * indexMap = new IndexMap(4,4); //Mapping 4x4 matrix

	cout << indexMap->getXFromNodes(3,1,3,2)<<endl; //Get index of x between 31 and 32 nodes

	NodePair pair = indexMap->getNodesFromX(22); //Get nodes terminals of x
	cout << "pair" << endl;
	pair.printPair();

	//Initializer list test
	//Matrix<int> B = {{1,2,3,4,5,6,7,8,9,10,11},{1,2,3,4,5,6,7,8,9,10,11},{1,2,3,4,5,6,7,8,9,10,11},{1,2,3,4,5,6,7,8,9,10,11}};
	//printMatrix(B);

	//PathFinder
	uploadImage->upload();
	PathFinder<float>  * pathFinder = new PathFinder<float>(0,0,2,2,uploadImage->getImage());
	Matrix<int> A = pathFinder->getA();
	printMatrix(A);



	return 0;
}
#else
int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	//Non- SIMD Matrix test
	Matrix<float> M = Matrix<float>(5, 5, float(4), Matrix<float>::Padded); //Padded ignored
	cout <<"dcols"<<M.dcols()<<endl; //dcols = cols
	cout <<"ecols"<<M.ecols()<<endl; //ecols = 0
	Matrix<float> M2 = Matrix<float>(5, 5, float(5), Matrix<float>::Padded); //Padded ignored

	Matrix<float> M3 = M2 - M;


	printMatrix(M3);


	return 0;
}
#endif
