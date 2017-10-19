/**
 * @file MatrixDescomposition.h
 * @brief Class that implements decomposition LU of matrix method to find equations solutions. And uses SIMD to optimize the method of lu.
 * @author Denporras
 * @date 28 de sep. de 2017
 */


#ifndef DESCOMPOSITION_MATRIXDESCOMPOSITION_H_
#define DESCOMPOSITION_MATRIXDESCOMPOSITION_H_

#include <iostream>	//Inputs and outputs of text
#include <cmath>	//Math operations
#include <limits>	//Limits precision
#include <vector>	//vector 
#include <stdexcept> //Common Exceptions

#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE)
#include <xmmintrin.h>
#include <emmintrin.h>
#endif

#include "Matrix.hpp"

using namespace std;

namespace anpi {

template<typename T>
class MatrixDescomposition {
public:
	MatrixDescomposition();
	void lu(const Matrix<T> &A, Matrix<T> &LU);
	bool solveLU(const Matrix<T> &A, vector<T> &x, const vector<T> &b);
	void inverse(const Matrix<T> &A ,Matrix<T> &Ai);
private:
	void inverse_aux(const Matrix<T> &B, Matrix<T> &X);
	bool inverseFlag;
	int n;
	vector<T> index;
	Matrix<T> luMatrix;
};


/**
 * @brief Constructor of the class that initialize some variables.
 */
#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE)
template<typename T>
inline MatrixDescomposition<T>::MatrixDescomposition() {
	this->n = 0;
	this->inverseFlag = false;
}

/**
 * @brief Method that take the references of the matrix and calculates LU descomposition.
 * @param A Matrix.
 * @param LU Decomposition of A.
 */
template<typename T>
void MatrixDescomposition<T>::lu(const Matrix<T>& A,
		Matrix<T>& LU) {
	this->n = A.rows();
		if(this->n != A.cols())
			throw runtime_error("'A' matrix is not square in method: void lu(const Matrix<T>& A, Matrix<T>& LU)");
		this->index.clear();
		LU = A;

		const T SMALL = 1.0e-40;
		int i, i_max, j, k;
		T big, tmp;
		vector<T> scaling;
		for(i = 0; i < this->n; i++){
			big = T(0);
			for(j = 0; j < this->n; j++){
				if((tmp = abs(LU(i, j))) > big)
					big = tmp;
			}
			if(abs(big) < numeric_limits<T>::epsilon()){
				throw runtime_error("Singular matrix in method: void lu(const Matrix<T>& A, Matrix<T>& LU)");
			}
			scaling.push_back(T(1)/big);
		}

		for(k = 0; k < this->n; k++){
			big = T(0);
			for(i = k; i < this->n; i++){
				tmp = scaling.at(i) * abs(LU(i,k));
				if(tmp > big){
					big = tmp;
					i_max = i;
				}
			}
			if(k != i_max){

				for(j = 0; j < this->n; j++){
					tmp = LU(i_max, j);
					LU(i_max, j) = LU(k, j);
					LU(k, j) = tmp;
				}
				scaling.at(i_max) = scaling.at(k);
			}
			this->index.push_back(i_max);
			if(abs(LU[k][k]) < numeric_limits<T>::epsilon())
				LU[k][k] = SMALL;
			for(i = k+1;i < this->n; i++){
				tmp = LU[i][k] /= LU[k][k];
				for(j= k+1; j < this->n; j++)
					LU[i][j] -= tmp*LU[k][j];
			}

	}
}

/**
 * @brief Lu decomposition with optimization for float
 * @param A
 * @param LU
 */
template<>
void MatrixDescomposition<float>::lu(const Matrix<float>& A,
		Matrix<float>& LU) {
	cout << "Optimized float" << endl;
	this->n = A.rows();
	if(this->n != A.cols())
		throw runtime_error("'A' matrix is not square in method: void lu(const Matrix<T>& A, Matrix<T>& LU)");
	this->index.clear();
	LU = A;

	const float SMALL = 1.0e-40;
	int i, i_max, j, k;
	float big, tmp;
	vector<float> scaling;
	for(i = 0; i < this->n; i++){
		big = float(0);
		for(j = 0; j < this->n; j++){
			if((tmp = abs(LU(i, j))) > big)
				big = tmp;
		}
		if(abs(big) < numeric_limits<float>::epsilon()){
			throw runtime_error("Singular matrix in method: void lu(const Matrix<T>& A, Matrix<T>& LU)");
		}
		scaling.push_back(float(1)/big);
	}

	for(k = 0; k < this->n; k++){
		big = float(0);
		for(i = k; i < this->n; i++){
			tmp = scaling.at(i) * abs(LU(i,k));
			if(tmp > big){
				big = tmp;
				i_max = i;
			}
		}
		if(k != i_max){
			__m128 * zeroArr = (__m128 *) calloc(4,  sizeof(float));
			float * tmpData = (float*)calloc(LU.dcols(),sizeof(float));
			memcpy(tmpData,LU[i_max],sizeof(float)*LU.dcols());
			__m128 * tmpSSE1 = (__m128 *) tmpData;
			__m128 * tmpSSE2 = (__m128 *) LU[i_max];
			__m128 * tmpSSE3 = (__m128 *) LU[k];
			for(j = 0; j < (LU.dcols()/4); j++){
				tmpSSE2[j] = _mm_add_ps(tmpSSE3[j],zeroArr[0]);
				tmpSSE3[j] = _mm_add_ps(tmpSSE1[j],zeroArr[0]);
			}
			free(zeroArr);
			free(tmpData);
			scaling.at(i_max) = scaling.at(k);
		}
		this->index.push_back(i_max);
		if(abs(LU[k][k]) < numeric_limits<float>::epsilon())
			LU[k][k] = SMALL;
		for(i = k+1;i < this->n; i++){
			tmp = LU[i][k] /= LU[k][k];
			int var = k+1;
			if(var%4 != 0){//Verified if the index can be optimized, if not iterates until it can
				for(j= k+1; j%4 != 0; j++)
					LU[i][j] -= tmp*LU[k][j];
				var = j;
			}
			if(var < this->n){
				__m128 * rowsSSE1 = (__m128 *) LU[i];
				__m128 * rowsSSE2 = (__m128 *) LU[k];
				float *tmpVec = (float*) calloc(4,  sizeof(float));
				tmpVec[0] = tmpVec[1] = tmpVec[2] = tmpVec[3] = tmp;
				__m128 * vecSSE = (__m128 *) tmpVec;
				for(int y = var; y < LU.dcols(); y+=4)
					rowsSSE1[y/4] = _mm_sub_ps (rowsSSE1[y/4], _mm_mul_ps(rowsSSE2[y/4],vecSSE[0]));
				free(tmpVec);
			}

		}
	}
}

/**
 * @brief Lu decomposition with optimization for double
 * @param A
 * @param LU
 */
template<>
void MatrixDescomposition<double>::lu(const Matrix<double>& A,
		Matrix<double>& LU) {
	this->n = A.rows();
	if(this->n != A.cols())
		throw runtime_error("'A' matrix is not square in method: void lu(const Matrix<T>& A, Matrix<T>& LU)");
	this->index.clear();
	LU = A;

	const double SMALL = 1.0e-40;
	int i, i_max, j, k;
	double big, tmp;
	vector<double> scaling;
	for(i = 0; i < this->n; i++){
		big = double(0);
		for(j = 0; j < this->n; j++){
			if((tmp = abs(LU(i, j))) > big)
				big = tmp;
		}
		if(abs(big) < numeric_limits<double>::epsilon()){
			throw runtime_error("Singular matrix in method: void lu(const Matrix<T>& A, Matrix<T>& LU)");
		}
		scaling.push_back(double(1)/big);
	}

	for(k = 0; k < this->n; k++){
		big = double(0);
		for(i = k; i < this->n; i++){
			tmp = scaling.at(i) * abs(LU(i,k));
			if(tmp > big){
				big = tmp;
				i_max = i;
			}
		}
		if(k != i_max){//Optimized pivoting
			__m128d * zeroArr = (__m128d *) calloc(2,  sizeof(double));
			double * tmpData = (double*)calloc(LU.dcols(),sizeof(double));
			memcpy(tmpData,LU[i_max],sizeof(double)*LU.dcols());
			__m128d * tmpSSE1 = (__m128d *) tmpData;
			__m128d * tmpSSE2 = (__m128d *) LU[i_max];
			__m128d * tmpSSE3 = (__m128d *) LU[k];
			for(j = 0; j < (LU.dcols()/2); j++){
				tmpSSE2[j] = _mm_add_pd(tmpSSE3[j],zeroArr[0]);
				tmpSSE3[j] = _mm_add_pd(tmpSSE1[j],zeroArr[0]);
			}
			free(zeroArr);
			free(tmpData);
			scaling.at(i_max) = scaling.at(k);
		}
		this->index.push_back(i_max);
		if(abs(LU[k][k]) < numeric_limits<double>::epsilon())
			LU[k][k] = SMALL;
		for(i = k+1;i < this->n; i++){
			tmp = LU[i][k] /= LU[k][k];
			int var = k+1;
			if(var%2 != 0){//Verifies if the index can be optimized, if not iterates until it can
				for(j= k+1; j%2 != 0; j++)
					LU[i][j] -= tmp*LU[k][j];
				var = j;
			}
			if(var < this->n){
				__m128d * rowsSSE1 = (__m128d *) LU[i];
				__m128d * rowsSSE2 = (__m128d *) LU[k];
				double *tmpVec = (double*) calloc(2,  sizeof(double));
				tmpVec[0] = tmpVec[1] = tmpVec[2] = tmpVec[3] = tmp;
				__m128d * vecSSE = (__m128d *) tmpVec;
				for(int y = var; y < LU.dcols(); y+=2)
					rowsSSE1[y/2] = _mm_sub_pd(rowsSSE1[y/2], _mm_mul_pd(rowsSSE2[y/2],vecSSE[0]));
				free(tmpVec);
			}

		}
	}
}

/**
 * @brief Solves an systems of linear equations by LU decomposition.
 * @param A Matrix.
 * @param x Vector variables.
 * @param b Right side vector.
 * @return true if the method could solve.
 */
template<typename T>
inline bool MatrixDescomposition<T>::solveLU(const Matrix<T>& A,
		vector<T>& x, const vector<T>& b) {
	bool result = true;
	int i, ip, j;
	int ii = 0;
	T sum;
	if(!this->inverseFlag)
		this->lu(A, this->luMatrix);

	if(b.size() != this->n)
		throw runtime_error("The rows dimension is not correct in void LUDescomposition<T>::solve(vector<T>& b, vector<T>& x)");

	x = b;
	for(i = 0; i < this->n; i++){
		ip = this->index.at(i);
		sum = x.at(ip);
		x.at(ip) = x.at(i);
		if(ii != 0){
			for(j = ii-1; j < i; j++){
				sum -= (this->luMatrix[i][j])*(x.at(j));
			}
		}else if(abs(sum) > numeric_limits<T>::epsilon())
			ii = i+1;
		x.at(i) = sum;
	}for(i = n-1; i >= 0; i--){
		sum = x.at(i);
		for(j = i+1; j < n; j++)
			sum -= this->luMatrix[i][j]*x.at(j);
		x.at(i) = sum/(this->luMatrix[i][i]);
	}
	return result;
}

/**
 * @brief Method that calculates the inverse of A Matrix.
 * @param A Matrix.
 * @param Ai Inverse of the A matrix.
 */
template<typename T>
inline void MatrixDescomposition<T>::inverse(const Matrix<T>& A,
		Matrix<T>& Ai) {
	this->lu(A,this->luMatrix);
	this->inverseFlag = true;
	Ai = anpi::Matrix<T>(this->n, this->n, T(0));
	for(int i = 0; i < this->n; i++){
		Ai[i][i] = 1;
	}
	this->inverse_aux(Ai,Ai);
	this->inverseFlag = false;
}

/**
 * @brief Auxiliary method that calculates the inverse by using LU descomposition.
 * @param B Identity matrix.
 * @param X Inverse of A.
 */
template<typename T>
inline void MatrixDescomposition<T>::inverse_aux(const Matrix<T>& B,
		Matrix<T>& X) {
	int i, j;
	int m = B.cols();
	if(B.rows() != this->n || X.rows() != this->n || B.cols() != X.cols()){
		throw runtime_error("Bad sizes for the matrix in void LUDescomposition<T>::inverse_aux(const anpi::Matrix<T>& B, anpi::Matrix<T>& X)");
	}
	vector<T> tmp;
	for(j = 0; j < m; j++){
		tmp.clear();
		for(i = 0; i < this->n; i++)
			tmp.push_back(B[i][j]);
		this->solveLU(B, tmp, tmp);
		for(i = 0; i < this->n; i++)
			X[i][j] = tmp.at(i);
	}
}
#else
template<typename T>
inline MatrixDescomposition<T>::MatrixDescomposition() {
	this->n = 0;
	this->inverseFlag = false;
}

/**
 * @brief Method that take the references of the matrix and calculates LU descomposition.
 * @param A Matrix.
 * @param LU Decomposition of A.
 */
template<typename T>
inline void MatrixDescomposition<T>::lu(const Matrix<T>& A,
		Matrix<T>& LU) {
	this->n = A.rows();
		if(this->n != A.cols())
			throw runtime_error("'A' matrix is not square in method: void lu(const Matrix<T>& A, Matrix<T>& LU)");
		this->index.clear();
		LU = A;

		const T SMALL = 1.0e-40;
		int i, i_max, j, k;
		T big, tmp;
		vector<T> scaling;
		for(i = 0; i < this->n; i++){
			big = T(0);
			for(j = 0; j < this->n; j++){
				if((tmp = abs(LU(i, j))) > big)
					big = tmp;
			}
			if(abs(big) < numeric_limits<T>::epsilon()){
				throw runtime_error("Singular matrix in method: void lu(const Matrix<T>& A, Matrix<T>& LU)");
			}
			scaling.push_back(T(1)/big);
		}

		for(k = 0; k < this->n; k++){
			big = T(0);
			for(i = k; i < this->n; i++){
				tmp = scaling.at(i) * abs(LU(i,k));
				if(tmp > big){
					big = tmp;
					i_max = i;
				}
			}
			if(k != i_max){

				for(j = 0; j < this->n; j++){
					tmp = LU(i_max, j);
					LU(i_max, j) = LU(k, j);
					LU(k, j) = tmp;
				}
				scaling.at(i_max) = scaling.at(k);
			}
			this->index.push_back(i_max);
			if(abs(LU[k][k]) < numeric_limits<T>::epsilon())
				LU[k][k] = SMALL;
			for(i = k+1;i < this->n; i++){
				tmp = LU[i][k] /= LU[k][k];
				for(j= k+1; j < this->n; j++)
					LU[i][j] -= tmp*LU[k][j];
			}

	}
}

/**
 * @brief Solves an systems of linear equations by LU decomposition.
 * @param A Matrix.
 * @param x Vector variables.
 * @param b Right side vector.
 * @return true if the method could solve.
 */
template<typename T>
inline bool MatrixDescomposition<T>::solveLU(const Matrix<T>& A,
		vector<T>& x, const vector<T>& b) {
	bool result = true;
	int i, ip, j;
	int ii = 0;
	T sum;
	if(!this->inverseFlag)
		this->lu(A, this->luMatrix);

	if(b.size() != this->n)
		throw runtime_error("The rows dimension is not correct in void LUDescomposition<T>::solve(vector<T>& b, vector<T>& x)");

	x = b;
	for(i = 0; i < this->n; i++){
		ip = this->index.at(i);
		sum = x.at(ip);
		x.at(ip) = x.at(i);
		if(ii != 0){
			for(j = ii-1; j < i; j++){
				sum -= (this->luMatrix[i][j])*(x.at(j));
			}
		}else if(abs(sum) > numeric_limits<T>::epsilon())
			ii = i+1;
		x.at(i) = sum;
	}for(i = n-1; i >= 0; i--){
		sum = x.at(i);
		for(j = i+1; j < n; j++)
			sum -= this->luMatrix[i][j]*x.at(j);
		x.at(i) = sum/(this->luMatrix[i][i]);
	}
	return result;
}

/**
 * @brief Method that calculates the inverse of A Matrix.
 * @param A Matrix.
 * @param Ai Inverse of the A matrix.
 */
template<typename T>
inline void MatrixDescomposition<T>::inverse(const Matrix<T>& A,
		Matrix<T>& Ai) {
	this->lu(A,this->luMatrix);
	this->inverseFlag = true;
	Ai = anpi::Matrix<T>(this->n, this->n, T(0));
	for(int i = 0; i < this->n; i++){
		Ai[i][i] = 1;
	}
	this->inverse_aux(Ai,Ai);
	this->inverseFlag = false;
}

/**
 * @brief Auxiliary method that calculates the inverse by using LU descomposition.
 * @param B Identity matrix.
 * @param X Inverse of A.
 */
template<typename T>
inline void MatrixDescomposition<T>::inverse_aux(const Matrix<T>& B,
		Matrix<T>& X) {
	int i, j;
	int m = B.cols();
	if(B.rows() != this->n || X.rows() != this->n || B.cols() != X.cols()){
		throw runtime_error("Bad sizes for the matrix in void LUDescomposition<T>::inverse_aux(const anpi::Matrix<T>& B, anpi::Matrix<T>& X)");
	}
	vector<T> tmp;
	for(j = 0; j < m; j++){
		tmp.clear();
		for(i = 0; i < this->n; i++)
			tmp.push_back(B[i][j]);
		this->solveLU(B, tmp, tmp);
		for(i = 0; i < this->n; i++)
			X[i][j] = tmp.at(i);
	}
}
#endif
} /* namespace anpi */
#endif /* DESCOMPOSITION_MATRIXDESCOMPOSITION_H_ */
