/*
 * PathFinder.h
 *
 *  Created on: 8 de oct. de 2017
 *      Author: kevin
 */

#ifndef PATHFINDER_PATHFINDER_H_
#define PATHFINDER_PATHFINDER_H_

#include "../Matrix/Matrix.hpp"
#include "../Matrix/MatrixDescomposition.h"
#include "IndexMap.h"
#include "NodePair.h"
#include <iostream>
#include <opencv2/opencv.hpp>
using namespace anpi;
using namespace cv;
using namespace std;

template<typename T>
class PathFinder{
private:

	//Ax = b system
	Matrix<T>  A;
	vector<T>  x; //Solution vector
	vector<T>  b;

	//Map of vector-nodes
	IndexMap * indexMap;

	//Matrix A properties
	int rows;
	int cols;

	//Image matrix
	Mat imageMatrix;
	int imgRows;
	int imgCols;


public:
	PathFinder(int initialRow, int initialCol, int finalRow, int finalCol, Mat map);

	const Matrix<T>& getA() const {
		return A;
	}

	const vector<T>& getB() const {
		return b;
	}

	const vector<T>& getX() const {
		return x;
	}

private:
	void getNodeEquations();
	void getMeshEquations();

};

template<typename T>
PathFinder<T>::PathFinder(int initialRow, int initialCol, int finalRow, int finalCol, Mat map){

	this->imageMatrix = map;
	cout << map << endl;
	imgRows = map.rows;
	imgCols = map.cols;


	indexMap = new IndexMap(imgRows,imgCols);

	cols = 2*imgRows*imgCols -(imgRows+imgCols); //Incognites number
	rows = 2*imgRows*imgCols -(imgRows+imgCols); //Equations number

	A = Matrix<T>(rows, cols, T(0), Matrix<T>::Padded);
	b = *(new vector<T>(rows));

	//Input current
	int initialPosition = (initialCol) + imgCols*initialRow;
	b.at(initialPosition) = 1;

	//Output current
	int finalPosition = (finalCol) + imgCols*finalRow;
	b.at(finalPosition) = -1;

	getNodeEquations();
	getMeshEquations();

	MatrixDescomposition<T> * solver = new MatrixDescomposition<T>();
	solver->solveLU(A, x, b);

}

template<typename T>
void PathFinder<T>::getNodeEquations(){

	int position = 0;
	int equation_row = 0;

	for(int i = 0; i < this->imgRows; i++){
		for(int j = 0; j < this->imgCols; j++){

			//Right current
			if(j < this->imgCols-1){
				position = indexMap->getXFromNodes(i,j,i,j+1);
				A(equation_row,position) = -1;
			}
			//Down current
			if(i < this->imgRows-1){
				position = indexMap->getXFromNodes(i,j,i+1,j);
				A(equation_row,position) = -1;
			}

			//Left current
			if(j > 0){
				position = indexMap->getXFromNodes(i,j,i,j-1);
				A(equation_row,position) = 1;
			}
			//Up current
			if(i > 0){
				position = indexMap->getXFromNodes(i,j,i-1,j);
				A(equation_row,position) = 1;
			}

			equation_row++;

		}
	}
}

template<typename T>
void PathFinder<T>::getMeshEquations(){
	int position = 0;
	int equation_row = this->imgRows * this->imgCols-1;
	T value = 0;
	for(int i = 0; i < this->imgRows-1; i++){
		for(int j = 0; j < this->imgCols-1; j++){
			if ((int)this->imageMatrix.template at<uchar>(i, j) > 250 and (int)this->imageMatrix.template at<uchar>(i, j+1) > 250){
				value = 1;
				position = indexMap->getXFromNodes(i,j,i,j+1);
				A(equation_row,position) = value;
			} else {
				value = 1000000;
				position = indexMap->getXFromNodes(i,j,i,j+1);
				A(equation_row,position) = value;
			}

			if ((int)this->imageMatrix.template at<uchar>(i, j) > 250 and (int)this->imageMatrix.template at<uchar>(i+1, j) > 250){
				value = 1;
				position = indexMap->getXFromNodes(i,j,i+1,j);
				A(equation_row,position) = -value;
			} else {
				value = 1000000;
				position = indexMap->getXFromNodes(i,j,i+1,j);
				A(equation_row,position) = -value;
			}

			if ((int)this->imageMatrix.template at<uchar>(i, j+1) > 250 and (int)this->imageMatrix.template at<uchar>(i+1, j+1) > 250){
				value = 1;
				position = indexMap->getXFromNodes(i,j+1,i+1,j+1);
				A(equation_row,position) = value;
			} else {
				value = 1000000;
				position = indexMap->getXFromNodes(i,j+1,i+1,j+1);
				A(equation_row,position) = value;
			}

			if ((int)this->imageMatrix.template at<uchar>(i+1, j) > 250 and (int)this->imageMatrix.template at<uchar>(i+1, j+1) > 250){
				value = 1;
				position = indexMap->getXFromNodes(i+1,j,i+1,j+1);
				A(equation_row,position) = -value;
			} else {
				value = 1000000;	
				position = indexMap->getXFromNodes(i+1,j,i+1,j+1);
				A(equation_row,position) = -value;
			}
			equation_row++;

		}
	}

}

#endif /* PATHFINDER_PATHFINDER_H_ */
