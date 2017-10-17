/*
 * PathFinder.h
 *
 *  Created on: 8 de oct. de 2017
 *      Author: kevin
 */

#ifndef PATHFINDER_PATHFINDER_H_
#define PATHFINDER_PATHFINDER_H_

#include "../Matrix/Matrix.hpp"
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
	Matrix<int>  A;
	Matrix<T>  x; //Solution vector
	Matrix<int>  b;

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

	const Matrix<int>& getA() const {
		return A;
	}

	const Matrix<int>& getB() const {
		return b;
	}

	const Matrix<T>& getX() const {
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
	rows = 2*imgRows*imgCols -(imgRows+imgCols)+1; //Equations number

	A = Matrix<int>(rows, cols, 0, Matrix<int>::Padded);
	x = Matrix<T>(rows, 1, T(0), Matrix<T>::Padded);
	b = Matrix<int>(rows, 1, 0, Matrix<int>::Padded);

	//Input current
	int initialPosition = (initialCol) + imgCols*initialRow;
	b(initialPosition,0) = 1;

	//Output current
	int finalPosition = (finalCol) + imgCols*finalRow;
	b(finalPosition,0) = -1;

	getNodeEquations();
	getMeshEquations();


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
	int value = 0;
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
