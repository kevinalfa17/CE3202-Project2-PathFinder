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
using namespace anpi;

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



public:
	PathFinder(int initialRow, int initialCol, int finalRow, int finalCol);

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
	void getNodeEquations(int imgRows, int imgCols);

};

template<typename T>
PathFinder<T>::PathFinder(int initialRow, int initialCol, int finalRow, int finalCol){

	//Just for test, replace with OpenCV image properties
	int imgRows = 3;
	int imgCols = 3;


	indexMap = new IndexMap(imgRows,imgCols);

	cols = 2*imgRows*imgCols -(imgRows+imgCols); //Incognites number
	rows = imgRows*imgCols; //Equations number

	A = Matrix<int>(rows, cols, 0, Matrix<int>::Padded);
	//x = Matrix<T>(rows, 1, T(0), Matrix<T>::Padded);
	//b = Matrix<int>(rows, 1, 0, Matrix<int>::Padded);

	//Input current
	int initialPosition = (initialCol) + imgCols*initialRow;
	//b(initialPosition,0) = 1;

	//Output current
	int finalPosition = (finalCol) + imgCols*finalRow;
	//b(finalPosition,0) = -1;

	getNodeEquations(imgRows, imgCols);


}

template<typename T>
void PathFinder<T>::getNodeEquations(int imgRows, int imgCols){

	int position = 0;
	int equation_row = 0;

	for(int i = 0; i < imgRows; i++){
		for(int j = 0; j < imgCols; j++){

			//Right current
			if(j < imgCols-1){
				position = indexMap->getXFromNodes(i,j,i,j+1);
				A(equation_row,position) = -1;
			}
			//Down current
			if(i < imgRows-1){
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


#endif /* PATHFINDER_PATHFINDER_H_ */
