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
#include <opencv2/core.hpp>
using namespace anpi;
using namespace cv;
using namespace std;

template <typename T>
class PathFinder
{
  private:
	//Ax = b system
	Matrix<T> A;
	vector<T> x; //Solution vector
	vector<T> b;

	//Map of vector-nodes
	IndexMap *indexMap;

	//Matrix A properties
	int rows;
	int cols;

	Matrix<T> x_axis;
	Matrix<T> y_axis;

	//Image matrix
	Mat imageMatrix;
	int imgRows;
	int imgCols;

  public:
	PathFinder(int initialRow, int initialCol, int finalRow, int finalCol, Mat map);

	const Matrix<T> &getA() const
	{
		return A;
	}

	const vector<T> &getB() const
	{
		return b;
	}

	const vector<T> &getX() const
	{
		return x;
	}

	const Matrix<T> &getXAxis() const
	{
		return x_axis;
	}

	const vector<Point> &getPathPoints();

  private:
	void getNodeEquations(int initialRow, int initialCol, int finalRow, int finalCol);
	void getMeshEquations();
	void getXAxisMatrix();
};

template <typename T>
PathFinder<T>::PathFinder(int initialRow, int initialCol, int finalRow, int finalCol, Mat map)
{

	this->imageMatrix = map;
	imgRows = map.rows;
	imgCols = map.cols;

	indexMap = new IndexMap(imgRows, imgCols);

	cols = 2 * imgRows * imgCols - (imgRows + imgCols); //Incognites number
	rows = cols;										//Equations number

	A = Matrix<T>(rows, cols, T(0), Matrix<T>::Padded);
	x_axis = Matrix<T>(imgRows, imgCols, Matrix<T>::Padded);
	y_axis = Matrix<T>(imgRows, imgCols, Matrix<T>::Padded);
	b = *(new vector<T>(rows));

	getNodeEquations(initialRow, initialCol, finalRow, finalCol);
	getMeshEquations();

	MatrixDescomposition<T> *solver = new MatrixDescomposition<T>();
	solver->solveLU(A, x, b);
	getXAxisMatrix();
}

template <typename T>
void PathFinder<T>::getNodeEquations(int initialRow, int initialCol, int finalRow, int finalCol)
{

	int position = 0;
	int equation_row = 0;
	int flag = 0;
	int startJ = 0;
	int endJ = 0;
	int initialPosition = (initialCol) + this->imgCols * initialRow;
	int finalPosition = (finalCol) + this->imgCols * finalRow;

	//Check if 0,0 edge is free
	if (!(initialRow == 0 && initialCol == 0) && !(finalRow == 0 && finalCol == 0))
	{
		flag = 0;

		//Input current
		initialPosition = (initialCol - 1) + this->imgCols * initialRow;
		//Output current
		finalPosition = (finalCol - 1) + this->imgCols * finalRow;
	}
	//Check if 0,cols-1 edge is free (Upper right)
	else if (!(initialRow == 0 && initialCol == this->imgCols - 1) && !(finalRow == 0 && finalCol == this->imgCols - 1))
	{
		flag = 1;
		if (!(initialRow == 0 && initialCol == 0))
		{
			//Input current
			initialPosition = (initialCol - 1) + this->imgCols * initialRow;
		}
		if (!(finalRow == 0 && finalCol == 0))
		{
			//Output current
			finalPosition = (finalCol - 1) + this->imgCols * finalRow;
		}
	}
	//rows-1,0 edge is free (Lower left)
	else
	{
		flag = 2;
		if (initialRow == this->imgRows - 1 && initialCol == this->imgCols - 1)
		{
			//Input current
			initialPosition = (initialCol - 1) + this->imgCols * initialRow;
		}
		if (finalRow == this->imgRows - 1 && finalCol == this->imgCols - 1)
		{
			//Output current
			finalPosition = (finalCol - 1) + this->imgCols * finalRow;
		}
	}

	cout <<"ip "<<initialPosition<<endl;
	cout <<"fp "<<finalPosition<<endl;

	b.at(initialPosition) = 1;
	b.at(finalPosition) = -1;

	for (int i = 0; i < this->imgRows; i++)
	{

		//Remove one free node
		if ((flag == 0 && i == 0) || (flag == 2 && i == this->imgRows - 1))
		{
			startJ = 1;
			endJ = this->imgCols;
		}
		else if (flag == 1 && i == 0)
		{
			startJ = 0;
			endJ = this->imgCols - 1;
		}
		else
		{
			startJ = 0;
			endJ = this->imgCols;
		}

		for (int j = startJ; j < endJ; j++)
		{

			//Right current
			if (j < this->imgCols - 1)
			{
				position = indexMap->getXFromNodes(i, j, i, j + 1);
				A(equation_row, position) = -1;
			}
			//Down current
			if (i < this->imgRows - 1)
			{
				position = indexMap->getXFromNodes(i, j, i + 1, j);
				A(equation_row, position) = -1;
			}

			//Left current
			if (j > 0)
			{
				position = indexMap->getXFromNodes(i, j, i, j - 1);
				A(equation_row, position) = 1;
			}
			//Up current
			if (i > 0)
			{
				position = indexMap->getXFromNodes(i, j, i - 1, j);
				A(equation_row, position) = 1;
			}

			equation_row++;
		}
	}
}

template <typename T>
void PathFinder<T>::getMeshEquations()
{
	int position = 0;
	int equation_row = this->imgRows * this->imgCols - 1;
	T value = 0;
	for (int i = 0; i < this->imgRows - 1; i++)
	{
		for (int j = 0; j < this->imgCols - 1; j++)
		{
			if ((int)this->imageMatrix.template at<uchar>(i, j) > 250 and (int) this->imageMatrix.template at<uchar>(i, j + 1) > 250)
			{
				value = 1;
				position = indexMap->getXFromNodes(i, j, i, j + 1);
				A(equation_row, position) = value;
			}
			else
			{
				value = 1000000;
				position = indexMap->getXFromNodes(i, j, i, j + 1);
				A(equation_row, position) = value;
			}

			if ((int)this->imageMatrix.template at<uchar>(i, j) > 250 and (int) this->imageMatrix.template at<uchar>(i + 1, j) > 250)
			{
				value = 1;
				position = indexMap->getXFromNodes(i, j, i + 1, j);
				A(equation_row, position) = -value;
			}
			else
			{
				value = 1000000;
				position = indexMap->getXFromNodes(i, j, i + 1, j);
				A(equation_row, position) = -value;
			}

			if ((int)this->imageMatrix.template at<uchar>(i, j + 1) > 250 and (int) this->imageMatrix.template at<uchar>(i + 1, j + 1) > 250)
			{
				value = 1;
				position = indexMap->getXFromNodes(i, j + 1, i + 1, j + 1);
				A(equation_row, position) = value;
			}
			else
			{
				value = 1000000;
				position = indexMap->getXFromNodes(i, j + 1, i + 1, j + 1);
				A(equation_row, position) = value;
			}

			if ((int)this->imageMatrix.template at<uchar>(i + 1, j) > 250 and (int) this->imageMatrix.template at<uchar>(i + 1, j + 1) > 250)
			{
				value = 1;
				position = indexMap->getXFromNodes(i + 1, j, i + 1, j + 1);
				A(equation_row, position) = -value;
			}
			else
			{
				value = 1000000;
				position = indexMap->getXFromNodes(i + 1, j, i + 1, j + 1);
				A(equation_row, position) = -value;
			}
			equation_row++;
		}
	}

}


template<typename T>
void PathFinder<T>::getXAxisMatrix(){
			int position;
			T xl,xr;
			for (int i = 0; i < this->imgRows; i++)
			{
				for (int j = 0; j < this->imgCols; j++)
				{
					position = 0;
					xl = xr = 0;
					
					cout << i << " " << j <<endl;
					//Left current
					if(j > 0){
						position = indexMap->getXFromNodes(i,j,i,j-1);
						xl = x[position];
					}

					//Right current
					if(j < this->imgCols-1){
						position = indexMap->getXFromNodes(i,j,i,j+1);
						xr = x[position];
					}

					cout << "Left "<<xl << " Rigth "<<xr <<endl;
					cout << "culo "<<xl +xr<<endl;
					this->x_axis[i][j] = xl + xr; 
					cout << x_axis(i,j) <<endl;

					/**
					//Down current
					if(i < this->imgRows-1){
						position = indexMap->getXFromNodes(i,j,i+1,j);
						A(equation_row,position) = -1;
					}

					//Up current
					if(i > 0){
						position = indexMap->getXFromNodes(i,j,i-1,j);
						A(equation_row,position) = 1;
					}
					*/
				}
			}
			printMatrixx(x_axis);
}



//const vector<Point> &getPathPoints()
	//{
	//}

#endif /* PATHFINDER_PATHFINDER_H_ */
