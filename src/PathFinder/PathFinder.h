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
#include <limits>
#include <cmath>
#include "../plot/plotpy.h"
#include "math.h"

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
	vector<T> xPosition;
	vector<T> yPosition;

	//Image matrix
	Mat imageMatrix;
	int imgRows;
	int imgCols;

	//Initial and final
	int initialRow; 
	int initialCol; 
	int finalRow;
	int finalCol;

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

	const vector<Point> * getPathPoints();

  private:
	void getNodeEquations();
	void getMeshEquations();
	void getXAxisMatrix();
	void getYAxisMatrix();
	void normalize();
	void getPathPositions(T alpha);
	T bilinearInterpolation(T q11, T q12, T q21, T q22, T x1, T x2, T y1, T y2, T x, T y);
};

template <typename T>
PathFinder<T>::PathFinder(int initialRow, int initialCol, int finalRow, int finalCol, Mat map)
{

	this->imageMatrix = map;
	imgRows = map.rows;
	imgCols = map.cols;

	this->initialRow = initialRow;
	this->initialCol = initialCol;
	this->finalRow =  finalRow;
	this->finalCol =  finalCol;

	indexMap = new IndexMap(imgRows, imgCols);

	cols = 2 * imgRows * imgCols - (imgRows + imgCols); //Incognites number
	rows = cols;										//Equations number

	A = Matrix<T>(rows, cols, T(0), Matrix<T>::Padded);
	x_axis = Matrix<T>(imgRows, imgCols, T(0), Matrix<T>::Padded);
	y_axis = Matrix<T>(imgRows, imgCols, T(0), Matrix<T>::Padded);
	b = *(new vector<T>(rows));

	getNodeEquations();
	getMeshEquations();

	MatrixDescomposition<T> *solver = new MatrixDescomposition<T>();
	solver->solveLU(A, x, b);
	getXAxisMatrix();
	getYAxisMatrix();
	normalize();
	getPathPositions(1);
	printMatrixx(x_axis);
	cout<<endl;
	printMatrixx(y_axis);

	plotpy::Plot2d<T> plt;

	plt.initialize(1);

	plt.settitle("Path found");

	plt.quiver(x_axis, y_axis);

	plt.showallplots();

}

template <typename T>
void PathFinder<T>::getNodeEquations()
{

	int position = 0;
	int equation_row = 0;
	int flag = 0;
	int startJ = 0;
	int endJ = 0;
	int initialPosition = (this->initialCol) + this->imgCols * this->initialRow;
	int finalPosition = (finalCol) + this->imgCols *this-> finalRow;

	

	//Check if 0,0 edge is free
	if (!(this->initialRow == 0 && this->initialCol == 0) && !(finalRow == 0 && this-> finalCol == 0))
	{
		flag = 0;

		//Input current
		initialPosition = (initialCol - 1) + this->imgCols * this-> initialRow;
		//Output current
		finalPosition = (finalCol - 1) + this->imgCols * this->finalRow;
	}
	//Check if 0,cols-1 edge is free (Upper right)
	else if (!(initialRow == 0 && this-> initialCol == this->imgCols - 1) && !(finalRow == 0 && this->finalCol == this->imgCols - 1))
	{
		flag = 1;
		if (!(initialRow == 0 && this-> initialCol < this->imgCols - 1 ))
		{
			//Input current
			initialPosition = (initialCol - 1) + this->imgCols *this-> initialRow;
		}
		if (!(finalRow == 0 && this->finalCol < this->imgCols - 1))
		{
			//Output current
			finalPosition = (finalCol - 1) + this->imgCols *this->finalRow;
		}
	}
	//rows-1,0 edge is free (Lower left)
	else
	{
		flag = 2;
		if (initialRow == this->imgRows - 1 && this->initialCol > 0 )
		{
			//Input current
			initialPosition = (initialCol - 1) + this->imgCols * this-> initialRow;
		}
		if (finalRow == this->imgRows - 1 && this->finalCol > 0)
		{
			//Output current
			finalPosition = (finalCol - 1) + this->imgCols * this->finalRow;
		}
	}


	this->b.at(initialPosition) = 1;
	this->b.at(finalPosition) = -1;


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
					this->x_axis[i][j] = -(xr + xl) ; 
				}
			}
}

template<typename T>
void PathFinder<T>::getYAxisMatrix(){
			int position;
			T yu,yb;
			for (int i = 0; i < this->imgRows; i++)
			{
				for (int j = 0; j < this->imgCols; j++)
				{
					position = 0;
					yu = yb = 0;
					//Down current
					if(i < this->imgRows-1){
						position = indexMap->getXFromNodes(i,j,i+1,j);
						yb = x[position];
					}
					//Up current
					if(i > 0){
						position = indexMap->getXFromNodes(i,j,i-1,j);
						yu = x[position];
					}
					this->y_axis[i][j] = -(yu + yb); 		
				}
			}
}


template<typename T>
void PathFinder<T>::normalize(){
	T big, tmp = T(0);
	for (int i = 0; i < imgRows; i++){
		for(int j = 0; j < imgCols; j++){
			tmp = sqrt(x_axis(i,j)*x_axis(i,j)+y_axis(i,j)*y_axis(i,j));
			if(tmp > big){
					big = tmp;
			}
		}
	}
	
	for(int i = 0; i < imgRows; i++){
		for(int j = 0; j < imgCols; j++){
			x_axis(i,j)=x_axis(i,j)/big;
			y_axis(i,j)=y_axis(i,j)/big;
		}
	}
}


template <typename T>
const vector<Point> * PathFinder<T>::getPathPoints()
{
	vector<Point> * points = new vector<Point>;
	vector<T> solutions(this->x);
	
	int actualRow = this->initialRow;
	int actualCol = this->initialCol;
	T maxCurrent = 0;
	T nextCurrent = 0;
	int position = 0;
	int nextX = 0;
	int nextY = 0;

	points->push_back(Point(actualRow,actualCol)); //Initial node
	
	
	while(!(actualRow == this->finalRow && actualCol == this->finalCol)){
		maxCurrent = 0;

		//Right current
		if (actualCol < this->imgCols - 1)
		{	
			position = indexMap->getXFromNodes(actualRow, actualCol , actualRow, actualCol  + 1);			
			nextCurrent = abs(solutions.at(position));
			if(nextCurrent > maxCurrent ){
				maxCurrent = nextCurrent;
				nextX = actualRow;
				nextY = actualCol + 1;
				solutions.at(position) = 0;
			}
			
		}
		//Down current
		if (actualRow < this->imgRows - 1)
		{
			position = indexMap->getXFromNodes(actualRow, actualCol , actualRow + 1, actualCol );
			nextCurrent = abs(solutions.at(position));
			if(nextCurrent > maxCurrent ){
				maxCurrent = nextCurrent;
				nextX = actualRow + 1;
				nextY = actualCol;
				solutions.at(position) = 0;
			}
		}
		//Left current
		if (actualCol > 0)
		{
			position = indexMap->getXFromNodes(actualRow, actualCol , actualRow, actualCol  - 1);
			nextCurrent = abs(solutions.at(position));
			if(nextCurrent > maxCurrent ){
				maxCurrent = nextCurrent;
				nextX = actualRow;
				nextY = actualCol - 1;
				solutions.at(position) = 0;
			}
		}
		//Up current
		if (actualRow > 0)
		{
			position = indexMap->getXFromNodes(actualRow, actualCol , actualRow - 1, actualCol);
			nextCurrent = abs(solutions.at(position));
			if(nextCurrent > maxCurrent ){
				maxCurrent = nextCurrent;
				nextX = actualRow - 1;
				nextY = actualCol;
				solutions.at(position) = 0;
			}
		}

		actualRow = nextX;
		actualCol = nextY;
		points->push_back(Point(nextX,nextY));

	}//End while
	return points;

}

template <typename T>
void  PathFinder<T>::getPathPositions(T alpha)
{
	T row = this->initialRow;
	T col = this->initialCol;
	T tmpRow, tmpCol;
	this->xPosition.push_back(row);
	this->yPosition.push_back(col);
	tmpRow = row + alpha * x_axis(row,col);
	tmpCol = col + alpha * y_axis(row,col);
		
	
	//while (row != this.finalRow and col != this.finalCol){
	for (int i = 0; i<5;i++){
		T q11x = x_axis(row,col);
		T q12x = x_axis(row,col+1);
		T q21x = x_axis(row+1,col);
		T q22x = x_axis(row+1,col+1);
		T dx = bilinearInterpolation(q11x,q12x,q21x,q22x,col,col+1,row,row+1,tmpCol,tmpRow);
		T q11y = y_axis(row,col);
		T q12y = y_axis(row,col+1);
		T q21y = y_axis(row+1,col);
		T q22y = y_axis(row+1,col+1);
		T dy = bilinearInterpolation(q11y,q12y,q21y,q22y,col,col+1,row,row+1,tmpCol,tmpRow);
		cout << dx<<" "<<dy <<endl;
		tmpRow = tmpRow + alpha * dx;
		tmpCol = tmpCol + alpha * dy;
		cout << tmpRow<<" saddas "<<tmpCol <<endl;
		this->xPosition.push_back(tmpRow);
		this->yPosition.push_back(tmpCol);
		row = floor(tmpRow);
		col = floor(tmpCol);
		//cout << "x " << x_axis(row,col) << " y " <<y_axis(row,col)<<endl;
		//cout << "row " <<row << " col "<<col<<endl;
	}
}

template <typename T>
T PathFinder<T>::bilinearInterpolation(T q11, T q12, T q21, T q22, T x1, T x2, T y1, T y2, T x, T y) 
{
    T x2x1, y2y1, x2x, y2y, yy1, xx1, x1x2,y1y2;
    x2x1 = x2 - x1;
    y2y1 = y2 - y1;
    x1x2 = -x2x1;
    y1y2 = -y2y1;
    x2x = x2 - x;
    y2y = y2 - y;
    yy1 = y - y1;
    xx1 = x - x1;
    return 1.0 / (x2x1 * y2y1) * (
        q11 * x2x * y2y +
        q21 * xx1 * y2y +
        q12 * x2x * yy1 +
        q22 * xx1 * yy1
    );
}


#endif /* PATHFINDER_PATHFINDER_H_ */
