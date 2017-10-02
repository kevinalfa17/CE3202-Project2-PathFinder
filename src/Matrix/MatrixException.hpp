/**
 * @file MatrixException.h
 * @brief Class that throw custom exception
 * @author Kevin Alfaro
 * @date 25 de sept. de 2017
 */

#ifndef MATRIX_MATRIXEXCEPTION_HPP_
#define MATRIX_MATRIXEXCEPTION_HPP_


#include <iostream>
#include <stdexcept>

class MatrixException : public std::exception
{
public:

	/**
	 * @brief Constructor by default.
	 */
	explicit MatrixException(){}

	/**
	 * @brief what() method with custom message
	 */
	virtual const char* what() const throw(){
		return "Number of columns of the first matrix is different to the number of rows of the second";
	}
};


#endif /* MATRIX_MATRIXEXCEPTION_HPP_ */
