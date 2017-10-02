/**
 * @file PaddingWarningException.h
 * @brief Class that throw custom warning exception
 * @author Kevin Alfaro
 * @date 25 de sept. de 2017
 */

#ifndef MATRIX_PADDINGWARNINGEXCEPTION_HPP_
#define MATRIX_PADDINGWARNINGEXCEPTION_HPP_


#include <iostream>
#include <stdexcept>

class PaddingWarningException : public std::exception
{
public:

	/**
	 * @brief Constructor by default.
	 */
	explicit PaddingWarningException(){}

	/**
	 * @brief what() method with custom message
	 */
	virtual const char* what() const throw(){
		return "Warning: your columns number is not divisible by 4, is recommended to use padding functions!";
	}
};


#endif /* MATRIX_PADDINGWARNINGEXCEPTION_HPP_ */
