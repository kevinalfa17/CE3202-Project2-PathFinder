/*
 * NodePair.h
 *
 *  Created on: 7 de oct. de 2017
 *  Author: kevin
 * 	
 */

#ifndef NODEPAIR_H_
#define NODEPAIR_H_

#include <cassert>
#include <string>
#include <iostream>
#include <boost/unordered_set.hpp>

class NodePair
{
  private:
	int firstNodeRow;
	int firstNodeCol;
	int secondNodeRow;
	int secondNodeCol;

  public:
	void setNodes(int node1row, int node1col, int node2row, int node2col);
	void getNextNode(int actualRow, int actualCol, int &nextRow, int &nextCol);
	NodePair(int node1row, int node1col, int node2row, int node2col);

	void printPair();

	//Getters of the nodes
	int getFirstNodeCol() const
	{
		return firstNodeCol;
	}

	int getFirstNodeRow() const
	{
		return firstNodeRow;
	}

	int getSecondNodeCol() const
	{
		return secondNodeCol;
	}

	int getSecondNodeRow() const
	{
		return secondNodeRow;
	}

	//Struct required for be used in hash tables
	struct Hash
	{
		size_t operator()(const NodePair &a) const
		{
			size_t seed = 0;
			boost::hash_combine(seed, a.getFirstNodeRow());
			boost::hash_combine(seed, a.getFirstNodeCol());
			boost::hash_combine(seed, a.getSecondNodeRow());
			boost::hash_combine(seed, a.getSecondNodeCol());
			return seed;
		}
	};
	//Struct required to determine if two node pairs are equal
	struct Equality
	{
		bool operator()(const NodePair &p1, const NodePair &p2) const
		{
			return (p1.getFirstNodeRow() == p2.getFirstNodeRow() && p1.getFirstNodeCol() == p2.getFirstNodeCol() &&
					p1.getSecondNodeRow() == p2.getSecondNodeRow() && p1.getSecondNodeCol() == p2.getSecondNodeCol());
		}
	};
};

/**
 * @brief Constructor by default
 * @param node1row first row
 * @param node1col first col
 * @param node2row second row
 * @param node2col second col
 */
NodePair::NodePair(int node1row, int node1col, int node2row, int node2col)
{
	firstNodeRow = 0;
	firstNodeCol = 0;
	secondNodeRow = 0;
	secondNodeCol = 0;
	setNodes(node1row, node1col, node2row, node2col);
}

/**
 * @brief This function determine the first and second node of the class based in the rows and cols
 * @param node1row first row
 * @param node1col first col
 * @param node2row second row
 * @param node2col second col
 */
void NodePair::setNodes(int node1row, int node1col, int node2row, int node2col)
{
	//assert(node1row != node2row && node1col != node2col);
	if (node1row < node2row)
	{
		firstNodeRow = node1row;
		firstNodeCol = node1col;
		secondNodeRow = node2row;
		secondNodeCol = node2col;
	}
	else if (node1row == node2row)
	{
		if (node1col < node2col)
		{
			firstNodeRow = node1row;
			firstNodeCol = node1col;
			secondNodeRow = node2row;
			secondNodeCol = node2col;
		}
		else
		{
			firstNodeRow = node2row;
			firstNodeCol = node2col;
			secondNodeRow = node1row;
			secondNodeCol = node1col;
		}
	}
	else
	{
		firstNodeRow = node2row;
		firstNodeCol = node2col;
		secondNodeRow = node1row;
		secondNodeCol = node1col;
	}
}

/**
 * @brief This function print the pair in the format 00,11
 */
void NodePair::printPair()
{

	std::cout << firstNodeRow << firstNodeCol << "," << secondNodeRow << secondNodeCol << std::endl;
}

/**
 * @brief This function get one node and determines the pair
 * @param actualRow actual node row
 * @param actualCol actual node col
 * @param nextRow Passed by reference
 * @param nextCol Passed by reference
 */
void NodePair::getNextNode(int actualRow, int actualCol, int &nextRow, int &nextCol)
{
	//Check if the pair of the actual node is the first or the second
	if (actualRow == firstNodeRow && actualCol == firstNodeCol)
	{
		nextRow = secondNodeRow;
		nextCol = secondNodeCol;
	}
	else
	{
		nextRow = firstNodeRow;
		nextCol = firstNodeCol;
	}
}

#endif /* NODEPAIR_H_ */
