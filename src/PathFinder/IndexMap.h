/*
 * 	IndexMap.h
 *	This class maps currents assigning index to them based on the related nodes, functions
 *	works in general using hash tables. 	
 * 
 *  Created on: 7 de oct. de 2017
 *  Author: kevin
 */

#ifndef INDEXMAP_H_
#define INDEXMAP_H_

#include "NodePair.h"
#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>

using namespace boost::bimaps;

class IndexMap{
	typedef boost::bimap<unordered_set_of<NodePair, NodePair::Hash, NodePair::Equality>,set_of<int>> xVectorMap;
	typedef xVectorMap::value_type relation;
private:
	xVectorMap map;
	int matrix_cols;
	int matrix_rows;
public:

	IndexMap(int rows, int cols);
	void setMap(int rows, int cols);
	NodePair getNodesFromX(int xIndex);
	void getNodesFromX(int & node1row, int & node1col, int & node2row, int & node2col, int xIndex);
	int getXFromNodes(int  node1row, int  node1col, int  node2row, int  node2col);
	int getXFromNodes(NodePair nodes);
};

/**
 * @brief Constructor by default
 * @param rows Number of rows of the image (height)
 * @param cols Number of cols of the image (weidth)
 */
IndexMap::IndexMap(int rows, int cols){

	setMap(rows,cols);

}

/**
 * @brief This function maps all the index of currents between nodes
 * @param rows Number of rows of the image (height)
 * @param cols Number of cols of the image (weidth)
 */
void IndexMap::setMap(int rows, int cols){

	int indexNumber = 0;

	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){

			//Last row (only right node pairs)
			if(i == (rows-1)){

				//Right down edge have no node pair
				if(j != (cols-1)){
					NodePair pair1 =  NodePair(i,j,i,j+1); //Right node pair

					map.insert(relation(pair1,indexNumber));

					indexNumber = indexNumber + 1;
				}

			}
			//Other rows (Have right and down node pair)
			else{

				//last column (Have only down pair)
				if(j == (cols-1)){
					NodePair pair1 = NodePair(i,j,i+1,j); //Down node pair

					map.insert(relation(pair1,indexNumber));

					indexNumber = indexNumber + 1;

				}
				//Other columns
				else{
					NodePair pair1 = NodePair(i,j,i,j+1);//Rigth node pair
					NodePair pair2 = NodePair(i,j,i+1,j);//Down node pair

					map.insert(relation(pair1,indexNumber));
					map.insert(relation(pair2,indexNumber+1));

					indexNumber = indexNumber + 2;
				}

			}//End other rows


		}//End for j

	}//End for i
}

/**
 * @brief This function find the nodes related with an specific current
 * @param xIndex Index of the current
 * @return NodePair related with the current
 */
NodePair IndexMap::getNodesFromX(int xIndex){

	NodePair pointFind = map.right.find(xIndex)->second;
	return pointFind;
}

/**
 * @brief This function find the nodes related with an specific current
 * @param node1row passed by reference
 * @param node1col passed by reference
 * @param node2row passed by reference
 * @param node2col passed by reference
 * @param xIndex Index of the current
 */
void IndexMap::getNodesFromX(int & node1row, int & node1col, int & node2row, int & node2col, int xIndex){
	NodePair pointFind = map.right.find(xIndex)->second;
	node1row = pointFind.getFirstNodeRow();
	node1col = pointFind.getFirstNodeCol();
	node2row = pointFind.getSecondNodeRow();
	node2col = pointFind.getSecondNodeCol();
}

/**
 * @brief This function find current index between two nodes
 * @param node1row first row
 * @param node1col first col
 * @param node2row second row
 * @param node2col second col
 * @return Index of the current between nodes
 */
int IndexMap::getXFromNodes(int  node1row, int  node1col, int  node2row, int  node2col){

	NodePair pointFind =  NodePair(node1row,node1col,node2row,node2col);
	int indexFind = map.left.find(pointFind)->second;

	return indexFind;

}

/**
 * @brief This function find current index between two nodes
 * @param nodes NodePair with the two nodes
 * @return Index of the current between nodes
 */
int IndexMap::getXFromNodes(NodePair nodes){
	int indexFind = map.left.find(nodes)->second;
	return indexFind;
}


#endif /* INDEXMAP_H_ */
