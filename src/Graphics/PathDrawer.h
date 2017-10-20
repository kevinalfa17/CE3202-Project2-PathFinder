/*
 *  PathDrawer.h
 *  This class is used to draw pixels on screen
 * 
 *  Created on: 18 de oct. de 2017
 *  Author: Kevin
 */

#ifndef PATHDRAWER_H_
#define PATHDRAWER_H_

#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>

class PathDrawer{
public:
    PathDrawer(Mat mat);
    void drawPath(const vector<Point> * points);

private:
    Mat imageMatrix;
};


/**
 * @brief Constructor by default
 * @param mat Image to draw path
 */
PathDrawer::PathDrawer(Mat mat){

    cvtColor( mat, this->imageMatrix, CV_GRAY2RGB);
}

/**
 * @brief Method that draws a previously calculated path over an image
 * @param points vector of points to draw
 */
void PathDrawer::drawPath(const vector<Point> * points){
    Vec3b redLine(0,0,255);
    for (int i=0; i < points->size(); i++){
        this->imageMatrix.at<Vec3b>(points->at(i).x,points->at(i).y) = redLine;        
    }

    namedWindow( "PathFinder Strategy 1", WINDOW_NORMAL );
    imshow( "PathFinder Strategy 1", this->imageMatrix );                   

    waitKey(0); 
}


#endif /* PATHDRAWER_H_ */