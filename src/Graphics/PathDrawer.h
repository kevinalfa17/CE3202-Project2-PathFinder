/*
 * NodePair.h
 *
 *  Created on: 18 de oct. de 2017
 *      Author: kevin
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

PathDrawer::PathDrawer(Mat mat){

    cvtColor( mat, this->imageMatrix, CV_GRAY2RGB);
}

void PathDrawer::drawPath(const vector<Point> * points){
    Vec3b redLine(0,0,255);
    for (int i=0; i < points->size(); i++){
        this->imageMatrix.at<Vec3b>(points->at(i).x,points->at(i).y) = redLine;        
    }

    namedWindow( "PathFinder Strategy 1", WINDOW_AUTOSIZE );
    imshow( "PathFinder Strategy 1", this->imageMatrix );                   

    waitKey(0); 
}


#endif /* PATHDRAWER_H_ */