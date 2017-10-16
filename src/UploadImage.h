/*
 * UploadImage.h
 *
 *  Created on: 15 de oct. de 2017
 *      Author: daedgomez
 */

#ifndef UPLOADIMAGE_H_
#define UPLOADIMAGE_H_

#include <opencv2/opencv.hpp>

#include <iostream>
#include <string>
using namespace cv;
using namespace std;

class UploadImage {
public:
	UploadImage();
	void upload();
};


UploadImage::UploadImage() {
}
void UploadImage::upload() {
	   Mat image = imread( "Ruta.png", IMREAD_COLOR ); // Read the file
			       namedWindow( "Display window", WINDOW_AUTOSIZE ); // Create a window for display.
			       imshow( "Display window", image );                // Show our image inside it.
			       waitKey(0); // Wait for a keystroke in the window

}



#endif /* UPLOADIMAGE_H_ */
