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
	Mat getImage();
	Mat image;
};


UploadImage::UploadImage() {
}
void UploadImage::upload() {

	   this->image = imread( "Ruta10px.png", 0 ); // Read the file

			      /**
			       namedWindow( "Display window", WINDOW_AUTOSIZE ); // Create a window for display.
			       cout << image<< endl;
			       cout << image.rows << endl;
			       cout << image.cols << endl;
			       imshow( "Display window", image );                // Show our image inside it.
			       waitKey(0); // Wait for a keystroke in the window
			      */

}

Mat UploadImage::getImage(){
	return this->image;
}


#endif /* UPLOADIMAGE_H_ */
