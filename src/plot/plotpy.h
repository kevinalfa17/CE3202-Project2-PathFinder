/**

ANPI
Metodos para realizar graficas x,y usando Matplotlib (Python)
@Autor: David Badilla S.
@Description: Clase que se utiliza para poder realizar graficas de puntos x y
	      usando como base la biblioteca Matplotlib de Python. La comuni-
	      se establece mediante Python.h, lo cual permite ejecutar codigo 
	      de Python desde c++.
@Contact: davbs94@gmail.com


*/

#ifndef PLOTPY_H
#define PLOTPY_H

#include <python2.7/Python.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "../Matrix/Matrix.hpp"



/**
El paquete plotpy contiene los metodos necesarios para realizar
graficas xy.
*/

namespace plotpy{



template<typename T>
class  Plot2d{

	private:
		//Titulo de la grafica.
		std::string _title;
		//Estiqueta del eje x
		std::string _xlabel;
		//Etiqueta del eje y
		std::string _ylabel;
		//Tamano de la cuadricula
		int _sizegrid;


	public:
		
		Plot2d();

		~Plot2d();

		void initialize(int id);

		void settitle(std::string title);

		void setxlabel(std::string label);

		void setylabel(std::string label);

		void setgridsize(int sizegrid);

		void setxrange(int xi,int xs);

		void setyrange(int yi,int ys);

		void plot(std::vector<T> datax, std::vector<T> datay,std::string label);

        void setVecRoute(anpi::Matrix<T> &X, anpi::Matrix<T> &Y, std::vector<T> &xv, std::vector<T> &yv, int id);

		void showallplots();


}; //class Plot2d

/*
*@brief Constructor
*/
template <typename T>
Plot2d<T>::Plot2d(){}

/*
*@brief Destructor
*/
template <typename T>
Plot2d<T>::~Plot2d(){}

/*
*@brief Inicializa la figura donde colocaremos la grafica(s). 
*Se le debe asigna un ID UNICO.
*@param id ID de la figura que contendra las grafica(s).
*/
template <typename T>
void Plot2d<T>::initialize(int id){
    Py_Initialize();
    PyRun_SimpleString("import matplotlib.pyplot as plt");
    PyRun_SimpleString("from matplotlib.path import Path");
    PyRun_SimpleString("import matplotlib.patches as patches");
    std::string tmp1 = "plt.figure("+std::to_string(id)+")"; 
    PyRun_SimpleString(tmp1.c_str());
    _title = "";
    _xlabel = "";
    _ylabel = "";
    _sizegrid = 0;
}

/*
*@brief Se asigna un titulo a la grafica.
*@param title Nombre de la grafica.
*/
template <typename T>
void Plot2d<T>::settitle(std::string title){
    _title = title;
    std::string strtitle = "plt.title('"+ _title +"')";
    PyRun_SimpleString(strtitle.c_str());
}

/*
*@brief Se asigna una etiqueta de nombre para el eje x.
*@param xlabel Nombre de la etiqueta.
*/
template <typename T>
void Plot2d<T>::setxlabel(std::string xlabel){
    _xlabel = xlabel;
    std::string strxlabel = "plt.xlabel('"+ _ylabel +"')";
    PyRun_SimpleString(strxlabel.c_str());
}

/*
*@brief Se asigna una etiqueta de nombre para el eje y.
*@param ylabel Nombre de la etiqueta.
*/
template <typename T>
void Plot2d<T>::setylabel(std::string ylabel){
    _ylabel = ylabel;
    std::string strylabel = "plt.ylabel('"+ _ylabel +"')";
    PyRun_SimpleString(strylabel.c_str());
}

/*
*@brief Se escoge el tamano del cuadriculado de la grafica
*Por default no hay cuadriculado
*@param sizegride Tamano de la cuadricula
*/
template <typename T>
void Plot2d<T>::setgridsize(int sizegrid){
    _sizegrid = sizegrid;
    std::string strgrid = "plt.grid("+std::to_string(_sizegrid) +")";
    PyRun_SimpleString(strgrid.c_str());
}

/*
*@brief Escoge el rango de valores que tomara el eje x.
*Si no se escoje rango, se ajusta un default.
*@param xi Minimo valor para el eje x
*@param xs Maximo valor para eje x.
*/
template <typename T>
void Plot2d<T>::setxrange(int xi, int xs){
    std::string strxlim = "plt.xlim("+  std::to_string(xi) + "," + std::to_string(xs) + ")";
    PyRun_SimpleString(strxlim.c_str());
}

/*
*@brief Escoge el rango de valores que tomara el eje y.
*@param yi Minimo valor para el eje y
*@param ys Maximo valor para eje y.
*/
template <typename T>
void Plot2d<T>::setyrange(int yi, int ys){
    std::string strylim = "plt.ylim("+  std::to_string(yi) + "," + std::to_string(ys) + ")";
    PyRun_SimpleString(strylim.c_str());
}


/*
*@brief Se encarga de hacer una grafica. Se pueden llamar varias veces, una seguida de otra
*para hacer varias graficas en una misma figura
*@param datax Vector de valores en el eje x.
*@param datay Vector de valores en el eje y.
*@param label Etiqueta para simbologia de la grafica. Se le puede enviar "" y no agrega la 
*simbologia, de momento lanza un warning, cuestion que se corregira luego.
*/
template <typename T>
void Plot2d<T>::plot(std::vector<T> datax, std::vector<T> datay, std::string label){
    std::string tmp1 = "datax = [";
    std::string tmp2 = "datay = [";
    std::string tmp3 = "plt.plot(datax,datay,label='" +label + "')";
    for(unsigned int i = 0; i < datax.size(); i++){
                if(i == datax.size() - 1){
                    tmp1.append(std::to_string(datax[i]) + "]");
            		tmp2.append(std::to_string(datay[i]) + "]");
            }
                else{
                    tmp1.append(std::to_string(datax[i]) + ", ");
            		tmp2.append(std::to_string(datay[i]) + ",");
        }
    }

    const char * a = tmp1.c_str();
    const char * b = tmp2.c_str();
    PyRun_SimpleString(a);
    PyRun_SimpleString(b);
    PyRun_SimpleString(tmp3.c_str());
    PyRun_SimpleString("plt.legend()");
}

template<typename T>
void Plot2d<T>::setVecRoute(anpi::Matrix<T> &X, anpi::Matrix<T> &Y, std::vector<T> &xv, std::vector<T> &yv, int id) {
	std::string tmp1 = "datax = [";
	std::string tmp2 = "datay = [";
    std::string tmp3 = "plt.quiver(datax,datay,units='width')";

	for (int i = 0; i < X.rows(); i++)
	{
		tmp1.append("[");
		tmp2.append("[");
		for(int j = 0; j < X.cols(); j++){
			if(j == X.cols()-1){
                if(i == X.rows()-1){
                    tmp1.append(std::to_string(X[i][j]) + "]");
                    tmp2.append(std::to_string(Y[i][j]) + "]");
                }else{
				    tmp1.append(std::to_string(X[i][j]) + "],");
				    tmp2.append(std::to_string(Y[i][j]) + "],");
                }
			}else{
				tmp1.append(std::to_string(X[i][j]) + ",");
				tmp2.append(std::to_string(Y[i][j]) + ",");
			}
		}
	}
    tmp1.append("]");
    tmp2.append("]");

    std::string tmp4 = "verts = [";
    std::string tmp5 = "path = Path(verts)";

    for(int i = 0; i < xv.size(); i++){
    	if(i == xv.size()-1)
    		tmp4.append("("+std::to_string(xv.at(i))+","+std::to_string(yv.at(i))+")");
    	else
    		tmp4.append("("+std::to_string(xv.at(i))+","+std::to_string(yv.at(i))+"),");
    }
    tmp4.append("]");

    std::string tmp6 = "ax = plt.figure("+std::to_string(id)+").add_subplot(111)";
    const char * a = tmp1.c_str();
    const char * b = tmp2.c_str();
    const char * c = tmp4.c_str();
    PyRun_SimpleString(a);
    PyRun_SimpleString(b);
    PyRun_SimpleString(c);
    PyRun_SimpleString(tmp3.c_str());
    PyRun_SimpleString(tmp5.c_str());
    PyRun_SimpleString(tmp6.c_str());
    PyRun_SimpleString("patch = patches.PathPatch(path, edgecolor='orange' , facecolor='none', lw=5)");
    PyRun_SimpleString("ax.add_patch(patch)");
    
}

/*
*@brief Metodo que despliega todas las figuras que se han creado antes de hacer este llamado
*
*
*
*/
template <typename T>
void Plot2d<T>::showallplots(){

    PyRun_SimpleString("plt.show()");
}




}

#endif // PLOTPY_H
