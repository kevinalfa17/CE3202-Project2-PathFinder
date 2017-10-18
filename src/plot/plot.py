import matplotlib.pyplot as plt

def plot2d(datax,datay,xd,xu,yd,yu,title,xlabel,ylabel):
	line = plt.plot(datax,datay,label='hola')
	plt.title(title)
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.ylim(yd,yu)
	plt.xlim(xd,xu)
	plt.legend()
	plt.show()

plot2d([1,2,4,5],[5,8,12,16],0,20,0,20,"example",'etiquetax','etiquetay')
