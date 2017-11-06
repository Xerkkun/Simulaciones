import os
import numpy as np
import matplotlib.pyplot as plt
import commands
from decimal import Decimal,getcontext

n=225 #particulas
k=10000 #pasos
dx=10 #intervalo
l=0

#plt.xlabel()
#plt.ylabel()
for j in range(10,k,dx):
    l=l+1
    archivo = str(j)+".dat"
    datos = np.loadtxt(archivo)
    x = []
    y = []

    plt.title("ISTEP"+" "+str(j))

    for i in range(0,n):
        x.append(datos[i][0])
        y.append(datos[i][1])

    plt.axis([-8.5, 8.5, -8.5, 8.5]) #[xmin, xmax, ymin, ymax]
#    plt.axis('off')
    plt.plot(x,y, marker='o',linestyle="None", label=str(j),color="b")
    salida = str(l).zfill(3)+".png"
    plt.savefig(salida)
    plt.clf()

plt.close()

commands.getoutput('convert -delay 10 -loop 0 *.png animated.gif')
#commands.getoutput('ffmpeg -framerate 1/2 -i img%04d.png -c:v libx264 -r 30 out.mp4')
#commands.getoutput('rm *.png')
#commands.getoutput('rm *.dat')
