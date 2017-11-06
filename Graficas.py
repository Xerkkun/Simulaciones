import os
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal,getcontext

nombre_archivo = input("Archivo de entrada	= ")
datos = np.loadtxt(nombre_archivo)
n = len(datos)
presion = []
densidad = []

for i in range(0,n):
    presion.append(datos[i][1])
    densidad.append(datos[i][2])

plt.plot(presion,densidad)
