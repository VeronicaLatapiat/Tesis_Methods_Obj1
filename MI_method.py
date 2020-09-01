#!/usr/bin/env python

# Metodo MI para cuantificar la fuerza de relacion entre dos variables continuas aleatorias

import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
#from scipy.special import gamma, psi
#from numpy import pi
from itertools import combinations


import matplotlib.pyplot as plt
import seaborn

# Leer matriz, GenExpressionMatrix
datos=pd.read_csv('topology_expression_data.txt', sep='\t', engine='python', index_col = 0)
datos

# Ingresar valor de k

#numero = input('Introduce el valor de k: ')
#print ('El valor de k que ingreso es:\n',numero)
k=5

samples= len(datos.columns) # numero de samples
print("Número de samples:",samples)

N = len(datos.index) # numero de genes, nodos
print("Número de genes:",N)

seaborn.distplot(datos.loc['A', :], bins = 50)
# This function combines the matplotlib hist function (with automatic calculation of a good default bin size)
# with the seaborn kdeplot() and rugplot() functions. It can also fit scipy.stats distributions and plot the
# estimated PDF over the data.

# Combinatoria para ver el total de veces que se realizara el calculo de MI
# dado por el nro de genes

#####
#lista(range(0, N))
# Todo eso se puede reemplazar por
combinations(res, 2) # esto es un generador

#obtener todas las combinaciones de numeros de largo 2

# tomo el primer gen(todas las muestras) y lo comparo con el segundo, con el tercero hasta n
# tomo el segundo gen(todas las muestras) y lo comparo con el tercero, el cuarto...
# la combinatoria es N(N-1)/2 veces
for a, b in combinations(res, 2):
    lista.append([a, b]) # fancy
len(lista)
print(lista)

#for d in lista:
#    print(datos.iloc[d[0],:]) #columna, todas las filas
#    print(datos.iloc[d[1],:])
    
# se puede comprimir
for genA, genB in combinations(res, 2):
    print(datos.iloc[[genA, genB], :]) # fila idxA, todas las columnas
    
    ## Primer paso para el calculo de la Informacion Mutua

from sklearn.neighbors import KernelDensity

# I(X:Y)=H(X)+H(Y)-H(X,Y) -> X e Y son diferentes genes, para eso la combinatoria
# Aqui solo calculamos H(X)

k=5

# nodos a comparar (para entender que hará el for sobre las combinaciones)
genA = 0
genB = 1

#Aqui solo estoy calculando para un gen, H(X) ya despues hare la combinatoria con una funcion...

x_prueba=datos.iloc[0:1,1:samples+1] # primera fila con info de 1 gen en n samples, menos el nombre del nodo
#x_prueba_2=datos.iloc[[genA, genB], :] # primera fila con info de 1 gen en n samples, menos el nombre del nodo
#x_prueba=datos.iloc[genA,:]
x_prueba=datos.iloc[genA]
#nombre_nodo=datos.iloc[0,0]
nombre_nodo = datos.index[genA] # por qué solo uno? la entropia se calcula en un gen y despues en otro y despues se calcula la MI
print("x_prueba",x_prueba) 
#print("x_prueba2",x_prueba_2)
print("nombre nodo",nombre_nodo)

# tenemos N muestras y calcularemos las k muestras mas cercanas en los valores de x(expresion de un gen a esta)

#antiguo
#datos_num = x_prueba.to_numpy()
#datos_num_transpose = datos_num.transpose()

#print(datos_num_transpose)

#print("x",x_prueba)
datos_num_transpose = x_prueba.T.values # hace lo mismo


print("t2",datos_num_transpose)

tama=len(datos_num_transpose)
print("tamaño",tama)

print(datos_num_transpose.shape)

datos_num_transpose_2 = datos_num_transpose.reshape(-1, 1)
print(datos_num_transpose_2)       # Produces: [[ 0.58  0.76]]
print(datos_num_transpose_2.shape) # Produces (1, 2) meaning 1 row, 2 cols

#print("Sample",samples)
muestra=datos.iloc[0,0] # primera fila con info de 1 gen en n samples, menos el nombre del nodo
print("Sample",muestra)
muestra2=datos.iloc[0,99] # primera fila con info de 1 gen en n samples, menos el nombre del nodo
print("Sample",muestra2)
#muestra3=datos.iloc[0,100] # primera fila con info de 1 gen en n samples, menos el nombre del nodo
#print("Sample",muestra3)

# Centrar datos

#z= (1/k+1)sumatoria de j hasta k de xi^j

listaY = []
#listaY2 = []

#print(samples)
   
for j in range(0,samples):#[-1:]: #deberian ser 100 samples, aca van desde el 0 al 99, estamos viendo solo un gen, con k vecinos mas cercanos de x en cada uno
    jj = np.int(j) #otro genes se veran con for a la combinatoria
    print("Sample",jj+1)
    print(k,"vecinos más cercanos centrados a x dentro de la muestra",jj+1)
    #print(xi[jj,0:k+1])# k-vecinos mas cercanos en diferentes samples...
    #sumar datos
    suma = 0
    
    for l in range(0,k+1):
        suma += xi[jj,l]
        #print("suma...... empieza de 1?",xi[jj,l])
        
        if(l==k):
        #if(l==k):
            print ("suma:",suma)
            #print ("z:",suma*(1/(k+1)))
            #z= (1/k+1)sumatoria de j hasta k de xi^j  
            z= ((k+1)**-1)*suma
            print("z:",z)
        else:
            continue
        #y^j= Xi^j-z
        
        #Aqui esta mal... se debe sacar xi desde 
        
        #print("xi[jj,l]")
        #print(xi[jj,l])
        #print("xi[jj,0:k+1]")
        print(xi[jj,0:k+1])# k-vecinos mas cercanos en j -> diferentes samples...
        
        for ki in range(0,k+1):
            
            y= xi[jj,ki] - z
                #y= xi[jj,l] - z
                #y= xi[jj,l] - z
            print("ki")
            print(xi[jj,ki])
            print ("Y?",y)
            yr=round(y,3) 
            print ("y-for*****",yr)# guardar y como matriz Y
            #LISTA
            #loop para agregar sucesores a lista hasta llegar a r2   !!?buscar forma mas eficiente   
            
            listaY.append(yr) 
            #print("listaY")
            #print("listaY",listaY)
        
print("listaY-fuera for2",listaY) #vector de expresión de gen centrado en las k-muestras cercanas, por las 100 muestras

print(listaY)
vector=np.array(listaY)
tam=len(vector)
print("tamaño",tam)
print(vector)
matrix = vector.reshape(samples,k+1) 
print(matrix)

Ymatrix = matrix[1,:] # en el caso del paper Ymatrix es de dimensión k+1*d
print(Ymatrix)

Ymatrix_tt = Ymatrix.reshape(1,-1)
#print(Ymatrix_tt)
Ymatrix_trans=Ymatrix_tt.transpose()
print(Ymatrix_trans)

print("----SVD----")

#Ymatrix_trans(k+1*d)
u, s, vt = np.linalg.svd(Ymatrix_trans)

print("Left Singular Vectors:")
print(u) # (k+1*k+1)

print("Singular Values:") 
print(np.diag(s)) # (k+1*d)

print("Right Singular Vectors:") 
print(vt) # (d*d)

#Since V is unitary, the singular vectors, vi(l)
#are of unit
#length and orthogonal. The first singular vector, vi(1)
#points in the direction in which the data is stretched the
#most, and each subsequent singular vector points in the
#direction which is orthogonal to all previous singular vectors and that accounts for the most stretching.
