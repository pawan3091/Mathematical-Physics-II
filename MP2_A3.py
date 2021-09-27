from sympy import *
from sympy import simplify
from scipy.interpolate import lagrange
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#2a
def lang(x_i,y_i):
    x= symbols('x')
    if len(x_i)==len(y_i):
       yp=0
       for i in range(len(x_i)):
           p=1
           for j in range(len(x_i)):
               if j!=i:
                  p*=(x-x_i[j])/(x_i[i]-x_i[j])
           yp+=y_i[i]*p
       p_x=simplify(yp)
       fx=lambdify(x,p_x,modules=['numpy'])
       return fx
    else:
        return("len(x_i) should be equal to len(y_i)")

#2b
def in_lang(x_i,y_i):
    return lang(y_i,x_i)


#2c
def inbuilt(x_i,y_i,x1):
    poly = lagrange(x_i, y_i) 
    L=poly(x1)
    return L
    
#Graph

def graph1(x_i,y_i,arr_x,arr_y_f,arr_y_in,pt_x,pt_y,xlab,ylab,titl):
    plt.scatter(x_i,y_i,marker='*',c="red",label="Given discrete points")
    plt.scatter(pt_x,pt_y,marker='o',c='black',label="Interpolated point")
    plt.plot(arr_x,arr_y_in,c='blue',label="Scipy's Inbuilt function",linestyle="-.")
    plt.plot(arr_x,arr_y_f,c="green",label="Interpolated langrange function")
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(titl)
    plt.grid(True)
    plt.legend()
    plt.show()
    
#3a 
a=[0.00,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0]
b=[1,0.99,0.96,0.91,0.85,0.76,0.67,0.57,0.46,0.34,0.22,0.11,0.00,-0.1,-0.18,-0.26]
g=lang(a,b)
k=in_lang(a,b)

arr_beta=np.linspace(0,3,1000)
arr_Jbeta_f=[]
arr_Jbeta_i=[]
for i in arr_beta:
    arr_Jbeta_f.append(g(i))
    arr_Jbeta_i.append(inbuilt(a,b,i))

arr_Jbeta=np.linspace(1,-0.26,1000)
arr_beta_f=[]
arr_beta_i=[]
for i in arr_Jbeta:
    arr_beta_f.append(k(i))
    arr_beta_i.append(inbuilt(b,a,i))
    
print("The value of bessel function for \u03B2 = 0.5 is ",g(2.3))
print("The value of \u03B2 for which the value of bessel function is 2.3 = ",k(0.5))   
graph1(a,b,arr_beta,arr_Jbeta_f,arr_Jbeta_i,2.3,g(2.3),"\u03B2","J0_\u03B2","3a. (i) Bessel Function") 
graph1(b,a,arr_Jbeta,arr_beta_f,arr_beta_i,0.5,k(0.5),"J0_\u03B2","\u03B2","3a. (ii) Inverse Bessel Function")
    
#3b
I=[2.81,3.24,3.80,4.30,4.37,5.29,6.03]
V=[0.5,1.2,2.1,2.9,3.6,4.5,5.7]
s=in_lang(I,V)
z=lang(I,V)

arr_I=np.linspace(2.81,6.03,1000)
arr_V_f=[]
arr_V_i=[]
for i in arr_I:
    arr_V_f.append(z(i))
    arr_V_i.append(inbuilt(I,V,i))

arr_V=np.linspace(0.5,5.7,1000)
arr_I_f=[]
arr_I_i=[]
for i in arr_V:
    arr_I_f.append(s(i))
    arr_I_i.append(inbuilt(V,I,i))

graph1(I,V,arr_I,arr_V_f,arr_V_i,3.79,z(3.79),"I","V","3b. (i) Photoelectric Effect") 
graph1(V,I,arr_V,arr_I_f,arr_I_i,2.4,s(2.4),"V","I","3b. (ii) Inverse Photoelectric effect")
print("The value of I for V= 2.4 is ",s(2.4))

#Comparison
j=["3a (i)","3a (ii)","3b"]
d=[g(2.3),k(0.5),s(2.4)]
c=[inbuilt(a,b,2.3),inbuilt(b,a,0.5),inbuilt(V,I,2.4)]
error=np.array(d)-np.array(c)
print("# Comparison Table")
Data={"Ques":j,"Scipy":d,"My function":c,"Error":error}
print(pd.DataFrame(Data))
