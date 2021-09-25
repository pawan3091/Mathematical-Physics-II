import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from numpy import trapz
from scipy.integrate import quad
x1= np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])  
y1=np.array([0.0,0.5,2.0,4.05,8.0,12.5,18.0,24.5,32.0,40.5,50.0])  
h=0.1 ; a=x1[0] ; b=x1[9]
n=int((b-a)/h) ; S=0.5*(y1[0]+y1[10]) ; S1 = y1[0]+y1[10] 
for i in range(1,10):
       if i%2 == 0:
           S1 = S1 + 2 * y1[i]
       else:
           S1 = S1 + 4 * y1[i]
       S+= y1[i]
Integral = S1 * h/3 ; Power = S * h
print("Power using trapezoidal and simpson",("%f"%Power,"%f"%Integral))
def trap(f,a,b,n):
    h=(b-a)/n
    S=0.5*(f(a)+f(b))
    for i in range(1,n):
        S+= f(a+i*h)
    Integral = S * h
    return Integral
def simpson(f,a,b,n):
    h=(b-a)/n
    S = f(a) + f(b)
    for i in range(1,n):       
        if i%2 == 0:
            S = S + 2 * f(a + i*h)
        else:
            S = S + 4 * f(a + i*h)
    Integral = S * h/3
    return Integral
     
def integration(f,a,b,n):
    print("a =",a,"; b =" ,b, " ; Number of intervals =",n)
    h=(b-a)/n
    x = np.linspace(a,b,n+1)
    y = f(x)
    inte , err = quad(f, a, b)   #absolute value taken for comparing my results 
    Integral_trap= trap(f,a,b,n)
    Integral_simp=simpson(f,a,b,n)
    print("Integral using scipy function (Quad) = ",inte) 
    print("Integral (using composite trapezoidal formula) = %f" %Integral_trap)
    print("Error in composite trapezoidal (by subtracting from the one i got from inbuilt quad)= ",abs(Integral_trap-inte))   
    print("Integral (using composite simpson formula) = %f" %Integral_simp)
    print("Error in composite simpson (by subtracting from the one i got from inbuilt quad) = ",abs(Integral_simp-inte))
    y_data_2,y_data,y_data_3=[],[],[]
    geo= np.array([10**i for i in range(4) ])
    n_array=np.arange(30,400,2)
    h_array=(b-a)/n_array
    for n in n_array:
        x_data=np.linspace(a,b,n+1)
        y_data.append(trap(f,a,b,n))
        y_data_2.append(simpson(f,a,b,n))
    d=[inte]*len(h_array)
    y_data_3=np.array(d)
    plt.plot(h_array,y_data,label="Trapezoidal rule")
    plt.scatter(h_array,y_data,label="Trapezoidal rule")
    plt.plot(h_array,y_data_2,label="Simpson's rule")
    plt.scatter(h_array,y_data_2,label="Trapezoidal rule")
    plt.plot(h_array,y_data_3,linestyle='--',label="scipy's simpson implementation")
    plt.xlabel("h")
    plt.ylabel("I(h)")
    plt.title("I(h) vs h plot [Convergence Test]")
    plt.legend()
    plt.grid()
    plt.xscale('log')
    plt.show()
f = lambda x : x*x
integration(f,0,10,100)
