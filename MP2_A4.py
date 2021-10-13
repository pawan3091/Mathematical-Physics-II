import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def euler_forward(f,t0,x0,b,h):
    n= int((b-t0)/h)
    t_arr=[t0]
    X_arr=[x0]
    for i in range(n):
        xn = x0 + h * f(t0,x0)
        x0 = xn
        t0 = t0+h
        X_arr.append(x0)
        t_arr.append(t0)
    return t_arr,X_arr

def RK2(f,t0,x0,b,h):
    n= int((b-t0)/h)
    t_arr=[t0];X_arr=[x0]
    for i in range(n):
        k1=h*f(t0,x0)
        k2=h*f(t0+h,x0+k1)
        k_RK2=(k1+k2)/2
        xn=x0+k_RK2
        x0 = xn
        t0 = t0+h
        X_arr.append(x0)
        t_arr.append(t0)
    return t_arr,X_arr
    
def graph(x1_1,y1_1,x1_2,y1_2,x1_3,y1_3,x1_lab,y1_lab,tit_1,x2_1,y2_1,x2_2,y2_2,x2_3,y2_3,x2_lab,y2_lab,tit_2,x3_1,y3_1,x3_2,y3_2,x3_3,y3_3,x3_lab,y3_lab,tit_3,x4,y4_1,y4_2,x4_lab,y4_lab,tit_4,x_an,y_an):
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(x1_1, y1_1,marker="*",label="EULER") 
    axs[0, 0].plot(x1_2, y1_2,marker="*",label="RK2")
    axs[0, 0].plot(x1_3, y1_3,marker="o",label="ANALYTIC")
    axs[0,0].set_xlabel(x1_lab)
    axs[0,0].set_ylabel(y1_lab)
    axs[0, 0].set_title(tit_1,c="red")
    axs[0, 1].plot(x_an, y_an,marker="|",label="ANALYTIC")
    axs[0, 1].plot(x2_1, y2_1,marker="*", label="h ={}".format(x2_1[1]-x2_1[0]))
    axs[0, 1].plot(x2_2, y2_2,marker="*", label="h ={}".format(x2_2[1]-x2_2[0]))
    axs[0, 1].plot(x2_3, y2_3,marker="*", label="h ={}".format(x2_3[1]-x2_3[0]))
    axs[0,1].set_xlabel(x2_lab)
    axs[0,1].set_ylabel(y2_lab)
    axs[0, 1].set_title(tit_2,c="red")
    axs[1,0].plot(x_an, y_an,marker="|",label="ANALYTIC")
    axs[1,0].plot(x3_1, y3_1,marker="*", label="h ={}".format(x3_1[1]-x3_1[0]))
    axs[1,0].plot(x3_2, y3_2,marker="*", label="h ={}".format(x3_2[1]-x3_2[0]))
    axs[1,0].plot(x3_3, y3_3,marker="*", label="h ={}".format(x3_3[1]-x3_3[0]))
    axs[1,0].set_xlabel(x3_lab)
    axs[1,0].set_ylabel(y3_lab)
    axs[1,0].set_title(tit_3,c="red")
    axs[1, 1].scatter(x4,y4_1 , label="EULER")
    axs[1, 1].scatter(x4,y4_2 , label="RK2")
    axs[1,1].set_xlabel(x4_lab)
    axs[1,1].set_ylabel(y4_lab)
    axs[1,1].set_title(tit_4,c="red")
    axs[0,0].grid()
    axs[1,0].grid()
    axs[0,1].grid()
    axs[1,1].grid()
    axs[0,0].legend()
    axs[1,0].legend()
    axs[0,1].legend()
    axs[1,1].legend()
    plt.show()

def Anal(x0,t0,b,h,tau):
    t_arr=[]
    X_arr=[]
    e=np.arange(t0,b+h,h)
    for t in e:
        D=x0*np.exp(-1*t/tau)
        X_arr.append(D)
        t_arr.append(t)
    return t_arr,X_arr

def q3a(t0,x0,t_h):
    b=5*t_h
    tou=t_h/0.693
    f= lambda x,N : -1*N/tou
    a=euler_forward(f,t0,x0,b,t_h/8)
    a1=euler_forward(f,t0,x0,b,t_h/2)
    a2=euler_forward(f,t0,x0,b,t_h/4)
    b1=RK2(f,t0,x0,b,t_h/8)
    b2=RK2(f,t0,x0,b,t_h/2)
    b3=RK2(f,t0,x0,b,t_h/4)
    c=Anal(x0,t0,b,t_h/8,tou)
    eu=[]
    rk=[]
    for i in range(len(a[1])):  
        eu.append((c[1][i]-a[1][i])/c[1][i])
        rk.append((c[1][i]-b1[1][i])/c[1][i])
    data = {"t":a[0],"N (Euler)":a[1],"N (RK2)":b1[1],"N (Analytic)":c[1],"Ab. Error in Euler":eu,"Ab. Error in RK2":rk}
    print(pd.DataFrame(data))
    graph(np.array(a[0])/t_h,a[1],np.array(b1[0])/t_h,b1[1],np.array(c[0])/t_h,c[1],"No. of half lives","No. of particles","Euler and RK2",a[0],a[1],a1[0],a1[1],a2[0],a2[1],"t","N","Euler method for various h",b1[0],b1[1],b2[0],b2[1],b3[0],b3[1],"t","N","RK2 method for various h",a[0],eu,rk,"Time(t)","Absolute Error","Error plot",c[0],c[1])

def q3b(t0,x0,R,C):
    b=5*R*C
    tou=R*C
    f= lambda x,V : -1*V/tou
    a=euler_forward(f,t0,x0,b,tou/10)
    a1=euler_forward(f,t0,x0,b,tou/5)
    a2=euler_forward(f,t0,x0,b,tou/15)
    b1=RK2(f,t0,x0,b,tou/10)
    b2=RK2(f,t0,x0,b,tou/5)
    b3=RK2(f,t0,x0,b,tou/15)
    c=Anal(x0,t0,b,tou/10,tou)
    eu=[]
    rk=[]
    for i in range(len(a[1])):  
        eu.append((c[1][i]-a[1][i])/c[1][i])
        rk.append((c[1][i]-b1[1][i])/c[1][i])
    data = {"t":a[0],"V (Euler)":a[1],"V (RK2)":b1[1],"N (Analytic)":c[1],"Ab. Error in Euler":eu,"Ab. Error in RK2":rk}
    print(pd.DataFrame(data))
    graph(a[0],np.array(a[1])/a[1][0],b1[0],np.array(b1[1])/b1[1][0],c[0],np.array(c[1])/c[1][0],"Time (t)","V/V0","Euler and RK2",a[0],a[1],a1[0],a1[1],a2[0],a2[1],"t","V","Euler method for various h",b1[0],b1[1],b2[0],b2[1],b3[0],b3[1],"t","V","RK2 method for various h",a[0],eu,rk,"Time(t)","Absolute Error","Error plot",c[0],c[1])
    
def q3c(t0,x0,eta,rad,m):
    tou=m/(np.pi*6*rad*eta)
    b=4*tou
    f= lambda x,v : -1*v/tou
    a=euler_forward(f,t0,x0,b,tou/10)
    a1=euler_forward(f,t0,x0,b,tou/5)
    a2=euler_forward(f,t0,x0,b,tou/15)
    b1=RK2(f,t0,x0,b,tou/10)
    b2=RK2(f,t0,x0,b,tou/5)
    b3=RK2(f,t0,x0,b,tou/15)
    c=Anal(x0,t0,b,tou/10,tou)
    eu=[]
    rk=[]
    for i in range(len(a[1])):  
        eu.append((c[1][i]-a[1][i])/c[1][i])
        rk.append((c[1][i]-b1[1][i])/c[1][i])
    data1 = {"t":a[0],"v (Euler)":a[1],"v (RK2)":b1[1],"v (Analytic)":c[1],"Ab. Error in Euler":eu,"Ab. Error in RK2":rk}
    print(pd.DataFrame(data1))
    graph(a[0],a[1],b1[0],b1[1],c[0],np.array(c[1]),"Time (t)","Terminal Velocity","Euler and RK2",a[0],a[1],a1[0],a1[1],a2[0],a2[1],"t","v","Euler method for various h",b1[0],b1[1],b2[0],b2[1],b3[0],b3[1],"t","v","RK2 method for various h",a[0],eu,rk,"Time(t)","Absolute Error","Error plot",c[0],c[1])
    
t0=0;x0=2e4;t_h=4     
print("----------------------------Q. 3a For h= ",t_h/8,"-----------------------------")
q3a(t0,x0,t_h)

t0=0;x0=10;R=1e3;C=1e-6
print("----------------------------Q. 3b For h= ",R*C/10,"-----------------------------")
q3b(t0,x0,R,C)


t0=0;x0=10;eta=10;rad=0.2;m=200
print("----------------------------Q. 3c For h= ",m/(np.pi*6*rad*eta)/10,"-----------------------------")
q3c(t0,x0,eta,rad,m)
