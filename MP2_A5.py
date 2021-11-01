import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def Eul(func, X0 ,tmin,tmax,h):
    N=int((tmax-tmin)/h)
    t = np.linspace(tmin,tmax,N)
    X  = np.zeros([N, len(X0)])
    X[0] = X0
    for i in range(N-1):
        X[i+1] = X[i] + func(t[i],X[i]) * h
    return X,t

def RK2(func, X0,tmin,tmax,h):
    N=int((tmax-tmin)/h)
    t = np.linspace(tmin,tmax,N)
    X  = np.zeros([N, len(X0)])
    X[0] = X0
    for i in range(N-1):
        k1 =h* func(t[i],X[i])
        k2 = h*func( t[i] + h,X[i] +  k1)
        X[i+1] = X[i] +  (k1 +k2 )/2
    return X,t
    
def RK4(func, X0,tmin,tmax,h):
    N=int((tmax-tmin)/h)
    t = np.linspace(tmin,tmax,N)
    X  = np.zeros([N, len(X0)])
    X[0] = X0
    for i in range(N-1):
        k1 = func(t[i],X[i])
        k2 = func( t[i] + h/2,X[i] + (h* k1)/2)
        k3 = func( t[i] + h/2,X[i] + h/2* k2)
        k4 = func(t[i] + h,X[i] + h   * k3)
        X[i+1] = X[i] + h / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
    return X,t

def func(t,X):
    x,y=X
    dx_dt=y+x-x**3
    dy_dt=-x
    return np.array([dx_dt, dy_dt])

IC=[[0,-1],[0,-2],[0,-3],[0,-4]];tmax=15;tmin=0;h=0.1
EA=[];R2A=[];R4A=[]
def main1():
    for d in range(0,4):
        res1=RK4(func,IC[d],tmin,tmax,h)
        w1,w2 = res1[0].T
        eu1=Eul(func,IC[d],tmin,tmax,h)
        f1,f2 = eu1[0].T
        re1=RK2(func,IC[d],tmin,tmax,h)
        k1,k2= re1[0].T
        EA.append(f1);EA.append(f2)
        R2A.append(k1);R2A.append(k2)
        R4A.append(w1);R4A.append(w2)
    t=res1[1]
    return EA,R2A,R4A,t
L=main1()
EA,R2A,R4A,t=L[0],L[1],L[2],L[3]
data1= {"t":t,"x (Eul)":EA[0],"x (RK2)":R2A[0],"x (RK4)":R4A[0],"y (Eul)":EA[1],"y (RK2)":R2A[1],"y(RK4)":R4A[1]}
data2= {"t":t,"x (Eul)":EA[2],"x (RK2)":R2A[2],"x (RK4)":R4A[2],"y (Eul)":EA[3],"y (RK2)":R2A[3],"y (RK4)":R4A[3]}
data3= {"t":t,"x (Eul)":EA[4],"x (RK2)":R2A[4],"x (RK4)":R4A[4],"y (Eul)":EA[5],"y (RK2)":R2A[5],"y (RK4)":R4A[5]}
data4= {"t":t,"x (Eul)":EA[6],"x (RK2)":R2A[6],"x (RK4)":R4A[6],"y (Eul)":EA[7],"y (RK2)":R2A[7],"y (RK4)":R4A[7]}

print("-----For x(0) = 0 and y(0) = -1-----")
print(pd.DataFrame(data1))
print("-----For x(0) = 0 and y(0) = -2-----")
print(pd.DataFrame(data2))
print("-----For x(0) = 0 and y(0) = -3-----")
print(pd.DataFrame(data3))
print("-----For x(0) = 0 and y(0) = -4-----")
print(pd.DataFrame(data4))

def xyt_g(t,x1_1,x1_2,x1_3,y1_1,y1_2,y1_3,x2_1,x2_2,x2_3,y2_1,y2_2,y2_3,x3_1,x3_2,x3_3,y3_1,y3_2,y3_3,x4_1,x4_2,x4_3,y4_1,y4_2,y4_3,tit1,tit2,tit3,tit4,titm):
     fig1, axs = plt.subplots(2, 2)
     axs[0, 0].plot(t, y1_1,marker="*",label="y (eul)") 
     axs[0, 0].plot(t,x1_1,marker="*",label="x(eul)") 
     axs[0, 0].plot(t,y1_2,marker="*",label="y(rk2") 
     axs[0, 0].plot(t,x1_2,marker="*",label="x(rk2)") 
     axs[0, 0].plot(t,y1_3,marker="*",label="y(rk4)") 
     axs[0, 0].plot(t,x1_3,marker="*",label="x(rk4)") 
     axs[0,0].set(xlabel = "time",title=tit1)
     axs[0,0].legend(loc = "best")
     axs[0,0].grid()
     axs[1, 0].plot(t,y2_1,marker="*",label="y (eul)") 
     axs[1, 0].plot(t,x2_1,marker="*",label="x(eul)") 
     axs[1, 0].plot(t,y2_2,marker="*",label="y(rk2") 
     axs[1, 0].plot(t,x2_2,marker="*",label="x(rk2)") 
     axs[1, 0].plot(t,y2_3,marker="*",label="y(rk4)") 
     axs[1, 0].plot(t,y2_3,marker="*",label="x(rk4)") 
     axs[1,0].set(xlabel = "time",title=tit2)
     axs[1,0].legend(loc = "best")
     axs[1,0].grid()
     axs[0, 1].plot(t, y3_1,marker="*",label="y (eul)") 
     axs[0, 1].plot(t,x3_1,marker="*",label="x(eul)")
     axs[0, 1].plot(t,y3_2,marker="*",label="y(rk2") 
     axs[0, 1].plot(t,x3_2,marker="*",label="x(rk2)") 
     axs[0, 1].plot(t,y3_3,marker="*",label="y(rk4)") 
     axs[0, 1].plot(t,x3_3,marker="*",label="x(rk4)") 
     axs[0,1].set(xlabel = "time",title=tit3)
     axs[0,1].legend(loc = "best")
     axs[0,1].grid()
     axs[1, 1].plot(t,y4_1,marker="*",label="y (eul)") 
     axs[1, 1].plot(t,x4_1,marker="*",label="x(eul)")
     axs[1, 1].plot(t,y4_2,marker="*",label="y(rk2") 
     axs[1, 1].plot(t,x4_2,marker="*",label="x(rk2)") 
     axs[1, 1].plot(t,y4_3,marker="*",label="y(rk4)") 
     axs[1, 1].plot(t,x4_3,marker="*",label="x(rk4)") 
     axs[1,1].set(xlabel = "time",title=tit4)
     axs[1,1].legend(loc = "best")
     axs[1,1].grid()
     fig1.suptitle(titm,c="brown")
     plt.show()
xyt_g(t,EA[0],R2A[0],R4A[0],EA[1],R2A[1],R4A[1],EA[2],R2A[2],R4A[2],EA[3],R2A[3],R4A[3],EA[4],R2A[4],R4A[4],EA[5],R2A[5],R4A[5],EA[6],R2A[6],R4A[6],EA[7],R2A[7],R4A[7],"x(0) = 0 and y(0) = -1","x(0) = 0 and y(0) = -2","x(0) = 0 and y(0) = -3","x(0) = 0 and y(0) = -4","x/y vs t plot")
def xy_g(x1_1,x1_2,x1_3,y1_1,y1_2,y1_3,x2_1,x2_2,x2_3,y2_1,y2_2,y2_3,x3_1,x3_2,x3_3,y3_1,y3_2,y3_3,x4_1,x4_2,x4_3,y4_1,y4_2,y4_3,tit1,tit2,tit3,tit4,main):
    fig2, axs = plt.subplots(2, 2)
    axs[0, 0].plot(x1_1, y1_1,marker=".",label="Eul") 
    axs[0, 0].plot(x1_2, y1_2,marker=".",label="RK2") 
    axs[0, 0].plot(x1_3, y1_3,marker=".",label="RK4") 
    axs[0,0].set(xlabel = "x",ylabel="y",title=tit1)
    axs[0,0].grid()
    axs[0,0].legend()
    axs[1, 0].plot(x2_1, y2_1,marker=".",label="Eul") 
    axs[1, 0].plot(x2_2, y2_2,marker=".",label="RK2") 
    axs[1, 0].plot(x2_3, y2_3,marker=".",label="RK4")
    axs[1,0].set(xlabel = "x",ylabel="y",title=tit2)
    axs[1,0].grid()
    axs[1,0].legend()
    axs[0, 1].plot(x3_1, y3_1,marker=".",label="Eul") 
    axs[0, 1].plot(x3_2, y3_2,marker=".",label="RK2") 
    axs[0, 1].plot(x3_3, y3_3,marker=".",label="RK4")
    axs[0,1].set(xlabel = "x",ylabel="y",title=tit3)
    axs[0,1].grid()
    axs[0,1].legend()
    axs[1, 1].plot(x4_1, y4_1,marker=".",label="Eul") 
    axs[1, 1].plot(x4_2, y4_2,marker=".",label="RK2") 
    axs[1, 1].plot(x4_3, y4_3,marker=".",label="RK4")
    axs[1,1].set(xlabel = "x",ylabel="y",title=tit4)
    axs[1,1].grid()
    axs[1,1].legend()
    fig2.suptitle(main,c="brown")
    plt.show()
xy_g(EA[0],R2A[0],R4A[0],EA[1],R2A[1],R4A[1],EA[2],R2A[2],R4A[2],EA[3],R2A[3],R4A[3],EA[4],R2A[4],R4A[4],EA[5],R2A[5],R4A[5],EA[6],R2A[6],R4A[6],EA[7],R2A[7],R4A[7],"x(0) = 0 and y(0) = -1","x(0) = 0 and y(0) = -2","x(0) = 0 and y(0) = -3","x(0) = 0 and y(0) = -4","y vs x plot")
def func2(t,X):
    x,y=X
    dx_dt=-4*x+4*y+12
    dy_dt=-1.6*x+1.2*y+4.8
    return np.array([dx_dt, dy_dt])
IC=[0,0];tmax=6;tmin=0;h=0.1
re1=Eul(func2,IC,tmin,tmax,h)
k1,k2= re1[0].T
r1=RK2(func2,IC,tmin,tmax,h)
g1,g2= r1[0].T
e1=RK4(func2,IC,tmin,tmax,h)
f1,f2= e1[0].T
N=int((tmax-tmin)/h)
t=np.linspace(0,6,N,endpoint=True)
A1=[ ];A2=[ ]
for i in t :
    I1=-8*np.exp(-2*i)+5*np.exp(-0.8*i)+3
    I2=-4*np.exp(-2*i)+4*np.exp(-0.8*i)
    A1.append(I1)
    A2.append(I2)
    
fig, (ax1, ax2) = plt.subplots(2)
fig.suptitle('Ques. 2 (Coupled equation in Linear Electric Circuit)',c="r")
ax1.plot(t,A1,marker="*",label=" I1 Ana")
ax1.plot(t,A2,marker="*",label="I2 Ana")
ax1.plot(t,f1,marker="*",label=" I1 RK4")
ax1.plot(t,f2,marker="*",label="I2 RK4")
ax1.plot(t,g1,marker="*",label=" I1 RK2")
ax1.plot(t,g2,marker="*",label="I2 RK2")
ax1.plot(t,k1,marker="*",label=" I1 Eul")
ax1.plot(t,k2,marker="*",label="I2 Eul")
ax1.set(xlabel="time (in sec)",ylabel="Current(in A)",title="Current vs Time plot")
ax1.legend()
ax1.grid()
ax2.plot(A1,A2,marker="*",label="Ana")
ax2.plot(k1,k2,marker="*",label="Eul")
ax2.plot(g1,g2,marker="*",label="RK2")
ax2.plot(f1,f2,marker="*",label="RK4")
ax2.set(xlabel="I1",ylabel="I2",title="I2 vs I1 plot")
ax2.legend()
ax2.grid()
plt.show()
data={"t":re1[1],"I1 (Euler)":k1,"I2(Euler)":k2,"I1 (RK2)":g1,"I2(RK2)":g2,"I1 (RK4)":f1,"I2(RK4)":f2,"I1 (Ana)":A1,"I2(Ana)":A2,}
print(pd.DataFrame(data))
