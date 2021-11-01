import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def RK2(func, X0,tmin,tmax,h,cons):
    N=int((tmax-tmin)/h)
    t = np.linspace(tmin,tmax,N)
    X  = np.zeros([N, len(X0)])
    X[0] = X0
    for i in range(N-1):
        k1 =h* func(t[i],X[i],cons)
        k2 = h*func( t[i] + h,X[i] +  k1,cons)
        X[i+1] = X[i] +  (k1 +k2 )/2
    return X,t

def f2(t,X,cons):
    k,m,b=cons
    x,z=X
    dx_dt=z
    dx2_dt2=-b*z/m-k*x/m
    return np.array([dx_dt, dx2_dt2])

def f3(t,X,cons):
    g,L=cons
    th,z=X
    dth_dt=z
    dth2_dt2=-g*th/L
    return np.array([dth_dt, dth2_dt2])

def f4(t,X,cons):
    k,m,g,L=cons
    w0=np.sqrt(g/L)
    xa,z1,xb,z2=X
    xb=np.radians(xb)
    xa=np.radians(xa)
    e1=z1
    e2=-(w0**2)*xa-(k/m)*(xa-xb)
    e3=z2
    e4=-(w0**2)*xb+(k/m)*(xa-xb)
    return np.array([e1,e2,e3,e4])

#3a
IC=[2,0];tmin=0;m=0.5;k=4;h=0.01;b=0#SI units
tp=2*np.pi*np.sqrt(m/k)
tmax=5*tp
cons=(k,m,b)
S=RK2(f2,IC,tmin,tmax,h,cons)
x,v= S[0].T
t=S[1]
t1=t/tp
fig, axs = plt.subplots(2)
fig.suptitle('SIMPLE HARMONIC OSCILLATOR',c="r")
axs[0].plot(t1,x,marker="*",c="orange")
axs[0].set(xlabel="time/time period ",ylabel="Dispacement")
axs[0].grid()
axs[1].plot(t1,v,marker="*",c="cyan")
axs[1].set(xlabel=" Time/Time Period ",ylabel="Velocity")
axs[1].grid()
plt.show()
print("----------------------------- SIMPLE HARMONIC OSCILLATOR---------------------------------")
d={"No. of time periods":t1,"Displacement":x,"Velocity":v}
print(pd.DataFrame(d))


#3b https://scipy-lectures.org/intro/scipy/auto_examples/plot_odeint_damped_spring_mass.html (ref)
IC=[1,0];tmin=0;h=0.1;m=0.5;k=4 #SI units
tmax=100
dis=[];vel=[];time=[]
ba=[0.2,np.sqrt(4*m*k),3]
for b in ba:
     cons=(m,k,b)
     S=RK2(f2,IC,tmin,tmax,h,cons)
     x,v= S[0].T
     t=S[1]
     dis.append(x)
     vel.append(v)
fig, ax = plt.subplots(3,2)
fig.suptitle('DAMPED HARMONIC OSCILLATOR',c="r")
ax[0,0].plot(t,dis[0],marker="*",c="r")
ax[0,0].set(xlabel="time ",ylabel="Dispacement",title="Underdamped")
ax[0,0].grid()
ax[0,1].plot(t,vel[0],marker="*",c="green")
ax[0,1].set(xlabel="time ",ylabel="Velocity",title="Underdamped")
ax[0,1].grid()
ax[1,0].plot(t,dis[1],marker="*",c="purple")
ax[1,0].set(xlabel="time ",ylabel="Dispacement",title="Critically Damped")
ax[1,0].grid()
ax[1,1].plot(t,vel[1],marker="*",c="darkblue")
ax[1,1].set(xlabel="time ",ylabel="Velocity",title="Critically Damped")
ax[1,1].grid()
ax[2,0].plot(t,dis[2],marker="*",c="brown")
ax[2,0].set(xlabel="time ",ylabel="Dispacement",title="Overdamped")
ax[2,0].grid()
ax[2,1].plot(t,vel[2],marker="*",c="violet")
ax[2,1].set(xlabel="time ",ylabel="Velocity",title="Overdamped")
ax[2,1].grid()
plt.show()
print("---------------------- DAMPED HARMONIC OSCILLATOR  (UNDERDAMPED)-----------------------")
d={"Time":t,"Displacement":dis[0],"Velocity":vel[0]}
print(pd.DataFrame(d))
print("---------------------- DAMPED HARMONIC OSCILLATOR  (CRITICALLY DAMPED)-----------------------")
d={"Time":t,"Displacement":dis[1],"Velocity":vel[1]}
print(pd.DataFrame(d))
print("---------------------- DAMPED HARMONIC OSCILLATOR  (OVERDAMPED)-----------------------")
d={"Time":t,"Displacement":dis[2],"Velocity":vel[2]}
print(pd.DataFrame(d))

#3c
IC=[1,0];tmin=0;g=9.8;L=2  #SI units
cons1=(g,L)
tp=2*np.pi*np.sqrt(L/g)
tmax=10*tp
h=tp/50
S=RK2(f3,IC,tmin,tmax,h,cons1)
x,v= S[0].T
t=S[1]
t1=t/tp
fig, axs = plt.subplots(2)
fig.suptitle('SIMPLE PENDULUM',c="r")
axs[0].plot(t1,x,marker="*",c="r")
axs[0].set(xlabel="Time/Time Period ",ylabel="Angular Dispacement")
axs[0].grid()
axs[1].plot(t1,v,marker="*",c="green")
axs[1].set(xlabel=" Time/Time Period ",ylabel="Angular Velocity")
axs[1].grid()
plt.show()
print("------------------------------ SIMPLE PENDULUM------------------------------")
d={"No. of time periods":t,"Angular Displacement":x,"Angular Velocity":v}
print(pd.DataFrame(d))


#3d
IC=[10,0,-10,0];tmax=80;tmin=0;h=0.01;m1=10;k1=90;g1=9.8;l1=10  #SI units
cons1=(k1,m1,g1,l1)
S=RK2(f4,IC,tmin,tmax,h,cons1)
x1,v1,x2,v2= S[0].T
t3=S[1]
fig, axs = plt.subplots(2, 2)
fig.suptitle('COUPLED SYSTEM',c="r")
axs[0, 0].plot(t3,x1,c="r")
axs[0,0].set(xlabel=" Time",ylabel="Displacement",title="Mass A")
axs[0,0].grid()
axs[0, 1].plot(t3,v1,c="green")
axs[0,1].set(xlabel=" Time",ylabel="Velocity",title="Mass A")
axs[0,1].grid()
axs[1, 0].plot(t3,x2,c="magenta")
axs[1,0].set(xlabel=" Time",ylabel="Displacement",title="Mass B")
axs[1,0].grid()
axs[1, 1].plot(t3,v2,c="green")
axs[1,1].set(xlabel=" Time",ylabel="Velocity",title="Mass B")
axs[1,1].grid()
plt.show()
print("---------------------- COUPLED PENDULUM-----------------------")
print("-------------------------------- MASS A ---------------------------------")
d={"Time":t3,"Displacement":x1,"Velocity":v1}
print(pd.DataFrame(d))
print("-------------------------------- MASS B ---------------------------------")
d={"Time":t3,"Displacement":x2,"Velocity":v2}
print(pd.DataFrame(d))