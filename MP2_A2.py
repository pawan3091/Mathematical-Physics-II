import numpy as np
from scipy.special import legendre
import matplotlib.pyplot as plt
from scipy.integrate import quad
import pandas as pd

# 1(a) Gamma Function
def gamma_(s):               #This works only for positive integers and positive half integers
    n= float(s)
    if n<0 or n%0.5!=0 :
       return ("Value must be either a positive integer or positive half integer.")
    else:
        if n==0:
           return "INFINITE"
        elif n==0.5:
           return np.sqrt(np.pi)
        elif n==1:
           return n
        else:
           return (n-1)*gamma_(n-1)    

# 1(b) Legendre 
def Lege(n,x):
    e=[]
    if n%2 == 0:       #even
        m=int(n/2)
    else :             #odd
        m = int((n-1)/2)   
    s = np.arange(0,m+1,1)
    for i in s:
        polyn = ( (-1)**i*gamma_(2*n-2*i+1)*(x**(n-2*i)))/(2**n*gamma_(i+1)*gamma_(n-i+1)*gamma_(n-2*i+1))
        e.append(polyn)
    return (sum(e))

# 1(c) Derivative of legendre 
def DLege(n,x):
    e=[]
    if n%2 == 0:
        m=int(n/2)
    else :
        m = int((n-1)/2)
    s = np.arange(0,m+1,1,dtype=float)
    for i in s:
        polyn= (((-1)**i*gamma_(2*n-2*i+1)*(n-2*i))*( x**(n-2*i-1)))/(2**n*gamma_(i+1)*gamma_(n-i+1)*gamma_(n-2*i+1))
        e.append(polyn)
    return sum(e)
    
# 1(d) Verifying with inbuilt functions    
def V_inbuilt(n,x):
    Pn = legendre(n)
    y = Pn(x)
    Pn_d= Pn.deriv()  
    evald = np.polyval(Pn_d,x)  
    return y,evald

def check(n,x):

    if x<-1 or x>1 or n<0 :
       print("Check your values . x should lie in [-1,1] and n>=0")
    else:
        print("Legendre polynomial & its derivative value at given x & n using my function = ",(Lege(n,x),DLege(n,x)))  
        print("Legendre polynomial & its derivative value at given x & n using inbuilt function =",V_inbuilt(n,x))  
           
# 2(c) Recurrence Relations

#1. n*Lege(n,x)= x*DLege(n,x)- DLege(n-1,x)
def vr1(n,x,P_nx,D_P_nx,D_P_n_1_x):
    q=np.array([1]*len(x))
    LHS=n*P_nx ; RHS= D_P_nx*x-D_P_n_1_x
    if np.allclose(LHS,RHS):
       print("Recurrence relation 1 is verified")
    else:
        print("Recurrence relation 1 is not verified")
    DataOut111 = np.column_stack((x,n*q,(n-1)*q,P_nx,D_P_nx,D_P_n_1_x,LHS,RHS))
    np.savetxt('data/leg02.dat', DataOut111,delimiter=',')
    Data={"x":x,"n":n*q,"n-1":(n-1)*q,"Pn(x)":P_nx,"P'n(x)":D_P_nx,"P'(n-1)(x)":D_P_n_1_x,"LHS":LHS,"RHS":RHS}
    print(pd.DataFrame(Data))
    return
  
#2. (2*n+1)*x*Lege(n,x)= (n+1)*Lege(n+1,x) + n*Lege(n-1,x)
def vr2(n,x,P_nx,P_n_1_x,P_n_2_x):
    q=np.array([1]*len(x))
    LHS=(2*n+1)*x*P_nx ; RHS= (n+1)*P_n_1_x + n*P_n_2_x
    if np.allclose(LHS,RHS):
       print("Recurrence relation 2 is verified")

    else:
        print("Recurrence relation 2 is not verified")
    DataOut111 = np.column_stack((x,n*q,(n-1)*q,(n+1)*q,P_nx,P_n_1_x,P_n_2_x,LHS,RHS))
    np.savetxt('data/leg03.dat', DataOut111,delimiter=',')
    Data={"x":x,"n":n*q,"n-1":(n-1)*q,"n+1":(n+1)*q,"Pn(x)":P_nx,"P(n+1)(x)":P_n_1_x,"P(n-1)(x)":P_n_2_x,"LHS":LHS,"RHS":RHS}
    print(pd.DataFrame(Data))
    return

#3. n*Lege(n,x)= (2*n-1)*x*Lege(n-1,x) - (n-1)*Lege(n-2,x)
def vr3(n,x,P_nx,P_n_1_x,P_n_2_x):
    q=np.array([1]*len(x))
    LHS=n*P_nx ; RHS= (2*n-1)*x*P_n_1_x - (n-1)*P_n_2_x
    if np.allclose(LHS,RHS):
       print("Recurrence relation 3 is verified")
    else:
        print("Recurrence relation 3 is not verified")
    DataOut111 = np.column_stack((x,n*q,(n-1)*q,(n+1)*q,P_nx,P_n_1_x,Lege(n+2,x),LHS,RHS))
    np.savetxt('data/leg04.dat', DataOut111,delimiter=',')
    Data={"x":x,"n":n*q,"n-1":(n-1)*q,"n+1":(n+1)*q,"Pn(x)":P_nx,"P(n-1)(x)":P_n_1_x,"P(n-2)(x)":Lege(n+2,x),"LHS":LHS,"RHS":RHS}
    print(pd.DataFrame(Data))
    return

#IMPLEMENTATION 

# 2(a)
p0=[];p1=[];p2=[];p3=[]
x=np.linspace(-1,1,100)
for i in x:
    p0.append(Lege(0,i))
    p1.append(Lege(1,i))
    p2.append(Lege(2,i))
    p3.append(Lege(3,i))
        
DataOut = np.column_stack((x,p0,p1,p2,p3))
Data1={"x":x,"P0(x)":p0,"P1(x)":p1,"P2(x)":p2,"P3(x)":p3}
print(pd.DataFrame(Data1))
np.savetxt('data/leg00.dat', DataOut,delimiter=',')
v=np.loadtxt('data/leg00.dat',unpack=True,delimiter=',',dtype='float')
#print(v)
for u in range(1,5):
    plt.plot(v[0],v[u],label= "$P_{0}$".format(u))   
plt.xlabel("x")
plt.ylabel("Pn(x)")
plt.grid()
plt.legend()
plt.title("Pn(x) vs x graph [Legendre Polynomial 2a]")
plt.show()

# 2(b)
d_p1=[];d_p2=[];d_p3=[]
for i in x:
    d_p1.append(DLege(1,i))
    d_p2.append(DLege(2,i))
    d_p3.append(DLege(3,i))

DataOut1 = np.column_stack((x,d_p1,d_p2,d_p3))
Data2={"x":x,"P'1(x)":d_p1,"P'2(x)":d_p2,"P'3(x)":d_p3}
print(pd.DataFrame(Data2))
np.savetxt('data/leg01.dat', DataOut1,delimiter=',')
r=np.loadtxt('data/leg01.dat',unpack=True,delimiter=',',dtype='float')
#print(r)
plt.plot(r[0],DLege(0,x),label="$P'_{0}$".format(0))
plt.plot(r[0],r[2],label="$P'_{0}$".format(2))
plt.plot(r[0],v[2],label="$P_{0}$".format(1))
plt.xlabel("x")
plt.ylabel("P'n(x)")
plt.grid()
plt.legend()
plt.title("P'n(x) vs x graph [2b]")
plt.show()

vr1(2,v[0],v[3],r[2],r[1]) # 2c 1.
vr2(2,v[0],v[3],v[4],v[2]) # 2c 2. 
vr3(3,v[0],v[4],v[3],v[2]) # 2c 3. 

LHS = np.ones((4,4))
RHS = np.identity(4)
for i in range(4):
     for j in range(4):
          LHS[i][j]=np.sum(v[i+1]*v[j+1]*np.diff(v[0])[0])
          RHS[i][j]*=2/(2*j+1)

l1,l2=[],[]
for n in range(3):
    for m in range(3):
        if n == m:
           l1.append(2/(2*n+1))
        else:
           l1.append(0)
        f=legendre(n)*legendre(m)
        inte , err = quad(f, -1, 1)
        l2.append(inte)
RHS = np.array(l2).reshape(3,3)
LHS = np.array(l1).reshape(3,3)
print("RHS = ")
print(RHS)
print("LHS = ")
print(LHS)
if np.allclose(LHS,RHS):
   print("Orthogonality relation verified")
else:
    print("Orthogonality relation not verified")
print()
print("Comparing our results for a given value of x and n with inbuilt one")
x1=float(input("Enter the value of x :"))
n=float(input("Enter the value of n :"))
check(n,x1)

