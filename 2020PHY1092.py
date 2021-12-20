#University Roll No. = 20068567038
#Name = Pawanpreet Kaur
#College Roll No. = 2020PHY1092
import numpy as np                                          #importing required libraries
from scipy.special import legendre
import matplotlib.pyplot as plt
import pandas as pd

# Gamma Function
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

# 1. (a) Legendre Polynomial
def Lege(n,x):
    e=[]
    if n%2 == 0:       #even
        m=int(n/2)
    else :             #odd
        m = int((n-1)/2)   
    s = np.arange(0,m+1,1)                 
    for i in s:
        polyn = ( (-1)**i*gamma_(2*n-2*i+1)*(x**(n-2*i)))/(2**n*gamma_(i+1)*gamma_(n-i+1)*gamma_(n-2*i+1))
        e.append(polyn)      #appending all terms in list
    return (sum(e))            #returning sum of all terms

# 1. (b) 1st Derivative of legendre 
def DLege(n,x):
    e=[]
    if n%2 == 0:   #even
        m=int(n/2)
    else :       
        m = int((n-1)/2)                #odd
    s = np.arange(0,m+1,1,dtype=float)
    for i in s:
        polyn= (((-1)**i*gamma_(2*n-2*i+1)*(n-2*i))*( x**(n-2*i-1)))/(2**n*gamma_(i+1)*gamma_(n-i+1)*gamma_(n-2*i+1))
        e.append(polyn)              #appending all terms in list
    return sum(e)                  #returning sum of all terms
    
#4.  Inbuilt functions    
def V_inbuilt(n,x):
    Pn = legendre(n)          #inbuilt legendre
    y = Pn(x)
    Pn_d= Pn.deriv()  
    evald = np.polyval(Pn_d,x)     #1st derivative of legendre (inbuilt)
    return y,evald                   #returns legendre,1st derivative of legendre

#2. Recurrence Relation - n*Lege(n,x)= (2*n-1)*x*Lege(n-1,x) - (n-1)*Lege(n-2,x)
def vr(n,x,P_nx,P_n_1_x,P_n_2_x):
    q=np.array([1]*len(x))
    LHS=n*P_nx ; RHS= (2*n-1)*x*P_n_1_x - (n-1)*P_n_2_x
    if np.allclose(LHS,RHS):
       print("Recurrence relation n*P_n(x)= (2*n-1)*x*P_n-1(x) - (n-1)*P_n-2(x) is verified")
    else:
        print("Recurrence relation n*P_n(x)= (2*n-1)*x*P_n-1(x) - (n-1)*P_n-2(x) is not verified")
    return

#1 Calling function for findinng P0(x) ,P1(x),P2(x),P3(x),P'1(x),P'2(x),P'3(x)
p0=[];p1=[];p2=[];p3=[]
d_p1=[];d_p2=[];d_p3=[]                   #creating empty lists
x=np.linspace(-0.999,0.999,100)            #taking range [-0.999,0.999]
for i in x:                    #calculating at all points in x array
    p0.append(Lege(0,i))                         
    p1.append(Lege(1,i))
    p2.append(Lege(2,i))
    p3.append(Lege(3,i))
    d_p1.append(DLege(1,i))
    d_p2.append(DLege(2,i)) 
    d_p3.append(DLege(3,i))

#1.Printing data    
Data={"x":x,"P0(x)":p0,"P'1(x)":d_p1,"P1(x)":p1,"P'2(x)":d_p2,"P2(x)":p2}
print(pd.DataFrame(Data))

#2.Verifying Recurrence Relation
x=np.array(x);p3=np.array(p3);p2=np.array(p2);p1=np.array(p1)
vr(3,x,p3,p2,p1)

#calling inbuit function to find P1(x),P2(x),P3(x),P'1(x),P'2(x),P'3(x)
z=V_inbuilt(3,x)         # for P3(x),P'3(x)
u=V_inbuilt(2,x)         # for P2(x),P'2(x)
v=V_inbuilt(1,x)          # for P1(x),P'1(x)

#PLOT1
plt.plot(x,p1,label="$P_{0}(x)$".format(1),marker='o',c='red')                                 # P1(x) vs x (using my function)
plt.plot(x,d_p2,label="$P'_{0}(x)$".format(2),marker='*',c='green')                      # P'2(x) vs x (using my function)
plt.plot(x,v[0],label="$P_{0}(x)[inbuilt]$".format(1),linestyle='--',c='black')                 # P1(x) vs x (using inbuilt function)
plt.plot(x,u[1],label="$P'_{0}(x)[inbuilt]$".format(2),linestyle='--',c='maroon')            # P'2(x) vs x (using inbuilt function)
plt.xlabel("x")
plt.ylabel("Pn(x),P'n(x)")
plt.grid()
plt.legend()
plt.title("Legendre Polynomial Properties (Plot 1)")
plt.show()

#PLOT2
plt.plot(x,p2,label="$P_{0}(x)$".format(2),marker='o',c='red')              # P1(x) vs x (using my function)
plt.plot(x,d_p3,label="$P'_{0}(x)$".format(3),marker='*',c='green')               # P'3(x) vs x (using my function)
plt.plot(x,u[0],label="$P_{0}(x)[inbuilt]$".format(2),linestyle='--',c='black')             # P1(x) vs x (using inbuilt function)
plt.plot(x,z[1],label="$P'_{0}(x)[inbuilt]$".format(3),linestyle='--',c='maroon')  # P'3(x) vs x (using inbuilt function)
plt.xlabel("x")
plt.ylabel("Pn(x),P'n(x)")
plt.grid()
plt.legend()
plt.title("Legendre Polynomial Properties (Plot 2)")
plt.show()





