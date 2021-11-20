import numpy as np
import pandas as pd
import sys

#2a
def Gauss_Elimination(A):
    n=len(A)
    for i in range(n):
        #if A[i][i] == 0:
           #return 0       
        for j in range(i+1, n):
             ratio = A[j][i]/A[i][i]
             for k in range(n+1):
                   A[j][k] = A[j][k] - ratio * A[i][k]
             print(np.array(A))
    return A

#2b
def Back_Subs(A):
    n=len(A)
    x=np.zeros([n])
    if A[n-1][n-1]==0:
        sys.exit("Divide by zero")
    x[n-1] = A[n-1][n]/A[n-1][n-1]
    for i in range(n-2,-1,-1):
        x[i] = A[i][n]
        for j in range(i+1,n):
            x[i] = x[i] - A[i][j]*x[j]
        x[i] = x[i]/A[i][i]
    return x

#3a      
A=[[-4,1,0,1],
   [-1,0,-1,2],
   [0,1,-3,1]]
print("Augmented matrix:")
print(np.array(A))
print("Intermediate Steps:")
j=Gauss_Elimination(A)
if j==0:
    print("System of equations is inconsistent")
else:
    print("Reduced Row Echelon Form :")    
    print(np.array(j))
    x = Back_Subs(j)
    print('Solution: ')
    print(' I1=','%0.6f'%x[0],' I2=','%0.6f'%x[1],'I3=' ,'%0.6f'%x[2] )
    
'''
#2c    
def Gauss_Seidal(f1,f2,f3,e):
      x0,y0,z0= 0,0,0
      t = 1
      x_a=[x0];y_a=[y0];z_a=[z0]
      condition = True
      while  condition :
          x1 = f1(x0,y0,z0)
          y1 = f2(x1,y0,z0)
          z1 = f3(x1,y1,z0)
          x_a.append(x1);y_a.append(y1);z_a.append(z1)
          e1 = abs(x0-x1); e2 = abs(y0-y1); e3 = abs(z0-z1)
          x0,y0,z0 = x1,y1,z1
          t += 1
          condition= e1>e and e2>e and e3>e                 
      data={"x":x_a,"y":y_a,"z":z_a}
      print(pd.DataFrame(data))
      print('Solution: x=','%0.6f'%x1,' y=','%0.6f'%y1,'z=' ,'%0.6f'%z1 )
'''
'''
f1 = lambda x,y,z: (y-1)/4
f2 = lambda x,y,z: (2+x+z)/5
f3 = lambda x,y,z: (y-1)/3
Gauss_Seidal(f1,f2,f3,1e-6)
'''
