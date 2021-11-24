import numpy as np

def check(A):
    dia=[];off =[]
    A=np.array(A)
    n=len(A)
    for i in range(n):
        temp=0
        for j in range(n):
            if i==j:
                dia.append(abs(A[i][j]))
            else:
                temp+=abs(A[i,j])
        off.append(temp)
    dia,off=np.array(dia),np.array(off)
    if np.all(dia > off):
       return True
    else:
         return False

def seidal(a,b,tol):
    n=len(a)
    x=np.zeros(n)
    it=0
    while True:
        x_old = np.array(x)
        for j in range(0, n):        
             d = b[j]                  
             for i in range(0, n):     
                  if j != i:
                     d-=a[j][i] * x[i]
             if a[j][j] ==0 :
                 sys.exit("divide by zero detected")
             else:
              x[j] = d / a[j][j]
        dif = abs( x- x_old)
        e=max(dif)
        it += 1
        if abs(e) <= tol:
            break
    return x,it

a=[[-4,1,0],[1,-5,1], [0,1,-3]]
b=[1,-2,1]

#a = [[-3,0,4],[1,8,6],[ 4,2, 8]]  #not dominant
#b = [4,7,3]

j,i=np.shape(a);k=len(b)
'''
if j==k and j==i :
    s=check(a)
    if s==False:
        print("Coefficient Matrix not diagonally dominant.")
    else :
         mat,itr=seidal(a,b,1e-12)
         print("Number of iterations = ",itr)
         print("Solution matrix = ",mat)
         print("Results using numpy.linalg.solve = ",np.linalg.solve(a,b))   
else:
   print("Number of equations should be equal to number of variables \nand number of rows in a must be equal to number of rows in b")
'''
if j==k and j==i :
         s=check(a)
         if s==False:
             print("Coefficient Matrix not diagonally dominant.")
         elif s==True:
             print("Coefficient Matrix diagonally dominant.")
         mat,itr=seidal(a,b,1e-5)
         mat,itr=seidal(a,b,1e-5)
         print("Number of iterations = ",itr)
         print("Solution matrix = ",mat)
         print("Results using numpy.linalg.solve = ",np.linalg.solve(a,b))   
else:
   print("Number of equations should be equal to number of variables \nand number of rows in a must be equal to number of rows in b")

