import numpy as np
import sys

def Swap(a,Row):
    n=len(a)
    r=Row
    a=np.array(a)
    a=a.astype(float)
    m=int(np.size(a)/n)
    num=1
    while (a[r][r]==0):
        if (r==n-1):
            a[[n-1, 0],:]=a[[0, n-1],:]
        else:
            a[[r+num, r],:]=a[[r, r+num],:]
        num=num+1
        if num==n:
            return a,0
    return a,1

def Rank_m(a):
    n=len(a)
    m=int(np.size(a)/n)
    a=np.array(a)
    no_zero_rows=len(np.where(~a.any(axis=1))[0])
    coeff_mat=np.zeros((n,m-1))
    for i in range(n):
        for j in range(m-1):
            coeff_mat[i][j]=a[i][j]
    no_zero_rows2=len(np.where(~coeff_mat.any(axis=1))[0])
    rank1=n-no_zero_rows
    rank2=n-no_zero_rows2
    return rank1,rank2

def Clmn_r(a,rmv):
    n=len(a)
    m=int((np.size(a))/n)
    new_a=np.zeros((n,m-1))
    for i in range(n):
        for j in range(1,m,1):
            new_a[i][j-1]=a[i][j]
    return new_a
    
def Gauss_Elimination(a):
    n=len(a) 
    m=int(np.size(a)/len(a))
    flag=1
    for i in range(n):
        if a[i][i]==0:
            print("row interchange")
            a,p=Swap(a,i)
            print(a)
            if p==0:
                print("removing redundant column")
                a=Clmn_r(a,i)
                m=m-1       
        for j in range(i+1,n): #elements lower than pivot element moving columnwise
            r=a[j][i]/a[i][i]  #ratio according to the columnn
            for k in range(m):  
                a[j][k]=round((a[j][k]-(r*a[i][k])),3) #row operation
            print("intermediate steps")
            print(np.array(a))
            rank1,rank2=Rank_m(a)
            if rank1<n:
                return a,rank1,rank2
    return a,rank1,rank2

def back_subs(a,n):
    x=np.zeros([n])
    x[n-1] =round(a[n-1][n]/a[n-1][n-1],4) #solving last row
    for i in range(n-2,-1,-1):     #moving upwards
        x[i] = a[i][n]             #assigning last column value to x
        for j in range(i+1,n):     #solving for each row 
            x[i] = x[i] - a[i][j]*x[j]  #subtracting the values from last column
        x[i] = round(x[i]/a[i][i],4)    #dividing it by appropriate pivot value 
    return x
    
Mat=[[-4,1,0,1],
         [1,-5,1,-2],
         [0,1,-3,1]]

#Mat=[[6,8,1,3],[3,7,9,1],[4,6,1,2],[3,3,8,1]]   #overdetermined
#Mat=[[5,7,2,3,6],[8,4,6,1,1,0],[7,4,5,6,7]]        #undetermined
    
n=len(Mat)
m=int(np.size(Mat)/n)
coeff_mat=np.zeros((n,m-1))

print("Coefficient Matrix:")
for i in range(n):
    for j in range(m-1):
        coeff_mat[i][j]=Mat[i][j]
print(np.array(coeff_mat))

print("Augmented Matrix")
print(np.array(Mat))

ans,r1,r2=Gauss_Elimination(Mat)

print("Reduced Row Echelon Form")
print(np.array(ans))
print("Rank (Augmented Matrix): ",r1)
print("Rank (Coefficient Matrix): ",r2)

if m-1==r1: #no. of variable is equal to rank
    x=back_subs(ans,r1)
    print("Unique Solution")
    for i in range(len(x)):
        print("x"+str(i+1),":",x[i])
elif(m-1<r1):
    print("Overdetermined")
else:
    print("Underdetermined")
