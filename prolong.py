import numpy as np
from sympy import *
from itertools import combinations_with_replacement

def jet(X,U):
    J = [U[i]+X[j] for i in range(len(U)) for j in range(len(X))]
    return J
                
def njet(X,U,n):
    J = U
    for i in range(n):
        J = jet(X,J)
    for i in range(len(J)):
        k = i+1
        for j in range(i+1,len(J)-1):
            if sorted(J[i]) == sorted(J[k]):
                J.remove(J[k])
                k = k-1
            k += 1
        if i >= len(J)-2:
            break            
    return J

def fulljet(X,U,n):
    J = X+U
    for i in range(n):
        J = J+njet(X,U,i+1)
    return J

def uJalpha(u,x,X,U):
    J = njet(X,U,len(str(u)))
    for i in range(len(J)):
        if sorted(str(u)+str(x)) == sorted(J[i]):
            return var(J[i])
            break
            
def divUJi(X,U,u,J,i):
    f = u
    for j in range(len(J)):
        f = uJalpha(f,X[J[j]],X,U)
    f = uJalpha(f,X[i],X,U)
    return f


def totDiv(X,U,P,i,n):
    expr = diff(P,X[i])
    J = filter(lambda v: v not in X, fulljet(X,U,n))
    for j in range(len(J)):
        expr = expr + diff(P,J[j])*uJalpha(J[j],X[i],X,U)
    return expr

def TotDiv(X,U,P,n,J):
    expr = totDiv(X,U,P,J[0],n)
    for i in range(1,len(J)):
        expr = totDiv(X,U,expr,J[i],n)
    return expr

def phiAlpha(X,U,v,J):
    phi = v[len(X):]
    xi = v[:len(X)]
    a = 0
    q = len(U)-1
    for i in range(len(X)):
        a = a+xi[i]*uJalpha(U[q],X[i],X,U)
    b = 0
    for i in range(len(X)):
        b = b+xi[i]*divUJi(X,U,U[q],J,i)
    c = TotDiv(X,U,phi[q]-a,len(J),J)+b
    return c

def diffOrd(p,n):
    b = []
    for i in range(1,n+1):
        b.append(list(combinations_with_replacement(range(p),i)))
    return b

def Prolong(X,U,v,n):
    a = v
    J = diffOrd(len(X),n)
    for i in range(n):
        for j in range(len(X)):
            a.append(phiAlpha(X,U,v,J[i][j]))
    return a
        

if __name__ == "__main__":
    
    #initialize symbols for sympy
    p = input("Input the number of independent variables:")
    q = input("Input the number of dependent variables:")
    X = []
    U = []
    for i in range(p):
        X.append(raw_input("What label would you like to assign to independent variable "+str(i+1)+"?"))
    for j in range(q):
        U.append(raw_input("What label would you like to assign to dependent variable "+str(j+1)+"?"))
    var(fulljet(X,U,3))
    P = filter(lambda v: v not in X, fulljet(X,U,3))
