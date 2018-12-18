#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np

A = np.array([[0.78523, 0.01605, 0.14915, -0.08940, 0.63189],
              [0.01605, 1.36564, -0.00184, 0.12625, 0.06519],
              [0.14915, -0.00184, 0.86805, -0.24137, 0.32770],
              [-0.08940, 0.12625, -0.24137, 0.88976, -0.05569],
              [0.63189, 0.06519, 0.32770, -0.05569, 1.19883]])


def Frobenius(A):
    F = A
    invS = np.eye(5)
    S = np.eye(5)
    for i in range(4, 0, -1):
        M = np.eye(5)
        invM = np.eye(5)
        M[i - 1] = -F[i] / F[i][i - 1]
        M[i - 1][i - 1] = 1 / F[i][i - 1]
        invM[i - 1] = F[i]
        S = np.matmul(S, M)
        invS = np.matmul(invM, invS)
        F = np.matmul(np.matmul(invM, F), M)
    return (F, invS, S)

def EigenValues(F):
    coeffs = ((-1) ** 5) * F[0][::-1]
    coeffs = np.append(coeffs, 1)
    eigenvalues = np.polynomial.polynomial.polyroots(coeffs)
    return eigenvalues

def EigenVectors(S, eigenvalues):
    eigenvectors = []
    for v in eigenvalues:
        y = np.array(np.zeros(np.shape(A[0])))
        for i in range(0, 5):
            y[4 - i] = v ** (i)
        eigenvector = np.matmul(S, y)
        eigenvectors.append(eigenvector / np.linalg.norm(eigenvector))
    return eigenvectors

if __name__ == "__main__":
    (F, invS, S) = Frobenius(A)
    eigenvalues = EigenValues(F)
    eigenvectors = EigenVectors(S, eigenvalues)
    i = 0
    for v in eigenvectors:
        print(v)
    for v in eigenvectors:
        print(np.matmul(A, v) - eigenvalues[i] * v)
        i = i + 1

