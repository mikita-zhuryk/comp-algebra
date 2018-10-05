#pragma once
#include <iostream>
#include "../Matrix/Matrix.h"
#include "../Vector/Vector.h"

using namespace std;

typedef unsigned int UINT;

class GaussMatrix {

	Matrix<double> A;
	Matrix<double> initial_A;
	Matrix<double> inverse;
	Vector<double> b;
	Vector<double> initial_b;
	Vector<double> x;
	double detA;
	size_t n;
	UINT swapCount;
	
	void setDimensions(size_t);

	void swapRows(size_t, size_t);

	size_t findMaxInRows(UINT);

public:

	GaussMatrix(size_t);

	~GaussMatrix();

	void makeUpperTriangular();

	void printUpperTriangular(ostream&);

	void getSolution();

	void printSolution(ostream&);

	void deficiency(ostream& out);

	double determinant() const;

	void printInverse(ostream& out);

	void inverseDeficiency(ostream& out);

	friend istream& operator>>(istream& in, GaussMatrix& obj);

	friend ostream& operator<<(ostream& out, GaussMatrix& obj);

};