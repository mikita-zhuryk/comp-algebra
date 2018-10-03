#pragma once
#include <iostream>
#include "../Matrix/Matrix.h"
#include "../Vector/Vector.h"
#define OUTPUT_WIDTH 15
#define DECIMAL_PRECISION 9

using namespace std;

typedef unsigned int UINT;

class GaussMatrix {

	Matrix<double> A;
	Vector<double> b;
	Vector<double> x;
	double det;
	size_t n;
	UINT swapCount;
	
	void setDimensions(size_t);

	void swapRows(size_t, size_t);

	size_t findMaxInRows(UINT);

public:

	GaussMatrix(size_t);

	~GaussMatrix();

	void makeUpperTriangular();

	void printUpperTriangular(ostream&) const;

	void getSolution();

	void printSolution(ostream&) const;

	void deficiency(ostream& out) const;

	double determinant() const;

	friend istream& operator>>(istream& in, GaussMatrix& obj);

	friend ostream& operator<<(ostream& out, const GaussMatrix& obj);

};