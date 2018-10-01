#pragma once
#include "IndexComparator.h"
#include <iostream>
#include <vector>
#define OUTPUT_WIDTH 15
#define DECIMAL_PRECISION 9

using namespace std;

typedef unsigned int UINT;
typedef vector<pair<size_t, double>> indexedVector;

class GaussMatrix {

	double** values;
	indexedVector b;
	indexedVector x;
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

	void calculate() const;

	friend istream& operator>>(istream& in, GaussMatrix& obj);

	friend ostream& operator<<(ostream& out, const GaussMatrix& obj);

};