#include "../Method/Method.h"
#include <fstream>
#pragma once

using namespace std;

typedef unsigned int UINT;

class GaussMatrix: public Method {

	GaussMatrix();
	UINT swapCount;

public:

	GaussMatrix(size_t);
	~GaussMatrix();
	void solve(ostream&) override;

private:

	void swapRows(size_t, size_t);
	size_t findMaxInRows(UINT);
	void makeUpperTriangular();
	void printUpperTriangular(ostream&);
	void getSolution();

public:

	double determinant() const override;
	friend istream& operator>>(istream& in, GaussMatrix& obj);
	friend ostream& operator<<(ostream& out, GaussMatrix& obj);

};