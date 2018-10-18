#include "../Method/Method.h"
#include <fstream>
#pragma once

using namespace std;

typedef unsigned int UINT;

class GaussMatrix: public Method {

	GaussMatrix();
	UINT swapCount;
	vector<double> pivot;

public:

	GaussMatrix(size_t);
	~GaussMatrix();
	void solve(ostream&) override;

private:

	void swapRows(size_t, size_t);
	size_t findMaxInRows(UINT);
	void makeUpperTriangular(Matrix<double>&, Vector<double>&, bool);
	void printUpperTriangular(ostream&);
	void getSolution(Matrix<double>&, Vector<double>&);
	double determinant() override;
	void findInverse() override;

};