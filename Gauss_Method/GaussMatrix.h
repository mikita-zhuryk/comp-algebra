#include "../General_Method/Method.h"
#include <fstream>
#pragma once

using namespace std;

typedef unsigned int UINT;

class GaussMatrix: public Method {

	GaussMatrix();

	UINT swapCount;

	void swapRows(size_t, size_t);

	size_t findMaxInRows(UINT);

public:

	GaussMatrix(size_t);

	~GaussMatrix();

	void makeUpperTriangular();

	void printUpperTriangular(ostream&);

	void getSolution();

	double determinant() const override;

	friend istream& operator>>(istream& in, GaussMatrix& obj);

	friend ostream& operator<<(ostream& out, GaussMatrix& obj);

};