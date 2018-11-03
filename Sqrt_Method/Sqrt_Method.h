#include "../Method/Method.h"
#include <fstream>
#pragma once

using namespace std;

class Sqrt_Method: public Method {

	Matrix<double> S;
	Vector<double> y;

public:

	Sqrt_Method(size_t);

private:

	void solve(ostream&) override;
	void calculateS();
	void solveForY();
	void solveForX();
	double determinant() override;
	void findInverse() override;
	void printInverse(ostream& out) override;
	void printInverseDeficiency(ostream&) override;

};