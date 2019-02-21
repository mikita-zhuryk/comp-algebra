#include "PowerMethod.h"
#include "../Danilevsky/Danilevsky.cpp"
#include <fstream>

using namespace CMA;

PowerMethod::PowerMethod(size_t dim, int acc) : EigenMethod(dim), accuracy(pow(10, acc)) {
	lambda = Vector<double>(n, 0);
}

void PowerMethod::find(ostream& out) {
	x = Vector<double>(n, 1);
	Vector<double> x1;
	Vector<double> tempLambda(n, 0);
	int k = 1;
	do {
		lambda = tempLambda;
		x1 = A * x;
		tempLambda = findLambda(x1, x);
		x = x1 / x1.norm();
		++k;
	} while ((tempLambda - lambda).firstNorm() > accuracy);
	out << "Number of iterations: " << k << endl;
	eigenValues = tempLambda;
	eigenVectors.push_back(x);
	CMA::Danilevsky dan(n);
	ifstream fIn("../EigenMatrix.txt");
	fIn >> dan;
	dan.run(cout);
	eigenPolynomial = dan.getPoly();
}

Vector<double> PowerMethod::findLambda(Vector<double> xk1, Vector<double> xk) {
	Vector<double> temp(n);
	for (size_t i = 0; i < n; ++i) {
		temp[i] = xk1[i] / xk[i];
	}
	return temp;
}