#include "GradientDescent.h"

using namespace CMA;

GradientDescent::GradientDescent(size_t dim, int acc) : Method(dim), accuracy(pow(10, acc)) {
}

void GradientDescent::solve(ostream& out) {
	b = A.transpose() * b;
	A = A.transpose() * A;
	x = b;
	int i = 1;
	Vector<double> def;
	double t;
	Vector<double> temp;
	do {
		if (i != 1) {
			x = temp;
		}
		def = calcDeficiency();
		t = def.dot(def) / (A * def).dot(def);
		temp = x - def * t;
		++i;
	} while ((temp - x).norm() > accuracy);
	out << "Number of iterations: " << i << endl << endl;
}

double GradientDescent::determinant() {
	return 0;
}