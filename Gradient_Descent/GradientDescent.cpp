#include "GradientDescent.h"

using namespace CMA;

GradientDescent::GradientDescent(size_t dim, int acc) : Method(dim), accuracy(pow(10, acc)) {
}

void GradientDescent::solve(ostream& out) {
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
		out << i << " iteration:\nDeficiency:\n" << def << endl;
		printSolution(out);
		++i;
	} while ((temp - x).norm() > accuracy);
}

double GradientDescent::determinant() {
	return 0;
}