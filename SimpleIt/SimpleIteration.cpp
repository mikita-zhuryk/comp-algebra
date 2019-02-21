#include "SimpleIteration.h"

using namespace CMA;

SimpleIteration::SimpleIteration(size_t dim, int acc) : Method(dim), accuracy(pow(10, acc)) {
	B = Matrix<double>(n, 0);
	g = Vector<double>(n, 0);
}

void SimpleIteration::solve(ostream& out) {
	buildBg();
	x = A.transpose() * b;
	int k = 1;
	Vector<double> temp(n, 0);
	out << "B: " << endl << B << endl << "B norm: " << B.norm() << endl << endl;
	do {
		if (k != 1) {
			x = temp;
		}
		temp = B * x + g;
		++k;
	} while ((temp - x).norm() > accuracy);
	out << "Number of iterations: " << k << endl << endl;
}

void SimpleIteration::buildBg() {
	auto ATA = A.transpose() * A;
	B = Matrix<double>::identity(n, 1) - ATA * (1.0 / ATA.norm());
	g = A.transpose() * b * (1.0 / ATA.norm());
}

double SimpleIteration::determinant() {
	return 0;
}