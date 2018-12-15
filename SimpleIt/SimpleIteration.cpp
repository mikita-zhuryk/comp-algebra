#include "SimpleIteration.h"

using namespace CMA;

SimpleIteration::SimpleIteration(size_t dim, int acc) : Method(dim), accuracy(pow(10, acc)) {
	B = Matrix<double>(n, 0);
	g = Vector<double>(n, 0);
}

void SimpleIteration::solve(ostream& out) {
	buildBg();
	x = b;
	int k = 1;
	Vector<double> temp(n, 0);
	out << "B: " << endl << B << endl << "B norm: " << endl << B.norm() << endl << endl;
	do {
		if (k != 1) {
			x = temp;
		}
		temp = B * x + g;
		out << k << " iteration:\nDeficiency:\n" << calcDeficiency() << endl;
		printSolution(out);
		++k;
	} while ((temp - x).norm() > accuracy);
}

void SimpleIteration::buildBg() {
	auto ATA = A.transpose() * A;
	B = Matrix<double>::identity(n, 1) - ATA * (1.0 / ATA.norm());
	g = A.transpose() * b * (1.0 / ATA.norm());
	/*B = A;
	g = b;
	for (size_t i = 0; i < n; ++i) {
		g[i] /= A[i][i];
		B[i] /= A[i][i];
		B[i][i] = 0;
	}
	B *= -1;*/
}

double SimpleIteration::determinant() {
	return 0;
}