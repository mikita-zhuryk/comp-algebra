#include "SimpleIteration.h"

using namespace CMA;

SimpleIteration::SimpleIteration(size_t dim) : Method(dim) {
	B = Matrix<double>(n, 0);
	g = Vector<double>(n, 0);
}

void SimpleIteration::solve(ostream& out) {
	buildBg();
	x = b;
	size_t i = 0;
	Vector<double> temp(B * x + g);
	out << "B: " << endl << B << endl << "B norm: " << endl << B.norm() << endl << endl;
	this->printSolution(out);
	double eps = pow(10, -6);
	while ((temp-x).norm() > eps) {
		x = temp;
		temp = B * x + g;
		out << i + 1 << " iteration " << endl;
		this->printSolution(out);
		++i;
	}
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