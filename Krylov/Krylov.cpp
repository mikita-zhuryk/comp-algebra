#include "Krylov.h"
#include "../Danilevsky/Danilevsky.cpp"
#include "../Reflections/Reflections.cpp"

using namespace CMA;

Krylov::Krylov(size_t dim) : EigenMethod(dim) {
	dimA = n;
	C = Vector<double>(n, 1);
	CMatrix = Matrix<double>(n, 0);
}

void Krylov::find(ostream& out) {
	buildPolynomial();
	findEigenVectors();
}

void Krylov::buildPolynomial() {
	Vector<double> temp(n);
	Matrix<double> powerA(Matrix<double>::identity(n, 1));
	for (uint16_t i = 0; i < dimA; ++i) {
		temp = powerA * C;
		CMatrix[n - 1 - i] = temp;
		powerA *= A;
	}
	temp = powerA * C;
	CMatrix = CMatrix.transpose();
	CMA::Reflections refl(CMatrix, temp);
	refl.run(cout);
	Vector<double> coeffs = refl.getX();
	eigenPolynomial[0] = 1;
	for (size_t i = 0; i < n; ++i) {
		eigenPolynomial[i + 1] = -coeffs[i];
	}
}

void Krylov::findEigenVectors() {
	for (size_t i = 0; i < n; ++i) {
		eigenVectors.push_back(findEigenVec(eigenValues[i]));
	}
}

Vector<double> Krylov::findEigenVec(double v) {
	Vector<double> beta(n, 0);
	beta[0] = 1;
	for (size_t i = 1; i < n; ++i) {
		beta[i] = v * beta[i - 1] + eigenPolynomial[i];
	}
	auto CT = CMatrix.transpose();
	Vector<double> eigen(n, 0);
	for (size_t i = 0; i < n; ++i) {
		eigen += CT[i] * beta[i];
	}
	return eigen / eigen.norm();
}

void Krylov::printDeficiency(ostream& out) {
	CMA::Danilevsky dan(n);
	ifstream fIn("../EigenMatrix.txt");
	fIn >> dan;
	dan.run(cout);
	auto temp = eigenPolynomial - dan.getPoly();
	out << "Difference between Krylov and Danilevsky polynomials:\n" << temp << endl;
	out << "Norm:\n" << temp.norm() << endl;
}