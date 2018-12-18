#pragma once
#include "../EigenMethod/EigenMethod.h"

namespace CMA {

	class Krylov : public EigenMethod {

		uint16_t dimA;
		Vector<double> C;
		Matrix<double> CMatrix;

	public:

		Krylov(size_t);

	private:

		void find(ostream&) override;
		void buildPolynomial();
		void findEigenVectors();
		Vector<double> findEigenVec(double);
		void printDeficiency(ostream&) override;

	};

}