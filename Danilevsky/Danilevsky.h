#pragma once
#include "../EigenMethod/EigenMethod.h"

namespace CMA {

	class Danilevsky : public EigenMethod {

		Matrix<double> S;
		Matrix<double> invS;

	public:
		
		Danilevsky(size_t);

	private:

		void solve(ostream&) override;
		void Frobenius();
		void collectPoly();
		void findEigenVectors();
		Vector<double> buildY(double);

	};

}