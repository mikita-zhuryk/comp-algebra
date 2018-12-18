#pragma once
#include "../EigenMethod/EigenMethod.h"

namespace CMA {

	class IterJacobi : public EigenMethod {

		double accuracy;
		Vector<double> calculatedEigenValues;
		Matrix<double> T;

	public:

		IterJacobi(size_t, int);

	private:

		pair<size_t, size_t> findMax();
		void buildDiag();

		void find(ostream&) override;

	};

}