#pragma once
#include "../Method/Method.h"
#include <fstream>

namespace CMA {

	class SimpleIteration : public Method {

		Matrix<double> B;
		Vector<double> g;
		const double accuracy;

	public:

		SimpleIteration(size_t, int);

	private:

		void solve(ostream&) override;
		void buildBg();
		double determinant() override;

	};

}