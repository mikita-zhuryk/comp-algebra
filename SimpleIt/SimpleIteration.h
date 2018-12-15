#pragma once
#include "../Method/Method.h"
#include <fstream>

namespace CMA {

	class SimpleIteration : public Method {

		Matrix<double> B;
		Vector<double> g;

	public:

		SimpleIteration(size_t);

		void solve(ostream&) override;
		void buildBg();
		double determinant() override;

	};

}