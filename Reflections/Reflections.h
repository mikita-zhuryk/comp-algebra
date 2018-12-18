#pragma once
#include "../Method/Method.h"
#include <fstream>

namespace CMA {

	class Reflections : public Method {

		Vector<double> omega;
		Vector<double> y;
		Vector<double> e;

	public:

		Reflections(size_t);
		Reflections(Matrix<double>, Vector<double>);

	private:

		void solve(ostream&) override;

		void calcOmega();
		void pickY(size_t);
		void getSolution();
		double determinant() override;

	};

}