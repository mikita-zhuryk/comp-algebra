#include "../Method/Method.h"
#include <fstream>
#pragma once

using namespace std;

namespace CMA {

	class PassMethod : public Method {

		vector<double> upper;
		vector<double> diag;
		vector<double> lower;
		vector<double> alpha;
		vector<double> beta;
		vector<double> ksi;
		vector<double> eta;
		vector<double> delta;
		size_t mid;
		bool transformed;

	public:

		PassMethod(size_t, bool = true);

	private:

		void buildTriDiag();
		void fillVectors();

		void solve(ostream&) override;

		void calcRight();
		void calcLeft();
		void printCoeffs(ostream&);
		void calcX();
		void printDelta(ostream& out);
		double determinant() override;

	};

}