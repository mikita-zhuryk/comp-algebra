#include "../Method/Method.h"
#include <fstream>
#pragma once

using namespace std;

typedef unsigned int UINT;

namespace CMA {

	class GaussMethod : public Method {

		GaussMethod();
		UINT swapCount;
		vector<double> pivot;

	public:

		GaussMethod(size_t);
		~GaussMethod();
		void solve(ostream&) override;

	private:

		void swapRows(size_t, size_t);
		size_t findMaxInRows(UINT);
		void makeUpperTriangular();
		void printUpperTriangular(ostream&);
		void getSolution(Vector<double>&);
		double determinant() override;
		void findInverse() override;

	};

}