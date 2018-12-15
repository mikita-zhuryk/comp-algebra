#pragma once
#include "../Method/Method.h"

namespace CMA {

	class GaussSeidel : public Method {

		const double accuracy;

	public:

		GaussSeidel(size_t, int);

	private:

		void solve(ostream&) override;
		double determinant() override;

	};

}