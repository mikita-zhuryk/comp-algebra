#pragma once
#include "../Method/Method.h"

namespace CMA {

	class Relaxation : public Method {

		const double accuracy;
		const double omega;

	public:

		Relaxation(size_t, double, int);

	private:

		void solve(ostream&) override;
		double determinant() override;

	};

}