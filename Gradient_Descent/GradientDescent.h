#pragma once
#include "../Method/Method.h"

namespace CMA {

	class GradientDescent : public Method {

		const double accuracy;

	public:

		GradientDescent(size_t, int);

	private:

		void solve(ostream&) override;
		double determinant() override;

	};

}