#pragma once
#include "../Method/Method.h"

namespace CMA {

	class GaussSeidel : public Method {

	public:

		GaussSeidel(size_t);
		void solve(ostream&) override;
		double determinant() override;

	};

}