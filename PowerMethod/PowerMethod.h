#pragma once
#include "../EigenMethod/EigenMethod.h"

namespace CMA {

	class PowerMethod : public EigenMethod {

		Vector<double> lambda;
		Vector<double> x;
		double accuracy;

	public:

		PowerMethod(size_t, int);

	private:

		void find(ostream&) override;

		Vector<double> findLambda(Vector<double>, Vector<double>);

	};

}