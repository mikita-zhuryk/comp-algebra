#include "../Method/Method.h"
#include <fstream>
#pragma once

using namespace std;

class Sqrt_Method: public Method {



public:

	Sqrt_Method(size_t);

	void solve(ostream&) override;

};