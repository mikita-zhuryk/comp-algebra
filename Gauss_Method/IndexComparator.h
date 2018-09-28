#pragma once
#include <vector>

using namespace std;

class IndexComparator {

public:

	bool operator()(pair<size_t, double>, pair<size_t, double>);

};