#include <iostream>
#include <vector>
#include "../Vector/Vector.h"
#pragma once

using namespace std;

template<class T>
class Matrix {

protected:

	vector<Vector<T>> _values;
	size_t _n;

public:

	Matrix() : _values({}), _n(0) {}
	Matrix(size_t);
	Matrix(size_t, const T&);
	Matrix(vector<Vector<T>>, size_t);
	Matrix(const Matrix&);
	~Matrix();
	Matrix& operator+=(const T&);
	Matrix& operator+=(const Matrix&);
	Matrix& operator-=(const T&);
	Matrix& operator-=(const Matrix&);
	Matrix& operator*=(const T&);
	Matrix& operator*=(Matrix&);
	Matrix operator+(const T&) const;
	Matrix operator+(const Matrix&) const;
	Matrix operator-(const T&) const;
	Matrix operator-(const Matrix&) const;
	Matrix operator*(const T&) const;
	Matrix operator*(Matrix&);
	Vector<T> operator*(Vector<T>&);
	Matrix transpose() const;
	bool operator==(const Matrix&) const;
	bool operator!=(const Matrix&) const;
	Vector<T>& operator[](size_t);
	Matrix& operator=(const Matrix&);

	double norm();

	static Matrix<T> identity(size_t, T);

	template<class T>
	friend ostream& operator<<(ostream&, const Matrix<T>&);

	template<class T>
	friend istream& operator>>(istream&, Matrix<T>&);

};

template<class T>
Matrix<T>::Matrix(size_t n) {
	_n = n;
	for (size_t i = 0; i < _n; ++i) {
		_values.push_back(Vector<T>(_n));
	}
}

template<class T>
Matrix<T>::Matrix(size_t n, const T& element) {
	_n = n;
	for (size_t i = 0; i < _n; ++i) {
		_values.push_back(Vector<T>(_n, element));
	}
}

template<class T>
Matrix<T>::Matrix(vector<Vector<T>> values, size_t n) {
	_n = n;
	for (size_t i = 0; i < _n; ++i) {
		_values[i] = values[i];
	}
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& obj) {
	_n = obj._n;
	_values = vector<Vector<T>>(obj._values);
}

template<class T>
Matrix<T>::~Matrix() {}

template<class T>
Matrix<T>& Matrix<T>::operator+=(const T& c) {
	for (size_t i = 0; i < _n; ++i) {
		_values[i][i] += c;
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& obj) {
	if (_n != obj._n) {
		throw invalid_argument("Matrix addition error. Different matrix sizes.");
	}
	for (size_t i = 0; i < _n; ++i) {
		_values[i] += obj._values[i];
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator-=(const T& c) {
	for (size_t i = 0; i < _n; ++i) {
		_values[i][i] -= c;
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& obj) {
	if (_n != obj._n) {
		throw invalid_argument("Matrix subtraction error. Different matrix sizes.");
	}
	for (size_t i = 0; i < _n; ++i) {
		_values[i] -= obj._values[i];
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator*=(const T& c) {
	for (size_t i = 0; i < _n; ++i) {
		_values[i] *= c;
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator*=(Matrix<T>& obj) {
	if (_n != obj._n) {
		throw invalid_argument("Matrix multiplication error. Different matrix sizes.");
	}
	Matrix<T> temp = obj.transpose();
	Vector<T> row;
	for (size_t i = 0; i < _n; ++i) {
		row = _values[i];
		for (size_t j = 0; j < _n; ++j) {
			_values[i][j] = row.dot(temp[j]);
		}
	}
	return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const T& c) const {
	Matrix<T> temp(*this);
	temp += c;
	return temp;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& obj) const {
	Matrix<T> temp(*this);
	temp += obj;
	return temp;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const T& c) const {
	Matrix<T> temp(*this);
	temp -= c;
	return temp;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& obj) const {
	Matrix<T> temp(*this);
	temp -= obj;
	return temp;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const T& c) const {
	Matrix<T> temp(*this);
	temp *= c;
	return temp;
}

template<class T>
Matrix<T> Matrix<T>::operator*(Matrix<T>& obj) {
	Matrix<T> temp(*this);
	temp *= obj;
	return temp;
}

template<class T>
Vector<T> Matrix<T>::operator*(Vector<T>& vec) {
	if (_n != vec.size()) {
		throw invalid_argument("Matrix multiplication error. Different matrix sizes.");
	}
	Vector<T> temp(_n, 0.0);
	for (size_t i = 0; i < _n; ++i) {
		for (size_t j = 0; j < _n; ++j) {
			temp[i] += _values[i][j] * vec[j];
		}
	}
	return temp;
}

template<class T>
Matrix<T> Matrix<T>::transpose() const {
	Matrix<T> temp(*this);
	for (size_t i = 0; i < _n; ++i) {
		for (size_t j = i + 1; j < _n; ++j) {
			swap(temp._values[i][j], temp._values[j][i]);
		}
	}
	return temp;
}

template<class T>
bool Matrix<T>::operator==(const Matrix<T>& obj) const {
	if (_n != obj._n) {
		return false;
	}
	for (size_t i = 0; i < _n; ++i) {
		if (_values[i] != obj._values[i]) {
			return false;
		}
	}
	return true;
}

template<class T>
bool Matrix<T>::operator!=(const Matrix<T>& obj) const {
	return !(*this == obj);
}

template<class T>
Vector<T>& Matrix<T>::operator[](size_t index) {
	if ((index < 0) || (index >= _n)) {
		throw out_of_range("Matrix row index out of range.");
	}
	return _values[index];
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix& obj) {
	if (this != &obj) {
		_n = obj._n;
		_values.resize(_n);
		for (size_t i = 0; i < _n; ++i) {
			_values[i] = obj._values[i];
		}
	}
	return *this;
}

template<class T>
double Matrix<T>::norm() {
	double sum = 0;
	for (size_t i = 0; i < _n; ++i) {
		for (size_t j = 0; j < _n; ++j) {
			sum += pow(_values[i][j], 2);
		}
	}
	return sqrt(sum);
}

template<class T>
Matrix<T> Matrix<T>::identity(size_t dim, T unit) {
	return Matrix<T>(dim, 0) + unit;
}

template<class T>
ostream& operator<<(ostream& out, const Matrix<T>& obj) {
	out << '[';
	for (size_t i = 0; i < obj._n; ++i) {
		out << obj._values[i];
		if (i != obj._n - 1) {
			out << ", " << endl;
		}
	}
	out << ']' << endl;
	return out;
}

template<class T>
istream& operator>>(istream& in, Matrix<T>& obj) {
	Vector<T> temp;
	while (in >> temp) {
		obj._values.push_back(temp);
	}
	return in;
}