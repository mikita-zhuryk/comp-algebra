#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#define OUTPUT_WIDTH 8

using namespace std;

template<class T>
class Matrix;

template<class T>
class Vector {

protected:
	
	vector<T> _vec;
	size_t _size;

public:

	Vector() : _vec({}), _size(0) {}
	Vector(size_t);
	Vector(size_t, const T&);
	Vector(const vector<T>&);
	Vector(T*, size_t);
	Vector(const Vector&);
	~Vector();
	Vector& operator+=(const T&);
	Vector& operator+=(const Vector&);
	Vector& operator-=(const T&);
	Vector& operator-=(const Vector&);
	Vector& operator*=(const T&);
	Vector& operator*=(const Matrix<T>&);
	Vector& operator/=(const T&);
	Vector operator+(const T&) const;
	Vector operator+(const Vector&) const;
	Vector operator-(const T&) const;
	Vector operator-(const Vector&) const;
	Vector operator*(const T&) const;
	Vector operator*(const Matrix<T>&) const;
	Vector operator/(const T&) const;
	T dot(const Vector&) const;
	Matrix<T> matrixMult(const Vector&) const;
	bool operator==(const Vector&) const;
	bool operator!=(const Vector&) const;
	T& operator[](size_t);
	Vector& operator=(const Vector&);
	size_t size();

	double norm() const;
	double firstNorm() const;

	template<class T>
	friend ostream& operator<<(ostream&, const Vector<T>&);

	template<class T>
	friend istream& operator>>(istream&, Vector<T>&);

};

template<class T>
Vector<T>::Vector(size_t size) {
	_size = size;
	_vec.resize(_size);
}

template<class T>
Vector<T>::Vector(size_t size, const T& obj) {
	_size = size;
	_vec.resize(_size);
	for (size_t i = 0; i < _size; i++) {
		_vec[i] = obj;
	}
}

template<class T>
Vector<T>::Vector(const vector<T>& vec) {
	_size = vec.size();
	_vec.resize(_size);
	for (size_t i = 0; i < _size; ++i) {
		_vec[i] = vec[i];
	}
}

template<class T>
Vector<T>::Vector(T* arr, size_t size) {
	_size = size;
	_vec.resize(_size);
	for (size_t i = 0; i < _size; i++) {
		_vec[i] = arr[i];
	}
}

template<class T>
Vector<T>::Vector(const Vector<T>& obj) {
	_size = obj._size;
	_vec.resize(_size);
	for (size_t i = 0; i < _size; i++) {
		_vec[i] = obj._vec[i];
	}
}

template<class T>
Vector<T>::~Vector() {}

template<class T>
Vector<T>& Vector<T>::operator+=(const T& c) {
	for (size_t i = 0; i < _size; i++) {
		_vec[i] += c;
	}
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& obj) {
	if (obj._size != _size) {
		throw invalid_argument("Vector addition error. Different vector sizes.");
	}
	for (size_t i = 0; i < _size; i++) {
		_vec[i] += obj._vec[i];
	}
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator-=(const T& c) {
	for (size_t i = 0; i < _size; i++) {
		_vec[i] -= c;
	}
	return *this;

}

template<class T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& obj) {
	if (obj._size != _size) {
		throw invalid_argument("Vector subtraction error. Different vector sizes.");
	}
	for (size_t i = 0; i < _size; i++) {
		_vec[i] -= obj._vec[i];
	}
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator*=(const T& c) {
	for (size_t i = 0; i < _size; ++i) {
		_vec[i] *= c;
	}
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator*=(const Matrix<T>& m) {
	auto temp = m.transpose();
	for (size_t i = 0; i < _size; ++i) {
		_vec[i] = (*this).dot(temp[i]);
	}
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator/=(const T& c) {
	for (size_t i = 0; i < _size; ++i) {
		_vec[i] /= c;
	}
	return *this;
}

template<class T>
Vector<T> Vector<T>::operator+(const T& c) const {
	Vector<T> temp(*this);
	temp += c;
	return temp;
}

template<class T>
Vector<T> Vector<T>::operator+(const Vector<T>& obj) const {
	Vector<T> temp(*this);
	temp += obj;
	return temp;
}

template<class T>
Vector<T> Vector<T>::operator-(const T& c) const {
	Vector<T> temp(*this);
	temp -= c;
	return temp;
}


template<class T>
Vector<T> Vector<T>::operator-(const Vector<T>& obj) const {
	Vector<T> temp(*this);
	temp -= obj;
	return temp;
}

template<class T>
Vector<T> Vector<T>::operator*(const T& c) const {
	Vector<T> temp(*this);
	temp *= c;
	return temp;
}

template<class T>
Vector<T> Vector<T>::operator*(const Matrix<T>& m) const {
	Vector<T> temp(*this);
	temp *= m;
	return temp;
}

template<class T>
Vector<T> Vector<T>::operator/(const T& c) const {
	Vector<T> temp(*this);
	temp /= c;
	return temp;
}

template<class T>
T Vector<T>::dot(const Vector<T>& obj) const {
	if (obj._size != _size) {
		throw invalid_argument("Vector multiplication error. Different vector sizes.");
	}
	T sum = 0;
	for (size_t i = 0; i < _size; ++i) {
		sum += _vec[i] * obj._vec[i];
	}
	return sum;
}

template<class T>
Matrix<T> Vector<T>::matrixMult(const Vector<T>& obj) const {
	if (obj._size != _size) {
		throw invalid_argument("Vector multiplication error. Different vector sizes.");
	}
	Matrix<T> result(_size);
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			result[i][j] = _vec[i] * obj._vec[j];
		}
	}
	return result;
}

template<class T>
bool Vector<T>::operator==(const Vector<T>& obj) const {
	if (_size != obj._size) {
		return false;
	}
	for (size_t i = 0; i < _size; i++) {
		if (_vec[i] != obj._vec[i]) {
			return false;
		}
	}
	return true;
}

template<class T>
bool Vector<T>::operator!=(const Vector<T>& obj) const {
	return !(*this == obj);
}

template<class T>
T& Vector<T>::operator[](size_t index) {
	if ((index < 0) || (index >= _size)) {
		throw out_of_range("Vector index is out of range");
	}
	return _vec[index];
}

template<class T>
Vector<T>& Vector<T>::operator=(const Vector<T>& obj) {
	if (this != &obj) {
		_size = obj._size;
		_vec.resize(_size);
		for (size_t i = 0; i < _size; i++) {
			_vec[i] = obj._vec[i];
		}
	}
	return *this;
}

template<class T>
size_t Vector<T>::size() {
	return _size;
}

template<class T>
double Vector<T>::norm() const {
	double sum = 0;
	for (size_t i = 0; i < _size; ++i) {
		sum += pow(_vec[i], 2);
	}
	return sqrt(sum);
}

template<class T>
double Vector<T>::firstNorm() const {
	double max = 0;
	for (size_t i = 0; i < _size; ++i) {
		if (_vec[i] > max) {
			max = _vec[i];
		}
	}
	return max;
}

template<class T>
ostream& operator<<(ostream& out, const Vector<T>& obj) {
	out << '[';
	for (size_t i = 0; i < obj._size; ++i) {
		out << setw(OUTPUT_WIDTH) << obj._vec[i];
		if (i != obj._size - 1) {
			out << ", ";
		}
	}
	out << ']';
	return out;
}

template<class T>
istream& operator>>(istream& in, Vector<T>& obj) {
	T temp;
	while (in >> temp) {
		obj._vec.push_back(temp);
	}
	return in;
}
