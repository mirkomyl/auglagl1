/*
 *  Created on: May 27, 2013
 *      Author: Mirko Myllykoski (mirko.myllykoski@gmail.com)
 */

#ifndef AUGLAGL1_Q_H
#define AUGLAGL1_Q_H

#include <stdexcept>
#include "common.h"
#include "R2.h"

namespace auglagL1 {

template <typename T>
class Q {
public:
	Q() {
		n1 = 0;
		n2 = 0;
		data = 0;
	}
	Q(std::size_t _n1, std::size_t _n2, T _h) {
		n1 = 0;
		n2 = 0;
		h = _h;
		data = 0;
		alloc(_n1, _n2);
	}
	Q(const Q<T>& old) {
		n1 = 0;
		n2 = 0;
		h = old.h;
		data = 0;
		alloc(old.n1, old.n2);

		#pragma omp parallel for
		for(std::size_t j = 1; j <= n2; j++)
			for(std::size_t i = 1; i <= n1; i++)
				(*this)(i,j) = old(i,j);

	}
	~Q() {
		if(data)
			delete [] data;
	}
	Q<T>& operator=(const Q<T>& a) {
		if(this == &a)
			return *this;

		alloc(a.n1, a.n2);

		#pragma omp parallel for
		for(std::size_t j = 1; j <= n2; j++)
			for(std::size_t i = 1; i <= n1; i++)
				(*this)(i,j) = a(i,j);

		h = a.h;

		return *this;
	}
	Q<T>& operator=(T a) {
		#pragma omp parallel for
		for(std::size_t j = 1; j <= n2; j++)
			for(std::size_t i = 1; i <= n1; i++)
				(*this)(i,j) = a;
		return *this;
	}
	R2<T>& operator()(std::size_t i, std::size_t j) {
#if DEVEL
		if(0 < i && 0 < j && i <= n1 && j <= n2)
			return data[(j-1)*n1+i-1];
		throw std::out_of_range("Q::operator() out of range.");
#else
		return data[(j-1)*n1+i-1];
#endif

	}
	const R2<T>& operator()(std::size_t i, std::size_t j) const {
#if DEVEL
		if(0 < i && 0 < j && i <= n1 && j <= n2)
			return data[(j-1)*n1+i-1];
		throw std::out_of_range("Q::operator() out of range.");
#else
		return data[(j-1)*n1+i-1];
#endif
	}

	Q<T>& operator+=(const Q<T>& a) {
		#pragma omp parallel for
		for(std::size_t j = 1; j <= n2; j++)
			for(std::size_t i = 1; i <= n1; i++)
				(*this)(i,j) += a(i,j);
		return *this;
	}

	Q<T>& operator-=(const Q<T>& a) {
		#pragma omp parallel for
		for(std::size_t j = 1; j <= n2; j++)
			for(std::size_t i = 1; i <= n1; i++)
				(*this)(i,j) -= a(i,j);
		return *this;
	}

	template <typename Y>
	Q<T>& operator*=(Y coef) {
		#pragma omp parallel for
		for(std::size_t j = 1; j <= n2; j++)
			for(std::size_t i = 1; i <= n1; i++)
				(*this)(i,j) *= coef;
		return *this;
	}

	template <typename Y>
	Q<T>& operator/=(Y coef) {
		#pragma omp parallel for
		for(std::size_t j = 1; j <= n2; j++)
			for(std::size_t i = 1; i <= n1; i++)
				(*this)(i,j) /= coef;
		return *this;
	}

	// (div q (P_{i,j}))
	T div(std::size_t i, std::size_t j) const {

#if DEVEL
		if(i < 1 || j < 1 || n1/2 <= i || n2 <= j)
			throw std::out_of_range("V::div() out of range.");
#endif

		T tmp = 0;

		tmp += dot((*this)(2*i,j+1), R2<T>(1.0, -1.0));
		tmp += dot((*this)(2*i+1,j+1), R2<T>(0.0, -1.0));
		tmp += dot((*this)(2*i+2,j+1), R2<T>(-1.0, 0.0));

		tmp += dot((*this)(2*i-1,j), R2<T>(1.0, 0.0));
		tmp += dot((*this)(2*i,j), R2<T>(0.0, 1.0));
		tmp += dot((*this)(2*i+1,j), R2<T>(-1.0, 1.0));

		return - tmp / (2.0 * h);
	}
	T getH() const { return h; }
	std::size_t getN1() const { return n1; }
	std::size_t getN2() const { return n2; }
private:
	void alloc(std::size_t _n1, std::size_t _n2) {
		if(0 < _n1 && 0 < _n2) {
			if(n1*n2 != _n1*_n2) {
				if(data)
					delete [] data;

				try {
					data = new R2<T>[_n1*_n2];
				} catch (...) {
					data = 0;
					n1 = 0;
					n2 = 0;
					throw;
				}
			}
			n1 = _n1;
			n2 = _n2;
		} else {
			if(data)
				delete [] data;
			data = 0;
			n1 = 0;
			n2 = 0;
		}
	}
	std::size_t n1;
	std::size_t n2;
	T h;
	R2<T> *data;
};

template <typename T, typename Y>
Q<T> operator*(Y coef, const Q<T>& d) {
	Q<T> tmp(d);
	std::size_t n1 = tmp.getN1();
	std::size_t n2 = tmp.getN2();
	#pragma omp parallel for
	for(std::size_t j = 1; j <= n2; j++)
		for(std::size_t i = 1; i <= n1; i++)
			tmp(i,j) *= coef;
	return tmp;
}

template <typename T, typename Y>
Q<T> operator/(const Q<T>& d, Y coef) {
	Q<T> tmp(d);
	std::size_t n1 = tmp.getN1();
	std::size_t n2 = tmp.getN2();
	#pragma omp parallel for
	for(std::size_t j = 1; j <= n2; j++)
		for(std::size_t i = 1; i <= n1; i++)
			tmp(i,j) /= coef;
	return tmp;
}

template <typename T>
Q<T> operator+(const Q<T>& a, const Q<T>& b) {
	Q<T> tmp(a);
	std::size_t n1 = tmp.getN1();
	std::size_t n2 = tmp.getN2();
	#pragma omp parallel for
	for(std::size_t j = 1; j <= n2; j++)
		for(std::size_t i = 1; i <= n1; i++)
			tmp(i,j) += b(i,j);
	return tmp;
}

template <typename T>
Q<T> operator-(const Q<T>& a, const Q<T>& b) {
	Q<T> tmp(a);
	std::size_t n1 = tmp.getN1();
	std::size_t n2 = tmp.getN2();
	#pragma omp parallel for
	for(std::size_t j = 1; j <= n2; j++)
		for(std::size_t i = 1; i <= n1; i++)
			tmp(i,j) -= b(i,j);
	return tmp;
}

// Dot-product for the vector representation
template <typename T>
T dot(const Q<T>& a, const Q<T>& b) {
	std::size_t n1 = a.getN1();
	std::size_t n2 = a.getN2();

	T tmp = 0;
	for(std::size_t j = 1; j <= n2; j++) {
		for(std::size_t i = 1; i <= n1; i++) {
			tmp += a(i,j).x1 * b(i,j).x1;
			tmp += a(i,j).x2 * b(i,j).x2;
		}
	}
	return tmp;
}

// Norm for the vector representation
template <typename T>
T vectorNorm2(const Q<T>& x) {
	return sqrt(dot(x, x));
}

template <typename T>
T vectorNormInf(const Q<T>& x) {
	std::size_t n1 = x.getN1();
	std::size_t n2 = x.getN2();

	T max = 0.0;

	for(std::size_t j = 1; j <= n2; j++)
		for(std::size_t i = 1; i <= n1; i++)
			max = std::max(max, ::fabs(x(i,j)));

	return max;
}

template <typename T>
void check(const Q<T>& x, int *nanCounter, int *infCounter) {
	*nanCounter = 0;
	*infCounter = 0;

	std::size_t n1 = x.getN1();
	std::size_t n2 = x.getN2();

	for(std::size_t j = 1; j <= n2; j++) {
		for(std::size_t i = 1; i <= n1; i++) {
			if(isNaN(x(i,j))) (*nanCounter)++;
			if(isInf(x(i,j))) (*infCounter)++;
		}
	}
}

template <typename T>
T avg1(const Q<T>& x) {
	T avg = 0;

	std::size_t n1 = x.getN1();
	std::size_t n2 = x.getN2();

	for(std::size_t j = 1; j <= n2; j++) {
		T tmp = 0;
		for(std::size_t i = 1; i <= n1; i++) {
			if(isNaN(x(i,j).x1))
				return NaN;
			tmp += x(i,j).x1;
		}
		avg += tmp/n1;
	}

	return avg/n2;
}

template <typename T>
T avg2(const Q<T>& x) {
	T avg = 0;

	std::size_t n1 = x.getN1();
	std::size_t n2 = x.getN2();

	for(std::size_t j = 1; j <= n2; j++) {
		T tmp = 0;
		for(std::size_t i = 1; i <= n1; i++) {
			if(isNaN(x(i,j).x2))
				return NaN;
			tmp += x(i,j).x2;
		}
		avg += tmp/n1;
	}

	return avg/n2;
}

template <typename T>
R2<T> avg(const Q<T>& x) {
	R2<T> avg = 0;

	std::size_t n1 = x.getN1();
	std::size_t n2 = x.getN2();

	for(std::size_t j = 1; j <= n2; j++) {
		R2<T> tmp = 0;
		for(std::size_t i = 1; i <= n1; i++) {
			if(isNaN(x(i,j).x1) || isNaN(x(i,j).x2))
				return NaN;
			tmp += x(i,j);
		}
		avg += tmp/(1.0*n1);
	}

	return avg/(1.0*n2);
}

}

#endif
