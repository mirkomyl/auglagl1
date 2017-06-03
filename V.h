/*
 *  Created on: May 27, 2013
 *      Author: Mirko Myllykoski (mirko.myllykoski@gmail.com)
 */

#ifndef AUGLAGL1_V_H
#define AUGLAGL1_V_H

#include <stdexcept>
#include "Q.h"

namespace auglagL1 {

// \in V
template <typename T>
class V {
public:
	V() {
		n1 = 0;
		n2 = 0;
		data = 0;
	}
	V(std::size_t _n1, std::size_t _n2, T _h) {
		n1 = 0;
		n2 = 0;
		h = _h;
		data = 0;
		malloc(_n1, _n2);
	}
	V(const V& old) {
		n1 = 0;
		n2 = 0;
		h = old.h;
		data = 0;
		malloc(old.n1, old.n2);

		if(0 < old.n1 && 0 < old.n2) {
			#pragma omp parallel for
			for(std::size_t j = 0; j < old.n2; j++)
				for(std::size_t i = 0; i < old.n1; i++)
					(*this)(i,j) = old(i,j);
		}
	}
	~V() {
		if(data)
			delete [] data;
	}
	V& operator=(const V& a) {
		if(this == &a)
			return *this;

		malloc(a.n1, a.n2);

		if(0 < a.n1 && 0 < a.n2) {
			#pragma omp parallel for
			for(std::size_t j = 0; j < a.n2; j++)
				for(std::size_t i = 0; i < a.n1; i++)
					(*this)(i,j) = a(i,j);
		}

		h = a.h;

		return *this;
	}
	V& operator=(T a) {

	#pragma omp parallel for
	for(std::size_t j = 0; j < n2; j++)
		for(std::size_t i = 0; i < n1; i++)
			(*this)(i,j) = a;

		return *this;
	}
	T& operator()(std::size_t i, std::size_t j) {
#if DEVEL
		if(0 <= i && 0 <= j && i < n1 && j < n2)
			return data[j*n1+i];
		throw std::out_of_range("V::operator() out of range.");
#else
		return data[j*n1+i];
#endif
	}
	T& operator()(std::size_t i, std::size_t j) const {
#if DEVEL
		if(0 <= i && 0 <= j && i < n1 && j < n2)
			return data[j*n1+i];
		throw std::out_of_range("V::operator() out of range.");
#else
		return data[j*n1+i];
#endif
	}
	V<T>& operator+=(const V<T>& a) {
		#pragma omp parallel for
		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1; i++)
				(*this)(i,j) += a(i,j);
		return *this;
	}

	V<T>& operator-=(const V<T>& a) {
		#pragma omp parallel for
		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1; i++)
				(*this)(i,j) -= a(i,j);
		return *this;
	}

	template <typename Y>
	V<T>& operator*=(Y coef) {
		#pragma omp parallel for
		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1; i++)
				(*this)(i,j) *= coef;
		return *this;
	}

	template <typename Y>
	V<T>& operator/=(Y coef) {
		#pragma omp parallel for
		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1; i++)
				(*this)(i,j) /= coef;
		return *this;
	}

	// (\grad u(\Omega^*_{i,j}))
	R2<T> grad(std::size_t i, std::size_t j) const {

		if(i & 0x1) {
			// Odd numbered triangle
			return (VER_ODD_UL(*this,i,j)*R2<T>(-1, 1) +
					VER_ODD_UR(*this,i,j)*R2<T>(1,0) +
					VER_ODD_LL(*this,i,j)*R2<T>(0,-1))/h;
		} else {
			// Even numbered triangle
			return (VER_EVEN_UR(*this,i,j)*R2<T>(0,1) +
					VER_EVEN_LL(*this,i,j)*R2<T>(-1,0) +
					VER_EVEN_LR(*this,i,j)*R2<T>(1,-1))/h;
		}
	}

	T getH() const { return h; }
	std::size_t getN1() const { return n1; }
	std::size_t getN2() const { return n2; }
private:
	void malloc(std::size_t _n1, std::size_t _n2) {
		if(2 <= _n1 && 2 <= _n2) {
			if(n1*n2 != _n1*_n2) {
				if(data)
					delete [] data;

				try {
					data = new T[_n1*_n2];
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
			this->data = 0;
			this->n1 = 0;
			this->n2 = 0;
		}
	}
	T* data;
	T h;
	std::size_t n1;
	std::size_t n2;
};

template <typename T, typename Y>
V<T> operator*(Y coef, const V<T>& d) {
	V<T> tmp(d);
	std::size_t n1 = tmp.getN1();
	std::size_t n2 = tmp.getN2();
	#pragma omp parallel for
	for(std::size_t j = 0; j < n2; j++)
		for(std::size_t i = 0; i < n1; i++)
			tmp(i,j) *= coef;
	return tmp;
}

template <typename T, typename Y>
V<T> operator*(const V<T>& d, Y coef) {
	V<T> tmp(d);
	std::size_t n1 = tmp.getN1();
	std::size_t n2 = tmp.getN2();
	#pragma omp parallel for
	for(std::size_t j = 0; j < n2; j++)
		for(std::size_t i = 0; i < n1; i++)
			tmp(i,j) *= coef;
	return tmp;
}


template <typename T, typename Y>
V<T> operator/(const V<T>& d, Y coef) {
	V<T> tmp(d);
	std::size_t n1 = tmp.getN1();
	std::size_t n2 = tmp.getN2();
	#pragma omp parallel for
	for(std::size_t j = 0; j < n2; j++)
		for(std::size_t i = 0; i < n1; i++)
			tmp(i,j) /= coef;
	return tmp;
}

template <typename T>
V<T> operator+(const Q<T>& a, const V<T>& b) {
	V<T> tmp(a);
	std::size_t n1 = tmp.getN1();
	std::size_t n2 = tmp.getN2();
	#pragma omp parallel for
	for(std::size_t j = 0; j < n2; j++)
		for(std::size_t i = 0; i < n1; i++)
			tmp(i,j) += b(i,j);
	return tmp;
}

template <typename T>
V<T> operator-(const V<T>& a, const V<T>& b) {
	V<T> tmp = a;
	std::size_t n1 = tmp.getN1();
	std::size_t n2 = tmp.getN2();
	#pragma omp parallel for
	for(std::size_t j = 0; j < n2; j++)
		for(std::size_t i = 0; i < n1; i++)
			tmp(i,j) -= b(i,j);
	return tmp;
}

template <typename T>
T dot(const V<T>& a, const V<T>& b) {
	std::size_t n1 = a.getN1();
	std::size_t n2 = a.getN2();
	T n = 0;
	for(std::size_t j = 0; j < n2; j++)
		for(std::size_t i = 0; i < n1; i++)
			n += a(i,j) * b(i,j);
	return n;
}

template <typename T>
T dot(const V<T>& a) {
	return dot(a, a);
}

template <typename T>
T vectorNorm2(const V<T>& a) {
	return sqrt(dot(a));
}

template <typename T>
T vectorNormInf(const V<T>& a) {
	std::size_t n1 = a.getN1();
	std::size_t n2 = a.getN2();

	T max = 0.0;

	for(std::size_t j = 0; j < n2; j++) {
		for(std::size_t i = 0; i < n1; i++) {
			max = std::max(max, ::fabs(a(i,j)));
		}
	}

	return max;
}

template <typename T>
void check(const V<T>& x, int *nanCounter, int *infCounter) {
	*nanCounter = 0;
	*infCounter = 0;

	std::size_t n1 = x.getN1();
	std::size_t n2 = x.getN2();

	for(std::size_t j = 0; j < n2; j++) {
		for(std::size_t i = 0; i < n1; i++) {
			if(isNaN(x(i,j))) (*nanCounter)++;
			if(isInf(x(i,j))) (*infCounter)++;
		}
	}
}

template <typename T>
T avg(const V<T>& x) {
	T avg = 0;

	std::size_t n1 = x.getN1();
	std::size_t n2 = x.getN2();

	for(std::size_t j = 0; j < n2; j++) {
		T tmp = 0.0;
		for(std::size_t i = 0; i < n1; i++) {
			if(isNaN(x(i,j)))
				return NaN;
			tmp += x(i,j);
		}
		avg += tmp/n1;
	}

	return avg/n2;
}

template <typename T>
V<T> fabs(const V<T>& x) {
	std::size_t n1 = x.getN1();
	std::size_t n2 = x.getN2();
	T h = x.getH();

	V<T> res(n1, n2, h);

	for(std::size_t j = 0; j < n2; j++)
		for(std::size_t i = 0; i < n1; i++)
			res(i,j) = ::fabs(x(i,j));

	return res;
}

}

#endif
