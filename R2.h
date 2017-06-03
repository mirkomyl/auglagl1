/*
 *  Created on: May 27, 2013
 *      Author: Mirko Myllykoski (mirko.myllykoski@gmail.com)
 */

#ifndef AUGLAGL1_R2_H
#define AUGLAGL1_R2_H

#include "common.h"

namespace auglagL1 {

// R^2-vector
template <typename T>
struct R2 {
	R2() {
		x1 = x2 = 0.0;
	}
	R2(T in) {
		x1 = in;
		x2 = in;
	}
	R2(T in1, T in2) {
		x1 = in1;
		x2 = in2;
	}
	R2& operator=(T in) {
		x1 = in;
		x2 = in;
		return *this;
	}
	R2& operator+=(const R2<T>& a) {
		x1 += a.x1;
		x2 += a.x2;
		return *this;
	}
	R2& operator-=(const R2<T>& a) {
		x1 -= a.x1;
		x2 -= a.x2;
		return *this;
	}
	template <typename Y>
	R2& operator*=(Y coef) {
		x1 *= coef;
		x2 *= coef;
		return *this;
	}
	template <typename Y>
	R2& operator/=(Y coef) {
		x1 /= coef;
		x2 /= coef;
		return *this;
	}
	bool operator==(const R2<T>& a) const {
		return x1 == a.x1 && x2 == a.x2;
	}
	bool operator!=(const R2<T>& a) const {
		return x1 != a.x1 || x2 != a.x2;
	}
	T x1;
	T x2;
};

template <typename T>
inline R2<T> operator*(T coef, const R2<T>& a) {
	return R2<T>(coef*a.x1, coef*a.x2);
}

template <typename T>
inline R2<T> operator*(const R2<T>& a, T coef) {
	return coef*a;
}

template <typename T>
inline R2<T> operator/(const R2<T>& a, T coef) {
	return R2<T>(a.x1/coef, a.x2/coef);
}

template <typename T>
inline R2<T> operator+(const R2<T>& a, const R2<T>& b) {
	return R2<T>(a.x1+b.x1, a.x2+b.x2);
}

template <typename T>
inline R2<T> operator-(const R2<T>& a, const R2<T>& b) {
	return R2<T>(a.x1-b.x1, a.x2-b.x2);
}

template <typename T>
inline T dot(const R2<T>& a, const R2<T>& b) {
	return a.x1 * b.x1 + a.x2 * b.x2;
}

template <typename T>
inline T dot(const R2<T>& a) {
	return dot(a,a);
}

template <typename T>
inline T norm(const R2<T>& a) {
	return sqrt(dot(a));
}

template<typename T>
inline bool isNaN(const R2<T>& value) {
	return isNaN(value.x1) || isNaN(value.x2);
}

template<typename T>
inline bool isInf(const R2<T>& value) {
	return isInf(value.x1) || isInf(value.x2);
}

template <typename T>
class ToStringHelper<R2<T> > {
public:
	static std::string to(const R2<T>& value) {
		return "(" + toString(value.x1) + ", " + toString(value.x2) + ")";
	}
};

}

#endif
