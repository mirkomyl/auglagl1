/*
 *  Created on: May 27, 2013
 *      Author: Mirko Myllykoski (mirko.myllykoski@gmail.com)
 */

#ifndef AUGLAGL1_COMMON_H
#define AUGLAGL1_COMMON_H

#include <limits>
#include <sstream>
#include <exception>
#include <sys/time.h>

#define POW4(a) (1<<(2*(a)))
#define LOG4(a) (log(a)/log(4))

#define INF (1.0/0.0)
#define NaN (0.0/0.0)

#define VER_ODD_UL(u,i,j) (u)(((i)-1)/2,j)
#define VER_ODD_UR(u,i,j) (u)(((i)-1)/2+1,j)
#define VER_ODD_LL(u,i,j) (u)(((i)-1)/2,j-1)

#define VER_EVEN_UR(u,i,j) (u)((i)/2,j)
#define VER_EVEN_LL(u,i,j) (u)((i)/2-1,j-1)
#define VER_EVEN_LR(u,i,j) (u)((i)/2,j-1)

namespace auglagL1 {

template <typename T>
inline bool isNaN(T value) {
	return value != value;
}

template <typename T>
inline bool isInf(T value) {
	return std::numeric_limits<T>::has_infinity &&
			value == std::numeric_limits<T>::infinity();
}

template <typename T>
inline T pow2(T x) { return x*x; }

template <typename T>
inline T pow3(T x) { return x*x*x; }

template <typename T>
inline T pow5(T x) { return x*x*x*x*x; }

template <typename T>
inline T pow3_2(T x) { return sqrt(pow3(x)); }

template <typename T>
inline T pow5_2(T x) { return sqrt(pow5(x)); }

template <typename T>
inline T sgn(T x) { return x < 0 ? -1 : 1; }

// Object to string functions and helpers

template <typename T>
class ToStringHelper {
public:
	static std::string to(const T& value) {
		std::stringstream ss;
		ss << value;
		return ss.str();
	}
};

template <typename T>
std::string toString(const T& value) {
	return ToStringHelper<T>::to(value);
}

class Timer {
public:
	Timer() {
		running = false;
		ready = false;
	}
	void begin() {
		gettimeofday(&begin_time, NULL);
		running = true;
		ready = false;
	}
	void end() {
		if(!running)
			throw std::exception();
		gettimeofday(&end_time, NULL);
		running = false;
		ready = true;
	}
	double getTime() const {
		return getEndTime() - getBeginTime();
	}
	double getBeginTime() const {
		if(!running && !ready)
			throw std::exception();
		return begin_time.tv_sec + begin_time.tv_usec*1.0E-6;
	}
	double getEndTime() const {
		if(!ready)
			throw std::exception();
		return end_time.tv_sec + end_time.tv_usec*1.0E-6;
	}
private:
	bool running ;
	bool ready;
	struct timeval begin_time;
	struct timeval end_time;
};

}

#endif
