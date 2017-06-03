/*
 *  Created on: May 27, 2013
 *      Author: Mirko Myllykoski (mirko.myllykoski@gmail.com)
 */

#ifndef IMAGE_H_
#define IMAGE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <stdexcept>
#include <math.h>
#include <limits>
#include <iostream>
#include <algorithm>

#include "io_png.h"

#define RANGE(x,a,b) std::max(std::min(x,b),a)

template <typename T>
bool isNaN(T value) {
	return value != value;
}

template <typename T>
bool isInf(T value) {
	return std::numeric_limits<T>::has_infinity &&
			value == std::numeric_limits<T>::infinity();
}

template <typename T>
T pow2(T x) { return x*x; }

template <typename T>
class Image {
public:
	explicit Image(const std::string& _fileName, T _low = 0.0, T _high = 1.0) {
		unsigned char *tmp  = io_png_read_u8_gray(_fileName.c_str(), &width, &height);

		if(tmp == NULL)
			throw std::exception();

		data = new T[width*height];

	    for(std::size_t j = 0; j < height; j++)
	    	for(std::size_t i = 0; i < width; i++)
	    		(*this)(i,j) = ((_high-_low)/255.0) * tmp[j*width+i] + _low;

	    delete [] tmp;

	    low = _low;
	    high = _high;
	}

	Image(std::size_t _width, std::size_t _height, T _low = 0.0, T _high = 1.0) {
		data = new T[_width*_height];
		width = _width;
		height = _height;
		low = _low;
		high = _high;
	}

	Image(const Image& _old) {
		height = _old.getHeight();
		width = _old.getWidth();

		data = new T[getHeight()*getWidth()];

		for(std::size_t j = 0; j < getHeight(); j++)
			for(std::size_t i = 0; i < getWidth(); i++)
				(*this)(i,j) = _old(i,j);

		low = _old.getLow();
		high = _old.getHigh();
	}
	~Image() {
		delete [] data;
	}
	Image& operator=(const Image& _a) {
		if(this == &_a)
			return *this;

		if(getHeight()*getWidth() != _a.getHeight()*_a.getWidth()) {
			delete [] data;
			data = new T[_a.getHeight()*_a.getWidth()];
		}

		width = _a.getWidth();
		height = _a.getHeight();

		for(std::size_t j = 0; j < getHeight(); j++)
			for(std::size_t i = 0; i < getWidth(); i++)
				(*this)(i,j) = _a(i,j);

		low = _a.getLow();
		high = _a.getHigh();

		return *this;
	}

	void save(const std::string& _fileName, T _eps = 0.0) const {
		save(_fileName, getLow(), getHigh(), _eps);
	}

	void save(const std::string& _fileName, T _low, T _high, T _eps = 0.0) const {

		unsigned char *tmp = new unsigned char[3*width*height];
		unsigned char *r = tmp;
		unsigned char *g = tmp +   width*height;
		unsigned char *b = tmp + 2*width*height;

	    T max = -1.0/0.0;
	    T min = 1.0/0.0;

	    for(std::size_t j = 0; j < getHeight(); j++) {
	    	for(std::size_t i = 0; i < getWidth(); i++) {
	    		min = std::min(min, (*this)(i,j));
	    		max = std::max(max, (*this)(i,j));
	    	}
	    }


	    for(std::size_t j = 0; j < getHeight(); j++) {
	    	unsigned char *tt_r = r + j*getWidth();
	    	unsigned char *tt_g = g + j*getWidth();
	    	unsigned char *tt_b = b + j*getWidth();
	    	for(std::size_t i = 0; i < getWidth(); i++) {
	    		//std::cout << "(" << i << "," << j << ") = " << (*this)(i,j) << std::endl;
	    		if(isNaN((*this)(i,j))) {
	    			tt_r[i] = 0;
	    			tt_g[i] = 255.0;
	    			tt_b[i] = 0;
	    		} else if(isInf((*this)(i,j))) {
	    			tt_r[i] = 255.0;
					tt_g[i] = 0;
					tt_b[i] = 0;
	    		} else if((*this)(i,j) < _low-_eps) {
	    			tt_r[i] = 0;
					tt_g[i] = 0;
					tt_b[i] = 192.0*((*this)(i,j)/min)+63;
	    		} else if(_high+_eps < (*this)(i,j)) {
	    			tt_r[i] = 192.0*((*this)(i,j)/max)+63;
					tt_g[i] = 0;
					tt_b[i] = 192.0*((*this)(i,j)/max)+63;
	    		} else {
	    			tt_r[i] = (255.0/(high-low))*RANGE((*this)(i,j), low, high);
	    			tt_g[i] = (255.0/(high-low))*RANGE((*this)(i,j), low, high);
	    			tt_b[i] = (255.0/(high-low))*RANGE((*this)(i,j), low, high);
	    		}
	    	}
	    }

	    int ret = io_png_write_u8(_fileName.c_str(), (const unsigned char *) tmp,
	    		getWidth(), getHeight(), 3);

	    if(ret != 0)
	    	throw std::exception();

		delete [] tmp;

	}
	T* getRaw() { return data; }
	std::size_t getWidth() const { return width; }
	std::size_t getHeight() const { return height; }
	T& operator()(std::size_t _i, std::size_t _j) const {
		if(0 <= _i && 0 <= _j && _i < getWidth() && _j < getHeight())
			return data[_j*getWidth()+_i];
		throw std::out_of_range("");

	}
	T& operator()(std::size_t _i, std::size_t _j) {
		if(0 <= _i && 0 <= _j && _i < getWidth() && _j < getHeight())
			return data[_j*getWidth()+_i];
		throw std::out_of_range("");
	}

	T getLow() const {
		return low;
	}
	T getHigh() const {
		return high;
	}
private:
	T* data;
	std::size_t width;
	std::size_t height;
	T low;
	T high;
};

template <typename T>
T diffNorm(const Image<T>& a, const Image<T>& b) {
	T tmp = 0.0;
	for(std::size_t j = 0; j < a.getHeight(); j++)
		for(std::size_t i = 0; i < a.getWidth(); i++)
			tmp += pow2(a(i,j)-b(i,j));
	return sqrt(tmp);
}

template <typename T>
void addNoise(Image<T>& image, T eps, bool threshold) {
	std::size_t width = image.getWidth();
	std::size_t height = image.getHeight();
	T low = image.getLow();
	T high = image.getHigh();

	if(threshold)
		for(std::size_t j = 1; j < height-1; j++)
			for(std::size_t i = 1; i < width-1; i++)
				image(i,j) = RANGE(image(i,j) + 2*eps*((T)rand()/(T)RAND_MAX) - eps, low, high);
	else
		for(std::size_t j = 1; j < height-1; j++)
			for(std::size_t i = 1; i < width-1; i++)
				image(i,j) = image(i,j) + 2*eps*((T)rand()/(T)RAND_MAX) - eps;
}

template <typename T>
void crossSection(const Image<T>& image, std::ostream& out, std::size_t x1, std::size_t x2, std::size_t y1, std::size_t y2, std::size_t count) {
	for(std::size_t i = 0; i < count; i++) {
		T t = 1.0*i/(count-1);
		out << t << " " << image((1-t)*x1+t*y1, (1-t)*x2+t*y2) << std::endl;
	}
}

#endif /* IMAGE_H_ */
