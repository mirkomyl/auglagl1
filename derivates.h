/*
 *  Created on: June 14, 2013
 *      Author: Mirko Myllykoski (mirko.myllykoski@gmail.com)
 */

#ifndef AUGLAGL1_DERIVATES_H
#define AUGLAGL1_DERIVATES_H

#include "V.h"
#include "Q.h"

namespace auglagL1 {

template <typename T>
Q<T> grad(const V<T>& v) {
	std::size_t n1 = 2*(v.getN1()-1);
	std::size_t n2 = v.getN2()-1;
	T h = v.getH();

	Q<T> tmp(n1, n2, h);

	for(std::size_t j = 1; j <= n2; j++)
		for(std::size_t i = 1; i <= n1; i++)
			tmp(i,j) = v.grad(i,j);

	return tmp;
}

template <typename T>
V<T> div(const Q<T> q) {
	std::size_t n1 = q.getN1()/2+1;
	std::size_t n2 = q.getN2()+1;
	T h = q.getH();

	V<T> tmp(n1, n2, h);

	for(std::size_t i = 0; i < n1; i++)
		tmp(i,0) = tmp(i,n2-1) = 0.0;

	for(std::size_t j = 1; j < n2-1; j++) {
		tmp(0,j) = tmp(0,n2-1) = 0.0;
		for(std::size_t i = 1; i < n1-1; i++)
			tmp(i,j) = q.div(i,j);
	}

	return tmp;
}

}

#endif
