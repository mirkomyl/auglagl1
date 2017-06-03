/*
 *  Created on: May 27, 2013
 *      Author: Mirko Myllykoski (mirko.myllykoski@gmail.com)
 */

#ifndef AUGLANGL1_H_
#define AUGLANGL1_H_

#include "Debugger.h"

namespace auglagL1 {


int augLagL1(const double *input, double *output, std::size_t width, std::size_t height, std::size_t ldfw, double r1, double r2, double r3, double eps, Phi phi, double beta, double h, double delta, int maxIt, const Debugger<double>& debugger);

//int augLagL1(const float *input, float *output, size_t width, size_t height, size_t ldfw, float r1, float r2, float r3, float eps, float delta, int maxIt);

}

#endif /* AUGLANGL1_H_ */
