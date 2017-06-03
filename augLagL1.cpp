/*
 *  Created on: May 27, 2013
 *      Author: Mirko Myllykoski (mirko.myllykoski@gmail.com)
 */

#include <math.h>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "augLagL1.h"
#include "pscrWrapper.h"

// Vertex numbering:
// (0,4)   (1,4)   (2,4)   (3,4)   (4,4)
//      +-------+-------+-------+-------+
//      |      /|      /|      /|      /|
//      |    /  |    /  |    /  |    /  |
//      |  /    |  /    |  /    |  /    |
// (0,3)|/ (1,3)|/ (2,3)|/ (3,3)|/ (4,3)|
//      |-------+-------+-------+-------|
//      |      /|      /|      /|******/|
//      |    /  |    /  |    /**|****/**|
//      |  /    |  /    |  /****|**/****| <- \Omega_{3,2}
// (0,2)|/ (1,2)|/ (2,2)|/*(3,2)|/*(4,2)|
//      |-------+-------+-------O-------| <- P_{3,1}
//      |      /|      /|******/|******/|
//      |    /  |    /  |****/**|****/  |
//      |  /    |  /    |**/****|**/    |
// (0,1)|/ (1,1)|/ (2,1)|/*(3,1)|/ (4,1)|
//      |-------+-------+-------+-------|
//      |      /|      /|      /|      /|
//      |    /  |    /  |    /  |    /  |
//      |  /    |  /    |  /    |  /    |
// (0,0)|/ (1,0)|/ (2,0)|/ (3,0)|/ (4,0)|
//      +-------+-------+-------+-------+

// ^
// | x2
// |
// |     x1
// +------->

// \Omega = (0, L1) X (0, L2)
// L1 = h * n1
// L2 = h * n2

// P_{i,j} = (i*h,j*h), i=0,...,n1, j=0,...,n2

// u \in V = \span{ v_{1,1}, ...  }
// u : \Omega -> R
// u = \sum_{i=0}^{n_1} \sum_{j=0}^{n_2} a_{i,j} v_{i,j}

// v_{i,j} : \Omega -> [0,1]
// v_{i,j}(P_{i,j}) = 1
// v_{i,j}(P_{a,b}) = 0, (a,b) != (i,j)
// v_{i,j}(\Omega_{i,j}) \in [0,1]
// v_{i,j}(\Omega \setminus \Omega_{i,j}) = 0

// x = (a_{1,0}; a_{2,0}; a_{3,0}; ... a_{1,1}; a_{2,1}; a_{3,1}; ...
// Ax = y

// Triangle numbering:
// +-------+-------+-------+-------+
// |(1,4) /|(3,4) /|(5,4) /|(7,4) /|
// |    /  |    /  |    /  |    /  |
// |  /    |  /    |  /    |  /    |
// |/ (2,4)|/ (4,4)|/ (6,4)|/ (8,4)|
// |-------+-------+-------+-------|
// |(1,3) /|(3,3) /|(5,3) /|(7,3) /|
// |    /  |    /  |    /  |    /  |
// |  /    |  /    |  /    |  /    |
// |/ (2,3)|/ (4,3)|/ (6,3)|/ (8,3)|
// |-------+-------+-------+-------|
// |(1,2) /|(3,2) /|(5,2) /|(7,2) /|
// |    /  |    /  |    /  |    /**|
// |  /    |  /    |  /    |  /****| <- \Omega^*_{8,2}
// |/ (2,2)|/ (4,2)|/ (6,2)|/*(8,2)|
// |-------+-------+-------+-------|
// |(1,1) /|(3,1) /|(5,1) /|(7,1) /|
// |    /  |    /  |    /  |    /  |
// |  /    |  /    |  /    |  /    |
// |/ (2,1)|/ (4,1)|/ (6,1)|/ (8,1)|
// +-------+-------+-------+-------+

// \Omega = \bigcup_{i=1}^{2n_1} \bigcup_{j=1}^{n_2} \Omega^*_{i,j}

// \Omega_{i,j} = \Omega^*_{2i,j+1} \cup \Omega^*_{2i+1,j+1} \cup \Omega^*_{2i+2,j+1} \cup
//                \Omega^*_{2i+1,j} \cup \Omega^*_{2i,j} \cup \Omega^*_{2i-1,j},

// \Omega^*_{i,j} = \emptyset if i < 1, j < 1, 2n1 < i or n2 < j

// p \in Q = \span{q_{1,1}_1, ... }
// p : \Omega -> R^2
// p = \sum_{i=1}^{2n_1} \sum_{j=1}^{n_2} (a_{i,j})_1 (q_{i,j})_1 + (a_{i,j})_2 (q_{i,j})_2

// q_{i,j} : \Omega -> {0,1}^2
// q_{i,j}_1(\Omega^*_{i,j}) = (1,0)
// q_{i,j}_1(\Omega \setminus \Omega^*_{i,j}) = 0
// q_{i,j}_2(\Omega^*_{i,j}) = (0,1)
// q_{i,j}_2(\Omega \setminus \Omega^*_{i,j}) = 0

// x = (a_{1,1}_1; a_{1,1}_2; a_{2,1}_1; a_{2,1}_2; a_{3,1}_1; a_{3,1}_2; ...
// Ax = y

#define DIVHQ_ODD1__UL 	-1.0
#define DIVHQ_ODD1__UR	1.0
#define DIVHQ_ODD2__UL 	1.0
#define DIVHQ_ODD2__LL 	-1.0

#define DIVHQ_EVEN1_LR	1.0
#define DIVHQ_EVEN1_LL	-1.0
#define DIVHQ_EVEN2_LR	-1.0
#define DIVHQ_EVEN2_UR	1.0

#define MASK_X1(x,i,j) ( \
	+1.0*(x)(2*(i),j+1).x1 \
	-1.0*(x)(2*(i)+2,j+1).x1 \
	+1.0*(x)(2*(i)-1,j).x1 \
	-1.0*(x)(2*(i)+1,j).x1 )

#define MASK_X2(x,i,j) ( \
	-1.0*(x)(2*(i),j+1).x2 \
	-1.0*(x)(2*(i)+1,j+1).x2 \
	+1.0*(x)(2*(i),j).x2 \
	+1.0*(x)(2*(i)+1,j).x2 )

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

namespace auglagL1 {

template <typename T>
class Matrix {
public:
	Matrix(T _r2, T _r3, T _h) {
		r2 = _r2;
		r3 = _r3;
		h = _h;
	}

	Q<T> operator*(const Q<T>& x) const {
		Q<T> y(x.getN1(), x.getN2(), x.getH());
		this->mul(y, x);
		return y;
	}

	// y = Ax
	void mul(Q<T>& y, const Q<T>& x) const {
		const std::size_t n1 = x.getN1();
		const std::size_t n2 = x.getN2();

		#pragma omp parallel for
		for(std::size_t j = 1; j <= n2; j++) {

			for(std::size_t i = 1; i <= n1; i+=2) {
				R2<T> odd, even;
				odd = even = 0;

				if(0 < (i-1)/2 && j < n2) {
					odd.x1 += DIVHQ_ODD1__UL*MASK_X1(x,(i-1)/2, j);
					odd.x2 += DIVHQ_ODD2__UL*MASK_X1(x,(i-1)/2, j);
					odd.x1 += DIVHQ_ODD1__UL*MASK_X2(x,(i-1)/2, j);
					odd.x2 += DIVHQ_ODD2__UL*MASK_X2(x,(i-1)/2, j);

				}

				if((i+1)/2 < n1/2 && j < n2) {
					odd.x1  += DIVHQ_ODD1__UR*MASK_X1(x,(i+1)/2, j);
					even.x2 += DIVHQ_EVEN2_UR*MASK_X1(x,(i+1)/2, j);
					odd.x1  += DIVHQ_ODD1__UR*MASK_X2(x,(i+1)/2, j);
					even.x2 += DIVHQ_EVEN2_UR*MASK_X2(x,(i+1)/2, j);
				}

				if(0 < (i-1)/2 && 0 < j-1) {
					odd.x2  += DIVHQ_ODD2__LL*MASK_X1(x,(i-1)/2, j-1);
					even.x1 += DIVHQ_EVEN1_LL*MASK_X1(x,(i-1)/2, j-1);
					odd.x2  += DIVHQ_ODD2__LL*MASK_X2(x,(i-1)/2, j-1);
					even.x1 += DIVHQ_EVEN1_LL*MASK_X2(x,(i-1)/2, j-1);
				}

				if((i+1)/2 < n1/2 && 0 < j-1) {
					even.x1 += DIVHQ_EVEN1_LR*MASK_X1(x,(i+1)/2, j-1);
					even.x2 += DIVHQ_EVEN2_LR*MASK_X1(x,(i+1)/2, j-1);
					even.x1 += DIVHQ_EVEN1_LR*MASK_X2(x,(i+1)/2, j-1);
					even.x2 += DIVHQ_EVEN2_LR*MASK_X2(x,(i+1)/2, j-1);
				}

				// Store and add mass
				y(i,j).x1   = r3 * odd.x1 / 4.0 + r2*h*h/2 * x(i,j).x1;
				y(i,j).x2   = r3 * odd.x2 / 4.0 + r2*h*h/2 * x(i,j).x2;
				y(i+1,j).x1 = r3 * even.x1 / 4.0 + r2*h*h/2 * x(i+1,j).x1;
				y(i+1,j).x2 = r3 * even.x2 / 4.0 + r2*h*h/2 * x(i+1,j).x2;
			}
		}
	}


	Q<T> precond(const Q<T>& x) const {
		Q<T> y(x.getN1(), x.getN2(), x.getH());
		this->precond(y, x);
		return y;
	}

	void precond(Q<T>& y, const Q<T>& x) const {
		std::size_t n1 = x.getN1();
		std::size_t n2 = x.getN2();
//y=x; return;
		// Bottom
		y(1,1).x1 = x(1,1).x1 / (0.25 * r3 * 1.0 + r2*h*h/2);
		y(1,1).x2 = x(1,1).x2 / (0.25 * r3 * 0.0 + r2*h*h/2);
		y(2,1).x1 = x(2,1).x1 / (0.25 * r3 * 0.0 + r2*h*h/2);
		y(2,1).x2 = x(2,1).x2 / (0.25 * r3 * 1.0 + r2*h*h/2);
		for(std::size_t i = 3; i <= n1-2; i+=2) {
			y(i,1).x1   = x(i,1).x1   / (0.25 * r3 * 2.0 + r2*h*h/2);
			y(i,1).x2   = x(i,1).x2   / (0.25 * r3 * 1.0 + r2*h*h/2);
			y(i+1,1).x1 = x(i+1,1).x1 / (0.25 * r3 * 0.0 + r2*h*h/2);
			y(i+1,1).x2 = x(i+1,1).x2 / (0.25 * r3 * 1.0 + r2*h*h/2);
		}
		y(n1-1,1).x1 = x(n1-1,1).x1 / (0.25 * r3 * 1.0 + r2*h*h/2);
		y(n1-1,1).x2 = x(n1-1,1).x2 / (0.25 * r3 * 1.0 + r2*h*h/2);
		y(n1,1).x1   = x(n1,1).x1   / (0.25 * r3 * 0.0 + r2*h*h/2);
		y(n1,1).x2   = x(n1,1).x2   / (0.25 * r3 * 0.0 + r2*h*h/2);

		// Middle
		for(std::size_t j = 2; j <= n2-1; j++) {
			y(1,j).x1 = x(1,j).x1 / (0.25 * r3 * 1.0 + r2*h*h/2);
			y(1,j).x2 = x(1,j).x2 / (0.25 * r3 * 0.0 + r2*h*h/2);
			y(2,j).x1 = x(2,j).x1 / (0.25 * r3 * 1.0 + r2*h*h/2);
			y(2,j).x2 = x(2,j).x2 / (0.25 * r3 * 2.0 + r2*h*h/2);
			for(std::size_t i = 3; i <= n1-2; i+=2) {
				y(i,j).x1   = x(i,j).x1   / (0.25 * r3 * 2.0 + r2*h*h/2);
				y(i,j).x2   = x(i,j).x2   / (0.25 * r3 * 2.0 + r2*h*h/2);
				y(i+1,j).x1 = x(i+1,j).x1 / (0.25 * r3 * 2.0 + r2*h*h/2);
				y(i+1,j).x2 = x(i+1,j).x2 / (0.25 * r3 * 2.0 + r2*h*h/2);
			}
			y(n1-1,j).x1 = x(n1-1,j).x1 / (0.25 * r3 * 1.0 + r2*h*h/2);
			y(n1-1,j).x2 = x(n1-1,j).x2 / (0.25 * r3 * 2.0 + r2*h*h/2);
			y(n1,j).x1   = x(n1,j).x1   / (0.25 * r3 * 1.0 + r2*h*h/2);
			y(n1,j).x2   = x(n1,j).x2   / (0.25 * r3 * 0.0 + r2*h*h/2);
		}

		// Top
		y(1,n2).x1 = x(1,n2).x1 / (0.25 * r3 * 0.0 + r2*h*h/2);
		y(1,n2).x2 = x(1,n2).x2 / (0.25 * r3 * 0.0 + r2*h*h/2);
		y(2,n2).x1 = x(2,n2).x1 / (0.25 * r3 * 1.0 + r2*h*h/2);
		y(2,n2).x2 = x(2,n2).x2 / (0.25 * r3 * 1.0 + r2*h*h/2);
		for(std::size_t i = 3; i <= n1-2; i+=2) {
			y(i,n2).x1   = x(i,n2).x1   / (0.25 * r3 * 0.0 + r2*h*h/2);
			y(i,n2).x2   = x(i,n2).x2   / (0.25 * r3 * 1.0 + r2*h*h/2);
			y(i+1,n2).x1 = x(i+1,n2).x1 / (0.25 * r3 * 2.0 + r2*h*h/2);
			y(i+1,n2).x2 = x(i+1,n2).x2 / (0.25 * r3 * 1.0 + r2*h*h/2);
		}
		y(n1-1,n2).x1 = x(n1-1,n2).x1 / (0.25 * r3 * 0.0 + r2*h*h/2);
		y(n1-1,n2).x2 = x(n1-1,n2).x2 / (0.25 * r3 * 1.0 + r2*h*h/2);
		y(n1,n2).x1   = x(n1,n2).x1   / (0.25 * r3 * 1.0 + r2*h*h/2);
		y(n1,n2).x2   = x(n1,n2).x2   / (0.25 * r3 * 0.0 + r2*h*h/2);

	}

private:
	T r2;
	T r3;
	T h;
};

template class Matrix<double>;

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


// template <typename T>
// int CG(const Matrix<T> &A, Q<T>& x, const Q<T>& b, int *max_iter, T *tol) {
// 
// 	std::size_t n1 = x.getN1();
// 	std::size_t n2 = x.getN2();
// 	T h = x.getH();
// 
// 	T resid;
// 	T alpha, beta, rho, rho_1;
// 	rho_1 = 0;
// 
// 	Q<T> p(n1, n2, h);
// 	Q<T> q(n1, n2, h);
// 	Q<T> r(n1, n2, h);
// 	Q<T> z(n1, n2, h);
// 
// 	T normb = vectorNorm2(b);
// 	A.mul(q, x);
// 	#pragma omp parallel for
// 	for(std::size_t j = 1; j <= n2; j++)
// 		for(std::size_t i = 1; i <= n1; i++)
// 			r(i,j) = b(i,j) - q(i,j);
// 
// 	if (normb == 0.0)
// 		normb = 1;
// 
// 	if ((resid = vectorNorm2(r) / normb) <= *tol) {
// 		*tol = resid;
// 		*max_iter = 0;
// 		return 0;
// 	}
// 
// 	for (int i = 1; i <= *max_iter; i++) {
// 		A.precond(z, r);
// 		rho = dot(r, z);
// 
// 		if (i == 1)
// 			p = z;
// 		else {
// 			beta = rho / rho_1;
// 			#pragma omp parallel for
// 			for(std::size_t j = 1; j <= n2; j++)
// 				for(std::size_t i = 1; i <= n1; i++)
// 					p(i,j) = z(i,j) + beta * p(i,j);
// 		}
// 
// 		A.mul(q, p);
// 		alpha = rho / dot(p, q);
// 
// 		#pragma omp parallel for
// 		for(std::size_t j = 1; j <= n2; j++) {
// 			for(std::size_t i = 1; i <= n1; i++) {
// 				x(i,j) += alpha * p(i,j);
// 				r(i,j) -= alpha * q(i,j);
// 			}
// 		}
// 
// 		if ((resid = vectorNorm2(r) / normb) <= *tol) {
// 			*tol = resid;
// 			*max_iter = i;
// 			return 0;
// 		}
// 
// #if DEBUG
// 		if(i % 100 == 0)
// 			std::cout << ">>> CG: " << i << " iterations, diff = " << resid << std::endl;
// #endif
// 
// 		rho_1 = rho;
// 	}
// 
// 	*tol = resid;
// 	return 1;
// }

template <typename T>
int CG(const Matrix<T> &A, Q<T>& x, const Q<T>& b, int *max_iter, T *tol) {

	std::size_t n1 = x.getN1();
	std::size_t n2 = x.getN2();
	T h = x.getH();

	T resid;
	T alpha, beta, rho;

	Q<T> r(n1, n2, h);
	Q<T> p(n1, n2, h);
	Q<T> Ap(n1, n2, h);

	T normb = vectorNorm2(b);
	A.mul(Ap, x);
	#pragma omp parallel for
	for(std::size_t j = 1; j <= n2; j++)
		for(std::size_t i = 1; i <= n1; i++)
			r(i,j) = b(i,j) - Ap(i,j);

	if (normb == 0.0)
		normb = 1;

	if ((resid = vectorNorm2(r) / normb) <= *tol) {
		*tol = resid;
		*max_iter = 0;
		return 0;
	}
	
	p = r;

	for (int i = 1; i <= *max_iter; i++) {
		rho = dot(r,r);
		
		A.mul(Ap, p);
		alpha = rho / dot(p,Ap);
	
		#pragma omp parallel for
		for(std::size_t j = 1; j <= n2; j++) {
			for(std::size_t i = 1; i <= n1; i++) {
				x(i,j) += alpha * p(i,j);
				r(i,j) -= alpha * Ap(i,j);
			}
		}
		
		if ((resid = vectorNorm2(r) / normb) <= *tol) {
			*tol = resid;
			*max_iter = i;
			return 0;
		}
		
		beta = dot(r,r) / rho;
			
		for(std::size_t j = 1; j <= n2; j++)
			for(std::size_t i = 1; i <= n1; i++)
				p(i,j) = r(i,j) + beta * p(i,j);

#if DEBUG
		if(i % 100 == 0)
			std::cout << ">>> CG: " << i << " iterations, diff = " << resid << std::endl;
#endif

	}

	*tol = resid;
	return 1;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

template <typename T>
T d_d0_g(T x, T b11, T b12, T b21, T b22, T r1, T r2, T beta) {
	return (pow2(x)/2)*(r1+r2/(beta+pow2(x)))-x*sqrt(pow2(b11+b21/sqrt(beta+pow2(x)))+pow2(b12+b22/sqrt(beta+pow2(x))));
}

template <typename T>
T d_d0_f(T x, T y, T b11, T b12, T b21, T b22, T r1, T r2, T beta) {
	return ((pow2(x)+pow2(y))/2)*(r1+r2/(beta+pow2(x)+pow2(y)))-((b11+b21/sqrt(beta+pow2(x)+pow2(y)))*x+(b12+b22/sqrt(beta+pow2(x)+pow2(y)))*y);
}

template <typename T>
T d_dx_f(T x, T y, T b11, T b12, T b21, T b22, T r1, T r2, T beta) {
	return x*(r1+r2/(pow2(x)+pow2(y)+beta))-(r2*x*(pow2(x)+pow2(y)))/pow2(pow2(x)+pow2(y)+beta)-b11+(b21*pow2(x))/pow3_2(pow2(x)+pow2(y)+beta)-b21/sqrt(pow2(x)+pow2(y)+beta)+(b22*x*y)/pow3_2(pow2(x)+pow2(y)+beta);
}

template <typename T>
T d_dy_f(T x, T y, T b11, T b12, T b21, T b22, T r1, T r2, T beta) {
	return y*(r1+r2/(pow2(x)+pow2(y)+beta))-(r2*y*(pow2(x)+pow2(y)))/pow2(pow2(x)+pow2(y)+beta)-b12+(b21*x*y)/pow3_2(pow2(x)+pow2(y)+beta)+(b22*pow2(y))/pow3_2(pow2(x)+pow2(y)+beta)-b22/sqrt(pow2(x)+pow2(y)+beta);
}

template <typename T>
T d_dxx_f(T x, T y, T b11, T b12, T b21, T b22, T r1, T r2, T beta) {
	return r1-(4*r2*pow2(x))/pow2(pow2(x)+pow2(y)+beta)+(4*r2*pow2(x)*(pow2(x)+pow2(y)))/pow3(pow2(x)+pow2(y)+beta)+r2/(pow2(x)+pow2(y)+beta)-(r2*(pow2(x)+pow2(y)))/pow2(pow2(x)+pow2(y)+beta)+(3*b21*x)/pow3_2(pow2(x)+pow2(y)+beta)-(3*b21*pow3(x))/pow5_2(pow2(x)+pow2(y)+beta)-(3*b22*pow2(x)*y)/pow5_2(pow2(x)+pow2(y)+beta)+(b22*y)/pow3_2(pow2(x)+pow2(y)+beta);
}

template <typename T>
T d_dyy_f(T x, T y, T b11, T b12, T b21, T b22, T r1, T r2, T beta) {
	return r1-(4*r2*pow2(y))/pow2(pow2(x)+pow2(y)+beta)+(4*r2*pow2(y)*(pow2(x)+pow2(y)))/pow3(pow2(x)+pow2(y)+beta)+r2/(pow2(x)+pow2(y)+beta)-(r2*(pow2(x)+pow2(y)))/pow2(pow2(x)+pow2(y)+beta)-(3*b21*x*pow2(y))/pow5_2(pow2(x)+pow2(y)+beta)+(b21*x)/pow3_2(pow2(x)+pow2(y)+beta)+(3*b22*y)/pow3_2(pow2(x)+pow2(y)+beta)-(3*b22*pow3(y))/pow5_2(pow2(x)+pow2(y)+beta);
}

template <typename T>
T d_dxy_f(T x, T y, T b11, T b12, T b21, T b22, T r1, T r2, T beta) {
	return -(4*r2*x*y)/pow2(pow2(x)+pow2(y)+beta)+(4*r2*x*y*(pow2(x)+pow2(y)))/pow3(pow2(x)+pow2(y)+beta)-(3*b21*pow2(x)*y)/pow5_2(pow2(x)+pow2(y)+beta)+(b21*y)/pow3_2(pow2(x)+pow2(y)+beta)+(b22*x)/pow3_2(pow2(x)+pow2(y)+beta)-(3*b22*x*pow2(y))/pow5_2(pow2(x)+pow2(y)+beta);
}

template <typename T>
T lambda1(T a, T b, T c, T d) {
	return 0.5 * (-sqrt(a*a-2*a*d+4*b*c+d*d)+a+d);
}

template <typename T>
T lambda2(T a, T b, T c, T d) {
	return 0.5 * (+sqrt(a*a-2*a*d+4*b*c+d*d)+a+d);
}

template <typename T>
R2<T> newton(const R2<T>& x0, const R2<T>& b1, const R2<T>& b2, T r1, T r2, T beta, int iterLimit, T tol, int* _k) {
	R2<T> x = x0;
	R2<T> old = x0;

	int k;
	for(k = 1; k <= iterLimit; k++) {
		T a = d_dxx_f(x.x1, x.x2, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);
		T b = d_dxy_f(x.x1, x.x2, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);
		T c = b;
		T d = d_dyy_f(x.x1, x.x2, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);

		T dx1 = d_dx_f(x.x1, x.x2, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);
		T dx2 = d_dy_f(x.x1, x.x2, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);

		if(sqrt(dx1*dx1+dx2*dx2) < tol && norm(x-old)/norm(x) < tol && 0 < lambda1(a,b,c,d) && 0 < lambda2(a,b,c,d))
			break;

		old = x;

		T det = a*d-b*c;
		x.x1 -= (1.0/det)*(d*dx1 - b*dx2);
		x.x2 -= (1.0/det)*(a*dx2 - c*dx1);
	}
	*_k = k;

	return x;
}

template <typename T>
T d_d1_g(T x, T b11, T b12, T b21, T b22, T r1, T r2, T beta) {
	return ( d_d0_g(x+1.0E-6, b11, b12, b21, b22, r1, r2, beta) -  d_d0_g(x, b11, b12, b21, b22, r1, r2, beta) ) / 1.0E-6;
}

template <typename T>
R2<T> secant(T _a_x, T _b_x, const R2<T>& b1, const R2<T>& b2, T r1, T r2, T beta, int iterLimit, T tol, int* _k) {
	T a_x = _a_x;
	T a = d_d1_g(a_x, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);

	T b_x = _b_x;
	T b = d_d1_g(b_x, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);
	
	int k;
	for(k = 1; k <= iterLimit; k++) {
		T c_x = std::max(0.0, b_x - (b_x - a_x)*b/(b - a));
		a_x = b_x;
		a = b;
		b_x = c_x;
		b = d_d1_g(c_x, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);
		
		if(::fabs(a_x-b_x)/::fabs(b_x) < tol && ::fabs((a+b)/2.0) < tol)
			break;
	}
	*_k = k;

	T rho = (b_x + a_x)/2.0;

	T lambda = rho / norm(b1+b2/sqrt(1.0+pow2(rho)));

	return lambda*(b1+b2/sqrt(1+pow2(rho)));
}

template <typename T>
R2<T> bisection(T _a_x, T _b_x, const R2<T>& b1, const R2<T>& b2, T r1, T r2, T beta, T tol, bool jump = false) {
	T a_x = _a_x;
	T a = d_d0_g(a_x, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);

	T b_x = _b_x;
	T b = d_d0_g(b_x, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);

	if(jump) {
		// Infinite looping is extremely rare but possible
		for(int i = 0; i < 10 && b < a; i++) {
			b_x *= 2;
			b = d_d0_g(b_x, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);
		}
	}

	while(tol*0.5*(b_x + a_x) < b_x - a_x) {
		T c1_x = a_x + (b_x - a_x)/2.0 - std::min(1.0E-2, (b_x - a_x)/10.0);
		T c1 = d_d0_g(c1_x, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);

		T c2_x = a_x + (b_x - a_x)/2.0 + std::min(1.0E-2, (b_x - a_x)/10.0);
		T c2 = d_d0_g(c2_x, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);

		if(c1 <= c2)
			b_x = c2_x;
		else
			a_x = c1_x;
	}

	T rho = (b_x + a_x)/2.0;

	T lambda = rho / norm(b1+b2/sqrt(1.0+pow2(rho)));

	return lambda*(b1+b2/sqrt(1+pow2(rho)));
}

template <typename T>
T checkDirection(const R2<T>& x, const R2<T>& b1, const R2<T>& b2, T beta) {
	T rho = norm(x);
	T lambda = rho / norm(b1+b2/sqrt(beta+pow2(rho)));
	R2<T> x0 = lambda*(b1+b2/sqrt(beta+pow2(rho)));
	return norm(x-x0)/norm(x0);
}

// Solve subproblem:
//
template <typename T>
int solveSub1(
		      Q<T>&      p1,
		      Q<T>&      p2,
		const V<T>& u,
		const Q<T>&      p3,
		const Q<T>&      l1,
		const Q<T>&      l2,
		      T               r1,
		      T               r2,
			  T          beta,
		      T               h,
		      T			tol,
		      bool       first = false) {

	const std::size_t n1 = p1.getN1();
	const std::size_t n2 = p1.getN2();

	const int iterLimit = 10;

	int maxIt = 0;
	int minIt = iterLimit+1;
	double avgIt = 0.0;
	int newtonCount = 0;
	int newtonSuccess = 0;
	int success = 0;

	#pragma omp parallel for
	for(std::size_t j = 1; j <= n2; j++) {
		int maxIt2 = 0;
		int minIt2 = iterLimit+1;
		double avgIt2 = 0.0;
		int newtonCount2 = 0;
		int newtonSuccess2 = 0;
		int success2 = 0;
		
		for(std::size_t i = 1; i <= n1; i++) {
			R2<T> b1 = r1*u.grad(i,j) + l1(i,j);
			R2<T> b2 = r2*p3(i,j) - l2(i,j);

			if(b1 == 0 && b2 == 0) {
				success2++;
				p1(i,j) = p2(i,j) = 0.0;
				continue;
			}

			newtonCount2++;

			T min1Val = INF;
			T min2Val = INF;
			T prevVal = d_d0_f(p1(i,j).x1, p1(i,j).x2, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);

			R2<T> min1, min2;

			int k;
			min1 = newton(p1(i,j), b1, b2, r1, r2, beta, iterLimit, tol, &k);
			min1Val = d_d0_f(min1.x1, min1.x2, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);
			
			if(!first && k <= iterLimit && checkDirection(min1, b1, b2, beta) < tol && min1Val <= prevVal) {
				maxIt2 = std::max(maxIt2, k);
				minIt2 = std::min(minIt2, k);
				avgIt2 += k;
				newtonSuccess2++;
			} else {
				T a = 0.0;
				T b = 1;
				R2<T> min = 0.0;
				T minVal = d_d0_f(min.x1, min.x2, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);
				do {
					R2<T> tmp = bisection(a, b, b1, b2, r1, r2, beta, 1.0E-3, false);
					T tmpVal = d_d0_f(tmp.x1, tmp.x2, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);
					if(tmpVal < minVal) {
 						min = tmp;
 						minVal = tmpVal;
					}
					a = b;
					b *= 4;
				} while(b < sqrt(2)/h);
				R2<T> tmp = bisection(a, b, b1, b2, r1, r2, beta, 1.0E-3, true);
				T tmpVal = d_d0_f(tmp.x1, tmp.x2, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);
				if(tmpVal < minVal)
					min = tmp;
				min2 = bisection(norm(min)*(1.0-1.0E-3), norm(min)*(1.0+1.0E-3), b1, b2, r1, r2, beta, tol, false);
				min2Val = d_d0_f(min2.x1, min2.x2, b1.x1, b1.x2, b2.x1, b2.x2, r1, r2, beta);
			}

			R2<T> min;
			T minVal;
			if(min1Val < min2Val) {
				min = min1;
				minVal = min1Val;
			} else {
				min = min2;
				minVal = min2Val;
			}

			if(minVal <= prevVal) {
				success2++;
				p1(i,j) = min;
				p2(i,j) = min/sqrt(beta+dot(min));
			}/* else {
#if DEBUG
				#pragma omp critical (print)
				std::cout << ">>> Failure at (" << i << ", " << j << "), b1 = " << toString(b1) << ", b2 = " << toString(b2) << std::endl;
#endif
			}*/
		}
		
		#pragma omp critical
		{
			maxIt = std::max(maxIt, maxIt2);
			minIt = std::min(minIt, minIt2);
			avgIt += avgIt2;
			newtonCount += newtonCount2;
			newtonSuccess += newtonSuccess2;
			success += success2;
		}
		
	}

#if DEBUG
	std::cout << ">>> Newton success rate = " << 1.0*newtonSuccess/newtonCount << ", total success rate = " << 1.0*success/(n1*n2) << std::endl;
	std::cout << ">>> Newton iterations: min = " << minIt << ", max = " << (maxIt <= iterLimit ? toString(maxIt) : "_max_") << ", avg = " << avgIt/newtonCount << std::endl;
#endif

	return 0;
}

template <typename T>
int solveSub2(
		      Q<T>&      p3,
		const V<T>& psi,
		const Q<T>&      p2,
		const Q<T>&      l2,
		const V<T>& l3,
		      T               r2,
		      T               r3,
		      T               h,
		      int             maxIter,
		      T               tol) {

	const std::size_t n1 = p2.getN1();
	const std::size_t n2 = p2.getN2();

	Q<T> b(n1,n2,h);

	#pragma omp parallel for
	for(std::size_t j = 1; j <= n2; j++) {
		for(std::size_t i = 1; i <= n1; i+=2) {
			b(i,j).x1 = (r2*p2(i,j).x1+l2(i,j).x1)*(h*h/2);

			if(1 < i && j < n2)
				b(i,j).x1 -= (1.0/2.0)*DIVHQ_ODD1__UL*h*(r3*VER_ODD_UL(psi,i,j)-VER_ODD_UL(l3,i,j));
			if(i < n1-1 && j < n2)
				b(i,j).x1 -= (1.0/2.0)*DIVHQ_ODD1__UR*h*(r3*VER_ODD_UR(psi,i,j)-VER_ODD_UR(l3,i,j));

			b(i,j).x2 = (r2*p2(i,j).x2+l2(i,j).x2)*(h*h/2);

			if(1 < i && j < n2)
				b(i,j).x2 -= (1.0/2.0)*DIVHQ_ODD2__UL*h*(r3*VER_ODD_UL(psi,i,j)-VER_ODD_UL(l3,i,j));
			if(1 < i && 1 < j)
				b(i,j).x2 -= (1.0/2.0)*DIVHQ_ODD2__LL*h*(r3*VER_ODD_LL(psi,i,j)-VER_ODD_LL(l3,i,j));

			b(i+1,j).x1 = (r2*p2(i+1,j).x1+l2(i+1,j).x1)*(h*h/2);

			if(1 < i && 1 < j)
				b(i+1,j).x1 -= (1.0/2.0)*DIVHQ_EVEN1_LL*h*(r3*VER_EVEN_LL(psi,i+1,j)-VER_EVEN_LL(l3,i+1,j));
			if(i < n1-1 && 1 < j)
				b(i+1,j).x1 -= (1.0/2.0)*DIVHQ_EVEN1_LR*h*(r3*VER_EVEN_LR(psi,i+1,j)-VER_EVEN_LR(l3,i+1,j));

			b(i+1,j).x2 = (r2*p2(i+1,j).x2+l2(i+1,j).x2)*(h*h/2);

			if(i < n1-1 && j < n2)
				b(i+1,j).x2 -= (1.0/2.0)*DIVHQ_EVEN2_UR*h*(r3*VER_EVEN_UR(psi,i+1,j)-VER_EVEN_UR(l3,i+1,j));
			if(i < n1-1 && 1 < j)
				b(i+1,j).x2 -= (1.0/2.0)*DIVHQ_EVEN2_LR*h*(r3*VER_EVEN_LR(psi,i+1,j)-VER_EVEN_LR(l3,i+1,j));
		}
	}

	int _maxIter = maxIter;
	T _tol = tol;

	Matrix<T> M(r2, r3, h);
	int ret = CG(M, p3, b, &_maxIter, &_tol);

#if DEBUG
	std::cout << ">>> CG: " << _maxIter << " iterations, tolerance = " << _tol << ", residual = " << vectorNorm2(M*p3-b)/vectorNorm2(p3) << std::endl;
#endif

	if(maxIter <= _maxIter)
		std::cout << ">>> Warning: CG reached the maximum amount of iterations." << std::endl;
	if(tol < _tol)
		std::cout << ">>> Warning: CG error residual is greater than the given tolerance." << std::endl;

	return ret;
}

template <typename T>
int solveSub3(V<T>& psi, const Q<T>& p3, const V<T>& l3, T r3, T eps, Phi phi, T h) {

	const std::size_t n1 = psi.getN1();
	const std::size_t n2 = psi.getN2();

	for(std::size_t i = 0; i < n1; i++) {
		psi(i,0) = psi(i,n2-1) = 0.0;
	}

	if(phi == abs) {
	
		#pragma omp parallel for
		for(std::size_t j = 1; j < n2-1; j++) {
			psi(0,j) = psi(n1-1,j) = 0.0;

			for(std::size_t i = 1; i < n1-1; i++) {
				T x = r3*p3.div(i,j) + l3(i,j);
				psi(i,j) = (1.0/r3)*sgn(x)*std::max(0.0, ::fabs(x)-eps);
			}
		}
	}
	
	if(phi == square) {
	
		#pragma omp parallel for
		for(std::size_t j = 1; j < n2-1; j++) {
			psi(0,j) = psi(n1-1,j) = 0.0;

			for(std::size_t i = 1; i < n1-1; i++) {
				T x = r3*p3.div(i,j) + l3(i,j);
				psi(i,j) = x/(2*eps+r3);
			}
		}
	}

	return 0;
}

template <typename T>
pscrWrapper::PSCR<T> initSub4(std::size_t n1, std::size_t n2, std::size_t ldf, T h, int *err) {
	T a1Diag[n2-2], a1OffDiag[n2-2], m1Diag[n2-2];
	T a2Diag[n1-2], a2OffDiag[n1-2], m2Diag[n1-2];
	for(std::size_t i = 0; i < n2-2; i++) {
		a1Diag[i] = 2.0;
		a1OffDiag[i] = -1.0;
		m1Diag[i] = 1;
	}

	for(std::size_t i = 0; i < n1-2; i++) {
		a2Diag[i] = 2.0;
		a2OffDiag[i] = -1.0;
		m2Diag[i] = 1;
	}

	return pscrWrapper::PSCR<T>(a1Diag, a1OffDiag, a2Diag,
			a2OffDiag, m1Diag, m2Diag,
			n2-2, n1-2, ldf, err);
}

template <typename T>
int solveSub4(V<T>& u, const pscrWrapper::PSCR<T>& solver, const V<T>& f, const Q<T>& p1, const Q<T>& l1, T r1, T h) {

	const std::size_t n1 = f.getN1();
	const std::size_t n2 = f.getN2();
	const std::size_t ldf = solver.getLdf();

	T *right = new T[ldf*(n2-2)];

	#pragma omp parallel for
	for(std::size_t j = 1; j < n2-1; j++) {
		for(std::size_t i = 1; i < n1-1; i++) {
			T tmp = 0;

			tmp += dot(r1*p1(2*i,j+1)-l1(2*i,j+1), R2<T>(1.0, -1.0));
			tmp += dot(r1*p1(2*i+1,j+1)-l1(2*i+1,j+1), R2<T>(0.0, -1.0));
			tmp += dot(r1*p1(2*i+2,j+1)-l1(2*i+2,j+1), R2<T>(-1.0, 0.0));

			tmp += dot(r1*p1(2*i-1,j)-l1(2*i-1,j), R2<T>(1.0, 0.0));
			tmp += dot(r1*p1(2*i,j)-l1(2*i,j), R2<T>(0.0, 1.0));
			tmp += dot(r1*p1(2*i+1,j)-l1(2*i+1,j), R2<T>(-1.0, 1.0));

			T g = 0;
			if(j == 1) g += f(i,0);
			if(i == 1) g += f(0,j);
			if(j == n2-2) g += f(i,n2-1);
			if(i == n1-2) g += f(n1-1,j);

			right[(j-1)*ldf+i-1] = g + 0.5*h*tmp/r1 + h*h*f(i,j)/r1;
		}
	}

	int err;
	err = solver.run(right, h*h/r1);

	for(std::size_t i = 0; i < n1; i++) {
		u(i,0) = f(i,0);
		u(i,n2-1) = f(i,n2-1);
	}

	#pragma omp parallel for
	for(std::size_t j = 1; j < n2-1; j++) {
		u(0,j) = f(0,j);
		u(n1-1,j) = f(n1-1,j);
		for(std::size_t i = 1; i < n1-1; i++)
			u(i,j) = right[(j-1)*ldf+i-1];
	}

	delete [] right;

	return err;
}

template <typename T>
int initHelper(V<T>& u, const pscrWrapper::PSCR<T>& solver, const V<T>& f, T h) {

	const std::size_t n1 = f.getN1();
	const std::size_t n2 = f.getN2();
	const std::size_t ldf = solver.getLdf();

	T *right = new T[ldf*(n2-2)];

	for(std::size_t j = 1; j < n2-1; j++) {
		for(std::size_t i = 1; i < n1-1; i++) {

			T g = 0;
			if(j == 1) g += f(i,0);
			if(i == 1) g += f(0,j);
			if(j == n2-2) g += f(i,n2-1);
			if(i == n1-2) g += f(n1-1,j);

			right[(j-1)*ldf+i-1] = g + f(i,j);
		}
	}

	int err;
	err = solver.run(right, 1.0);

	for(std::size_t i = 0; i < n1; i++) {
		u(i,0) = f(i,0);
		u(i,n2-1) = f(i,n2-1);
	}

	for(std::size_t j = 1; j < n2-1; j++) {
		u(0,j) = f(0,j);
		u(n1-1,j) = f(n1-1,j);
		for(std::size_t i = 1; i < n1-1; i++)
			u(i,j) = right[(j-1)*ldf+i-1];
	}

	delete [] right;

	return err;

/*
	const std::size_t n1 = f.getN1();
	const std::size_t n2 = f.getN2();
	const std::size_t ldf = solver.getLdf();

	for(std::size_t i = 0; i < n1; i++) {
		u(i,0) = f(i,0);
		u(i,n2-1) = f(i,n2-1);
	}

	for(std::size_t j = 1; j < n2-1; j++) {
		u(0,j) = f(0,j);
		u(n1-1,j) = f(n1-1,j);
		for(std::size_t i = 1; i < n1-1; i++) {
			T tmp[9];

			tmp[0] = f(i-1,j-1);
			tmp[1] = f(i-1,j);
			tmp[2] = f(i-1,j+1);

			tmp[3] = f(i,j-1);
			tmp[4] = f(i,j);
			tmp[5] = f(i,j+1);

			tmp[6] = f(i+1,j-1);
			tmp[7] = f(i+1,j);
			tmp[8] = f(i+1,j+1);

			std::sort(tmp, tmp+9);

			u(i,j) = tmp[4];
		}
	}

	return 0;*/
}

template <typename T>
T objectFunction(const DebugHelper<T>& helper) {
	return objectFunction(helper.u, helper.f, helper.eps, helper.phi, helper.beta);
}

template <typename T>
T objectFunction(const V<T>& v, const V<T> f, T eps, Phi phi, T beta) {
	T ret = 0.0;

	std::size_t n1 = v.getN1();
	std::size_t n2 = v.getN2();
	T h = v.getH();

	Q<T> vv(2*(n1-1), n2-1, h);

	#pragma omp parallel for
	for(std::size_t j = 1; j <= n2-1; j++)
		for(std::size_t i = 1; i <= 2*(n1-1); i++)
			vv(i,j) = v.grad(i,j)/sqrt((beta+dot(v.grad(i,j))));

	V<T> div_vv(n1,n2,h);

	for(std::size_t i = 0; i < n1; i++)
		div_vv(i,0) = div_vv(i,n2-1) = 0.0;

	#pragma omp parallel for
	for(std::size_t j = 1; j < n2-1; j++) {
		div_vv(0,j) = div_vv(n1-1,j) = 0.0;
		for(std::size_t i = 1; i < n1-1; i++) {
			div_vv(i,j) = vv.div(i,j);
		}
	}

	if(phi == abs) {
	
		for(std::size_t j = 1; j <= n2-1; j++) {
			for(std::size_t i = 1; i <= 2*(n1-1); i+=2) {
				ret += eps * (
						::fabs(VER_ODD_UL(div_vv, i, j)) +
						::fabs(VER_ODD_UR(div_vv, i, j)) +
						::fabs(VER_ODD_LL(div_vv, i, j)));

				ret += eps * (
						::fabs(VER_EVEN_UR(div_vv, i+1, j)) +
						::fabs(VER_EVEN_LL(div_vv, i+1, j)) +
						::fabs(VER_EVEN_LR(div_vv, i+1, j)));
			}
		}
	}
	
	if(phi == square) {
	
		for(std::size_t j = 1; j <= n2-1; j++) {
			for(std::size_t i = 1; i <= 2*(n1-1); i+=2) {
				ret += eps * (
						pow2(VER_ODD_UL(div_vv, i, j)) +
						pow2(VER_ODD_UR(div_vv, i, j)) +
						pow2(VER_ODD_LL(div_vv, i, j)));

				ret += eps * (
						pow2(VER_EVEN_UR(div_vv, i+1, j)) +
						pow2(VER_EVEN_LL(div_vv, i+1, j)) +
						pow2(VER_EVEN_LR(div_vv, i+1, j)));
			}
		}
	}
	
	for(std::size_t j = 1; j <= n2-1; j++) {
		for(std::size_t i = 1; i <= 2*(n1-1); i+=2) {

			ret += 0.5 * (
					pow2(VER_ODD_UL(f, i, j)-VER_ODD_UL(v, i, j)) +
					pow2(VER_ODD_UR(f, i, j)-VER_ODD_UR(v, i, j)) +
					pow2(VER_ODD_LL(f, i, j)-VER_ODD_LL(v, i, j)));

			ret += 0.5 * (
					pow2(VER_EVEN_UR(f, i+1, j)-VER_EVEN_UR(v, i+1, j)) +
					pow2(VER_EVEN_LL(f, i+1, j)-VER_EVEN_LL(v, i+1, j)) +
					pow2(VER_EVEN_LR(f, i+1, j)-VER_EVEN_LR(v, i+1, j)));
		}
	}

	return h*h/6.0 * ret;
}

template <typename T>
T augObjectFunction(const DebugHelper<T>& helper) {
	return augObjectFunction(helper.psi, helper.f, helper.u, helper.p1, helper.l1, helper.p2,
			helper.l2, helper.p3, helper.l3, helper.eps, helper.r1, helper.r2, helper.r3, helper.phi, helper.beta);
}

template <typename T>
T augObjectFunction(
		const V<T>& psi,
		const V<T>& f,
		const V<T>& v,
		const Q<T>& p1,
		const Q<T>& l1,
		const Q<T>& p2,
		const Q<T>& l2,
		const Q<T>& p3,
		const V<T>& l3,
		T eps, T r1, T r2, T r3, Phi phi, T beta) {

	T ret = 0.0;

	std::size_t n1 = v.getN1();
	std::size_t n2 = v.getN2();
	T h = v.getH();

	if(phi == abs) {
	
		for(std::size_t j = 1; j < n2-1; j++)
			for(std::size_t i = 1; i < n1-1; i++)
				ret += eps * ::fabs(psi(i,j)) * h*h;
		
	}
	
	if(phi == square) {
	
		for(std::size_t j = 1; j < n2-1; j++)
			for(std::size_t i = 1; i < n1-1; i++)
				ret += eps * pow2(psi(i,j)) * h*h;
		
	}

	ret += 0.5 * pow2(f(0,0)-v(0,0)) * h*h/3.0;
	for(std::size_t i = 1; i < n1-1; i++) {
		ret += 0.5 * pow2(f(i,0)-v(i,0)) * h*h/2.0;
	}
	ret += 0.5 * pow2(f(n1-1,0)-v(n1-1,0)) * h*h/6.0;

	for(std::size_t j = 1; j < n2-1; j++) {
		ret += pow2(f(0,j)-v(0,j)) * h*h/3.0;
		for(std::size_t i = 1; i < n1-1; i++) {
			ret += 0.5 * pow2(f(i,j)-v(i,j)) * h*h;
			ret += 0.5 * r3 * pow2(p3.div(i,j) - psi(i,j)) * h*h;
			ret += l3(i,j) * (p3.div(i,j)-psi(i,j)) * h*h;
		}
		ret += pow2(f(n1-1,j)-v(n1-1,j)) * h*h/3.0;
	}

	ret += 0.5 * pow2(f(0,n2-1)-v(0,n2-1)) * h*h/6.0;
	for(std::size_t i = 1; i < n1-1; i++) {
		ret += 0.5 * pow2(f(i,n2-1)-v(i,n2-1)) * h*h/2.0;
	}
	ret += 0.5 * pow2(f(n1-1,n2-1)-v(n1-1,n2-1)) * h*h/3.0;

	for(std::size_t j = 1; j <= n2-1; j++) {
		for(std::size_t i = 1; i <= 2*(n1-1); i++) {
			ret += 0.5 * r1 * dot(v.grad(i,j)-p1(i,j)) * h*h/2.0;
			ret += dot(l1(i,j), v.grad(i,j) - p1(i,j)) * h*h/2.0;

			ret += 0.5 * r2 * dot(p2(i,j)-p3(i,j)) * h*h/2.0;
			ret += dot(l2(i,j), p2(i,j)-p3(i,j)) * h*h/2.0;
		}
	}

	return ret;

}

template <typename T>
int augLangL1Helper(const T *input, T *output, std::size_t width, std::size_t height, std::size_t ldfw, T r1, T r2, T r3, T eps, Phi phi, double beta, T h, T delta, int maxIter, Debugger<T> debugger) {

	std::size_t s4n1 = 2*(width-1);
	std::size_t s4n2 = height-1;
	int err;

	Timer initTimer;
	initTimer.begin();
	
	Q<T> p1(s4n1, s4n2, h), p2(s4n1, s4n2, h), p3(s4n1, s4n2, h), l1(s4n1, s4n2, h), l2(s4n1, s4n2, h);
	V<T> u(width, height, h), f(width, height, h), psi(width, height, h), l3(width, height, h);

	// Initialization

#ifdef INFO
	std::cout << "augLangL1Helper(*, *, " << width << ", " << height << ", " << ldfw << ", " <<
			r1 << ", " << r2 << ", " << r3 << ", " << eps << ", " << maxIter << ")" << std::endl;

	std::cout << ">>> Initializing..." << std::endl;
#endif

	DebugHelper<T> helper(u, f, p1, p2, p3, psi, l1, l2, l3, eps, r1, r2, r3, phi, beta, 0);

#if DEBUG
	debugger.overallTests.add(0, PrintFValue<T>());

	debugger.subStep1Tests.add(0, PrintSubProblem1Values<T>());
	debugger.subStep2Tests.add(0, PrintSubProblem2Values<T>());
	debugger.subStep3Tests.add(0, PrintSubProblem3Values<T>());
	debugger.subStep4Tests.add(0, PrintSubProblem4Values<T>());
	debugger.lambdaTests.add(0, PrintLambdaValues<T>());

	debugger.stepTests.add(1, PrintStepProgress<T>());
	debugger.subStep1Tests.add(1, PrintSubProblemProgress<T>("Augmented Lagrangian change (1. subproblem)"));
	debugger.subStep2Tests.add(1, PrintSubProblemProgress<T>("Augmented Lagrangian change (2. subproblem)"));
	debugger.subStep3Tests.add(1, PrintSubProblemProgress<T>("Augmented Lagrangian change (3. subprobleme)"));
	debugger.subStep4Tests.add(1, PrintSubProblemProgress<T>("Augmented Lagrangian change (4. subproblem)"));
	debugger.lambdaTests.add(1, PrintSubProblemProgress<T>("Augmented Lagrangian change (lambda update)"));

	debugger.stepTests.add(0, PrintConstraints<T>());
#endif

#ifdef FULL_DEBUG
	debugger.subStep1Tests.add(-1, PrintSubProblem1Errors<T>());
	debugger.subStep2Tests.add(-1, PrintSubProblem2Errors<T>());
	debugger.subStep3Tests.add(-1, PrintSubProblem3Errors<T>());
	debugger.subStep4Tests.add(-1, PrintSubProblem4Errors<T>());
#endif

	pscrWrapper::PSCR<T> pscr = initSub4(width, height, width, h, &err);
	if(err) {
		std::cerr << ">>> Error while initializing the PSCR solver. Error code = " << err << std::endl;
		return err;
	}

	for(std::size_t j = 0; j < height; j++)
		for(std::size_t i = 0; i < width; i++)
			f(i,j) = input[j*ldfw+i];

	if(!debugger.getInitHelper().initializeAll(u, f, p1, p2, p3, psi, l1, l2, l3, helper)) {

		debugger.getInitHelper().initializeF(f, helper);

		if(!debugger.getInitHelper().initializeU(u, helper))
			u = 0;

		if(!debugger.getInitHelper().initializeP1(p1, helper))
			p1 = 0;

		if(!debugger.getInitHelper().initializeP2(p2, helper))
			p2 = 0;

		if(!debugger.getInitHelper().initializeP3(p3, helper))
			p3 = 0;

		if(!debugger.getInitHelper().initializePsi(psi, helper))
			psi = 0;

		if(!debugger.getInitHelper().initializeL1(l1, helper))
			l1 = 0.0;

		if(!debugger.getInitHelper().initializeL2(l2, helper))
			l2 = 0.0;

		if(!debugger.getInitHelper().initializeL3(l3, helper))
			l3 = 0.0;
	}
	
	initTimer.end();
	
	Timer solutionTimer;
	solutionTimer.begin();

	debugger.overallTests.runFirst(std::cout, helper);
	debugger.subStep1Tests.runFirst(std::cout, helper);
	debugger.subStep2Tests.runFirst(std::cout, helper);
	debugger.subStep3Tests.runFirst(std::cout, helper);
	debugger.subStep4Tests.runFirst(std::cout, helper);
	debugger.lambdaTests.runFirst(std::cout, helper);
	debugger.stepTests.runFirst(std::cout, helper);

	T last = augObjectFunction(psi, f, u, p1, l1, p2, l2, p3, l3, eps, r1, r2, r3, phi, beta);
	//V<T> lastU = u;

	double sub1TotalTime = 0.0;
	double sub2TotalTime = 0.0;
	double sub3TotalTime = 0.0;
	double sub4TotalTime = 0.0;
	double miscTotalTime = 0.0;
	
	int iter;
	for(iter = 1; iter <= maxIter; iter++) {
		helper.setIter(iter);

#if DEBUG
		std::cout << "========================================" << std::endl;
#endif
#ifdef INFO
		std::cout << ">>> Iteration " << iter << " / " << maxIter << std::endl;
#endif

		debugger.stepTests.runBefore(std::cout, helper);

		//
		// Subproblem 1
		//

		debugger.subStep1Tests.runBefore(std::cout, helper);

		Timer sub1Timer;
		sub1Timer.begin();
		
		err = solveSub1(p1, p2, u, p3, l1, l2, r1, r2, beta, h, 1.0E-8, iter <= 2);
		if(err) {
			std::cerr << ">>> Error while solving subproblem no. 1" << std::endl;
			return err;
		}
		
		sub1Timer.end();
		sub1TotalTime += sub1Timer.getTime();

		debugger.subStep1Tests.runAfter(std::cout, helper);

		//
		// Subproblem 2
		//

		debugger.subStep2Tests.runBefore(std::cout, helper);

		Timer sub2Timer;
		sub2Timer.begin();
		
		err = solveSub2(p3, psi, p2, l2, l3, r2, r3, h, 30, (T) 1.0E-3);
		if(err) {
			std::cerr << "Error while solving subproblem no. 2" << std::endl;
			return err;
		}

		sub2Timer.end();
		sub2TotalTime += sub2Timer.getTime();
		
		debugger.subStep2Tests.runAfter(std::cout, helper);

		//
		// Subproblem 3
		//

		debugger.subStep3Tests.runBefore(std::cout, helper);

		Timer sub3Timer;
		sub3Timer.begin();
		
		err = solveSub3(psi, p3, l3, r3, eps, phi, h);
		if(err) {
			std::cerr << ">>> Error while solving subproblem no. 3" << std::endl;
			return err;
		}

		sub3Timer.end();
		sub3TotalTime += sub3Timer.getTime();
		
		debugger.subStep3Tests.runAfter(std::cout, helper);

		//
		// Subproblem 4
		//

		debugger.subStep4Tests.runBefore(std::cout, helper);

		Timer sub4Timer;
		sub4Timer.begin();
		
		err = solveSub4(u, pscr, f, p1, l1, r1, h);
		if(err) {
			std::cerr << ">>> Error while solving subproblem no. 4" << std::endl;
			return err;
		}

		sub4Timer.end();
		sub4TotalTime += sub4Timer.getTime();
		
		debugger.subStep4Tests.runAfter(std::cout, helper);

		debugger.lambdaTests.runBefore(std::cout, helper);
		
		Timer miscTimer;
		miscTimer.begin();

		#pragma omp parallel for
		for(std::size_t j = 1; j <= s4n2; j++) {
			for(std::size_t i = 1; i <= s4n1; i++) {
				l1(i,j) += r1*(u.grad(i,j)-p1(i,j));
				l2(i,j) += r2*(p2(i,j)-p3(i,j));
			}
		}

		#pragma omp parallel for
		for(std::size_t j = 1; j < height-1; j++)
			for(std::size_t i = 1; i < width-1; i++)
				l3(i,j) += r3*(p3.div(i,j)-psi(i,j));
			
		miscTimer.end();
		miscTotalTime += miscTimer.getTime();

		debugger.lambdaTests.runAfter(std::cout, helper);

		debugger.stepTests.runAfter(std::cout, helper);

		miscTimer.begin();
		
		T current = augObjectFunction(psi, f, u, p1, l1, p2, l2, p3, l3, eps, r1, r2, r3, phi, beta);
		
		miscTimer.end();
		miscTotalTime += miscTimer.getTime();
		
		T ratio = ::fabs((last-current)/last);
		//T ratio = vectorNormInf(lastU-u)/vectorNormInf(lastU);

#if INFO
		std::cout << ">>> |last-current|/|last| = " << ratio << std::endl;
#endif

		if(ratio < delta && 5 <= iter) {
#ifdef INFO
			std::cout << ">>> Success." << std::endl;
#endif
			break;
		}

		last = current;
		//lastU = u;
	}

	debugger.overallTests.runAfter(std::cout, helper);

	debugger.overallTests.runLast(std::cout, helper);
	debugger.subStep1Tests.runLast(std::cout, helper);
	debugger.subStep2Tests.runLast(std::cout, helper);
	debugger.subStep3Tests.runLast(std::cout, helper);
	debugger.subStep4Tests.runLast(std::cout, helper);
	debugger.lambdaTests.runLast(std::cout, helper);
	debugger.stepTests.runLast(std::cout, helper);

	for(std::size_t j = 0; j < height; j++) {
		for(std::size_t i = 0; i < width; i++) {
			 output[j*ldfw+i] = u(i,j);
		}
	}
	
	solutionTimer.end();

#if DEBUG
	std::cout << ">>> Initialization time:             " << initTimer.getTime() << " s" << std::endl;
	std::cout << ">>> Full solution time:              " << solutionTimer.getTime() << " s" << std::endl;
	std::cout << ">>> Device solution time:            " << solutionTimer.getTime() << " s" << std::endl;
	std::cout << ">>> Average step time:               " << solutionTimer.getTime()/iter << " s" << std::endl;
	
	std::cout << ">>> Subproblem #1 total time:        " << sub1TotalTime << " s" << std::endl;
	std::cout << ">>> Subproblem #2 total time:        " << sub2TotalTime << " s" << std::endl;
	std::cout << ">>> Subproblem #3 total time:        " << sub3TotalTime << " s" << std::endl;
	std::cout << ">>> Subproblem #4 total time:        " << sub4TotalTime << " s" << std::endl;
	std::cout << ">>> Misc total time:                 " << miscTotalTime << " s" << std::endl;

	
	std::cout << ">>> Subproblem #1 average step time: " << sub1TotalTime/iter << " s" << std::endl;
	std::cout << ">>> Subproblem #2 average step time: " << sub2TotalTime/iter << " s" << std::endl;
	std::cout << ">>> Subproblem #3 average step time: " << sub3TotalTime/iter << " s" << std::endl;
	std::cout << ">>> Subproblem #4 average step time: " << sub4TotalTime/iter << " s" << std::endl;
	std::cout << ">>> Misc average step time:          " << miscTotalTime/iter << " s" << std::endl;
#endif
	
	return 0;
}

int augLagL1(const double *input, double *output, std::size_t width, std::size_t height, std::size_t ldfw, double r1, double r2, double r3, double eps, Phi phi, double beta, double h, double delta, int maxIt, const Debugger<double>& debugger) {
	return augLangL1Helper(input, output, width, height, ldfw, r1, r2, r3, eps, phi, beta, h, delta, maxIt, debugger);
}

//int augLagL1(const float *input, float *output, std::size_t width, std::size_t height, std::size_t ldfw, float r1, float r2, float r3, float eps, int maxIt) {
//
//	return augLangL1Helper(input, output, width, height, ldfw, r1, r2, r3, eps, maxIt);
//}

}
