/*
 *  Created on: May 27, 2013
 *      Author: Mirko Myllykoski (mirko.myllykoski@gmail.com)
 */

#ifndef PSCR_WRAPPER_H
#define PSCR_WRAPPER_H

#include <cstddef>

namespace pscrWrapper {

int pow2(int a) {
	return 1<<a;
}

int pow4(int a) {
	return 1<<(2*a);
}

int log4(int a) {
	return log(a)/log(4.0);
}

template <typename T>
T max(T a,T b) {
	return (a) > (b) ? (a) : (b);
}

template <typename T>
T min(T a,T b) {
	return (a) < (b) ? (a) : (b);
}

// 2D Fortran PSCR solver (dc2d.f)
extern"C" {
void dc2d_(int *n1, int *n2, double *f, int *ldf, double *a1, double *b1,
		double *d1, double *a2, double *b2, double *d2, double *ch, double *dw,
		int *ldw, int *iw, int *liw, int *init, int *ierr);
}

template <typename T>
class PSCR {  };

template <>
class PSCR<double> {
public:
	PSCR() {
		floatWork = 0;
		intWork = 0;
		floatWorkSize = 0;
		intWorkSize = 0;
		a1Diag = a1OffDiag = m1Diag = 0;
		a2Diag = a2OffDiag = m2Diag = 0;
		n1 = n2 = ldf = 0;
		initialized = false;
	}
	PSCR(
			const double *_a1Diag, const double *_a1OffDiag,
			const double *_a2Diag, const double *_a2OffDiag,
			const double *_m1Diag, const double *_m2Diag,
			int _n1, int _n2, int _ldf, int *err) {

		if(_n1 < 3 || _n2 < 3 || _ldf < _n2)
			throw std::exception();

		floatWork = 0;
		intWork = 0;
		floatWorkSize = 0;
		intWorkSize = 0;
		a1Diag = a1OffDiag = m1Diag = 0;
		a2Diag = a2OffDiag = m2Diag = 0;
		n1 = n2 = ldf = 0;
		initialized = false;

		malloc(_n1, _n2, _ldf);

		for(int i = 0; i < n1; i++) {
			a1Diag[i] = _a1Diag[i];
			a1OffDiag[i] = _a1OffDiag[i];
			m1Diag[i] = _m1Diag[i];
		}

		for(int i = 0; i < n2; i++) {
			a2Diag[i] = _a2Diag[i];
			a2OffDiag[i] = _a2OffDiag[i];
			m2Diag[i] = _m2Diag[i];
		}

		*err = dc2d(n1, n2, 0, ldf,
			a1OffDiag, a1Diag, m1Diag,
			a2OffDiag, a2Diag, m2Diag,
			0, floatWork, floatWorkSize, intWork, intWorkSize, true);

		if(*err) {
			free();
		} else {
			initialized = true;
		}
	}
	PSCR(const PSCR& old) {
		floatWork = 0;
		intWork = 0;
		floatWorkSize = 0;
		intWorkSize = 0;
		a1Diag = a1OffDiag = m1Diag = 0;
		a2Diag = a2OffDiag = m2Diag = 0;
		n1 = n2 = ldf = 0;
		initialized = false;

		if(old.initialized) {
			malloc(old.n1, old.n2, old.ldf);
			copy(old);
			this->initialized = true;
		}
	}
	PSCR& operator=(const PSCR& a) {
		if(this == &a)
			return *this;

		if(a.initialized) {
			malloc(a.n1, a.n2, a.ldf);
			copy(a);
			initialized = true;
		} else {
			free();
		}

		return *this;
	}
	~PSCR() {
		free();
	}
	int run(double *f, double ch) const {
		if(!initialized)
			throw std::exception();

		return dc2d(n1, n2, f, ldf,
			a1OffDiag, a1Diag, m1Diag,
			a2OffDiag, a2Diag, m2Diag,
			ch, floatWork, floatWorkSize, intWork, intWorkSize, false);
	}
	int getN1() const { return n1; }
	int getN2() const { return n2; }
	int getLdf() const { return ldf; }
private:
	void free() {
		if(a1Diag)
			delete [] a1Diag;
		if(a1OffDiag)
			delete [] a1OffDiag;
		if(m1Diag)
			delete [] m1Diag;
		if(a2Diag)
			delete [] a2Diag;
		if(a2OffDiag)
			delete [] a2OffDiag;
		if(m2Diag)
			delete [] m2Diag;
		if(floatWork)
			delete [] floatWork;
		if(intWork)
			delete [] intWork;

		floatWork = 0;
		intWork = 0;
		a1Diag = a1OffDiag = m1Diag = 0;
		a2Diag = a2OffDiag = m2Diag = 0;
		n1 = n2 = ldf = 0;
		initialized = false;
	}
	void malloc(int _n1, int _n2, int _ldf) {
		int nl = 1 + max((int)(log4(_n1)),0);
		int newFloatWorkSize = 6*nl*_n1 + max(9*_n1, 10*_n2);
		int newIntWorkSize = (pow4(nl) - 1)/3 + 2*nl + 5*_n1 + 4;

		if(floatWorkSize != newFloatWorkSize) {
			if(floatWork)
				delete [] floatWork;

			try {
				floatWork = new double[newFloatWorkSize];
			} catch (...) {
				free();
				throw;
			}

			floatWorkSize = newFloatWorkSize;
		}

		if(intWorkSize != newIntWorkSize) {
			if(intWork)
				delete [] intWork;

			try {
				intWork = new int[newIntWorkSize];
			} catch (...) {
				free();
				throw;
			}

			intWorkSize = newIntWorkSize;
		}

		if(n1 != _n1) {
			if(a1Diag)
				delete [] a1Diag;
			if(a1OffDiag)
				delete [] a1OffDiag;
			if(m1Diag)
				delete [] m1Diag;

			a1Diag = a1OffDiag = m1Diag = 0;

			try {
				a1Diag = new double[_n1];
				a1OffDiag = new double[_n1];
				m1Diag = new double[_n1];
			} catch (...) {
				free();
				throw;
			}

			n1 = _n1;
		}

		if(n2 != _n2) {
			if(a2Diag)
				delete [] a2Diag;
			if(a2OffDiag)
				delete [] a2OffDiag;
			if(m2Diag)
				delete [] m2Diag;

			a2Diag = a2OffDiag = m2Diag = 0;

			try {
				a2Diag = new double[_n2];
				a2OffDiag = new double[_n2];
				m2Diag = new double[_n2];
			} catch (...) {
				free();
				throw;
			}

			n2 = _n2;
		}

		ldf = _ldf;
	}
	void copy(const PSCR& a) {
		for(int i = 0; i < a.floatWorkSize; i++)
			floatWork[i] = a.floatWork[i];

		for(int i = 0; i < a.intWorkSize; i++)
			intWork[i] = a.intWork[i];

		for(int i = 0; i < n1; i++) {
			a1Diag[i] = a.a1Diag[i];
			a1OffDiag[i] = a.a1OffDiag[i];
			m1Diag[i] = a.m1Diag[i];
		}

		for(int i = 0; i < n2; i++) {
			a2Diag[i] = a.a2Diag[i];
			a2OffDiag[i] = a.a2OffDiag[i];
			m2Diag[i] = a.m2Diag[i];
		}
	}
	int dc2d(int n1, int n2, double *f, int ldf, const double *a1, const double *b1,
		const double *d1, const double *a2, const double *b2, const double *d2, double ch, const double *dw,
		int ldw, const int *iw, int liw, int init) const {

		if(init)
			throw std::exception();

		int err;

		dc2d_(&n1, &n2, const_cast<double*>(f), &ldf, const_cast<double*>(a1), const_cast<double*>(b1),
				const_cast<double*>(d1), const_cast<double*>(a2), const_cast<double*>(b2), const_cast<double*>(d2), &ch, const_cast<double*>(dw),
				&ldw, const_cast<int*>(iw), &liw, &init, &err);

		return err;
	}
	int dc2d(int n1, int n2, double *f, int ldf, const double *a1, const double *b1,
		const double *d1, const double *a2, const double *b2, const double *d2, double ch, double *dw,
		int ldw, int *iw, int liw, int init) {

		int err;
		dc2d_(&n1, &n2, const_cast<double*>(f), &ldf, const_cast<double*>(a1), const_cast<double*>(b1),
						const_cast<double*>(d1), const_cast<double*>(a2), const_cast<double*>(b2), const_cast<double*>(d2), &ch, const_cast<double*>(dw),
						&ldw, const_cast<int*>(iw), &liw, &init, &err);

		return err;
	}

	bool initialized;

	int n1;
	int n2;
	int ldf;

	int floatWorkSize;
	int intWorkSize;
	double *floatWork;
	int *intWork;

	double *a1Diag;
	double *a1OffDiag;
	double *a2Diag;
	double *a2OffDiag;
	double *m1Diag;
	double *m2Diag;
};

template <>
class PSCR<float> {
public:
	PSCR() { }
	PSCR(
			const float *a1Diag, const float *a1OffDiag,
			const float *a2Diag, const float *a2OffDiag,
			const float *m1Diag, const float *m2Diag,
			int n1, int n2, int ldf, int *err) {

		if(n1 < 3 || n2 < 3 || ldf < n2)
			throw std::exception();

		double *tmpA1Diag, *tmpA1OffDiag, *tmpM1Diag;
		double *tmpA2Diag, *tmpA2OffDiag, *tmpM2Diag;

		tmpA1Diag = tmpA1OffDiag = tmpM1Diag = 0;
		tmpA2Diag = tmpA2OffDiag = tmpM2Diag = 0;

		try {
			tmpA1Diag = new double[n1];
			tmpA1OffDiag = new double[n1];
			tmpM1Diag = new double[n1];
			tmpA2Diag = new double[n2];
			tmpA2OffDiag = new double[n2];
			tmpM2Diag = new double[n2];
		} catch (...) {
			if(tmpA1Diag)
				delete [] tmpA1Diag;
			if(tmpA1OffDiag)
				delete [] tmpA1OffDiag;
			if(tmpM1Diag)
				delete [] tmpM1Diag;
			if(tmpA2Diag)
				delete [] tmpA2Diag;
			if(tmpA2OffDiag)
				delete [] tmpA2OffDiag;
			if(tmpM2Diag)
				delete [] tmpM2Diag;
			throw;
		}

		for(int i = 0; i < n1; i++) {
			tmpA1Diag[i] = a1Diag[i];
			tmpA1OffDiag[i] = a1OffDiag[i];
			tmpM1Diag[i] = m1Diag[i];
		}

		for(int i = 0; i < n2; i++) {
			tmpA2Diag[i] = a2Diag[i];
			tmpA2OffDiag[i] = a2OffDiag[i];
			tmpM2Diag[i] = m2Diag[i];
		}

		handle = PSCR<double>(
				tmpA1Diag, tmpA1OffDiag, tmpA2Diag,
				tmpA2OffDiag, tmpM1Diag, tmpM2Diag,
				n1, n2, ldf, err);

		delete [] tmpA1Diag;
		delete [] tmpA1OffDiag;
		delete [] tmpM1Diag;
		delete [] tmpA2Diag;
		delete [] tmpA2OffDiag;
		delete [] tmpM2Diag;
	}
	int run(float *f, float ch) const {
		std::size_t n1 = handle.getN1();
		std::size_t n2 = handle.getN2();
		std::size_t ldf = handle.getLdf();

		double *tmpF = new double[n1*ldf];

		for(std::size_t i = 0; i < n1; i++)
			for(std::size_t j = 0; j < n2; j++)
				tmpF[i*ldf+j] = f[i*ldf+j];

		int err = handle.run(tmpF, (double)ch);

		for(std::size_t i = 0; i < n1; i++)
			for(std::size_t j = 0; j < n2; j++)
				tmpF[i*ldf+j] = f[i*ldf+j];

		delete [] tmpF;

		return err;
	}
	int getN1() const { return handle.getN1(); }
	int getN2() const { return handle.getN2(); }
	int getLdf() const { return handle.getLdf(); }
private:
	PSCR<double> handle;
};

}

#endif
