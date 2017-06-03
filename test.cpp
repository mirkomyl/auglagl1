/*
 *  Created on: May 27, 2013
 *      Author: Mirko Myllykoski (mirko.myllykoski@gmail.com)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <stdexcept>
#include <math.h>
#include <limits>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>

#include "Image.h"
#include "augLagL1.h"

using namespace std;

template <typename T>
class ConstantU : public auglagL1::AbstractInitHelper<T> {
public:
	ConstantU(T c) {
		this->c = c;
	}
	virtual bool initializeU(auglagL1::V<T>& u, const auglagL1::DebugHelper<T>& helper) const {
		u = this->c;
		return true;
	}
	virtual auglagL1::AbstractInitHelper<T>* clone() const {
		return new ConstantU<T>(*this);
	}
private:
	T c;
};

template <typename T>
class RandU : public auglagL1::AbstractInitHelper<T> {
public:
	RandU(T c, T eps) {
		this->c = c;
		this->eps = eps;
	}
	virtual bool initializeU(auglagL1::V<T>& u, const auglagL1::DebugHelper<T>& helper) const {
		size_t n1 = u.getN1();
		size_t n2 = u.getN2();

		for(size_t j = 0; j < n2; j++)
			for(size_t i = 0; i < n1; i++)
				u(i,j) = c + 2*eps*((T)rand()/(T)RAND_MAX) - eps;
		return true;
	}
	virtual auglagL1::AbstractInitHelper<T>* clone() const {
		return new RandU<T>(*this);
	}
private:
	T c;
	T eps;
};

template <typename T>
class LinU : public auglagL1::AbstractInitHelper<T> {
public:
	virtual bool initializeU(auglagL1::V<T>& u, const auglagL1::DebugHelper<T>& helper) const {
		size_t n1 = u.getN1();
		size_t n2 = u.getN2();
		T h = u.getH();

		for(size_t j = 0; j < n2; j++)
			for(size_t i = 0; i < n1; i++)
				u(i,j) = (1.0*i*h + 1.0*j*h)/2.0;
		return true;
	}
	virtual auglagL1::AbstractInitHelper<T>* clone() const {
		return new LinU<T>(*this);
	}
};

template <typename T>
class UequalF : public auglagL1::AbstractInitHelper<T> {
public:
	virtual bool initializeU(auglagL1::V<T>& u, const auglagL1::DebugHelper<T>& helper) const {
		u=helper.f;
		return true;
	}
	virtual auglagL1::AbstractInitHelper<T>* clone() const {
		return new UequalF<T>(*this);
	}
};

template <typename T>
class NaturalInit : public auglagL1::AbstractInitHelper<T> {
public:
	virtual bool initializeAll(auglagL1::V<T>& u, auglagL1::V<T>& f, auglagL1::Q<T>& p1, auglagL1::Q<T>& p2, auglagL1::Q<T>& p3, auglagL1::V<T>& psi, auglagL1::Q<T>& l1, auglagL1::Q<T>& l2, auglagL1::V<T>& l3, const auglagL1::DebugHelper<T>& helper) const {
		u = helper.f;
		p1 = grad(u);

		std::size_t n1 = p1.getN1();
		std::size_t n2 = p1.getN2();

		for(std::size_t j = 1; j <= n2; j++)
			for(std::size_t i = 1; i <= n1; i++)
				p2(i,j) = p1(i,j)/sqrt(1.0+auglagL1::dot(p1(i,j)));

		p3 = p2;
		psi = div(p3);

		return true;
	}
	virtual auglagL1::AbstractInitHelper<T>* clone() const {
		return new NaturalInit(*this);
	}
};

template <typename T>
void special(Image<T>& image) {
	size_t width = image.getWidth();
	size_t height = image.getHeight();
	T h = 1.0/max(width, height);

	for(size_t j = 0; j < height; j++) {
		for(size_t i = 0; i < width; i++) {
			T x = 1.0*i*h;
			T y = 1.0*j*h;
			T R = pow2(x+1.0) + pow2(y+1.0);
			image(i,j) = 0.36*log(R+sqrt(R*R-0.36*0.36));
		}
	}
}

template <typename T>
class P1DiffToAnalyt : public auglagL1::AbstractTest<T> {
public:
	P1DiffToAnalyt() {
		firstValue = -1.0;
		lastValue = -1.0;
	}
	virtual void after(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		T res = calc(helper);
		if(firstValue < 0.0) firstValue = res;
		if(lastValue < 0.0) lastValue = res;
		out << ">>> | ex_p1 - p1 | = " << sqrt(res) << " (diff = " << res - lastValue << ", total = " << res - firstValue << ")" << std::endl;
		lastValue = res;
	}
	virtual auglagL1::AbstractTest<T>* clone() const {
		return new P1DiffToAnalyt<T>(*this);
	}
private:
	T calc(const auglagL1::DebugHelper<T>& helper) const {
		std::size_t n1 = helper.p1.getN1();
		std::size_t n2 = helper.p1.getN2();
		T h = helper.p1.getH();

		T res = 0;

		for(std::size_t j = 1; j <= n2; j++) {
			for(std::size_t i = 1; i <= n1; i+=2) {
				T x1 = ((1.0*i+1)/2-(2.0/3.0))*h;
				T y1 = (1.0*j-(1.0/3.0))*h;
				T R1 = pow2(x1+1) + pow2(y1+1);
				auglagL1::R2<T> real1 = 0.36*auglagL1::R2<T>(2*x1+2, 2*y1+2)/sqrt(R1*R1-0.36*0.36);
				res += dot(helper.p1(i,j)-real1);

				T x2 = ((1.0*i+1)/2-(1.0/3.0))*h;
				T y2 = (1.0*j-(2.0/3.0))*h;
				T R2 = pow2(x2+1) + pow2(y2+1);
				auglagL1::R2<T> real2 = 0.36*auglagL1::R2<T>(2*x2+2, 2*y2+2)/sqrt(R2*R2-0.36*0.36);
				res += dot(helper.p1(i+1,j)-real2);
			}
		}
		return sqrt(res)/norm(helper.p1);
	}
	T firstValue;
	T lastValue;
};

template <typename T>
class ChangeTracker : public auglagL1::AbstractTest<T> {
public:
	ChangeTracker(Image<T>& _img) : image(_img) {}
	virtual void first(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		std::size_t n1 = helper.u.getN1();
		std::size_t n2 = helper.u.getN2();
		T h = helper.u.getH();

		change = auglagL1::V<T>(n1, n2, h);
		change = 0.0;
	}
	virtual void before(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		befor = helper.u;
	}
	virtual void after(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		change += fabs(befor - helper.u);
	}
	virtual void last(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		std::size_t n1 = change.getN1();
		std::size_t n2 = change.getN2();

		T max = 0.0;
		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1; i++)
				max = std::max(max, change(i,j));

		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1; i++)
				image(i,j) = change(i,j) / max;
	}
	virtual auglagL1::AbstractTest<T>* clone() const {
		return new ChangeTracker<T>(*this);
	}
private:
	auglagL1::V<T> change;
	auglagL1::V<T> befor;
	Image<T>& image;
};

template <typename T>
class PsiTracker : public auglagL1::AbstractTest<T> {
public:
	PsiTracker(Image<T>& _img) : image(_img) {}
	virtual void last(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		std::size_t n1 = helper.u.getN1();
		std::size_t n2 = helper.u.getN2();

		T max = 0.0;
		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1; i++)
				max = std::max(max, ::fabs(helper.psi(i,j)));

		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1; i++)
				image(i,j) = helper.psi(i,j) / max;
	}
	virtual auglagL1::AbstractTest<T>* clone() const {
		return new PsiTracker<T>(*this);
	}
private:
	Image<T>& image;
};

template <typename T>
class P1Tracker : public auglagL1::AbstractTest<T> {
public:
	P1Tracker(Image<T>& _img) : image(_img) {}
	virtual void last(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		std::size_t n1 = helper.p1.getN1();
		std::size_t n2 = helper.p1.getN2();

		T max = 0.0;
		for(std::size_t j = 1; j <= n2; j++)
			for(std::size_t i = 1; i <= n1; i++)
				max = std::max(max, norm(helper.p1(i,j)));

		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1/2; i++)
				image(i,j) = (norm(helper.p1(2*(i+1),j+1)) + norm(helper.p1(2*(i+1)+1,j+1)))/(2*max);
	}
	virtual auglagL1::AbstractTest<T>* clone() const {
		return new P1Tracker<T>(*this);
	}
private:
	Image<T>& image;
};

template <typename T>
class P2Tracker : public auglagL1::AbstractTest<T> {
public:
	P2Tracker(Image<T>& _img) : image(_img) {}
	virtual void last(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		std::size_t n1 = helper.p2.getN1();
		std::size_t n2 = helper.p2.getN2();

		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1/2; i++)
				image(i,j) = (norm(helper.p2(2*(i+1),j+1)) + norm(helper.p2(2*(i+1)+1,j+1)))/2;
	}
	virtual auglagL1::AbstractTest<T>* clone() const {
		return new P2Tracker<T>(*this);
	}
private:
	Image<T>& image;
};

template <typename T>
class P3Tracker : public auglagL1::AbstractTest<T> {
public:
	P3Tracker(Image<T>& _img) : image(_img) {}
	virtual void last(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		std::size_t n1 = helper.p3.getN1();
		std::size_t n2 = helper.p3.getN2();

		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1/2; i++)
				image(i,j) = (norm(helper.p3(2*(i+1),j+1)) + norm(helper.p3(2*(i+1)+1,j+1)))/2;
	}
	virtual auglagL1::AbstractTest<T>* clone() const {
		return new P3Tracker<T>(*this);
	}
private:
	Image<T>& image;
};
/*
template <typename T>
class P1VSP2Tracker : public auglagL1::AbstractTest<T> {
public:
	P1VSP2Tracker(Image<T>& _img) : image(_img) {}
	virtual void last(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		std::size_t n1 = helper.p1.getN1();
		std::size_t n2 = helper.p1.getN2();

		T max = 0.0;
		for(std::size_t j = 1; j <= n2; j++)
			for(std::size_t i = 1; i <= n1; i++)
				max = std::max(max, norm(helper.p1(i,j)));

		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1/2; i++)
				image(i,j) = (norm(helper.p1(2*(i+1),j+1)) + norm(helper.p1(2*(i+1)+1,j+1)))/(2*max);
	}
	virtual auglagL1::AbstractTest<T>* clone() const {
		return new P1VSP2Tracker<T>(*this);
	}
private:
	Image<T>& image;
};
*/
template <typename T>
class ContributionTracker : public auglagL1::AbstractTest<T> {
public:
	ContributionTracker(Image<T>& _org, Image<T>& _orgA, Image<T>& _orgB,
			Image<T>& _img, Image<T>& _imgA, Image<T>& _imgB) :
				original(_org), originalA(_orgA), originalB(_orgB),
				result(_img), resultA(_imgA), resultB(_imgB) {}
	virtual void last(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		calc(original, helper.f, helper.f, helper.eps, 1.0, 1.0);
		calc(result, helper.u, helper.f, helper.eps, 1.0, 1.0);

		calc(originalA, helper.f, helper.f, helper.eps, 1.0, 0.0);
		calc(resultA, helper.u, helper.f, helper.eps, 1.0, 0.0);

		calc(originalB, helper.f, helper.f, helper.eps, 0.0, 1.0);
		calc(resultB, helper.u, helper.f, helper.eps, 0.0, 1.0);
	}
	virtual auglagL1::AbstractTest<T>* clone() const {
		return new ContributionTracker<T>(*this);
	}
private:
	void calc(Image<T>& img, const auglagL1::V<T>& u, const auglagL1::V<T>& f, T eps, T w1, T w2) {
		std::size_t n1 = u.getN1();
		std::size_t n2 = u.getN2();
		T h = u.getH();

		auglagL1::Q<T> vv(2*(n1-1), n2-1, h);

		for(std::size_t j = 1; j <= n2-1; j++)
			for(std::size_t i = 1; i <= 2*(n1-1); i++)
				vv(i,j) = u.grad(i,j)/sqrt((1.0+dot(u.grad(i,j))));

		auglagL1::V<T> div_vv(n1,n2,h);

		for(std::size_t i = 0; i < n1; i++)
			div_vv(i,0) = div_vv(i,n2-1) = 0.0;

		for(std::size_t j = 1; j < n2-1; j++) {
			div_vv(0,j) = div_vv(n1-1,j) = 0.0;
			for(std::size_t i = 1; i < n1-1; i++) {
				div_vv(i,j) = vv.div(i,j);
			}
		}

		T max = 0.0;
		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1; i++)
				max = std::max(max, w1 * eps * ::fabs(div_vv(i,j)) + w2 * (1.0/2.0) * pow2(u(i,j)-f(i,j)));

		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1; i++)
				img(i,j) = (w1 * eps * ::fabs(div_vv(i,j)) + (1.0/2.0) * w2 * pow2(u(i,j)-f(i,j)))/max;

	}
	Image<T>& original;
	Image<T>& originalA;
	Image<T>& originalB;
	Image<T>& result;
	Image<T>& resultA;
	Image<T>& resultB;
};

template <typename T>
class CompareToOriginal : public auglagL1::AbstractTest<T> {
public:
	CompareToOriginal(const Image<T>& _org, T h) {
		std::size_t n1 = _org.getWidth();
		std::size_t n2 = _org.getHeight();

		u = auglagL1::V<T>(n1, n2, h);
		for(std::size_t j = 0; j < n2; j++)
			for(std::size_t i = 0; i < n1; i++)
				u(i,j) = _org(i,j);

		p1 = auglagL1::Q<T>(2*(n1-1), n2-1, h);
		for(std::size_t j = 1; j <= n2-1; j++)
			for(std::size_t i = 1; i <= 2*(n1-1); i++)
				p1(i,j) = u.grad(i,j);

		p3 = auglagL1::Q<T>(2*(n1-1), n2-1, h);
		for(std::size_t j = 1; j <= n2-1; j++)
			for(std::size_t i = 1; i <= 2*(n1-1); i++)
				p3(i,j) = p1(i,j)/sqrt(1.0+dot(p1(i,j)));

		psi = auglagL1::V<T>(n1, n2, h);
		for(std::size_t i = 0; i < n1; i++)
			psi(i,0) = psi(i,n2-1) = 0.0;

		for(std::size_t j = 1; j < n2-1; j++) {
			psi(0,j) = psi(n1-1,j) = 0.0;
			for(std::size_t i = 1; i < n1-1; i++) {
				psi(i,j) = p3.div(i,j);
			}
		}
	}
	virtual void first(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		firstU = vectorNorm2(helper.u-u)/vectorNorm2(u);
		firstP1 = vectorNorm2(helper.p1-p1)/vectorNorm2(p1);
		firstP3 =  vectorNorm2(helper.p3-p3)/vectorNorm2(p3);
		firstPsi = vectorNorm2(helper.psi-psi)/vectorNorm2(psi);
	}
	virtual void after(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		T newU = vectorNorm2(helper.u-u)/vectorNorm2(u);
		T newP1 = vectorNorm2(helper.p1-p1)/vectorNorm2(p1);
		T newP3 =  vectorNorm2(helper.p3-p3)/vectorNorm2(p3);
		T newPsi = vectorNorm2(helper.psi-psi)/vectorNorm2(psi);

		out << ">>> | u - _u | / |_u| = " << newU << " (" << newU - lastU << ", total = " << newU - firstU <<  ")" << std::endl;
		out << ">>> | p1 - _p1 | / |_p1| = " << newP1 << " (" << newP1 - lastP1 << ", total = " << newP1 - firstP1 <<  ")" << std::endl;
		out << ">>> | p3 - _p3 | / |_p3| = " << newP3 << " (" << newP3 - lastP3 << ", total = " << newP3 - firstP3 <<  ")" << std::endl;
		out << ">>> | psi - _psi | / |_psi| = " << newPsi << " (" << newPsi - lastPsi << ", total = " << newPsi - firstPsi <<  ")" << std::endl;

		lastU = newU;
		lastP1 = newP1;
		lastP3 = newP3;
		lastPsi = newPsi;
	}
	virtual auglagL1::AbstractTest<T>* clone() const {
		return new CompareToOriginal<T>(*this);
	}
private:
	auglagL1::V<T> u;
	auglagL1::Q<T> p1;
	auglagL1::Q<T> p3;
	auglagL1::V<T> psi;
	T lastU;
	T lastP1;
	T lastP3;
	T lastPsi;
	T firstU;
	T firstP1;
	T firstP3;
	T firstPsi;
};

template <typename T>
class Stages : public auglagL1::AbstractTest<T> {
public:
	Stages(std::vector<Image<T> > &_stages) : stages(_stages) { }
	virtual ~Stages() {}
	virtual void first(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		add(helper);
	}	
	virtual void after(std::ostream& out, const auglagL1::DebugHelper<T>& helper) {
		add(helper);
	}
	virtual Stages* clone() const {
		return new Stages(*this);
	}
private:
	void add(const auglagL1::DebugHelper<T>& helper) {
		stages.push_back(Image<T>(helper.u.getN1(), helper.u.getN2()));
		for(size_t j = 0; j < helper.u.getN1(); j++)
			for(size_t i = 0; i < helper.u.getN2(); i++)
				stages.back()(j,i) = helper.u(i+1,j+1);
	}
	std::vector<Image<T> > &stages;
};

void imageTest(std::string name, std::string _org, std::string _noise, double r0, double delta, int maxIt) {
	
	Image<double> org(_org);
	Image<double> noise(_noise);
	
	Image<double> noisy = org;
	
	for(size_t j = 0; j < noisy.getHeight(); j++)
		for(size_t i = 0; i < noisy.getWidth(); i++)
			noisy(i,j) = org(i,j) + noise(i,j) - 0.5;
	
	Image<double> result(noisy.getWidth(), noisy.getHeight());
	Image<double> diff(noisy.getWidth(), noisy.getHeight());

// 	Image<double> change(noisy.getWidth(), noisy.getHeight());

	//Image<double> contr_org(noisy.getWidth(), noisy.getHeight());
	//Image<double> contr_orgA(noisy.getWidth(), noisy.getHeight());
	//Image<double> contr_orgB(noisy.getWidth(), noisy.getHeight());

	//Image<double> contr(noisy.getWidth(), noisy.getHeight());
	//Image<double> contrA(noisy.getWidth(), noisy.getHeight());
	//Image<double> contrB(noisy.getWidth(), noisy.getHeight());

	Image<double> p1(noisy.getWidth()-1, noisy.getHeight()-1);
	Image<double> p2(noisy.getWidth()-1, noisy.getHeight()-1);
	Image<double> p3(noisy.getWidth()-1, noisy.getHeight()-1);
	Image<double> psi(noisy.getWidth(), noisy.getHeight());

	double h = 0.005;//1.0/max(noisy.getWidth()-1, noisy.getHeight()-1);
	
	auglagL1::Phi phi = auglagL1::abs;
	//auglagL1::Phi phi = auglagL1::square;
	double beta = BETA;
	
	//  5% noise => r0 = 0.001
	// 10% noise => r0 = 0.002
	// 20% noise => r0 = 0.004 - 0.005
	// 40% noise => r0 = 0.020

// abs, beta = 1
// 	double eps = 1.0 * r0 * h;
// 	double r1 = 5.0 * r0 * h;
// 	double r2 = 5.0 * r0;
// 	double r3 = 5.0 * r0 * h*h;

// abs, beta = 0.00001, DOES NOW WORK!
	double eps = 1.0 * r0 * h;
	double r1 = 50.0 * r0 * h;
	double r2 = 5.0 * r0;
	double r3 = 5.0 * r0 * h*h;
	
// square, beta = 0.00001
// 	double eps = 1.0 * r0 * h*h;
// 	double r1 = 50.0 * r0 * h;
// 	double r2 = 5.0 * r0;
// 	double r3 = 5.0 * r0 * h*h;

	std::cout << ">>> File: " << _org << " (" << org.getHeight() << " x " << org.getWidth() << ")" << std::endl;
	std::cout << ">>> r0 = " << r0 << std::endl;

	auglagL1::Debugger<double> debugger;
	
// 	std::vector<Image<double> > stages;
// 	Stages<double> stageTracker(stages);
// 	debugger.stepTests.add(100, stageTracker);

	//debugger.setInitHelper(NaturalInit<double>());

// 	debugger.stepTests.add(0, ChangeTracker<double>(change));
//	debugger.stepTests.add(3, CompareToOriginal<double>(org, h));
	debugger.overallTests.add(0, PsiTracker<double>(psi));
	debugger.overallTests.add(0, P1Tracker<double>(p1));
	debugger.overallTests.add(0, P2Tracker<double>(p2));
	debugger.overallTests.add(0, P3Tracker<double>(p3));
	//debugger.overallTests.add(0, ContributionTracker<double>(contr_org, contr_orgA, contr_orgB,
	//		contr, contrA, contrB));

	auglagL1::augLagL1(noisy.getRaw(), result.getRaw(), noisy.getWidth(), noisy.getHeight(), noisy.getWidth(), r1, r2, r3, eps, phi, beta, h, delta, maxIt, debugger);

	std::cout << ">>> Original norm = " << diffNorm(noisy, org) << std::endl;
	std::cout << ">>> New norm      = " << diffNorm(result, org) << std::endl;

	double maxDiff = 0.0;
	for(size_t j = 0; j < noisy.getHeight(); j++)
		for(size_t i = 0; i < noisy.getWidth(); i++)
			maxDiff = max(maxDiff, fabs(result(i,j) - noisy(i,j)));
	
	for(size_t j = 0; j < noisy.getHeight(); j++)
		for(size_t i = 0; i < noisy.getWidth(); i++)
			diff(i,j) =  fabs(result(i,j) - noisy(i,j)) / maxDiff;

	p1.save(name + "_p1.png", 0.0, 1.0, 0.05);
	p2.save(name + "_p2.png", 0.0, 1.0, 0.05);
	p3.save(name + "_p3.png");
	psi.save(name + "_psi.png", 0.0, 0.0, 0.0);

	//contr_org.save("contr_org.png");
	//contr_orgA.save("contr_orgA.png");
	//contr_orgB.save("contr_orgB.png");
	//contr.save("contr.png");
	//contrA.save("contrA.png");
	//contrB.save("contrB.png");

	diff.save(name + "_diff.png");
	org.save(name + "_original.png");
	noisy.save(name + "_noise.png", 0.5);
	result.save(name + "_output.png", 0.5);
// 	change.save(name + "_change.png");

	std::size_t sCount = 10;

	for(std::size_t i = 0; i < sCount; i++) {
		std::ofstream ofs;
		std::string fileName = name + "_section_org" + auglagL1::toString(i) + ".txt";
		ofs.open(fileName.c_str(), std::ofstream::out);
		crossSection(org, ofs, 1, ((1.0/(sCount+1))*(i+1))*result.getHeight(), result.getWidth()-1, ((1.0/(sCount+1))*(i+1))*result.getHeight(), result.getWidth()-2);
		ofs.close();
	}

	for(std::size_t i = 0; i < sCount; i++) {
		std::ofstream ofs;
		std::string fileName = name + "_section_noise" + auglagL1::toString(i) + ".txt";
		ofs.open(fileName.c_str(), std::ofstream::out);
		crossSection(noisy, ofs, 1, ((1.0/(sCount+1))*(i+1))*result.getHeight(), result.getWidth()-1, ((1.0/(sCount+1))*(i+1))*result.getHeight(), result.getWidth()-2);
		ofs.close();
	}

	for(std::size_t i = 0; i < sCount; i++) {
		std::ofstream ofs;
		std::string fileName = name + "_section" + auglagL1::toString(i) + ".txt";
		ofs.open(fileName.c_str(), std::ofstream::out);
		crossSection(result, ofs, 1, ((1.0/(sCount+1))*(i+1))*result.getHeight(), result.getWidth()-1, ((1.0/(sCount+1))*(i+1))*result.getHeight(), result.getWidth()-2);
		ofs.close();
	}
	
// 	for(int i = 0; i < (int) stages.size(); i++) {
// 		std::string idx;
// 		if(i < 10)
// 			idx = "00" + auglagL1::toString(i);
// 		else if(i < 100)
// 			idx = "0" + auglagL1::toString(i);
// 		else
// 			idx = auglagL1::toString(i);
// 		stages.at(i).save(name + "_step_" + idx + ".png", 0.5);
// 	}
}

int main(void) {
	srand (time(NULL));

	imageTest(NAME, ORG, NOISE, R0, 0.0, 300);

	return EXIT_SUCCESS;
}
