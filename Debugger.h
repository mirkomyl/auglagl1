/*
 *  Created on: May 27, 2013
 *      Author: Mirko Myllykoski (mirko.myllykoski@gmail.com)
 */

#ifndef AUGLAGL1_DEBUGGER_H
#define AUGLAGL1_DEBUGGER_H

#include <ostream>
#include <map>
#include "R2.h"
#include "Q.h"
#include "V.h"
#include "derivates.h"

namespace auglagL1 {
	
enum Phi { abs, square };

template <typename T>
class DebugHelper {
public:
	DebugHelper(const V<T>& _u, const V<T>& _f, const Q<T>& _p1, const Q<T>& _p2, const Q<T>& _p3, const V<T>& _psi, const Q<T>& _l1, const Q<T>& _l2, const V<T>& _l3, T _eps, T _r1, T _r2, T _r3, Phi _phi, T _beta, int _iter) :
		u(_u), f(_f), p1(_p1), p2(_p2), p3(_p3), psi(_psi), l1(_l1), l2(_l2), l3(_l3), eps(_eps), r1(_r1), r2(_r2), r3(_r3), phi(_phi), beta(_beta) {
		iter = _iter;
	}
	void setIter(int _iter) {
		iter = _iter;
	}
	const V<T>& u;
	const V<T>& f;
	const Q<T>& p1;
	const Q<T>& p2;
	const Q<T>& p3;
	const V<T>& psi;
	const Q<T>& l1;
	const Q<T>& l2;
	const V<T>& l3;
	const T eps;
	const T r1;
	const T r2;
	const T r3;
	const Phi phi;
	const T beta;
	int iter;
private:
	DebugHelper(const DebugHelper &old) {}
};

template <typename T>
class AbstractInitHelper {
public:
	virtual ~AbstractInitHelper() {}
	virtual bool initializeU(V<T>& u, const DebugHelper<T>& helper) const { return false; }
	virtual bool initializeF(V<T>& f, const DebugHelper<T>& helper) const { return false; }
	virtual bool initializeP1(Q<T>& p1, const DebugHelper<T>& helper) const { return false; }
	virtual bool initializeP2(Q<T>& p2, const DebugHelper<T>& helper) const { return false; }
	virtual bool initializeP3(Q<T>& p3, const DebugHelper<T>& helper) const { return false; }
	virtual bool initializePsi(V<T>& psi, const DebugHelper<T>& helper) const { return false; }
	virtual bool initializeL1(Q<T>& l1, const DebugHelper<T>& helper) const { return false; }
	virtual bool initializeL2(Q<T>& l2, const DebugHelper<T>& helper) const { return false; }
	virtual bool initializeL3(V<T>& l3, const DebugHelper<T>& helper) const { return false; }
	virtual bool initializeAll(V<T>& u, V<T>& f, Q<T>& p1, Q<T>& p2, Q<T>& p3, V<T>& psi, Q<T>& l1, Q<T>& l2, V<T>& l3, const DebugHelper<T>& helper) const { return false; }
	virtual AbstractInitHelper<T>* clone() const = 0;
};

template <typename T>
class NullInitHelper : public AbstractInitHelper<T> {
public:
	virtual AbstractInitHelper<T>* clone() const {
		return new NullInitHelper<T>(*this);
	}
};

template <typename T>
class AbstractTest {
public:
	virtual ~AbstractTest() {}
	virtual void first(std::ostream& out, const DebugHelper<T>& helper) {}
	virtual void before(std::ostream& out, const DebugHelper<T>& helper) {}
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {}
	virtual void last(std::ostream& out, const DebugHelper<T>& helper) {}
	virtual AbstractTest* clone() const = 0;
};

template <typename T>
class Debugger {
public:
	class Container {
	public:
		Container() {}
		Container(const Container& old) {
			const_iterator it;
			for(it = old._t.begin(); it != old._t.end(); it++)
				_t.insert(std::pair<int, AbstractTest<T>*>(it->first, it->second->clone()));
		}
		Container& operator=(const Container& a) {
			if(&a == this)
				return *this;

			iterator it;
			for(it = _t.begin(); it != _t.end(); it++)
				delete it->second;

			_t.clear();

			const_iterator it2;
			for(it2 = a._t.begin(); it2 != a._t.end(); it2++)
				_t.insert(std::pair<int, AbstractTest<T>*>(it2->first, it2->second->clone()));

			return *this;
		}
		~Container() {
			iterator it;
			for(it = _t.begin(); it != _t.end(); it++)
				delete it->second;

			_t.clear();
		}

		void add(int weight, const AbstractTest<T>& test) {
			_t.insert(std::pair<int, AbstractTest<T>*>(weight, test.clone()));
		}

		void runFirst(std::ostream& out, const DebugHelper<T>& helper) {
			iterator it;
			for(it = _t.begin(); it != _t.end(); it++)
				(it->second)->first(out, helper);
		}

		void runBefore(std::ostream& out, const DebugHelper<T>& helper) {
			iterator it;
			for(it = _t.begin(); it != _t.end(); it++)
				(it->second)->before(out, helper);
		}

		void runAfter(std::ostream& out, const DebugHelper<T>& helper) {
			iterator it;
			for(it = _t.begin(); it != _t.end(); it++)
				(it->second)->after(out, helper);
		}
		void runLast(std::ostream& out, const DebugHelper<T>& helper) {
			iterator it;
			for(it = _t.begin(); it != _t.end(); it++)
				(it->second)->last(out, helper);
		}
	private:
		typedef typename std::multimap<int, AbstractTest<T>* >::iterator iterator;
		typedef typename std::multimap<int, AbstractTest<T>* >::const_iterator const_iterator;
		std::multimap<int, AbstractTest<T>*> _t;
	};
	Debugger() {
		this->initHelper = new NullInitHelper<T>();
	}
	Debugger(const Debugger<T>& old) {
		
		this->initHelper = old.getInitHelper().clone();
		this->overallTests = old.overallTests;
		this->stepTests = old.stepTests;
		this->subStep1Tests = old.subStep1Tests;
		this->subStep2Tests = old.subStep2Tests;
		this->subStep3Tests = old.subStep3Tests;
		this->subStep4Tests = old.subStep4Tests;
		this->lambdaTests = old.lambdaTests;
	}
	Debugger<T>& operator=(const Debugger<T>& a) {
		if(&a == this)
			return *this;

		this->setInitHelper(a.getInitHelper());

		this->overallTests = a.overallTests;
		this->stepTests = a.stepTests;
		this->subStep1Tests = a.subStep1Tests;
		this->subStep2Tests = a.subStep2Tests;
		this->subStep3Tests = a.subStep3Tests;
		this->subStep4Tests = a.subStep4Tests;
		this->lambdaTests = a.lambdaTests;

		return *this;
	}
	~Debugger() {
		delete this->initHelper;
	}

	void setInitHelper(const AbstractInitHelper<T>& helper) {
		delete this->initHelper;
		this->initHelper = helper.clone();
	}
	const AbstractInitHelper<T>& getInitHelper() const {
		return *initHelper;
	}
	
	Container overallTests;
	Container stepTests;
	Container subStep1Tests;
	Container subStep2Tests;
	Container subStep3Tests;
	Container subStep4Tests;
	Container lambdaTests;
private:
	AbstractInitHelper<T>* initHelper;
};

template <typename T>
class PrintFValue : public AbstractTest<T> {
public:
	virtual void first(std::ostream& out, const DebugHelper<T>& helper) {
		out <<
				">>> Average value f = " << avg(helper.f) << std::endl;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintFValue(*this);
	}
};

template <typename T>
class PrintSubProblem1Values : public AbstractTest<T> {
public:
	virtual void first(std::ostream& out, const DebugHelper<T>& helper) {
		out <<
				">>> Average value p1 = " << toString(avg(helper.p1)) << std::endl <<
				">>> Average value p2 = " << toString(avg(helper.p2)) << std::endl;
	}
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {
		out <<
				">>> Average value p1 = " << toString(avg(helper.p1)) << std::endl <<
				">>> Average value p2 = " << toString(avg(helper.p2)) << std::endl;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintSubProblem1Values(*this);
	}
};

template <typename T>
class PrintSubProblem2Values : public AbstractTest<T> {
public:
	virtual void first(std::ostream& out, const DebugHelper<T>& helper) {
		out <<
				">>> Average value p3 = " << toString(avg(helper.p3)) << std::endl;
	}
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {
		out <<
				">>> Average value p3 = " << toString(avg(helper.p3)) << std::endl;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintSubProblem2Values(*this);
	}
};

template <typename T>
class PrintSubProblem3Values : public AbstractTest<T> {
public:
	virtual void first(std::ostream& out, const DebugHelper<T>& helper) {
		out <<
				">>> Average value psi = " << avg(helper.psi) << std::endl;
	}
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {
		out <<
				">>> Average value psi = " << avg(helper.psi) << std::endl;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintSubProblem3Values(*this);
	}
};

template <typename T>
class PrintSubProblem4Values : public AbstractTest<T> {
public:
	virtual void first(std::ostream& out, const DebugHelper<T>& helper) {
		out <<
				">>> Average value u = " << avg(helper.u) << std::endl;
	}
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {
		out <<
				">>> Average value u = " << avg(helper.u) << std::endl;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintSubProblem4Values(*this);
	}
};

template <typename T>
class PrintLambdaValues : public AbstractTest<T> {
public:
	virtual void first(std::ostream& out, const DebugHelper<T>& helper) {
		out <<
				">>> Average value l1 = " << toString(avg(helper.l1)) << std::endl <<
				">>> Average value l2 = " << toString(avg(helper.l2)) << std::endl <<
				">>> Average value l3 = " << toString(avg(helper.l3)) << std::endl;
	}
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {
		out <<
				">>> Average value l1 = " << toString(avg(helper.l1)) << std::endl <<
				">>> Average value l2 = " << toString(avg(helper.l2)) << std::endl <<
				">>> Average value l3 = " << toString(avg(helper.l3)) << std::endl;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintLambdaValues(*this);
	}
};

template <typename T>
class PrintStepProgress : public AbstractTest<T> {
public:
	virtual void first(std::ostream& out, const DebugHelper<T>& helper) {
		noisyValue = objectFunction(helper.f, helper.f, helper.eps, helper.phi, helper.beta);
		firstValue = objectFunction(helper);
		firstAugmentedValue = augObjectFunction(helper);
		lastValue = firstValue;
		lastAugmentedValue = firstAugmentedValue;
		out <<
				">>> Object function (u_0) = " << firstValue << std::endl <<
				">>> Object function (noisy) = " << noisyValue << std::endl <<
				">>> Augmented object function = " << firstAugmentedValue << std::endl;
	}
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {
		T newValue = objectFunction(helper);
		T newAugmentedValue = augObjectFunction(helper);

		out << ">>> Object function = " << newValue <<
				" (diff = " << newValue - lastValue << ")" << std::endl;
		out << ">>> Total change = " << newValue - firstValue << " (diff ratio = " << newValue / firstValue << ")" << std::endl;
		out << ">>> Total diff to noisy = " << newValue - noisyValue << " (diff ratio = " << newValue / noisyValue << ")" << std::endl;

		out << ">>> Augmented object function = " << newAugmentedValue <<
				" (diff = " << newAugmentedValue - lastAugmentedValue << ")" << std::endl;

		out << ">>> Augmented total change = " << newAugmentedValue - firstAugmentedValue << " (diff ratio = " << newAugmentedValue / firstAugmentedValue << ")" << std::endl;

		lastValue = newValue;
		lastAugmentedValue = newAugmentedValue;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintStepProgress(*this);
	}
private:
	T noisyValue;
	T firstValue;
	T firstAugmentedValue;
	T lastValue;
	T lastAugmentedValue;
};

template <typename T>
class PrintSubProblemProgress : public AbstractTest<T> {
public:
	explicit PrintSubProblemProgress(const std::string& _msg) : msg(_msg) {}
	virtual void before(std::ostream& out, const DebugHelper<T>& helper) {
		lastAugmentedValue = augObjectFunction(helper);
	}
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {
		T newAugmentedValue = augObjectFunction(helper);

		out << ">>> " << msg << " = " <<
				newAugmentedValue - lastAugmentedValue <<
				" (diff ratio = " << (newAugmentedValue - lastAugmentedValue) / lastAugmentedValue << ")" << std::endl;

		lastAugmentedValue = newAugmentedValue;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintSubProblemProgress(*this);
	}
private:
	std::string msg;
	T lastAugmentedValue;
};

template <typename T>
class PrintConstraints : public AbstractTest<T> {
public:
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {
		std::size_t s4n1 = helper.p1.getN1();
		std::size_t s4n2 = helper.p1.getN2();

		Q<T> gradU = grad(helper.u);
		out << ">>> (1/2) | p1 - grad u | / max(|p1|,|grad u|) = ";
		out << 0.5 * vectorNorm2(helper.p1-gradU) / std::max(vectorNorm2(helper.p1),vectorNorm2(gradU)) << std::endl;

		out << ">>> (1/2) | p2 - p1/sqrt(beta+dot(p1)) | = ";
		T tmp = 0.0;
		for(std::size_t j = 1; j <= s4n2; j++) {
			for(std::size_t i = 1; i <= s4n1; i++) {
				tmp += dot(helper.p2(i,j)-helper.p1(i,j)/sqrt(helper.beta+dot(helper.p1(i,j))));
			}
		}
		out << 0.5 * sqrt(tmp) << std::endl;

		out << ">>> (1/2) | p3 - p2 | / max(|p3|,|p2|)= ";
		out << 0.5 * vectorNorm2(helper.p3-helper.p2) / std::max(vectorNorm2(helper.p3),vectorNorm2(helper.p2)) << std::endl;

		V<T> divP3 = div(helper.p3);
		out << ">>> (1/2) | psi - div(p3) | / max(|psi|,|div(p3)|) = ";
		out << 0.5 * vectorNorm2(helper.psi-divP3) / std::max(vectorNorm2(helper.psi),vectorNorm2(divP3)) << std::endl;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintConstraints(*this);
	}
};

template <typename T>
class PrintSubProblem1Errors : public AbstractTest<T> {
public:
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {
		std::size_t s4n1 = helper.p1.getN1();
		std::size_t s4n2 = helper.p1.getN2();

		int nanCounter1 = 0;
		int nanCounter2 = 0;
		int infCounter1 = 0;
		int infCounter2 = 0;
		check(helper.p1, &nanCounter1, &infCounter1);
		check(helper.p2, &nanCounter2, &infCounter2);
		if(0 < nanCounter1)
			out << ">>> Warning: Subproblem solver no. 1 generates " <<
				100.0*nanCounter1/((s4n1)*(s4n2)) << "% (p1) NaNs" << std::endl;
		if(0 < nanCounter2)
			out << ">>> Warning: Subproblem solver no. 1 generates " <<
				100.0*nanCounter2/((s4n1)*(s4n2)) << "% (p2) NaNs" << std::endl;
		if(0 < infCounter1)
			out << ">>> Warning: Subproblem solver no. 1 generates " <<
				100.0*infCounter1/((s4n1)*(s4n2)) << "% (p1) Infs" << std::endl;
		if(0 < infCounter2)
			out << ">>> Warning: Subproblem solver no. 1 generates " <<
				100.0*infCounter2/((s4n1)*(s4n2)) << "% (p2) Infs" << std::endl;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintSubProblem1Errors(*this);
	}
};

template <typename T>
class PrintSubProblem2Errors : public AbstractTest<T> {
public:
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {
		std::size_t s4n1 = helper.p3.getN1();
		std::size_t s4n2 = helper.p3.getN2();

		int nanCounter1 = 0;
		int infCounter1 = 0;

		check(helper.p3, &nanCounter1, &infCounter1);
		if(0 < nanCounter1)
			out << ">>> Warning: Subproblem solver no. 2 generates " <<
				100.0*nanCounter1/((s4n1)*(s4n2)) << "% (p3) NaNs" << std::endl;
		if(0 < infCounter1)
			out << ">>> Warning: Subproblem solver no. 2 generates " <<
				100.0*infCounter1/((s4n1)*(s4n2)) << "% (p3) Infs" << std::endl;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintSubProblem2Errors(*this);
	}
};

template <typename T>
class PrintSubProblem3Errors : public AbstractTest<T> {
public:
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {
		std::size_t width = helper.psi.getN1();
		std::size_t height = helper.psi.getN2();

		int nanCounter1 = 0;
		int infCounter1 = 0;

		check(helper.psi, &nanCounter1, &infCounter1);
		if(0 < nanCounter1)
			out << ">>> Warning: Subproblem solver no. 3 generates " <<
				100.0*nanCounter1/(width*height) << "% (psi) NaNs" << std::endl;
		if(0 < infCounter1)
			out << ">>> Warning: Subproblem solver no. 3 generates " <<
				100.0*infCounter1/(width*height) << "% (psi) Infs" << std::endl;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintSubProblem3Errors(*this);
	}
};

template <typename T>
class PrintSubProblem4Errors : public AbstractTest<T> {
public:
	virtual void after(std::ostream& out, const DebugHelper<T>& helper) {
		std::size_t width = helper.u.getN1();
		std::size_t height = helper.u.getN2();

		int nanCounter1 = 0;
		int infCounter1 = 0;

		check(helper.u, &nanCounter1, &infCounter1);
		if(0 < nanCounter1)
			out << ">>> Warning: Subproblem solver no. 4 generates " <<
				100.0*nanCounter1/(width*height) << "% (u) NaNs" << std::endl;
		if(0 < infCounter1)
			out << ">>> Warning: Subproblem solver no. 4 generates " <<
				100.0*infCounter1/(width*height) << "% (u) Infs" << std::endl;
	}
	virtual AbstractTest<T>* clone() const {
		return new PrintSubProblem4Errors(*this);
	}
};

}

#endif
