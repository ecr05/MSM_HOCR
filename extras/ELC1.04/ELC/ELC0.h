/*//////////////////////////////////////////////////////////////////////////////////////////////////
///  ELC0.h   Clique Reduction by Excludable Local Configuration: helper classes and functions
///  Version 1.04         September 12th, 2014
////////////////////////////////////////////////////////////////////////////////////////////////////

Copyright 2014 Hiroshi Ishikawa. All rights reserved.
This software can be used for research purposes only.
This software or its derivatives must not be publicly distributed
without a prior consent from the author (Hiroshi Ishikawa).

THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

For the latest version, check: http://www.f.waseda.jp/hfs/indexE.html

////////////////////////////////////////////////////////////////////////////////////////////////////

ELC0.h
This file contains the internal classes and functions used in ELC.h.
For the usage of the Excludable Local Configuration software, see ELC.h.

//////////////////////////////////////////////////////////////////////////////////////////////////*/


#ifndef __ELC0_H__
#define __ELC0_H__

#if defined(_DEBUG)
template<typename IT>
void index_ordered(IT i1, IT iend)
{
	for (IT i2 = i1 + 1; i2 != iend; i1++, i2++)
		assert(*i1 < *i2);
}
#else
template<typename IT> void index_ordered(IT i1, IT iend) {}
#endif



// Container for monomials of a fixed degree.
template <typename REAL> // REAL: real number
class TermsD
{
public:
	TermsD(int d) : degree(d) {}
	void clear() {vars.clear(); coefs.clear();}
	int size() const {return (int)coefs.size();}
	VID maxID() const {return vars.empty() ? -1 : *std::max_element(vars.begin(), vars.end());}
	int add(REAL c, VID* vs) {index_ordered(vs, vs + degree); vars.insert(vars.end(), vs, vs + degree); coefs.push_back(c); return (int)coefs.size() - 1;}
	int add1(REAL c, VID v) {assert(degree == 1); vars.push_back(v); coefs.push_back(c); return (int)coefs.size() - 1;}
	int add2(REAL c, VID v1, VID v2) {assert(degree == 2); assert(v1 < v2); vars.push_back(v1); vars.push_back(v2); coefs.push_back(c); return (int)coefs.size() - 1;}
	// find the term, if it exists, add coefficient and return true
	bool addCoef(int idx, VID* vs, REAL c) {VVec::const_iterator i = vars.begin() + idx * degree; if (!std::equal(i , i + degree, vs)) return false; coefs[idx] += c; return true;}
	bool addCoef1(int idx, VID v, REAL c) {assert(degree == 1); if (v != vars[idx]) return false; coefs[idx] += c; return true;}
	bool addCoef2(int idx, VID v1, VID v2, REAL c) {assert(degree == 2); assert(v1 < v2); if (v1 != vars[2 * idx] || v2 != vars[2 * idx + 1]) return false; coefs[idx] += c; return true;}

	void shrink()
	{
		int i1 = 0, i2, i3 = 0;
		while (i3 < coefs.size())
		{
			for (i2 = i3; i2 < coefs.size() && !coefs[i2]; i2++);
			for (i3 = i2; i3 < coefs.size() && coefs[i3]; i3++);
			if (i1 < i2)
			{
				std::copy(coefs.begin() + i2, coefs.begin() + i3, coefs.begin() + i1);
				std::copy(vars.begin() + i2 * degree, vars.begin() + i3 * degree, vars.begin() + i1 * degree);
			}
			i1 += i3 - i2;
		}
		if (i1 < i3)
		{
			coefs.erase(coefs.begin() + i1, coefs.begin() + i3);
			vars.erase(vars.begin() + i1 * degree, vars.begin() + i3 * degree);
		}
	}
	VVecIt getVIt(int i) {return vars.begin() + degree * i;}
	REAL coef(int i) {return coefs[i];}
private:
	VVec vars;
	std::vector<REAL> coefs;
	int degree;
};

// Container for monomials of all degrees.
template <typename REAL>
class Terms
{
public:
	Terms(int maxVars) : termsD(3), termIdx(maxVars) 
	{
		termsD[0] = 0;
		for (int d = 1; d <= 2; d++)
			termsD[d] = new TermsD<REAL>(d);
	}
	~Terms() {for (auto it = termsD.begin(); it != termsD.end(); it++) if (*it) delete *it;}
	void clear() {termIdx.clear(); for (auto it = termsD.begin(); it != termsD.end(); it++) if (*it) delete *it; termsD.clear();}
	void shrink()
	{
		for (auto it = termsD.begin(); it != termsD.end(); it++)
			if (*it) 
				(*it)->shrink(); 
		rebuildTermIdx();
	}
	struct TermPointer
	{
		TermPointer(int degree, int index) : d(degree), idx(index) {}
		int degree() const {return d;}
		bool operator==(const TermPointer& i) const {return d == i.d && idx == i.idx;}
		bool operator!=(const TermPointer& i) const {return !operator==(i);}
		bool operator<(const TermPointer& i) const {if (d == i.d) return idx < i.idx; else return d < i.d;}
	private:
		friend class Terms;
		int d;
		int idx;
	};
	typedef std::vector<TermPointer> TPVec;
	struct iterator 
	{
		iterator(Terms const * p, int degree, int index) : pObject(p), pTerm(degree, index) {}
		VVecIt vars() const {return pObject->getVIt(pTerm);}
		int degree() const {return pTerm.d;}
		REAL coef() const {return pObject->coef(pTerm);}
		void operator++()
		{
			pTerm.idx++;
			if (pTerm.idx >= pObject->termsD[pTerm.d]->size()) 
			{
				pTerm.idx = 0; 
				pTerm.d--;
				while (pTerm.d > 0 && (!pObject->termsD[pTerm.d] || !pObject->termsD[pTerm.d]->size()))
					pTerm.d--;
			}
		}
		bool operator==(const iterator& i) const {return pObject == i.pObject && pTerm == i.pTerm;}
		bool operator!=(const iterator& i) const {return !operator==(i);}
	private:
		friend class Terms;
		Terms const * pObject;
		TermPointer pTerm;
	};
	iterator begin() const
	{
		int d = (int)termsD.size() - 1;
		while (d > 0 && (!termsD[d] || !termsD[d]->size()))
			d--;
		return iterator(this, d, 0);
	}
	iterator end() const {return iterator(this, 0, 0);}

	struct variables 
	{
		variables(const iterator& it) : sz(it.degree()), z(it.vars()+sz) {assert(sz > 0);}
		variables(const Terms& terms, const TermPointer& pTerm) : sz(pTerm.degree()), z(terms.getVIt(pTerm)+sz) {assert(sz > 0);}
		int size() const {return sz;}
		VVecIt begin() const {return z - sz;}
		VVecIt end() const {return z;}
		template<typename C> int copy(C v) const {std::copy(z - sz, z, v); return sz;}
		bool includes(const variables& vi) const
		{
			if (sz < vi.sz || *(z-1) < *(vi.z-vi.sz) || *(vi.z-1) < *(z-sz))
				return false;
			VVecIt i1 = z - sz, i2 = vi.z - vi.sz; 
			while (true)
			{
				if (*i2 <*i1)
					return false;
				if (*i1 == *i2)
				{
					if (++i2 == vi.z)
						return true;
				}
				else if (++i1 == z)
					return false;
			}
		}
	protected:
		int sz;
		VVecIt z;
	};

	VVecIt getVIt(const TermPointer& pTerm) const {return termsD[pTerm.d]->getVIt(pTerm.idx);}
	REAL coef(const TermPointer& pTerm) const {return termsD[pTerm.d]->coef(pTerm.idx);}
	VID maxID() const {VID m = -1; for (auto it = termsD.begin(); it != termsD.end(); it++) if (*it) m = std::max(m, (*it)->maxID()); return m;}

	bool addCoef(REAL c, int degree, VID* vars) // find the term, if it exists, add coefficient and return true
	{
		if (termIdx.size() <= vars[0])
			return false;
		TPVec& terms = termIdx[vars[0]]; // all the terms containing the variable vars[0]
		for (auto it = terms.begin(); it != terms.end(); it++)
			if ((*it).d == degree && termsD[degree]->addCoef((*it).idx, vars, c))
				return true;
		return false;
	}

	bool addCoef1(REAL c, VID v) // find the term, if it exists, add coefficient and return true
	{
		if (termIdx.size() <= v)
			return false;
		TPVec& terms = termIdx[v]; // all the terms containing the variable v
		for (auto it = terms.begin(); it != terms.end(); it++)
			if ((*it).d == 1 && termsD[1]->addCoef1((*it).idx, v, c))
				return true;
		return false;
	}

	bool addCoef2(REAL c, VID v1, VID v2) // find the term, if it exists, add coefficient and return true
	{
		if (termIdx.size() <= v1)
			return false;
		TPVec& terms = termIdx[v1]; // all the terms containing the variable v1
		for (auto it = terms.begin(); it != terms.end(); it++)
			if ((*it).d == 2 && termsD[2]->addCoef2((*it).idx, v1, v2, c))
				return true;
		return false;
	}

	void resizeTermIdx(int degree, VID* vars)
	{
		int maxvar = *std::max_element(vars, vars + degree);
		if (termIdx.size() <= maxvar)
			termIdx.resize(maxvar + 1);
	}

	void add(REAL c, int degree, VID* vars) 
	{
		if (termsD.size() <= degree)
			termsD.insert(termsD.end(), degree - termsD.size() + 1, 0);
		if (termsD[degree] == 0)
			termsD[degree] = new TermsD<REAL>(degree);
		else if (addCoef(c, degree, vars))
			return;
		int idx = termsD[degree]->add(c, vars);
		resizeTermIdx(degree, vars);
		for (int i = 0; i < degree; i++)
			termIdx[vars[i]].push_back(TermPointer(degree, idx));
	}

	void add1(REAL c, VID v) {if (addCoef1(c, v)) return; int idx = termsD[1]->add1(c, v); if (termIdx.size() <= v) termIdx.resize(v + 1); termIdx[v].push_back(TermPointer(1, idx));}
	void add2(REAL c, VID v1, VID v2) {if (addCoef2(c, v1, v2)) return; int idx = termsD[2]->add2(c, v1, v2); VID vmax = std::max(v1,v2); if (termIdx.size() <= vmax) termIdx.resize(vmax + 1); termIdx[v1].push_back(TermPointer(2, idx)); termIdx[v2].push_back(TermPointer(2, idx));}
	int maxD() const {int i = (int)termsD.size() - 1; while (!termsD[i]) i--; return i;}
	TPVec& getTermsContaining(VID v) {return termIdx[v];}

protected:
	void rebuildTermIdx()
	{
		termIdx.resize(maxID() + 1);
		fill(termIdx.begin(), termIdx.end(), TPVec());
		for (iterator it = begin(); it != end(); ++it)
		{
			variables va = it;
			for (VVecIt vi = va.begin(); vi != va.end(); ++vi)
				termIdx[*vi].push_back(it.pTerm);
		}
	}

	std::vector<TermsD<REAL>*> termsD; // termsD[d] : degree d terms
	// for each variable id i, termIdx[i] lists all the terms that contains the variable
	std::vector<TPVec> termIdx;
};


struct evenodd
{
	evenodd() : degree(0) {}
	void create(int d)
	{
		if (degree >= d)
			return;
		degree = d;
		int s = 1 << (d-1);
		ebits.resize(s+1);
		obits.resize(s+1);
		obits[0] = 1;
		ebits[0] = 0;
		int k = 1;
		for (int i = 1; i <= s; i++)
		{
			if (i >= 2 * k)
				k *= 2;
			obits[i] = ebits[i % k] | (2 * k);
			ebits[i] = obits[i % k] | (2 * k);
		}
	}
	int degree;
	std::vector<BITMASK> ebits, obits;
};

// Pseudo Boolean functions with up to MASKSIZE variables.
template<typename REAL>
class smallPBF
{
	typedef typename Terms<REAL>::variables variables;
public:
	smallPBF() : constant(0), numvars(0), MAX_ELC(std::numeric_limits<REAL>::max()), MIN_ELC(std::numeric_limits<REAL>::min()) {}

	void setVars(int n, VVecIt vs) 		// vs must be sorted
	{
		assert(n < MASKSIZE);
		numvars = n;
		std::copy(vs, vs + n, vars);
	}

	void setVars(const variables& vs) {numvars = vs.copy(vars);}

	void setVars(smallPBF& pbf)
	{
		numvars = pbf.numvars;
		std::copy(pbf.vars, pbf.vars + pbf.numvars, vars);
	}

	void add(const variables& vi, REAL c) {terms.push_back(std::make_pair(varvector(vi), c));}
	void addConst(REAL c)	 {constant += c;}

	REAL eval(BITMASK x) const
	{
		REAL rv = constant;
		for (auto it = terms.begin(); it != terms.end(); it++)
			if (!(~x & (*it).first))
				rv += (*it).second;
		return rv;
	}

	// odd: if false, find the even ELC (with even number of 1s), if true, find the odd ELC
	// q: if the gap is smaller than this value, opposite parity is fine, too.
	int getELC(const variables& core, bool odd, REAL& q)
	{
		int numcore = core.size();
		BITMASK tn = 1 << numcore;
		REAL maxmin1 = MIN_ELC, maxmin2 = MIN_ELC, minmax = MAX_ELC;
		eo.create(numcore);
		std::vector<BITMASK>& eobits = odd ? eo.obits : eo.ebits;
		BITMASK corex= varvector(core), i = 0;
		int maxminIdx1 = -1, maxminIdx2 = -1;
		smallPBF collapsed;
		for (BITMASK ix = 0; ix < tn; ix++)
		{
			getCollapsed(collapsed, numcore, corex, ix);
			REAL maxv, minv;
			collapsed.getBruteForceMinMax(maxv, minv);
			if (maxv < minmax)
				minmax = maxv;
			if (eobits[i] == ix)
			{
				if (minv > maxmin1)
				{
					maxmin1 = minv;
					maxminIdx1 = ix;
				}
				i++;
			}
			else if (minv > maxmin2)
			{
				maxmin2 = minv;
				maxminIdx2 = ix;
			}
			if (minmax < maxmin1)
				return maxminIdx1;
			else if (minmax + q < maxmin2) // (maxmin2 - minmax) is the gap
			{
				q *= -1;
				return maxminIdx2;
			}
		}
		return -1;
	}

	int getMax(bool odd) // odd: if false, find the even max (with even number of 1s), if true, find the odd max
	{
		int maxIdx = -1, tn = 1 << (numvars - 1);
		REAL maxv = MIN_ELC;
		eo.create(numvars);
		for (int i = 0; i < tn; i++)
		{
			BITMASK ix = odd ? eo.obits[i] : eo.ebits[i];
			REAL val = eval(ix);
			if (val > maxv)
			{
				maxv = val;
				maxIdx = ix;
			}
		}
		return maxIdx;
	}

private:
	BITMASK varvector(const variables& va)
	{
		BITMASK b = 1 << (numvars - 1), x = 0;
		VID* j = vars;
		for (VVecIt vi = va.begin(); vi != va.end(); ++vi)
		{
			while (*vi != *j)
			{
				j++;
				b >>= 1;
			}
			x |= b;
		}
		return x;
	}

	void getBruteForceMinMax(REAL &maxv, REAL &minv)
	{
		if (!numvars || terms.empty())
		{
			maxv = minv = constant;
			return;
		}
		BITMASK tn = 1 << numvars;
		REAL vals[1 << MAXVARELC];
		std::fill(vals, vals + tn, constant);
		for (auto it = terms.begin(); it != terms.end(); it++)
			for (BITMASK ix = 0; ix < tn; ix++)
				if (!(~ix & (*it).first))
					vals[ix] += (*it).second;
		maxv = MIN_ELC;
		minv = MAX_ELC;
		for (BITMASK ix = 0; ix < tn; ix++)
		{
			REAL val = vals[ix];
			if (maxv < val)
				maxv = val;
			if (minv > val)
				minv = val;
		}
	}

	BITMASK idxCore(int numcore, BITMASK corex, int ix)
	{
		assert(corex != 0);
		assert(numcore < numvars);
		BITMASK ixcore = 0, b1 = 1 << (numcore - 1), b2 = 1 << (numvars - 1);
		for (; b1 > 0; b1 >>= 1, b2 >>= 1)
		{
			while (!(corex & b2))
				b2 >>= 1;
			if (ix & b1)
				ixcore |= b2;
		}
		return ixcore;
	}

	void getCollapsed(smallPBF& collapsed, int numcore, BITMASK corex, BITMASK ix)
	{
		BITMASK ixcore = idxCore(numcore, corex, ix);
		collapsed.constant = constant;
		collapsed.terms.clear();
		BITMASK mustbezero = corex & ~ixcore, z = 0;
		for (auto it = terms.begin(); it != terms.end(); it++)
			if (!(mustbezero & (*it).first))
			{
				if (~corex & (*it).first)
					z |= ~corex & (*it).first;
				else
					collapsed.constant += (*it).second;
			}
		collapsed.numvars = 0;
		if (!z)
			return;
		int j, k = 0;
		BITMASK b2 = 1 << (numvars - 1), m[MASKSIZE];
		for (j = 0; j < numvars; j++, b2 >>= 1)
			if (z & b2)
			{
				collapsed.vars[k] = vars[j];
				m[k++] = b2;
			}
		collapsed.numvars = k;
		for (auto it = terms.begin(); it != terms.end(); it++)
			if (!(mustbezero & (*it).first))
			{
				if (~corex & (*it).first)
				{
					BITMASK y = 0, b1 = 1 << (k - 1);
					for (int j = 0; j < k; j++, b1 >>= 1)
						if ((*it).first & m[j])
							y |= b1;
					assert(y);
					collapsed.terms.push_back(std::make_pair(y, (*it).second));
				}
			}
	}

	VID vars[MASKSIZE];
	REAL constant;
	int numvars;
	std::vector<std::pair<BITMASK, REAL> > terms;
	static evenodd eo;
	REAL MAX_ELC, MIN_ELC;
};

template<typename REAL> evenodd smallPBF<REAL>::eo;

inline int matrixElement(int m, int n)
{
	if (~m & n)
		return 0;
	for (m ^= n, n = 0; m; m >>= 1)
		n ^= m & 1;
	return n ? -1 : 1;
}

// Gets the coefficient for the i'th monomial
// From an energy E that has n values
template<typename REAL>
inline REAL getCoef(int i, int n, REAL E[])
{
	REAL rv = 0;
	for (int j = 0; j < n; j++)
		rv += E[j] * matrixElement(i, j);
	return rv;
}



#endif // __ELC0_H__

