//  
//  Declarations/template-bodies for sparse matrix class SpMat
//
//  SpMat.h
//
//  Implements bare-bones sparse matrix class. 
//  Main considerations has been efficiency when constructing
//  from Compressed Column format, when multiplying with vector,
//  transposing and multiplying with a vector and when concatenating.
//  Other operations which have not been prioritised such as 
//  for example inserting elements in a random order may be
//  a bit slow.  
//
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2007 University of Oxford 
//

#ifndef SpMat_h
#define SpMat_h

#include <vector>
#include <fstream>
#include <iomanip>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "cg.h"
#include "bicg.h"
#include "miscmaths.h"

namespace MISCMATHS {

class SpMatException: public std::exception
{
private:
  std::string m_msg;
public:
  SpMatException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("SpMat::" + m_msg).c_str();
  }

  ~SpMatException() throw() {}
};

enum MatrixType {UNKNOWN, ASYM, SYM, SYM_POSDEF};

template<class T>
class Preconditioner;

template<class T>
class Accumulator;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class SpMat:
// Interface includes:
// Multiplication with scalar:  A*=s, B=s*A, B=A*s, A and B SpMat
// Multiplication with vector:  b=A*x, A SpMat, b and x ColumnVector
// Transpose and mul with vector: b=A.trans_mult(x), A SpMat, b and x ColumnVector
// Multiplication with sparse matrix: C=A*B, A, B and C SpMat
// Addition with sparse matrix: A+=B, C=A+B, A, B and C SpMat
// Horisontal concatenation: A|=B, C=A|B, A, B and C SpMat
// Vertical concatenation: A&=B, C=A&B, A, B and C SpMat
//
// Multiplications and addition with NEWMAT matrices are
// accomplished through type-conversions. For example
// A = B*SpMat(C), A and B SpMat, C NEWMAT
// A = B.AsNewmat()*C, B SpMat, A and C NEWMAT
//
// Important implementation detail:
// _nz or .NZ() isn't strictly speaking the # of non-zero elements,
// but rather the number of elements that has an explicit 
// representation, where that representation may in principle
// be 0. This is in contrast to e.g. Matlab which chooses
// not to represent an element when its value is zero. I have
// chosen this variant because of my main use of the class where 
// it is very convenient if e.g. my Hessian and the Gibbs form
// of membrane energy has the same sparsity pattern.
// For most users this is of no consequence and they will 
// never explicitly represent a zero.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template<class T>
class SpMat
{
public:
  SpMat() : _m(0), _n(0), _nz(0), _ri(0), _val(0), _pw(false) {}
  SpMat(unsigned int m, unsigned int n) : _m(m), _n(n), _nz(0), _ri(n), _val(n), _pw(false) {}
  SpMat(unsigned int m, unsigned int n, const unsigned int *irp, const unsigned int *jcp, const double *sp);
  SpMat(const NEWMAT::GeneralMatrix& M);
  SpMat(const std::string& fname);
  ~SpMat() {}

  unsigned int Nrows() const {return(_m);}
  unsigned int Ncols() const {return(_n);}
  unsigned int NZ() const {return(_nz);}

  NEWMAT::ReturnMatrix AsNEWMAT() const;
  void Save(const std::string&   fname,
            unsigned int         precision) const;
  void Save(const std::string&   fname) const {Save(fname,8);}
  void Print(const std::string&  fname,
             unsigned int        precision) const;
  void Print(const std::string&  fname) const {Print(fname,8);}
  void Print(unsigned int        precision) const {Print(std::string(""),precision);}
  void Print() const {Print(8);}
  void WarningsOn() {_pw=true;}
  void WarningsOff() {_pw=false;}


  T Peek(unsigned int r, unsigned int c) const;
  T operator()(unsigned int r, unsigned int c) const {return(Peek(r,c));}    // Read-only

  void Set(unsigned int r, unsigned int c, const T& v) {here(r,c) = v;}             // Set a single value
  void SetColumn(unsigned int c, const NEWMAT::ColumnVector& col, double eps=0.0);  // Set a whole column (obliterating what was there before) 
  void AddTo(unsigned int r, unsigned int c, const T& v) {here(r,c) += v;}          // Add value to a single (possibly existing) value

  SpMat<T>& operator+=(const SpMat& M) 
  {
    if (same_sparsity(M)) return(add_same_sparsity_mat_to_me(M,1));
    else return(add_diff_sparsity_mat_to_me(M,1));
  }
  SpMat<T>& operator-=(const SpMat& M) 
  {
    if (same_sparsity(M)) return(add_same_sparsity_mat_to_me(M,-1)); 
    else return(add_diff_sparsity_mat_to_me(M,-1));
  }
  const NEWMAT::ReturnMatrix operator*(const NEWMAT::ColumnVector& x) const;       // Multiplication with column vector
  const NEWMAT::ReturnMatrix trans_mult(const NEWMAT::ColumnVector& x) const;      // Multiplication of transpose with column vector
  const NEWMAT::ReturnMatrix TransMult(const NEWMAT::ColumnVector& x) const {
    return(trans_mult(x));                                                         // Duplication for compatibility with IML++
  }
  const SpMat<T> TransMult(const SpMat<T>& B) const;                               // Multiplication of transpose(*this) with sparse matrix B
  SpMat<T>& operator*=(double s);                                                  // Multiplication of self with scalar
  SpMat<T> operator-(const SpMat<T>& M) const {return(SpMat<T>(M) *= -1.0);}       // Unary minus
  SpMat<T>& operator|=(const SpMat<T>& rh);                                        // Hor concat to right
  SpMat<T>& operator&=(const SpMat<T>& bh);                                        // Vert concat below

  const SpMat<T> TransMultSelf() const {return(TransMult(*this));}                 // Returns transpose(*this)*(*this)

  const SpMat<T> t() const;                                                        // Returns transpose(*this). Avoid, if at all possible.
  
  friend class Accumulator<T>; 

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //
  // Class ColumnIterator
  //
  // This implements a const forward iterator for a given column
  //
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // template<class TT>
  class ColumnIterator
  {
  public:
    ColumnIterator(const SpMat<T>& mat, unsigned int col, bool end=false) {
      if (end) { _ri_it = mat._ri[col-1].end(); _val_it = mat._val[col-1].end(); }
      else { _ri_it = mat._ri[col-1].begin(); _val_it = mat._val[col-1].begin(); }
    }
    ~ColumnIterator() {}
    T operator*() const { return(*_val_it); }
    unsigned int Row() const { return((*_ri_it)+1); }
    bool operator==(const ColumnIterator& rhs) const { return(_val_it == rhs._val_it); }
    bool operator!=(const ColumnIterator& rhs) const { return(!(*this == rhs)); }
    // Prefix increment. Use whenever possible
    ColumnIterator& operator++() { ++_ri_it; ++_val_it; return(*this); }
    // Postfix increment. Avoid.
    ColumnIterator& operator++(int dummy) { ColumnIterator clone(*this); ++_ri_it; ++_val_it; return(clone); }
  private:
    typename std::vector<T>::const_iterator    _val_it;
    std::vector<unsigned int>::const_iterator  _ri_it;
  };

  ColumnIterator begin(unsigned int col) const { if (col>_n) throw SpMatException("ColumnIterator: col out of range"); return(ColumnIterator(*this,col)); }
  ColumnIterator end(unsigned int col) const { if (col>_n) throw SpMatException("ColumnIterator: col out of range"); return(ColumnIterator(*this,col,true)); }

  template<class TT>
  friend const SpMat<TT> operator*(const SpMat<TT>& lh, const SpMat<TT>& rh);      // Multiplication of two sparse matrices

  NEWMAT::ReturnMatrix SolveForx(const NEWMAT::ColumnVector&           b,          // Solve for x in b=(*this)*x
                                 MatrixType                            type = UNKNOWN,
                                 double                                tol = 1e-4, 
                                 unsigned int                          miter = 200,
				 boost::shared_ptr<Preconditioner<T> > C = boost::shared_ptr<Preconditioner<T> >()) const;

  NEWMAT::ReturnMatrix SolveForx(const NEWMAT::ColumnVector&            b,
			         MatrixType                             type,
			         double                                 tol,
			         unsigned int                           miter,
				 const NEWMAT::ColumnVector&            x_init) const;

  NEWMAT::ReturnMatrix SolveForx(const NEWMAT::ColumnVector&            b,
			         MatrixType                             type,
			         double                                 tol,
			         unsigned int                           miter,
				 boost::shared_ptr<Preconditioner<T> >  C,
				 const NEWMAT::ColumnVector&            x_init) const;

protected:
  const std::vector<unsigned int>& get_ri(unsigned int i) const;
  const std::vector<T>&            get_val(unsigned int i) const;
  T& here(unsigned int r, unsigned int c);

private:
  unsigned int                              _m;
  unsigned int                              _n;
  unsigned long                             _nz;
  std::vector<std::vector<unsigned int> >   _ri;
  std::vector<std::vector<T> >              _val;
  bool                                      _pw;   // Print Warnings

  bool found(const std::vector<unsigned int>&  ri, unsigned int key, int& pos) const;  
  void insert(std::vector<unsigned int>& vec, int indx, unsigned int val);
  void insert(std::vector<T>& vec, int indx, const T& val);
  bool same_sparsity(const SpMat<T>& M) const;
  SpMat<T>& add_same_sparsity_mat_to_me(const SpMat<T>& M, double s);
  SpMat<T>& add_diff_sparsity_mat_to_me(const SpMat<T>& M, double s);
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class Preconditioner:
//
// I haven't used conditioner for close to 20 years now, so writing
// this class was a special treat for me. A preconditioner is used
// to render the coefficient-matrix corresponding to some set of
// linear equations better conditioned. A concrete example would be
// when some set of columns/rows have a different scale than the 
// others, resulting in poor convergence of for example a conjugate
// gradient search. The simplest form of preconditioner might then 
// be inv(diag(A)), where A is the original matrix. It simply scales
// the columns of A with the inverse of the diagonal elements. This
// simple conditioning works fine when A is diagonal domninant, which
// i typically the case with e.g. Hessians in spatial normalisation.
// If not, a more sophisticated version like incomplete Cholesky
// decomposition might be needed.
// As of yet only diagonal preconditioners have been implemented.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template<class T>
class Preconditioner
{
public:
  Preconditioner(const SpMat<T>& M) : _m(M.Nrows())
  {
    if (M.Nrows() != M.Ncols()) throw SpMatException("Preconditioner: Matrix to condition must be square");
  }
  virtual ~Preconditioner() {}
  unsigned int Nrows() const {return(_m);}
  virtual NEWMAT::ReturnMatrix solve(const NEWMAT::ColumnVector& x) const = 0;
  virtual NEWMAT::ReturnMatrix trans_solve(const NEWMAT::ColumnVector& x) const = 0;
private:
  unsigned int     _m;
};

template<class T>
class DiagPrecond: public Preconditioner<T>
{
public:
  DiagPrecond(const SpMat<T>& M) : Preconditioner<T>(M), _diag(M.Nrows())
  {
    for (unsigned int i=0; i<Preconditioner<T>::Nrows(); i++) {
      _diag[i] = M(i+1,i+1);
      if (_diag[i] == 0.0) throw SpMatException("DiagPrecond: Cannot condition singular matrix");
    }
  }
  ~DiagPrecond() {}
  NEWMAT::ReturnMatrix solve(const NEWMAT::ColumnVector& x) const
  {
    if (x.Nrows() != int(Preconditioner<T>::Nrows())) throw SpMatException("DiagPrecond::solve: Vector x has incompatible size");

    NEWMAT::ColumnVector  b(Preconditioner<T>::Nrows());
    double                *bptr = static_cast<double *>(b.Store());
    double                *xptr = static_cast<double *>(x.Store());
    for (unsigned int i=0; i<Preconditioner<T>::Nrows(); i++) bptr[i] = xptr[i]/static_cast<double>(_diag[i]);
    b.Release();

    return(b);
  }
  NEWMAT::ReturnMatrix trans_solve(const NEWMAT::ColumnVector& x) const {return(solve(x));}

private:
  std::vector<T>   _diag;
};  
  

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class Accumulator:
//
// The concept of an accumulator was "borrowed" from Gilbert et al. 
// 92. It is intended as a helper class for SpMat and is used to
// hold the content of one column of a matrix. This column can then
// be accessed both by indexing a certain element, and also by indexing
// only non-zero elements.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template<class T>
class Accumulator
{
public:
  Accumulator(unsigned int sz) : _no(0), _sz(sz), _sorted(true), _occ(new bool [sz]), _val(new T [sz]), _occi(new unsigned int [sz]) 
  {
    for (unsigned int i=0; i<_sz; i++) {_occ[i]=false; _val[i]=static_cast<T>(0.0);}
  }
  ~Accumulator() {delete [] _occ; delete [] _val; delete [] _occi;}
  void Reset() {for (unsigned int i=0; i<_no; i++) {_occ[_occi[i]]=false; _val[_occi[i]]=static_cast<T>(0.0);} _no=0;}
  T& operator()(unsigned int i);
  unsigned int NO() const {return(_no);}
  unsigned int ri(unsigned int i) {                             // Index of i'th non-zero value.
    if (!_sorted) {sort(_occi,&(_occi[_no])); _sorted=true;}
    return(_occi[i]);
  }
  const T& val(unsigned int i) {                                // i'th non-zero value. Call ri(i) to find what index that corresponds to
    if (!_sorted) {sort(_occi,&(_occi[_no])); _sorted=true;}    
    return(_val[_occi[i]]);
  }
  const T& val_at(unsigned int i) const {return(_val[i]);}      // Value for index i (or i+1)
  const bool& occ_at(unsigned int i) const {return(_occ[i]);}   // Is value for index i non-zero
  
  const Accumulator<T>& ExtractCol(const SpMat<T>& M, unsigned int c);
  
private:
  unsigned int                   _no;       // Number of occupied positions
  unsigned int                   _sz;       // Max size of accumulated vector
  bool                           _sorted;   // True if _occi is ordered
  bool                           *_occ;     // True if position is "occupied"  
  T                              *_val;     // "Value" in position
  unsigned int                   *_occi;    // Unordered list of occupied indicies
};

/////////////////////////////////////////////////////////////////////
//
// Constructs sparse matrix from Compressed Column Storage representation
//
/////////////////////////////////////////////////////////////////////

template<class T>
SpMat<T>::SpMat(unsigned int m, unsigned int n, const unsigned int *irp, const unsigned int *jcp, const double *sp)
  : _m(m), _n(n), _nz(0), _ri(n), _val(n), _pw(false)
{
  _nz = jcp[n];
  unsigned long nz = 0;
  for (unsigned int c=0; c<_n; c++) {
    if (int len = jcp[c+1]-jcp[c]) {
      std::vector<unsigned int>&  ri = _ri[c];
      std::vector<T>&             val = _val[c];
      const unsigned int          *iptr = &(irp[jcp[c]]);
      const double                *vptr = &(sp[jcp[c]]);
      ri.resize(len);
      val.resize(len);
      for (int i=0; i<len; i++) {
        ri[i] = iptr[i];
        val[i] = static_cast<T>(vptr[i]);
        nz++;
      }
    }
  }
  if (nz != _nz) throw SpMatException("SpMat: Compressed column input not self consistent");
}

/////////////////////////////////////////////////////////////////////
//
// Constructs sparse matrix from NEWMAT Matrix or Vector
//
/////////////////////////////////////////////////////////////////////

template<class T>
SpMat<T>::SpMat(const NEWMAT::GeneralMatrix& M)
  : _m(M.Nrows()), _n(M.Ncols()), _nz(0), _ri(M.Ncols()), _val(M.Ncols()), _pw(false)
{
  double *m = static_cast<double *>(M.Store());

  for (unsigned int c=0; c<_n; c++) {
    // First find # of non-zeros elements in column
    unsigned int cnz = 0;
    for (unsigned int i=0; i<_m; i++) {
      if (m[i*_n+c]) cnz++;
    }
    if (cnz) {
      std::vector<unsigned int>&   ri = _ri[c];
      std::vector<T>&              val = _val[c];
      ri.resize(cnz);
      val.resize(cnz);
      for (unsigned int rii=0, i=0; i<_m; i++) {
        if (double v = m[i*_n+c]) {
          ri[rii] = i;
          val[rii] = v;
          rii++;
	}
      }
      _nz += cnz;
    }
  }  
}

/////////////////////////////////////////////////////////////////////
//
// Constructs matrix from row col val format produced by 
// Save/Print below.
//
/////////////////////////////////////////////////////////////////////

template<class T>
SpMat<T>::SpMat(const std::string&  fname)
: _m(0), _n(0), _nz(0), _ri(0), _val(0), _pw(false)
{
  // First read data into (nz+1)x3 NEWMAT matrix
  NEWMAT::Matrix rcv;
  try {
    rcv = read_ascii_matrix(fname);
  }
  catch(...) {
    throw SpMatException("SpMat::SpMat(string& fname): cannot read file given by fname");
  }
  // Then interpret it
  if (rcv(rcv.Nrows(),3)) throw SpMatException("SpMat::SpMat(string& fname): Last row must have zero value and indicate matrix size");
  _m = static_cast<unsigned int>(rcv(rcv.Nrows(),1)+0.5);
  _n = static_cast<unsigned int>(rcv(rcv.Nrows(),2)+0.5);
  // cout << "rcv = " << endl << rcv << endl << "_n = " << _n << endl;
  _ri.resize(_n);
  _val.resize(_n);
  // First pass to see how many elements in each colum
  std::vector<unsigned int> col_count(_n,0);
  unsigned int col = static_cast<unsigned int>(rcv(1,2)+0.5);
  for (unsigned int indx=1; indx<static_cast<unsigned int>(rcv.Nrows()); indx++) {
    if (static_cast<unsigned int>(rcv(indx,2)+0.5) != col) {
      if (static_cast<unsigned int>(rcv(indx,2)+0.5) < col) throw SpMatException("SpMat::SpMat(string& fname): Column index must be monotonously increasing");
      else col = static_cast<unsigned int>(rcv(indx,2)+0.5);
      if (col > _n) throw SpMatException("SpMat::SpMat(string& fname): File internally inconsistent");
    }
    // cout << "col = " << col << endl;
    col_count[col-1]++;
  }
  // Second pass to allocate and fill vectors
  unsigned int indx=1;
  for (col=0; col<_n; col++) {
    std::vector<unsigned int>& ri = _ri[col];
    std::vector<T>&            val = _val[col];
    ri.resize(col_count[col]);
    val.resize(col_count[col]);
    for (unsigned int i=0; i<col_count[col]; i++, indx++) {
      if (i && ri[i] >= static_cast<unsigned int>(rcv(indx,1)+0.5)) throw SpMatException("SpMat::SpMat(string& fname): Row index must be monotonously increasing");
      if (static_cast<unsigned int>(rcv(indx,1)+0.5) < 1 || static_cast<unsigned int>(rcv(indx,1)+0.5) > _m) {
        throw SpMatException("SpMat::SpMat(string& fname): Row index outside 1 -- -m range");
      } 
      ri[i] = static_cast<unsigned int>(rcv(indx,1)+0.5) - 1;
      val[i] = rcv(indx,3);
      _nz++;
    }
  }
}


/////////////////////////////////////////////////////////////////////
//
// Returns matrix in NEWMAT matrix format. Useful for debugging
//
/////////////////////////////////////////////////////////////////////

template<class T>
NEWMAT::ReturnMatrix SpMat<T>::AsNEWMAT() const
{
  NEWMAT::Matrix M(_m,_n);
  M = 0.0;
  for (unsigned int c=0; c<_n; c++) {
    if (_ri[c].size()) {
      const std::vector<unsigned int>&    ri = _ri[c];
      const std::vector<T>&               val = _val[c];
      for (unsigned int i=0; i<ri.size(); i++) {
        M(ri[i]+1,c+1) = static_cast<double>(val[i]);
      }
    }
  }
  M.Release();
  return(M);
}

/////////////////////////////////////////////////////////////////////
//
// Saves matrix in a row col val format that is useful for
// exporting it to Matlab (use Matlab function spconvert).
// Is really the same as Print below, but only writes to
// file as opposed to Print that optionally prints to the
// screen.
//
/////////////////////////////////////////////////////////////////////

template<class T>
void SpMat<T>::Save(const std::string&  fname,
                    unsigned int        precision) const
{
  if (!fname.length()) throw SpMatException("SpMat::Save: Must specify filename");
  else Print(fname,precision);
}

/////////////////////////////////////////////////////////////////////
//
// Prints matrix in a row col val format that is useful for
// exporting it to Matlab (use Matlab function spconvert).
//
/////////////////////////////////////////////////////////////////////

template<class T>
void SpMat<T>::Print(const std::string&  fname,
                     unsigned int        precision) const
{
  ostream  *sptr=0;

  if (!fname.length()) {
    sptr = &cout;
  }
  else {
    try {
      sptr = new ofstream(fname.c_str());
    }
    catch(...) {
      std::string  errmsg("BFMatrix::print: Failed to write to file " + fname);
      throw SpMatException(errmsg);
    }
  }
  (*sptr) << setprecision(precision);

  for (unsigned int c=0; c<_n; c++) {
    for (unsigned int i=0; i<_ri[c].size(); i++) {
      if (_val[c][i]) (*sptr) << _ri[c][i]+1 << "  " << c+1 << "  " << _val[c][i] << endl;
    }
  }
  (*sptr) << _m << "  " << _n << "  " << 0 << endl;
  
  if (fname.length()) delete sptr;
}

/////////////////////////////////////////////////////////////////////
//
// Solves for x in expression b=(*this)*x. Uses the IML++ templates
// to obtain an iterative solution. It is presently a little stupid
// when a matrix of UNKNOWN type is passed. It will then assume worst
// case (asymmetric) rather than testing for symmetry and positive
// definiteness. That really should be changed, but at the moment
// I don't have the time.
//
/////////////////////////////////////////////////////////////////////

template<class T>
NEWMAT::ReturnMatrix SpMat<T>::SolveForx(const NEWMAT::ColumnVector&            b,
			                 MatrixType                             type,
			                 double                                 tol,
			                 unsigned int                           miter,
					 const NEWMAT::ColumnVector&            x_init) const
{
  return this->SolveForx(b,type,tol,miter,boost::shared_ptr<Preconditioner<T> >(),x_init);
}

template<class T>
NEWMAT::ReturnMatrix SpMat<T>::SolveForx(const NEWMAT::ColumnVector&            b,
			                 MatrixType                             type,
			                 double                                 tol,
			                 unsigned int                           miter,
					 boost::shared_ptr<Preconditioner<T> >  C) const
{
  NEWMAT::ColumnVector x_init;
  return this->SolveForx(b,type,tol,miter,C,x_init);
}

template<class T>
NEWMAT::ReturnMatrix SpMat<T>::SolveForx(const NEWMAT::ColumnVector&            b,
			                 MatrixType                             type,
			                 double                                 tol,
			                 unsigned int                           miter,
					 boost::shared_ptr<Preconditioner<T> >  C,
					 const NEWMAT::ColumnVector&           x_init) const
{
  if (_m != _n) throw SpMatException("SolveForx: Matrix must be square");
  if (int(_m) != b.Nrows()) throw SpMatException("SolveForx: Mismatch between matrix and vector");

  NEWMAT::ColumnVector  x(_n);
  if (x.Nrows() == x_init.Nrows()) {
    x = x_init;
  } else {
    if (x_init.Nrows()>0) {
      throw SpMatException("SolveForx: initialisation vector has incorrect size");
    } else {
      x = 0.0;
    }
  }
  int                   status = 0;
  int                   liter = int(miter);
  double                ltol = tol;
  // Use diagonal conditioner if no user-specified one
  boost::shared_ptr<Preconditioner<T> > M = boost::shared_ptr<Preconditioner<T> >();
  if (!C) M = boost::shared_ptr<Preconditioner<T> >(new DiagPrecond<T>(*this));
  else M = C;

  switch (type) {
  case SYM_POSDEF:
    status = CG(*this,x,b,*M,liter,tol);
    break;
  case SYM:
  case ASYM:
  case UNKNOWN:
    status = BiCG(*this,x,b,*M,liter,tol);
    break;
  default:
    throw SpMatException("SolveForx: No idea how you got here. But you shouldn't be here, punk.");
  }

  if (status && _pw) {
    cout << "SpMat::SolveForx: Warning requested tolerence not obtained." << endl;
    cout << "Requested tolerance was " << ltol << ", and achieved tolerance was " << tol << endl;
    cout << "This may or may not be a problem in your application, but you should look into it" << endl;
  }
  
  x.Release();
  return(x);
}

/////////////////////////////////////////////////////////////////////
//
// Returns a sparse matrix that is the transpose of *this
//
/////////////////////////////////////////////////////////////////////

template<class T>
const SpMat<T> SpMat<T>::t() const
{
  SpMat<T>         t_mat(_n,_m);
  Accumulator<T>   t_col(_n);
  for (unsigned int new_col=0; new_col<_m; new_col++) {    // For all columns of transposed matrix
    t_col.Reset();
    for (unsigned int old_col=0; old_col<_n; old_col++) {  // Search old colums for row-index corresponding to new_col
      int pos = 0;
      if (found(_ri[old_col],new_col,pos)) {
	t_col(old_col) = _val[old_col][pos];
      }
    }
    t_mat._ri[new_col].resize(t_col.NO());  
    t_mat._val[new_col].resize(t_col.NO());
    std::vector<unsigned int>&   t_mat_ri = t_mat._ri[new_col];
    std::vector<T>&              t_mat_val = t_mat._val[new_col];
    for (unsigned int i=0; i<t_col.NO(); i++) {
      t_mat_ri[i] = t_col.ri(i);
      t_mat_val[i] = t_col.val(i);
    }
    t_mat._nz += t_col.NO();
  }
  return(t_mat);  
}

/////////////////////////////////////////////////////////////////////
//
// Sets the values of an entire column, destroying any previous content.
//
/////////////////////////////////////////////////////////////////////

template<class T>
void SpMat<T>::SetColumn(unsigned int                 c,      // Column #
			 const NEWMAT::ColumnVector&  col,    // The values in that column
			 double                       eps)    // Any value <= eps is treated as a zero
{
  if (c < 1 || c > _n) throw SpMatException("SetColumn: column index out of range");
  if (static_cast<unsigned int>(col.Nrows()) != _m) throw SpMatException("SetColumn: column size mismatch");

  Accumulator<T>    acc(_m);
  double            *colp = col.Store();
  
  for (unsigned int i=0; i<_m; i++) {
    if (colp[i] > eps) acc(i) = static_cast<T>(colp[i]);
  }
  
  std::vector<unsigned int>&   ri = _ri[c-1];
  std::vector<T>&              val = _val[c-1];
  unsigned int                 old_sz = ri.size();
  if (old_sz) {
    ri = std::vector<unsigned int>(acc.NO());
    val = std::vector<T>(acc.NO());
  }
  else {
    ri.resize(acc.NO()); 
    val.resize(acc.NO());
  }
  for (unsigned int i=0; i<acc.NO(); i++) {
    ri[i] = acc.ri(i);
    val[i] = acc.val(i);
  }
  _nz += (acc.NO() - old_sz);
}

/////////////////////////////////////////////////////////////////////
//
// Returns value at position i,j (one offset)
//
/////////////////////////////////////////////////////////////////////

template<class T>
T SpMat<T>::Peek(unsigned int r, unsigned int c) const
{
  if (r<1 || r>_m || c<1 || c>_n) throw SpMatException("Peek: index out of range");

  int i=0;
  if (found(_ri[c-1],r-1,i)) return(_val[c-1][i]);

  return(static_cast<T>(0.0));
}

/////////////////////////////////////////////////////////////////////
//
// Multiply with vector x returning vector b (b = A*x)
//
/////////////////////////////////////////////////////////////////////

template<class T>
const NEWMAT::ReturnMatrix SpMat<T>::operator*(const NEWMAT::ColumnVector& x) const
{
  if (_n != static_cast<unsigned int>(x.Nrows())) throw SpMatException("operator*: # of rows in vector must match # of columns in matrix");
  NEWMAT::ColumnVector b(_m);
  b = 0.0;
  const double *xp = static_cast<double *>(x.Store());
  double       *bp = static_cast<double *>(b.Store());

  for (unsigned int c=0; c<_n; c++) {
    if (_ri[c].size()) {
      double                            wgt = xp[c];
      const std::vector<unsigned int>&  ri = _ri[c];
      const std::vector<T>&             val = _val[c];
      for (unsigned int i=0; i<ri.size(); i++) {
        bp[ri[i]] += static_cast<double>(wgt*val[i]);
      }
    }
  }
  b.Release();
  return(b);
}  

/////////////////////////////////////////////////////////////////////
//
// Multiply transpose with sparse matrix B returning matrix C (C = A'*B)
//
/////////////////////////////////////////////////////////////////////

template<class T>
const SpMat<T> SpMat<T>::TransMult(const SpMat<T>& B) const
{
  if (_m != B._m) throw SpMatException("TransMult(SpMat& ): Left hand matrix must have same # of rows as right hand");
  
  SpMat<T>        C(_n,B._n);
  Accumulator<T>  outacc(_n);
  Accumulator<T>  Bcol(B._m);

  for (unsigned int Bc=0; Bc<B._n; Bc++) {
    outacc.Reset();
    Bcol.Reset();
    Bcol.ExtractCol(B,Bc);
    for (unsigned int Ac=0; Ac<_n; Ac++) {
      const std::vector<unsigned int>&   ri = _ri[Ac];
      const std::vector<T>&              val = _val[Ac];
      T tmp = static_cast<T>(0);
      for (unsigned int i=0; i<ri.size(); i++) {
        if (Bcol.occ_at(ri[i])) {
          tmp += val[i] * Bcol.val_at(ri[i]);
	}
      }
      if (tmp) outacc(Ac) += tmp;
    }
    C._ri[Bc].resize(outacc.NO());
    C._val[Bc].resize(outacc.NO());
    std::vector<unsigned int>&   Cri = C._ri[Bc];
    std::vector<T>&              Cval = C._val[Bc];
    for (unsigned int i=0; i<outacc.NO(); i++) {
      Cri[i] = outacc.ri(i);
      Cval[i] = outacc.val(i);
    }
    C._nz += outacc.NO();
  }
 
  return(C);
}

/////////////////////////////////////////////////////////////////////
//
// Multiply transpose with vector x returning vector b (b = A'*x)
//
/////////////////////////////////////////////////////////////////////

template<class T>
const NEWMAT::ReturnMatrix SpMat<T>::trans_mult(const NEWMAT::ColumnVector& x) const
{
  if (_m != static_cast<unsigned int>(x.Nrows())) throw SpMatException("trans_mult: # of rows in vector must match # of columns in transpose of matrix");
  NEWMAT::ColumnVector b(_n);
  const double *xp = static_cast<double *>(x.Store());
  double       *bp = static_cast<double *>(b.Store());

  for (unsigned int c=0; c<_n; c++) {
    double                            res = 0.0;
    if (_ri[c].size()) {
      const std::vector<unsigned int>&  ri = _ri[c];
      const std::vector<T>&             val = _val[c];
      for (unsigned int i=0; i<ri.size(); i++) {
        res += val[i]*xp[ri[i]];
      }
    }
    bp[c] = res;
  }
  b.Release();
  return(b);
}  

/////////////////////////////////////////////////////////////////////
//
// Multiplication of self with scalar
//
/////////////////////////////////////////////////////////////////////

template<class T>
SpMat<T>& SpMat<T>::operator*=(double s)
{
  for (unsigned int c=0; c<_n; c++) {
    if (_val[c].size()) {
      std::vector<T>&             val = _val[c];
      for (unsigned int i=0; i<val.size(); i++) val[i] *= s;
    }
  }
  return(*this);  
}
      
/////////////////////////////////////////////////////////////////////
//
// Concatenates rh to right of *this
//
/////////////////////////////////////////////////////////////////////

template<class T>
SpMat<T>& SpMat<T>::operator|=(const SpMat<T>& rh)
{
  if (_m != rh._m) throw SpMatException("operator|=: Matrices must have same # of rows");
  
  _ri.resize(_n+rh._n);
  _val.resize(_n+rh._n);
  for (unsigned int c=0; c<rh._n; c++) {
    _ri[_n+c] = rh._ri[c];
    _val[_n+c] = rh._val[c];
  }
  _n += rh._n;
  _nz += rh._nz;
  
  return(*this);    
}

/////////////////////////////////////////////////////////////////////
//
// Concatenates bh below *this
//
/////////////////////////////////////////////////////////////////////

template<class T>
SpMat<T>& SpMat<T>::operator&=(const SpMat<T>& bh)
{
  if (_n != bh._n) throw SpMatException("operator&=: Matrices must have same # of columns");

  for (unsigned int c=0; c<_n; c++) {
    if ((bh._ri[c]).size()) {
      std::vector<unsigned int>&        ri = _ri[c];
      const std::vector<unsigned int>&  bhri = bh._ri[c];
      std::vector<T>&                   val = _val[c];
      const std::vector<T>&             bhval = bh._val[c];
      unsigned int                      os = ri.size();
      unsigned int                      len = bhri.size();
      ri.resize(os+len);
      val.resize(os+len);
      for (unsigned int i=0; i<len; i++) {
        ri[os+i] = _m + bhri[i];
        val[os+i] = bhval[i];
      }
    }
  }
  _m += bh._m;
  _nz += bh._nz;

  return(*this);
}

/*###################################################################
##
## Here starts global functions
##
###################################################################*/

/////////////////////////////////////////////////////////////////////
//
// Global function for multiplication of two SpMat matrices
//
/////////////////////////////////////////////////////////////////////

template<class T>
const SpMat<T> operator*(const SpMat<T>& lh, const SpMat<T>& rh)
{
  if (lh._n != rh._m) throw SpMatException("operator*: Left hand matrix must have same # of columns as right hand has rows");

  SpMat<T>              out(lh._m,rh._n);
  Accumulator<T>        acc(lh._m);
  
  for (unsigned int cr=0; cr<rh._n; cr++) {
    acc.Reset();
    if (rh._ri[cr].size()) {
      const std::vector<unsigned int>&   rri = rh._ri[cr];
      const std::vector<T>&             rval = rh._val[cr];
      for (unsigned int i=0; i<rri.size(); i++) {
        if (lh._ri[rri[i]].size()) {
          double wgt = rval[i];
          const std::vector<unsigned int>&   lri = lh._ri[rri[i]];
          const std::vector<T>&             lval = lh._val[rri[i]];
          for (unsigned int j=0; j<lri.size(); j++) {
            acc(lri[j]) += wgt*lval[j];
	  }
	}
      }
    }
    out._ri[cr].resize(acc.NO());
    out._val[cr].resize(acc.NO());
    std::vector<unsigned int>&  ori = out._ri[cr];
    std::vector<T>&            oval = out._val[cr];
    for (unsigned int i=0; i<acc.NO(); i++) {
      ori[i] = acc.ri(i);
      oval[i] = acc.val(i);
    }
    out._nz += acc.NO();
  }
               
  return(out);     
}

/////////////////////////////////////////////////////////////////////
//
// Global functions for left and right multiplication with scalar
//
/////////////////////////////////////////////////////////////////////

template<class T>
const SpMat<T> operator*(double s, const SpMat<T>& rh)
{
  return(SpMat<T>(rh) *= s);
}

template<class T>
const SpMat<T> operator*(const SpMat<T>& lh, double s)
{
  return(SpMat<T>(lh) *= s);
}

/////////////////////////////////////////////////////////////////////
//
// Global function for adding two sparse matrices
//
/////////////////////////////////////////////////////////////////////

template<class T>
const SpMat<T> operator+(const SpMat<T>& lh, const SpMat<T>& rh)
{
  return(SpMat<T>(lh) += rh);
}

/////////////////////////////////////////////////////////////////////
//
// Global function for subtracting sparse from sparse matrix
//
/////////////////////////////////////////////////////////////////////

template<class T>
const SpMat<T> operator-(const SpMat<T>& lh, const SpMat<T>& rh)
{
  return(SpMat<T>(lh) -= rh);
}

/////////////////////////////////////////////////////////////////////
//
// Global functions for horisontally concatenating sparse-sparse,
// full-sparse, sparse-full
//
/////////////////////////////////////////////////////////////////////

template<class T>
const SpMat<T> operator|(const SpMat<T>& lh, const SpMat<T>& rh)
{
  return(SpMat<T>(lh) |= rh);
}

template<class T>
const SpMat<T> operator|(const NEWMAT::GeneralMatrix& lh, const SpMat<T>& rh)
{
  return(SpMat<T>(lh) |= rh);
}

template<class T>
const SpMat<T> operator|(const SpMat<T>& lh, const NEWMAT::GeneralMatrix& rh)
{
  return(SpMat<T>(lh) |= SpMat<T>(rh));
}

/////////////////////////////////////////////////////////////////////
//
// Global function for vertically concatenating sparse-sparse,
// full-sparse and sparse-full
//
/////////////////////////////////////////////////////////////////////

template<class T>
const SpMat<T> operator&(const SpMat<T>& th, const SpMat<T>& bh)
{
  return(SpMat<T>(th) &= bh);
}

template<class T>
const SpMat<T> operator&(const NEWMAT::GeneralMatrix& th, const SpMat<T>& bh)
{
  return(SpMat<T>(th) &= bh);
}

template<class T>
const SpMat<T> operator&(const SpMat<T>& th, const NEWMAT::GeneralMatrix& bh)
{
  return(SpMat<T>(th) &= SpMat<T>(bh));
}

/*###################################################################
##
## Here starts protected functions
##
###################################################################*/

/////////////////////////////////////////////////////////////////////
//
// The following two functions give read access to _ri and _val
// vectors (corresponding to one column).
//
/////////////////////////////////////////////////////////////////////

template<class T>
const std::vector<unsigned int>& SpMat<T>::get_ri(unsigned int i) const
{
  if (i >= _n) throw SpMatException("SpMat::get_ri: Index out of range");
  return(_ri[i]);
}

template<class T>
const std::vector<T>& SpMat<T>::get_val(unsigned int i) const
{
  if (i >= _n) throw SpMatException("SpMat::get_val: Index out of range");
  return(_val[i]);
}

/*###################################################################
##
## Here starts hidden functions
##
###################################################################*/

/////////////////////////////////////////////////////////////////////
//
// Binary search. Returns true if key already exists. pos contains
// current position of key, or position to insert it in if key does
// not already exist.
//
/////////////////////////////////////////////////////////////////////

template<class T>
bool SpMat<T>::found(const std::vector<unsigned int>&  ri, unsigned int key, int& pos) const
{
  if (!ri.size() || key<ri[0]) {pos=0; return(false);}
  else if (key>ri.back()) {pos=ri.size(); return(false);}
  else {
    int mp=0;
    int ll=-1;       
    pos=int(ri.size());
    while ((pos-ll) > 1) {
      mp = (pos+ll) >> 1;   // Possibly faster than /2. Bit geeky though.
      if (key > ri[mp]) ll = mp;
      else pos = mp;
    }
  }
  if (ri[pos] == key) return(true);
  return(false);
} 

/////////////////////////////////////////////////////////////////////
//
// Return read/write reference to position i,j (one offset)
// N.B. should _not_ be used for read-only referencing since
// it will insert a value (0.0) at position i,j
//
/////////////////////////////////////////////////////////////////////

template<class T>
T& SpMat<T>::here(unsigned int r, unsigned int c)
{
  if (r<1 || r>_m || c<1 || c>_n) throw SpMatException("here: index out of range");
  
  int i = 0;
  if (!found(_ri[c-1],r-1,i)) {
    insert(_ri[c-1],i,r-1);
    insert(_val[c-1],i,static_cast<T>(0.0));
    _nz++;
  }
  return(_val[c-1][i]);
}

/////////////////////////////////////////////////////////////////////
//
// Open gap in vec at indx and fill with val. 
// Should have been templated, but I couldn't figure out how
// to, and still hide it inside SpMat
//
/////////////////////////////////////////////////////////////////////

template<class T>
void SpMat<T>::insert(std::vector<unsigned int>& vec, int indx, unsigned int val)
{
  vec.resize(vec.size()+1);
  for (int j=vec.size()-1; j>indx; j--) {
    vec[j] = vec[j-1];
  }
  vec[indx] = val;
}

template<class T>
void SpMat<T>::insert(std::vector<T>& vec, int indx, const T& val)
{
  vec.resize(vec.size()+1);
  for (int j=vec.size()-1; j>indx; j--) {
    vec[j] = vec[j-1];
  }
  vec[indx] = val;
}

/////////////////////////////////////////////////////////////////////
//
// Returns true if M has the same sparsity pattern as *this
//
/////////////////////////////////////////////////////////////////////

template<class T>
bool SpMat<T>::same_sparsity(const SpMat<T>& M) const
{
  if (_m != M._m || _n != M._n) return(false);
  for (unsigned int c=0; c<_n; c++) {
    if (_ri[c].size() != M._ri[c].size()) return(false);
  }
  for (unsigned int c=0; c<_n; c++) {
    const std::vector<unsigned int>& ri = _ri[c];
    const std::vector<unsigned int>& Mri = M._ri[c];
    for (unsigned int i=0; i<ri.size(); i++) {
      if (ri[i] != Mri[i]) return(false);
    }
  }
  return(true);
}

/////////////////////////////////////////////////////////////////////
//
// Adds a matrix to *this assuming identical sparsity patterns
//
/////////////////////////////////////////////////////////////////////

template<class T>
SpMat<T>& SpMat<T>::add_same_sparsity_mat_to_me(const SpMat<T>& M, double s)
{
  for (unsigned int c=0; c<_n; c++) {
    if (_val[c].size()) {
      std::vector<T>&          val = _val[c];
      const std::vector<T>&    Mval = M._val[c];
      for (unsigned int i=0; i<val.size(); i++) {
        val[i] += s*Mval[i];
      }
    }
  }
  return(*this);
}

/////////////////////////////////////////////////////////////////////
//
// Adds a matrix to *this assuming non-identical sparsity patterns
//
/////////////////////////////////////////////////////////////////////

template<class T>
SpMat<T>& SpMat<T>::add_diff_sparsity_mat_to_me(const SpMat<T>& M, double s)
{
  if (_m != M._m || _n != M._n) throw SpMatException("add_diff_sparsity_mat_to_me: Size mismatch between matrices");

  Accumulator<T> acc(_m);

  _nz = 0;
  for (unsigned int c=0; c<_n; c++) {
    acc.Reset();
    if (M._ri[c].size()) {
      const std::vector<unsigned int>&        Mri = M._ri[c];
      const std::vector<T>&                   Mval = M._val[c];
      for (unsigned int i=0; i<Mri.size(); i++) {
        acc(Mri[i]) += s*Mval[i];
      }
      std::vector<unsigned int>&        ri = _ri[c];
      std::vector<T>&                   val = _val[c];
      for (unsigned int i=0; i<ri.size(); i++) {
        acc(ri[i]) += s*val[i];
      }
      ri.resize(acc.NO());
      val.resize(acc.NO());
      for (unsigned int i=0; i<acc.NO(); i++) {
        ri[i] = acc.ri(i);
        val[i] = acc.val(i);
      }
      _nz += acc.NO();
    }
  }
  return(*this);
}

/*
template<class T>
SpMat<T>& SpMat<T>::add_diff_sparsity_mat_to_me(const SpMat<T>& M, double s)
{
  if (_m != M._m || _n != M._n) throw SpMatException("add_diff_sparsity_mat_to_me: Size mismatch between matrices");

  for (unsigned int c=0; c<_n; c++) {
    if (M._ri[c].size()) {
      const std::vector<unsigned int>&        Mri = M._ri[c];
      const std::vector<T>&                   Mval = M._val[c];
      for (unsigned int i=0; i<Mri.size(); i++) {
        AddTo(Mri[i]+1,c+1,s*Mval[i]);
      }
    }
  }
  return(*this);
}
*/

/*###################################################################
##
## Here starts functions for helper class Accumulator
##
###################################################################*/

template<class T>
T& Accumulator<T>::operator()(unsigned int i)
{
  if (!_occ[i]) {
    if (_sorted && _no && i < _occi[_no-1]) _sorted = false;
    _occ[i] = true;
    _occi[_no++] = i;
  }
  return(_val[i]);
}

template<class T>
const Accumulator<T>& Accumulator<T>::ExtractCol(const SpMat<T>& M, unsigned int c)
{
  if (_sz != M._m) throw ;
  if (c>(M._n-1)) throw ;
  if (_no) Reset();
  const std::vector<unsigned int>&      ri = M._ri[c];
  const std::vector<T>&                 val = M._val[c];
  for (unsigned int i=0; i<ri.size(); i++) {
    _occ[ri[i]] = true;
    _val[ri[i]] = val[i];
    _occi[_no++] = ri[i];
  }
  _sorted = true;  // Assuming M is sorted (should be)

  return(*this);
}

} // End namespace MISCMATHS

#endif // End #ifndef SpMat_h
