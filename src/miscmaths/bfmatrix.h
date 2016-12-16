// Declarations for class BFMatrix.
//
// The purpose of class BFmatrix is to have a class from which
// to derive 2 other classes; FullBFMatrix and SparseBFMatrix.
// The reason for this is that the two classes SplineField and
// DCTField will return Hessian matrices that are either Sparse
// (SplineField) or full (DCTField). By defining a pure virtual
// class BFMatrix with a minimal (only what is needed for non-
// linear reg.) functionality I will be able to write code that
// is independent of type of matrix returned, and hence of type
// field.
//
// The syntax for the (little) functionality is sort of a mixture
// of Newmat and SparseMatrix. Mostly SparseMatrix actually.
// I hope this will not complicate the use of the nonlin package
// for those who are only interested in the full (normal) case.
//
// At one point SparseMatrix was replaced by SpMat as the underlying
// sparse matrix representation in SparseBFMatrix. SpMat was written
// with an API that largely mimicks that of NEWMAT. This is the 
// "historical" reason why a wrapper class was written, rather than
// using templatisation which would have been possible given the
// similarities in API between SpMat and NEWMAT.
//
/*    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#ifndef BFMatrix_h
#define BFMatrix_h

#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "SpMat.h"
#include "cg.h"
#include "bicg.h"

namespace MISCMATHS {


class BFMatrixException: public std::exception
{
private:
  std::string m_msg;
public:
  BFMatrixException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("BFMatrix::" + m_msg).c_str();
  }

  ~BFMatrixException() throw() {}
};

enum BFMatrixPrecisionType {BFMatrixDoublePrecision, BFMatrixFloatPrecision};

class BFMatrixColumnIterator;

class BFMatrix
{
protected:

public:
  // Constructors, destructors and stuff
  BFMatrix() {}
  BFMatrix(unsigned int m, unsigned int n) {}
  virtual ~BFMatrix() {}

  friend class BFMatrixColumnIterator;

  BFMatrixColumnIterator begin(unsigned int col) const;
  BFMatrixColumnIterator end(unsigned int col) const;

  // Access as NEWMAT::Matrix
  virtual NEWMAT::ReturnMatrix AsMatrix() const = 0;
  
  // Basic properties
  virtual unsigned int Nrows() const = 0;
  virtual unsigned int Ncols() const = 0;

  // Print matrix (for debugging)
  virtual void Print(const std::string fname=std::string("")) const = 0;

  // Setting, deleting or resizing the whole sparse matrix.
  // virtual void SetMatrix(const MISCMATHS::SpMat<double>& M) = 0;
  // virtual void SetMatrix(const MISCMATHS::SpMat<float>& M) = 0;
  virtual void SetMatrix(const NEWMAT::Matrix& M) = 0;
  virtual void Clear() = 0;
  virtual void Resize(unsigned int m, unsigned int n) = 0;

  // Accessing
  inline double operator()(unsigned int r, unsigned int c) const {return(Peek(r,c));}
  virtual double Peek(unsigned int r, unsigned int c) const = 0;
  NEWMAT::Matrix SubMatrix(unsigned int fr, unsigned int lr, unsigned int fc, unsigned int lc) const;

  // Assigning
  virtual void Set(unsigned int x, unsigned int y, double val) = 0;
  virtual void Insert(unsigned int x, unsigned int y, double val) = 0;
  virtual void AddTo(unsigned int x, unsigned int y, double val) = 0;  

  // Transpose
  virtual boost::shared_ptr<BFMatrix> Transpose() const = 0;

  // Concatenation. Note that desired polymorphism prevents us from using BFMatrix->NEWMAT::Matrix conversion
  // Concatenate two matrices yielding a third
  // AB = [*this B] in Matlab lingo
  virtual void HorConcat(const BFMatrix& B, BFMatrix& AB) const = 0;
  virtual void HorConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const = 0;
  // AB = [*this; B] in Matlab lingo
  virtual void VertConcat(const BFMatrix& B, BFMatrix& AB) const = 0;
  virtual void VertConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const = 0;
  // Concatenate another matrix to *this
  virtual void HorConcat2MyRight(const BFMatrix& B) = 0;
  virtual void HorConcat2MyRight(const NEWMAT::Matrix& B) = 0;
  virtual void VertConcatBelowMe(const BFMatrix& B) = 0;
  virtual void VertConcatBelowMe(const NEWMAT::Matrix& B) = 0;

  // Multiply by scalar
  virtual void MulMeByScalar(double s) = 0;
  // Multiply by vector
  virtual NEWMAT::ReturnMatrix MulByVec(const NEWMAT::ColumnVector& v) const = 0;
  // Add another matrix to this one 
  virtual void AddToMe(const BFMatrix& m, double s=1.0) = 0;
  // Given A*x=b, solve for x.
  virtual NEWMAT::ReturnMatrix SolveForx(const NEWMAT::ColumnVector& b,
					 MISCMATHS::MatrixType       type=SYM_POSDEF,
					 double                      tol=1e-6,
                                         int                         miter=200) const = 0;  
};

template<class T>
class SparseBFMatrix : public BFMatrix
{
private:
  boost::shared_ptr<MISCMATHS::SpMat<T> >    mp;

public:
  // Constructors, destructor and assignment
  SparseBFMatrix() 
  : mp(boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>())) {}
  SparseBFMatrix(unsigned int m, unsigned int n) 
  : mp(boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(m,n))) {}
  SparseBFMatrix(unsigned int m, unsigned int n, const unsigned int *irp, const unsigned int *jcp, const double *sp)
  : mp(boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(m,n,irp,jcp,sp))) {}
  SparseBFMatrix(const MISCMATHS::SpMat<T>& M) 
  : mp(boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(M))) {}
  SparseBFMatrix(const NEWMAT::Matrix& M) 
  : mp(boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(M))) {}
  virtual ~SparseBFMatrix() {}
  virtual const SparseBFMatrix& operator=(const SparseBFMatrix<T>& M) {
    mp = boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(*(M.mp))); return(*this);
  }

  friend class BFMatrixColumnIterator;

  // Access as NEWMAT::Matrix
  virtual NEWMAT::ReturnMatrix AsMatrix() const {NEWMAT::Matrix ret; ret = mp->AsNEWMAT(); ret.Release(); return(ret);}
  
  // Basic properties
  virtual unsigned int Nrows() const {return(mp->Nrows());}
  virtual unsigned int Ncols() const {return(mp->Ncols());}

  // Print matrix (for debugging)
  virtual void Print(const std::string fname=std::string("")) const {mp->Print(fname);}

  // Setting, deleting or resizing the whole sparse matrix.
  virtual void SetMatrix(const MISCMATHS::SpMat<T>& M) {mp = boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(M));}
  // virtual void SetMatrix(const MISCMATHS::SpMat<float>& M) {mp = boost::shared_ptr<MISCMATHS::SpMat<float> >(new MISCMATHS::SpMat<float>(M));}
  virtual void SetMatrix(const NEWMAT::Matrix& M) {mp = boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(M));}
  virtual void SetMatrixPtr(boost::shared_ptr<MISCMATHS::SpMat<T> >& mptr) {mp = mptr;}
  virtual void Clear() {mp = boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>());}
  virtual void Resize(unsigned int m, unsigned int n) {mp = boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(m,n));}

  // Accessing values
  virtual double Peek(unsigned int r, unsigned int c) const {return(mp->Peek(r,c));}

  // Setting and inserting values
  virtual void Set(unsigned int x, unsigned int y, double val) {mp->Set(x,y,val);}
  virtual void Insert(unsigned int x, unsigned int y, double val) {mp->Set(x,y,val);}
  virtual void AddTo(unsigned int x, unsigned int y, double val) {mp->AddTo(x,y,val);}

  // Transpose. 
  virtual boost::shared_ptr<BFMatrix> Transpose() const;
  
  // Concatenation of two matrices returning a third
  // AB = [*this B] in Matlab lingo
  virtual void HorConcat(const BFMatrix& B, BFMatrix& AB) const;
  virtual void HorConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const;
  // AB = [*this; B] in Matlab lingo
  virtual void VertConcat(const BFMatrix& B, BFMatrix& AB) const;
  virtual void VertConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const;

  // Concatenation of another matrix to *this
  virtual void HorConcat2MyRight(const BFMatrix& B);
  virtual void HorConcat2MyRight(const NEWMAT::Matrix& B);
  virtual void VertConcatBelowMe(const BFMatrix& B);
  virtual void VertConcatBelowMe(const NEWMAT::Matrix& B);

  // Multiply by scalar
  virtual void MulMeByScalar(double s) {(*mp)*=s;}

  // Multiply by vector
  virtual NEWMAT::ReturnMatrix MulByVec(const NEWMAT::ColumnVector& invec) const;

  // Add another matrix to this one
  virtual void AddToMe(const BFMatrix& m, double s=1.0);

  // Given A*x=b, solve for x
  virtual NEWMAT::ReturnMatrix SolveForx(const NEWMAT::ColumnVector& b,
					 MISCMATHS::MatrixType       type,
					 double                      tol,
                                         int                         miter) const;
};

class FullBFMatrix : public BFMatrix
{
private:
  boost::shared_ptr<NEWMAT::Matrix>    mp;
public:
  // Constructors, destructor and assignment
  FullBFMatrix() {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix());}
  FullBFMatrix(unsigned int m, unsigned int n) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(m,n));}
  FullBFMatrix(const MISCMATHS::SpMat<double>& M) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(M.AsNEWMAT()));}
  FullBFMatrix(const NEWMAT::Matrix& M) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(M));}
  virtual ~FullBFMatrix() {}
  virtual const FullBFMatrix& operator=(const FullBFMatrix& M) {
    mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(*(M.mp))); return(*this);
  }

  friend class BFMatrixColumnIterator;

  virtual NEWMAT::ReturnMatrix AsMatrix() const {NEWMAT::Matrix ret; ret = *mp; ret.Release(); return(ret);}
  virtual const NEWMAT::Matrix& ReadAsMatrix() const {return(*mp);} 

  // Basic properties
  virtual unsigned int Nrows() const {return(mp->Nrows());}
  virtual unsigned int Ncols() const {return(mp->Ncols());}

  // Print matrix (for debugging)
  virtual void Print(const std::string fname=std::string("")) const;

  // Setting, deleting or resizing the whole matrix.
  virtual void SetMatrix(const MISCMATHS::SpMat<double>& M) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(M.AsNEWMAT()));} 
  virtual void SetMatrix(const MISCMATHS::SpMat<float>& M) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(M.AsNEWMAT()));} 
  virtual void SetMatrix(const NEWMAT::Matrix& M) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(M));}
  virtual void SetMatrixPtr(boost::shared_ptr<NEWMAT::Matrix>& mptr) {mp = mptr;}
  virtual void Clear() {mp->ReSize(0,0);}
  virtual void Resize(unsigned int m, unsigned int n) {mp->ReSize(m,n);}

  // Accessing values
  virtual double Peek(unsigned int r, unsigned int c) const {return((*mp)(r,c));}

  // Setting and inserting values.
  virtual void Set(unsigned int x, unsigned int y, double val) {(*mp)(x,y)=val;}
  virtual void Insert(unsigned int x, unsigned int y, double val) {(*mp)(x,y)=val;}
  virtual void AddTo(unsigned int x, unsigned int y, double val) {(*mp)(x,y)+=val;}

  // Transpose. 
  virtual boost::shared_ptr<BFMatrix> Transpose() const;
  
  // Concatenation of two matrices returning a third
  virtual void HorConcat(const BFMatrix& B, BFMatrix& AB) const;
  virtual void HorConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const;
  virtual void VertConcat(const BFMatrix& B, BFMatrix& AB) const;
  virtual void VertConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const;

  // Concatenation of another matrix to *this
  virtual void HorConcat2MyRight(const BFMatrix& B);
  virtual void HorConcat2MyRight(const NEWMAT::Matrix& B);
  virtual void VertConcatBelowMe(const BFMatrix& B);
  virtual void VertConcatBelowMe(const NEWMAT::Matrix& B);

  // Multiply by scalar
  virtual void MulMeByScalar(double s);

  // Multiply by vector
  virtual NEWMAT::ReturnMatrix MulByVec(const NEWMAT::ColumnVector& invec) const;

  // Add another matrix to this one
  virtual void AddToMe(const BFMatrix& m, double s);

  // Given A*x=b, solve for x
  virtual NEWMAT::ReturnMatrix SolveForx(const NEWMAT::ColumnVector& b,
					 MISCMATHS::MatrixType       type,
					 double                      tol,
                                         int                         miter) const;
    
};

class BFMatrixColumnIterator {
public:
  BFMatrixColumnIterator(const BFMatrix& mat, unsigned int col, bool end=false) : _mat(mat), _col(col)
  {
    if (col > mat.Ncols()) throw BFMatrixException("BFMatrixColumnIterator: col out of range");
    const FullBFMatrix   *fp = dynamic_cast<const FullBFMatrix *>(&(_mat));
    if (fp) {
      if (end) _row=_mat.Nrows()+1;
      else _row=1;
      _is_sparse=false;
      _is_double=true;
    }
    else {
      const SparseBFMatrix<float> *sfp = dynamic_cast<const SparseBFMatrix<float> *>(&(_mat));
      if (sfp) {
	if (end) _sfi = new SpMat<float>::ColumnIterator(sfp->mp->end(_col));
	else _sfi = new SpMat<float>::ColumnIterator(sfp->mp->begin(_col));
        _is_sparse = true;
        _is_double = false;
      }
      else {
	const SparseBFMatrix<double> *sdp = dynamic_cast<const SparseBFMatrix<double> *>(&(_mat));
        if (sdp) {
	  if (end) _sdi = new SpMat<double>::ColumnIterator(sdp->mp->end(_col));
	  else _sdi = new SpMat<double>::ColumnIterator(sdp->mp->begin(_col));
	  _is_sparse = true;
	  _is_double = true;
	}
	else throw BFMatrixException("BFMatrixColumnIterator: No matching type for mat");
      }
    }
  }
  BFMatrixColumnIterator(const BFMatrixColumnIterator& rhs) : _mat(rhs._mat), _col(rhs._col), _is_sparse(rhs._is_sparse), _is_double(rhs._is_double) {
    if (_is_sparse) {
      if (_is_double) _sdi = new SpMat<double>::ColumnIterator(*(rhs._sdi));
      else _sfi = new SpMat<float>::ColumnIterator(*(rhs._sfi));
    }
    else _row = rhs._row; 
  }
  ~BFMatrixColumnIterator() { if (_is_sparse) { if (_is_double) free(_sdi); else free(_sfi); } }

  // Prefix case. Use this if at all possible.
  BFMatrixColumnIterator& operator++() {
    if (_is_sparse) { if (_is_double) ++(*_sdi); else ++(*_sfi); }
    else _row++;
    return(*this);
  }
  // Postfix case. Avoid.
  BFMatrixColumnIterator operator++(int dummy) {
    BFMatrixColumnIterator clone(*this);
    if (_is_sparse) { if (_is_double) ++(*_sdi); else ++(*_sfi); }
    else _row++;
    return(clone);    
  }
  bool operator==(const BFMatrixColumnIterator& rhs) const {
    if (_is_sparse!=rhs._is_sparse || _is_double!=rhs._is_double) return(false);
    if (_is_sparse) { if (_is_double) return(*_sdi==*(rhs._sdi)); else return(*_sfi==*(rhs._sfi)); }
    else {
      if (_col!=rhs._col || &_mat!=&(rhs._mat)) return(false);
      else return(_row==rhs._row);
    }
  }
  bool operator!=(const BFMatrixColumnIterator& rhs) const { return(!(*this==rhs)); }
  double operator*() const {
    if (_is_sparse) { if (_is_double) return(*(*_sdi)); else return(double(*(*_sfi))); }
    else return(_mat.Peek(_row,_col));
  }
  unsigned int Row() const { 
    if (_is_sparse) { if (_is_double) return(_sdi->Row()); else return(double(_sfi->Row())); }
    else return(_row);
  }
private:
  SpMat<double>::ColumnIterator                *_sdi; // ptr to Sparse Double Iterator
  SpMat<float>::ColumnIterator                 *_sfi; // ptr to Sparse Float Iterator
  const BFMatrix&                              _mat;
  unsigned int                                 _col;
  unsigned int                                 _row;
  bool                                         _is_sparse;
  bool                                         _is_double;
};

//
// Here comes member functions for SparseBFMatrix. Since it is templated
// these need to go here rather than in bfmatrix.cpp.
//

//
// Member functions for SparseBFMatrix
//

//
// Transpose
//

template<class T>
boost::shared_ptr<BFMatrix> SparseBFMatrix<T>::Transpose() const
{
  boost::shared_ptr<SparseBFMatrix<T> >   tm(new SparseBFMatrix<T>(mp->t()));
  return(tm);
}

//
// Concatenation of two matrices returning a third
//
template<class T>
void SparseBFMatrix<T>::HorConcat(const BFMatrix& B, BFMatrix& AB) const
{
  if (B.Nrows() && Nrows() != B.Nrows()) {throw BFMatrixException("SparseBFMatrix::HorConcat: Matrices must have same # of rows");}

  SparseBFMatrix<T> *pAB = dynamic_cast<SparseBFMatrix<T> *>(&AB);
  if (pAB) { // Means that output is sparse of type T
    *pAB = *this;
    pAB->HorConcat2MyRight(B);
  }
  else {
    FullBFMatrix *fpAB = dynamic_cast<FullBFMatrix *>(&AB);        
    if (fpAB) { // Means that output is full 
      *fpAB = FullBFMatrix(this->AsMatrix());
      fpAB->HorConcat2MyRight(B);
    }
    else throw BFMatrixException("SparseBFMatrix::HorConcat: dynamic cast error"); 
  }
}

template<class T>
void SparseBFMatrix<T>::HorConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const
{
  if (B.Nrows() && int(Nrows()) != B.Nrows()) {throw BFMatrixException("SparseBFMatrix::HorConcat: Matrices must have same # of rows");}

  SparseBFMatrix<T> *pAB = dynamic_cast<SparseBFMatrix<T> *>(&AB);   
  if (pAB) { // Means that output is sparse
    *pAB = *this;
    pAB->HorConcat2MyRight(B);
  }
  else {
    FullBFMatrix *fpAB = dynamic_cast<FullBFMatrix *>(&AB);     
    if (fpAB) {// Means that output is full
      *fpAB = FullBFMatrix(this->AsMatrix());
      fpAB->HorConcat2MyRight(B);
    }
    else throw BFMatrixException("SparseBFMatrix::HorConcat: dynamic cast error"); 
  }
}

template<class T>
void SparseBFMatrix<T>::VertConcat(const BFMatrix& B, BFMatrix& AB) const
{
  if (B.Ncols() && Ncols() != B.Ncols()) {throw BFMatrixException("SparseBFMatrix::VertConcat: Matrices must have same # of columns");}

  SparseBFMatrix<T> *pAB = dynamic_cast<SparseBFMatrix<T> *>(&AB);      
  if (pAB) { // Means that output is sparse
    *pAB = *this;
    pAB->VertConcatBelowMe(B);
  }
  else {
    FullBFMatrix *fpAB = dynamic_cast<FullBFMatrix *>(&AB);        
    if (fpAB) { // Means that output is full
      *fpAB = FullBFMatrix(this->AsMatrix());
      fpAB->VertConcatBelowMe(B);
    }
    else throw BFMatrixException("SparseBFMatrix::VertConcat: dynamic cast error"); 
  }
}

template<class T>
void SparseBFMatrix<T>::VertConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const
{
  if (B.Ncols() && int(Ncols()) != B.Ncols()) {throw BFMatrixException("SparseBFMatrix::VertConcat: Matrices must have same # of columns");}

  SparseBFMatrix<T> *pAB = dynamic_cast<SparseBFMatrix<T> *>(&AB);      
  if (pAB) { // Means that output is sparse
    *pAB = *this;
    pAB->VertConcatBelowMe(B);
  }
  else {
    FullBFMatrix *fpAB = dynamic_cast<FullBFMatrix *>(&AB);        
    if (fpAB) { // Means that output is full
      *fpAB = FullBFMatrix(this->AsMatrix());
      fpAB->VertConcatBelowMe(B);
    }
    else throw BFMatrixException("SparseBFMatrix::VertConcat: dynamic cast error"); 
  }
}

//
// Concatenate another matrix to *this
//
template<class T>
void SparseBFMatrix<T>::HorConcat2MyRight(const BFMatrix& B)
{
  if (!B.Nrows()) return;

  if (Nrows() != B.Nrows()) {throw BFMatrixException("SparseBFMatrix::HorConcat2MyRight: Matrices must have same # of rows");}

  const SparseBFMatrix<T> *pB = dynamic_cast<const SparseBFMatrix<T> *>(&B);   
  if (pB) { // Means that we want to concatenate a sparse matrix
    *mp |= *(pB->mp);
  }
  else {
    const FullBFMatrix *fpB = dynamic_cast<const FullBFMatrix *>(&B);     
    if (fpB) { // Means that we want to concatenate a full
      this->HorConcat2MyRight(fpB->AsMatrix());
    }
    else throw BFMatrixException("SparseBFMatrix::HorConcat2MyRight: dynamic cast error");
  }
}

template<class T>
void SparseBFMatrix<T>::HorConcat2MyRight(const NEWMAT::Matrix& B)
{
  if (!B.Nrows()) return;

  if (int(Nrows()) != B.Nrows()) {throw BFMatrixException("SparseBFMatrix::HorConcat2MyRight: Matrices must have same # of rows");}
  *mp |= B;
}

template<class T>
void SparseBFMatrix<T>::VertConcatBelowMe(const BFMatrix& B)
{
  if (!B.Ncols()) return;

  if (Ncols() != B.Ncols()) {throw BFMatrixException("SparseBFMatrix::VertConcatBelowMe: Matrices must have same # of columns");}

  const SparseBFMatrix<T> *pB = dynamic_cast<const SparseBFMatrix<T> *>(&B);   
  if (pB) { // Means that we want to concatenate a sparse matrix
    *mp &= *(pB->mp);
  }
  else {
    const FullBFMatrix *fpB = dynamic_cast<const FullBFMatrix *>(&B);     
    if (fpB) { // Means that we want to concatenate a full
      this->VertConcatBelowMe(fpB->AsMatrix());
    }
    else throw BFMatrixException("SparseBFMatrix::VertConcatBelowMe: dynamic cast error");
  }
}

template<class T>
void SparseBFMatrix<T>::VertConcatBelowMe(const NEWMAT::Matrix& B)
{
  if (!B.Ncols()) return;

  if (int(Ncols()) != B.Ncols()) {throw BFMatrixException("SparseBFMatrix::VertConcatBelowMe: Matrices must have same # of columns");}
  *mp &= B;
}

// Multiply by vector
template<class T>
NEWMAT::ReturnMatrix SparseBFMatrix<T>::MulByVec(const NEWMAT::ColumnVector& invec) const
{
  if (invec.Nrows() != int(Ncols())) {throw BFMatrixException("Matrix-vector size mismatch");}
  NEWMAT::ColumnVector   outvec = *mp * invec;
  outvec.Release();
  return(outvec);
}

// Add another matrix to this one
template<class T>
void SparseBFMatrix<T>::AddToMe(const BFMatrix& M, double s)
{
  if (Ncols() != M.Ncols() || Nrows() != M.Nrows()) {
    throw BFMatrixException("SparseBFMatrix::AddToMe: Matrix size mismatch");
  }

  const SparseBFMatrix<T> *pM = dynamic_cast<const SparseBFMatrix<T> *>(&M);  
  if (pM) { // Add sparse matrix to this sparse matrix
    if (s == 1.0) *mp += *(pM->mp);
    else *mp += s * *(pM->mp);
  }
  else {
    const FullBFMatrix *fpM = dynamic_cast<const FullBFMatrix *>(&M);      
    if (fpM) { // Add full matrix to this sparse matrix
      if (s == 1.0) *mp += SpMat<T>(fpM->ReadAsMatrix());
      else *mp += s * SpMat<T>(fpM->ReadAsMatrix());
    }
    else throw BFMatrixException("SparseBFMatrix::AddToMe: dynamic cast error"); 
  }
}

// Given A*x=b, solve for x
template<class T>
NEWMAT::ReturnMatrix SparseBFMatrix<T>::SolveForx(const NEWMAT::ColumnVector& b,
					          MISCMATHS::MatrixType       type,
					          double                      tol,
                                                  int                         miter) const
{
  if (b.Nrows() != int(Nrows())) {
    throw BFMatrixException("SparseBFMatrix::SolveForx: Matrix-vector size mismatch");
  }
  NEWMAT::ColumnVector  x = mp->SolveForx(b,type,tol,miter);
  x.Release();
  return(x);
}

} // End namespace MISCMATHS

#endif // End #ifndef BFMatrix_h
