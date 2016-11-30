/*  similarities.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
/*  class for calculating similarity between rows of connectivity matrices */
/* SIM KERNEL IS COLUMN ORIENTATED TO BE CONSISTENT WITH SPARSE MAT COLUMN COMPRESSION - MESH REG ASSUMES EACH VERTEX IS REPRESENTED BY A COLUMN *////


#if !defined(similarities_h)
#define similarities_h

#include <fstream>
#include <stdio.h>
#include "newmesh/featurespace.h"


#define MAXREAL 1e6
#define MAXSIM 1e6
#define bhatmin 10e-10

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;


/* histogram class- used for information theoretic similarity measures */
/* 1d histogram class exists in miscmaths - */
namespace DISCRETEOPT{ 
class histogram2D{ 
  
    
   int m_nbinsx;
   int m_nbinsy;
   int _nsamps;
   double m_min_x;
   double m_min_y;
   double m_max_x;
   double m_max_y;
   double m_width_x;
   double m_width_y;
   Matrix _bins;
   Matrix _weights;

  public:

    //Constructor
   histogram2D(){};
   
    // Destructor
   ~histogram2D();
   inline double operator()(int i, int j){return _bins(i,j);};

   void Initialize(int, int, double, double, double, double);
   
   /// this can be used to remove bin matrix if memory overheads are a consideration (i.e. when using a large number of separate histograms as in Costfunction.cc)  
   void Reset(){_bins.ReSize(0,0);_nsamps=0;};
   void Zero();
   void ReSize(int m, int n){_bins.ReSize(m,n);_bins=0;_nsamps=0;};
    
   // add count for sample x,y
   void AddSample(double, double);
   void AddSample(double,double, double); // for weighted NMI

   // delete count for sample x,y
   void DeleteSample(double, double);

   /// Entropy calculations   
   double MarginalEntropyX();
   double MarginalEntropyY();
   double JointEntropy();

   // Normalised mutual Information
   double normalisedmutualinformation();


  };



  class  DISCRETEOPTHOCRException : public std::exception{
		
  public:
    const char* errmesg;
   
    DISCRETEOPTHOCRException(const char* msg)
      {
	errmesg=msg;
      }
    
    virtual const char* what() const throw()
    {
      std::cout<<errmesg<<std::endl;
      return errmesg;
    }
    ~DISCRETEOPTHOCRException() throw() {}
    //	private:
    
  };
	
  /// THIS CLASS ALOWS CONSTRUCTION OF A FULL NODE BY NODE SIMILARITY KERNEL OR A SPARSE SIMILARITY KERNEL, POPULATION BY NIEGHBOURHOOD CONNECTIONS ONLY (saves memory during gradient driven registrations)
  //// also holds functions for simple calculations of similarity of vectors
  class simkernel{
   
  protected:
    
    // DATA
    boost::shared_ptr<BFMatrix > m_A;
    boost::shared_ptr<BFMatrix > m_B;
        
    RowVector _rmeanA; //ROW-WISE MEAN
    RowVector _rmeanB;
        
    double _meanA;
    double _meanB;
    int _sim;
    double _thr;
    double maxx;
    double maxy;
    double minx;
    double miny;
  
    histogram2D hist; /// for NMI

    boost::shared_ptr<RELATIONS > _rel; /// remembers nearest neighbours of each input mesh vertex

    bool _issparse; // true if BFMATRIX is sparse;
    bool _initialised;
  public:

    //Constructors - options for sparse and regular matrices
    
    simkernel(){ maxx=0;maxy=0;minx=MAXREAL;miny=MAXREAL; _sim=1; _thr=0.0; _issparse=false;_initialised=false;};
    simkernel(unsigned int m, unsigned int n) { maxx=0;maxy=0;minx=MAXREAL;miny=MAXREAL;_sim=1;  _thr=0.0;_issparse=false;_initialised=false;};
    
    // Destructor
    virtual ~simkernel(){  } ;
    
    // Access as NEWMAT::Matrix
    virtual NEWMAT::ReturnMatrix AsMatrix() const = 0;
  
    // Basic properties
    virtual unsigned int Nrows() const = 0;
    virtual unsigned int Ncols() const = 0;
  
    virtual void Clear() = 0;
    virtual void Resize(unsigned int m, unsigned int n) = 0;

    // Print matrix (for debugging)
    virtual void Print(const std::string fname=std::string("")) const = 0;
    
    // Accessing
    inline double operator()(unsigned int r, unsigned int c) const {return(Peek(r,c));}
    
    virtual double Peek(unsigned int r, unsigned int c) const = 0;
    
    // Assigning
    virtual void Set(unsigned int x, unsigned int y, double val) = 0;
    void set_simval(const int & val){_sim=val;}
    inline void Set_thr(double T){_thr=T;};
    //// INITIALIZE /////////////////////////
    
    void set_input(BFMatrix &in ){  FullBFMatrix *pin = dynamic_cast<FullBFMatrix *>(&in);  if(pin) m_A = boost::shared_ptr<BFMatrix >(new FullBFMatrix (*pin)); else { SparseBFMatrix<double> *psdin = dynamic_cast<SparseBFMatrix<double> *>(&in);  m_A = boost::shared_ptr<BFMatrix >(new SparseBFMatrix<double> (*psdin)); _issparse=true;}};  
    
    void set_reference(BFMatrix &ref){ FullBFMatrix *pref = dynamic_cast<FullBFMatrix *>(&ref);  if(pref) m_B = boost::shared_ptr<BFMatrix >(new FullBFMatrix (*pref)); else { SparseBFMatrix<double> *psdref = dynamic_cast<SparseBFMatrix<double> *>(&ref);  m_B = boost::shared_ptr<BFMatrix >(new SparseBFMatrix<double> (*psdref));_issparse=true;}};
    
    void set_input(boost::shared_ptr<BFMatrix > in ){ m_A = in;};   /// for use with featurespace class - featurespace creates the data and initialises the pointer
    
    void set_reference(boost::shared_ptr<BFMatrix > ref){ m_B =ref;};// for use with featurespace class - featurespace creates the data and initialises the pointer

    void set_input(string);
  
    void set_reference(string);
    
    void set_relations(boost::shared_ptr<RELATIONS > R) {_rel= R;};

    inline  void release_data(){m_A.reset(); m_B.reset();};

    void Initialize(int);

    void reset_histogram(){hist.Reset();};  // for minimising memory overheads
    void zero_histogram(){hist.Zero();};
    void resize_histogram(int m,int n){hist.ReSize(m,n);};
   
    //////// SIMILARITY MEASURES  /////////////
    RowVector meanvector(const BFMatrix &); // for correlation measure 
    
    void calc_range(boost::shared_ptr<BFMatrix > ,double&, double&); // for histogram based measures
    void calc_range(const vector<double> &, double &, double &);
    
    /// pearsons correlation 
    double corr(int, int);
    // sum of square differences
    double SSD(int, int);
        
    /// normalised mutual information
    double NMI(int, int);
    
    ////////////// FOR VECTOR FORMAT /////////////////
    void Initialize(int,const vector<double> &, const vector<double> &,const vector<double> & weights= vector<double>()); // weighted version

    double corr(const vector<double> &, const vector<double> &, const vector<double> & weights= vector<double>());
    double corr(const map<int,float> &,const map<int,float> &);

    double corrdebug(const int &i, Matrix &, const vector<double> &, const vector<double> &, const vector<double> & weights= vector<double>());

    double SSD(const vector<double> &,const vector<double> &,const vector<double> & weights= vector<double>());

    double nSSD(const vector<double> &,const vector<double> &,const vector<double> & weights= vector<double>() );

    double NMI(const vector<double> &, const vector<double> &, const vector<double> & weights= vector<double>());
   

    vector<double> alpha_entropy(vector<vector<double> > &, vector<double> &, int &, const double &,vector<vector<int> > &);
    

    /// calculate similiarity using similarity function chosen during initialisation i.e. SSD, corr etc
    double get_sim_for_min(const vector<double> &,const vector<double> &,const vector<double>& weights= vector<double>());
    double get_sim_for_mindebug(const int &i, Matrix &, const vector<double> &,const vector<double> &,const vector<double>& weights= vector<double>());
    
    double get_sim(int, int);
    double get_sim_for_min(int , int );
    
    /////////////////////////////////////////////
    virtual void calculate_sim_column(int)  =0;  /// one input vertex at a time
  
    virtual void calculate_full_sim_kernel()  =0;
    virtual void calculate_NN_sim_kernel()=0;

  };

  
  ///// full and sparse simkernel are not currently used for discrete optimisation
  template<class T>
    class sparsesimkernel: public simkernel
  {
  
  private:
    boost::shared_ptr<MISCMATHS::SpMat<T> >    mp;
    
  public:
    
    // Constructors, destructor and assignment - copied from bfmatrix simpler way? 
  sparsesimkernel() 
    : mp(boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>())) {};
  sparsesimkernel(unsigned int m, unsigned int n) 
    : mp(boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(m,n))){};
  sparsesimkernel(unsigned int m, unsigned int n, const unsigned int *irp, const unsigned int *jcp, const double *sp)
    : mp(boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(m,n,irp,jcp,sp))){};
  sparsesimkernel(const MISCMATHS::SpMat<T>& M) 
    : mp(boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(M))){};
  sparsesimkernel(const NEWMAT::Matrix& M) 
    : mp(boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(M))) {};
    virtual ~sparsesimkernel() {};
    virtual const sparsesimkernel& operator=(const sparsesimkernel<T>& M) {
      mp = boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(*(M.mp))); return(*this);
    };

    // Access as NEWMAT::Matrix
    virtual NEWMAT::ReturnMatrix AsMatrix() const {NEWMAT::Matrix ret; ret = mp->AsNEWMAT(); ret.Release(); return(ret);};
  
    // Basic properties
    virtual unsigned int Nrows() const {return(mp->Nrows());};
    virtual unsigned int Ncols() const {return(mp->Ncols());};
    
    virtual void Clear() {mp = boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>());}
    virtual void Resize(unsigned int m, unsigned int n) {mp = boost::shared_ptr<MISCMATHS::SpMat<T> >(new MISCMATHS::SpMat<T>(m,n));}
    
    // Accessing values
    virtual double Peek(unsigned int r, unsigned int c) const {return(mp->Peek(r,c));}

    // Setting and inserting values
    virtual void Set(unsigned int x, unsigned int y, double val) {mp->Set(x,y,val);}

    // Print matrix (for debugging)
    virtual void Print(const std::string fname=std::string("")) const {mp->Print(fname);};
    // connectivity matrices will require transposing due to column compressed sparse mat format - taking path rather than object to save duplication of memory.
    

    virtual void calculate_sim_column(int);
  
    virtual void calculate_full_sim_kernel() ;
    void calculate_NN_sim_kernel();
  };

  class fullsimkernel : public simkernel
  {
  private:
    boost::shared_ptr<NEWMAT::Matrix>    mp;
  public:
    // Constructors, destructor and assignment
    fullsimkernel() {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix());};
    fullsimkernel(unsigned int m, unsigned int n) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(m,n));};
    fullsimkernel(const MISCMATHS::SpMat<double>& M) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(M.AsNEWMAT()));};
    fullsimkernel(const NEWMAT::Matrix& M) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(M));};
    virtual ~fullsimkernel() {};
    virtual const fullsimkernel& operator=(const fullsimkernel& M) {
      mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(*(M.mp))); return(*this);
    };

    virtual NEWMAT::ReturnMatrix AsMatrix() const {NEWMAT::Matrix ret; ret = *mp; ret.Release(); return(ret);};
    virtual const NEWMAT::Matrix& ReadAsMatrix() const {return(*mp);} ;

    // Basic properties
    virtual unsigned int Nrows() const {return(mp->Nrows());};
    virtual unsigned int Ncols() const {return(mp->Ncols());};
    virtual void Clear() {mp->ReSize(0,0);}
    virtual void Resize(unsigned int m, unsigned int n) {mp->ReSize(m,n);};
    
    // Print matrix (for debugging)
    virtual void Print(const std::string fname=std::string("")) const;
    
    // Accessing values
    virtual double Peek(unsigned int r, unsigned int c) const {return((*mp)(r,c));};
    
 // Setting and inserting values.
    virtual void Set(unsigned int x, unsigned int y, double val) {(*mp)(x,y)=val;}
    
    // functions for calculating full similarity kernel for whichever sim measure
    virtual void calculate_sim_column(int); 
    
    virtual void calculate_full_sim_kernel() ;
    void calculate_NN_sim_kernel();
  };
  

  template<class T>
    void sparsesimkernel<T>::calculate_sim_column(int ind){
    RowVector A, B;
    
    if(!_rel.get()){ throw   DISCRETEOPTHOCRException("SIMILARITIES:: calculate_sim_column must set RELATIONS matrix");}

    

    for (int j=1; j <= _rel->Nrows(ind); j++){
      if((*_rel)(j,ind)){
	switch (_sim){
	case 1:
	  mp->Set((int)(*_rel)(j,ind),ind,-SSD(ind,(int)(*_rel)(j,ind)));

	  break;
	case 2:	  
	  mp->Set((int)(*_rel)(j,ind),ind,corr(ind,(int)(*_rel)(j,ind)));
	  break;
	case 3:	 
	  mp->Set((int)(*_rel)(j,ind),ind,NMI(ind,(int)(*_rel)(j,ind)));	  
	  break;
	}
	
      }
    }
   
  }
 
  template<class T> 
  void sparsesimkernel<T>::calculate_NN_sim_kernel(){

 
    if(!_rel.get()){ throw   DISCRETEOPTHOCRException("SIMILARITIES:: calculate_sim_column must set RELATIONS matrix");}
	
 
    if((unsigned int) _rel->Ncols() != (unsigned int) m_A->Ncols()){
      throw   DISCRETEOPTHOCRException("SIMILARITIES::connectivity matrix and mesh have incompatible dimensions ");}
 

    for (unsigned int i=1; i <= m_A->Ncols(); i++){
  
      if((*_rel)(1,i)){
	  switch (_sim){
	  case 1:
	     mp->Set((int)(*_rel)(1,i),i,SSD(i,(int) (*_rel)(1,i)));
	    break;
	  case 2:
	    mp->Set((int)(*_rel)(1,i),i,corr(i,(int) (*_rel)(1,i)));
	    break;
	  case 3:
	    mp->Set((int)(*_rel)(1,i),i,NMI(i,(int) (*_rel)(1,i)));
	    break;
	  }
      }
    }
  }
  
  



template<class T> 
void sparsesimkernel<T>::calculate_full_sim_kernel(){

 
  if(!_rel.get()){ throw   DISCRETEOPTHOCRException("SIMILARITIES:: calculate_sim_column must set RELATIONS matrix");}
	
 
  if((unsigned int) _rel->Ncols() != (unsigned int) m_A->Ncols()){
    throw   DISCRETEOPTHOCRException("SIMILARITIES::connectivity matrix and mesh have incompatible dimensions ");}
 


  for (unsigned int i=1; i <= m_A->Ncols(); i++){
  
    for (int j=1; j <= _rel->Nrows(i); j++){
      if((*_rel)(j,i)){
	
	  switch (_sim){
	  case 1:
	     mp->Set((int)(*_rel)(j,i),i,SSD(i,(int) (*_rel)(j,i)));
	    break;
	  case 2:
	    mp->Set((int)(*_rel)(j,i),i,corr(i,(int) (*_rel)(j,i)));
	    break;
	  case 3:
	    mp->Set((int)(*_rel)(j,i),i,NMI(i,(int) (*_rel)(j,i)));
	    break;
	  }
	}
    }
  }
  
  
 }

}
#endif
