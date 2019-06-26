/* newmesh.h

    Emma Robinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
/*  THIS CLASS CONTAINS ALL FUNCTIONS FOR READING, WRITING AND STORING SURFACES 
    AS WELL AS DATA HELD ON THE SURFACE*/
#if !defined (NEWMESH_H)
#define NEWMESH_H


#define NEW_ASCII 5
#define NEW_VTK   6
#define NEW_GIFTI 7
#define NEW_MATRIX 8
#define NEW_DPV 9

#include <stdio.h>
#include <iomanip>
#include "miscmaths/miscmaths.h"
#include "fslsurface/fslsurfaceio.h"
#include "fslsurface/fslsurface.h"
#include "point.h"
#include <boost/shared_ptr.hpp>
#include "giftiInterface.h"

using namespace fslsurface_name;
using namespace std;
using namespace MISCMATHS;

// ER Currently code does not fully implement tvalues ///////////////////////////
#define EPSILON    (1.0E-8)

namespace NEWMESH{

int  meshFileType(const string& filename);
bool meshExists(const string& filename);

 class NEWMESHException : public std::exception{
		
 public:
  const char* errmesg;
  NEWMESHException(const char* msg)
    {
      errmesg=msg;
    }
  virtual const char* what() const throw()
  {
    std::cout<<errmesg<<std::endl;
    return errmesg;
  }
  //	private:
		
 };

 class newmesh;

 /* MPOINT STORES THE VERTEX INFORMATION, THE INDICES OF
    THEIR NEIGHBOURS AND ADJOINING FACES */
 class Mpoint {  

 public:

   Mpoint(){ _coord=Pt(0, 0, 0); _no=0;};
   Mpoint(double x, double y, double z, int counter):_no(counter){ _coord=Pt(x, y, z); }
   Mpoint(const Pt p, int counter,float val=0):_no(counter){ _coord=p;}
   ~Mpoint(){}
   Mpoint(const Mpoint &p):_coord(p._coord),_no(p._no),_trID(p._trID),_nID(p._nID){*this=p;}

   Mpoint& operator=(const Mpoint& p){
    _coord=p._coord;
    _no=p._no;
    _trID=p._trID;
    _nID=p._nID;
    return(*this);
  }

  friend class newmesh;

  const Pt&    get_coord() const         {return _coord;}
  void         set_coord(const Pt& coord){ _coord=coord;}
  const int&   get_no() const            {return _no;}
 
  void         push_triangle(const int& i){_trID.push_back(i);}
  void         push_neighbour(const int& i){_nID.push_back(i);}

  int          ntriangles()const{return (int)_trID.size();}
  int          nneighbours()const{return (int)_nID.size();}

  const int&   get_trID(const int i)const{return _trID[i];}
  const int&   get_nID(const int i) const {return _nID[i];}

  inline bool is_neighbour(const int& );  // ER added, used as part of push_triangle
  bool is_triangle(const int& );

 private:

  Pt                  _coord;  /// vertex coordinate
  int                 _no;    /// counting index of vertex (where it appears in the _points array of mesh)
  vector<int>         _trID;  //// faces adjacent to vertex
  vector<int>         _nID;  /// neighbours of this vertex

  void         normalize(){_coord.normalize();}  //  used in make mesh_mesh_from_icosa() - made private to prevent accidental use in place of Pt::normalize elsewhere in the code
  
};

 class Triangle {
   /* TRIANGLE STORES THE FACE INFORMATION: ARRAY OF ALL ADJOINING VERTEX MPOINTS */
 private:
   vector<boost::shared_ptr<Mpoint> >    _vertice; 
   int               _no;
   double            _area;

 Triangle(boost::shared_ptr<Mpoint> p1, boost::shared_ptr<Mpoint> p2, boost::shared_ptr<Mpoint> p3,int no):_no(no)
   { _vertice.clear(); _vertice.push_back(p1); _vertice.push_back(p2); _vertice.push_back(p3);
     _area=area();
   }  /// newmesh should initialise the Mpoints these point to, want to avoid any chance that 2 sets of vertices are created - so made constructor private also, and newmesh a friend class.
   

 public:
   
   Triangle(){}
   
   ~Triangle(){}
   
   friend class newmesh;

 Triangle(const Triangle &t):_vertice(t._vertice),_no(t._no),_area(t._area) { *this=t;}   
   
 Triangle(Pt p1,Pt p2, Pt p3,int no):_no(no)
   { 
     boost::shared_ptr<Mpoint>  m1=boost::shared_ptr<Mpoint> ( new Mpoint (p1,0));
     boost::shared_ptr<Mpoint>  m2=boost::shared_ptr<Mpoint> ( new Mpoint (p2,1));
     boost::shared_ptr<Mpoint>  m3=boost::shared_ptr<Mpoint> ( new Mpoint (p3,2));
     _vertice.clear(); _vertice.push_back(m1); _vertice.push_back(m2); _vertice.push_back(m3);
   
     _area=area();

   }

   Triangle& operator=(const Triangle& t){
    _vertice=t._vertice;_no=t._no;return *this;
   }
 
   Triangle copy() const;
   
   void set(Pt p1,Pt p2, Pt p3,int no)
   { 
     boost::shared_ptr<Mpoint>  m1=boost::shared_ptr<Mpoint> ( new Mpoint (p1,0));
     boost::shared_ptr<Mpoint>  m2=boost::shared_ptr<Mpoint> ( new Mpoint (p2,1));
     boost::shared_ptr<Mpoint>  m3=boost::shared_ptr<Mpoint> ( new Mpoint (p3,2));
     _vertice.clear(); _vertice.push_back(m1); _vertice.push_back(m2); _vertice.push_back(m3);
     _no=no;
     _area=area();

   }

   Pt centroid() const{  
     Pt p ((_vertice[0]->get_coord().X +_vertice[1]->get_coord().X +_vertice[2]->get_coord().X)/3,
	   (_vertice[0]->get_coord().Y +_vertice[1]->get_coord().Y +_vertice[2]->get_coord().Y)/3,
	   (_vertice[0]->get_coord().Z +_vertice[1]->get_coord().Z +_vertice[2]->get_coord().Z)/3);
     
     return p;
   }
   
   Pt normal() const{ 
     Pt result=(_vertice[2]->get_coord()-_vertice[0]->get_coord())*(_vertice[1]->get_coord()-_vertice[0]->get_coord());
     result.normalize();
    return result;     
   }
  
  
   double area() const{ 
     Pt result=(_vertice[2]->get_coord()-_vertice[0]->get_coord())*(_vertice[1]->get_coord()-_vertice[0]->get_coord());
     
     return 0.5*result.norm();     
   }

   inline void swap(){   //changes triangle orientation -added from mesh as part of make_mesh_from_icosa() 
     boost::shared_ptr<Mpoint> p = _vertice[1];
    _vertice[1] = _vertice[2];
    _vertice[2] = p;      	  
  }

  bool isinside(const Pt& x)const;  
  double dist_to_point(const Pt& x0)const; 

  const double get_area() const{return _area;}

  const Mpoint& get_vertice(const int& i) const{return *_vertice[i];}

  void set_vertice(const int& i, const Pt & p) {boost::shared_ptr<Mpoint>  m=boost::shared_ptr<Mpoint> ( new Mpoint (p,0)); 
    if(i>=3) throw NEWMESHException("index exceeds triangle  dimensions");
    if(_vertice.size()==0) _vertice.resize(3);
    _vertice[i]=m;}

  const Pt& get_vertex_coord(const int& i) const{return   _vertice[i]->get_coord();}  // ER added
  const bool intersect(const vector<Pt> & p)const;          // checks if a segment intersects the triangle
  const bool intersect(const vector<Pt> & p,int& ind)const; // checks if a segment intersects the triangle+gives the index of the closest vertex
  const int get_no()const{return _no;}
  const vector<double> get_angles() const;
  const int get_vertex_no(const int& i)const{return _vertice[i]->get_no();} // ER added
};


 const bool operator ==(const NEWMESH::Mpoint &p2, const NEWMESH::Mpoint &p1);

 const bool operator ==(const NEWMESH::Mpoint &p2, const NEWMESH::Pt &p1);

 const Pt operator -(const NEWMESH::Mpoint &p1, const NEWMESH::Mpoint &p2);

 const Pt operator -(const NEWMESH::Pt &p1, const NEWMESH::Mpoint &p2);

 const bool operator ==(const NEWMESH::newmesh &M1, const NEWMESH::newmesh &M2);

class newmesh {

 private: 

  vector<boost::shared_ptr<Mpoint> >   _points;
  vector<NEWMESH::Triangle> _triangles;
  vector<Pt> _normals;

  vector<vector<float> >       _pvalues;  // enables multivariate indexed as _pvalues[feat][vertex]
  vector<vector<float> >      _tvalues;

  ////////////// META DATA /////////////////////////  note might benefit from allowing user to define some attributes, especially when loading an ascii and outputting a surf.gii?
  std::vector<GIFTImeta> global_metaData;  
  std::vector<GIFTImeta> global_Attributes; //Any other attributes that the tag has
  std::map<int,GIFTIlabel> global_GIFTIlabels;
  std::vector<GIFTIcoordinateSystem> global_defaultcoord;

 
  void push_triangle(const Triangle& t); // made private to try to safeguared initialisation of the vertex coordinates to points
  void push_point(boost::shared_ptr<Mpoint> m){_points.push_back(m);}

 public:


  newmesh();

  newmesh(const newmesh &m);

  newmesh& operator=(const newmesh& m);

  void copy_meta(const newmesh& m);

  void clear(){
    _points.clear();_triangles.clear();_pvalues.clear();_tvalues.clear();
  }

  void clear_data(){
    _pvalues.clear();_tvalues.clear();
  }

  void make_mesh_from_icosa(int n);  /// functions for generatating a regularly spaced icosahedron of order n 
  int get_ico_resolution() const ;

  std::vector<unsigned int> cluster(const double threshold,const int field=-1);
  void retessellate();  
  void retessellate(vector<vector<int> > &);
   //------------ACCESS-------------------------------------//

  int   nvertices() const{return (int)_points.size();}  
  int   ntriangles() const{return (int)_triangles.size();}
  int   npvalues(int dim=0) const{return (int)_pvalues[dim].size();}

  const NEWMESH::Mpoint&    get_point(const int& n)const{return *_points[n];}
  const NEWMESH::Pt&    get_coord(const int& n)const{if(n>=(int) _points.size() || (int) _points.size()==0){throw NEWMESHException("get_coord: index exceeds data dimensions");} return  _points[n]->get_coord();}
  const NEWMESH::Pt&    get_normal(const int& n)const {if((int)_normals.size()<(int) _points.size()){throw NEWMESHException("get_normal: normals have not been calculated, apply estimate_normals() first");}    return _normals[n];
}
  const NEWMESH::Triangle&  get_triangle(const int& n)const{if(n>= (int) _triangles.size() || (int) _triangles.size()==0){throw NEWMESHException("get_triangle: index exceeds face dimensions");}else return  _triangles[n];}
  const double  get_triangle_area(const int& n)const{if(n>= (int) _triangles.size() ||(int) _triangles.size()==0){throw NEWMESHException("get_triangle: index exceeds face dimensions");} else return  _triangles[n].get_area();}
  const Pt  get_triangle_normal(const int& n)const{if(n>=(int) _triangles.size() || (int) _triangles.size()==0){throw NEWMESHException("get_triangle: index exceeds face dimensions");} else return  _triangles[n].normal();}

  const int get_triangleID_from_vertexIDs(const int & n0, const int & n1,const int & n2)const;

  const NEWMESH::Triangle& get_triangle_from_vertex(const int & n, const int & ID)const{if(n>=(int) _triangles.size() || (int) _triangles.size()==0){throw NEWMESHException("get_triangle: index exceeds face dimensions");} else return _triangles[_points[n]->get_trID(ID)];}
  const NEWMESH::Pt&  get_triangle_vertex(const int& n, const int& i)const{if(n>=(int) _triangles.size() || (int) _triangles.size()==0){throw NEWMESHException("get_triangle: index exceeds face dimensions");} else return  _triangles[n].get_vertex_coord(i);}
  const int get_neighbour(const int & n, const int & ID)const{if(n>=(int) _points.size() || (int) _points.size()==0){throw NEWMESHException("get_coord: index exceeds data dimensions");} return _points[n]->get_nID(ID);}
  int  get_triangle_vertexID(const int& n, const int& i)const{if(n>=(int) _triangles.size() || (int) _triangles.size()==0){throw NEWMESHException("get_triangle: index exceeds face dimensions");} else return  _triangles[n].get_vertex_no(i);}

  bool is_triangle(const int &n, const int & ID)const{if(n>=(int) _points.size() || (int) _points.size()==0){throw NEWMESHException("get_coord: index exceeds data dimensions");} return _points[n]->is_triangle(ID);}

  float get_pvalue(const int& i, const int dim =0)const{if((int)_pvalues.size()<=dim || (int)_pvalues[dim].size()<i){throw NEWMESHException("get_pvalue: index exceeds data dimensions");} return _pvalues[dim][i];}
  float get_tvalue(const int& i, const int dim =0)const{if((int)_tvalues.size()<=dim || (int) _tvalues[dim].size()<i){throw NEWMESHException("get_tvalue: index exceeds data dimensions");} return _tvalues[dim][i];}
  int get_dimension()const{ return _pvalues.size();}

  int get_total_neighbours(const int& i)const{if(i>=(int) _points.size() || (int) _points.size()==0){throw NEWMESHException("get_coord: index exceeds data dimensions");} return _points[i]->nneighbours();}
  int get_total_triangles(const int& i)const{if(i>= (int) _points.size() || (int) _points.size()==0){throw NEWMESHException("get_coord: index exceeds data dimensions");} return _points[i]->ntriangles();}

  Matrix get_pvalues()const;
  vector<boost::shared_ptr<Mpoint> > get_points()const{return _points; };


  Pt local_normal(const int& pt)const;

  const vector<vector<double> > get_face_angles() const;

  vector< float > getPointsAsVectors()const;  /// accessing as vectors is used during write 
  vector< float > getValuesAsVectors(const int &d=0)const;

  vector< vector<unsigned int> > getTrianglesAsVectors()const;
  vector<int> getTrianglesAsVector()const;

  ///////////////// for identfying which face (and therefore vertices) surround a given point //////////////////////////////////////////////
  ReturnMatrix closest_triangle(int , const NEWMESH::Pt&) const ;  
  int closest_triangle_ID(int , const NEWMESH::Pt&) const ; 
 
  bool return_closest_points(int, NEWMESH::Pt&, NEWMESH::Pt&,NEWMESH::Pt&,const NEWMESH::Pt &, int &, int &, int &) const ;
   

  //-----------FOR ACCESS USING ITERATORS ----------------//   

  vector<boost::shared_ptr<Mpoint> >::const_iterator vbegin() const { return _points.begin();};   /// const to prevent reset of pointers - note this means _points can be modified
  vector<boost::shared_ptr<Mpoint> >::const_iterator vend() const {return _points.end();};

  vector<NEWMESH::Triangle>::iterator tbegin(){ return _triangles.begin();};
  vector<NEWMESH::Triangle>::const_iterator tbegin() const{ return _triangles.begin();};
  vector<NEWMESH::Triangle>::iterator tend(){ return _triangles.end();};
  vector<NEWMESH::Triangle>::const_iterator tend() const { return _triangles.end();};

  vector<int>::const_iterator nbegin(const int &i) const{ return _points[i]->_nID.begin();}; // neighbours of a particular vertex, access only
  vector<int>::const_iterator nend(const int &i) const { return _points[i]->_nID.end();};

  vector<int>::const_iterator tIDbegin(const int &i) const{ return _points[i]->_trID.begin();}; //faces associated with a particular vertex
  vector<int>::const_iterator tIDend(const int &i) const { return _points[i]->_trID.end();};

  /////-----------ASSIGNMENT --------------------------------//
  inline void set_coord(const int& i,const Pt& p){_points[i]->set_coord(p);}
  void set_pvalue(const unsigned int& i,const float& val, const int dim=0);
  void set_tvalue(const unsigned int& i,const float& val, const int dim=0);
  void set_pvalues(const vector<int>& ids,const float& val,const int dim=0);
  void set_pvalues(const Matrix& M,const bool appendFieldData=false);


  void estimate_normals();

  void initialize_pvalues(const int dim=1,const bool appendFieldData=false); // initialize new data array
  void initialize_tvalues(const int dim=1,const bool appendFieldData=false); // initialize new tvalue data array

  int step_sign(const int& vertind,const Pt& step)const;

  Pt estimate_origin();
  ////////////////  LOAD FUNCTIONS /////////////////
  ////////////////  appendFieldData refers only to pvalue/tvalue fields. As such if you are reading a .func or a .shape following a .surf appendFieldData=false. ////////

  void load(const string& filename, const bool loadSurfaceData=true, const bool appendFieldData=false); /////////// ASSUME EITHER SURFACE FILE OR DATA FILE NOT BOTH //////////////
  void load_gifti(const string& filename, const bool loadSurfaceData=true, const bool appendFieldData=false); 
  void load_ascii(const string& filename, const bool loadSurfaceData=true, const bool appendFieldData=false);  /// loads either surface or data
  void load_ascii_file(const string& filename); /// loads surface & data 
  void load_vtk(const string& filename); // cannot provide field data - surface only
  void load_matrix(const string& filename,const int& type) ;// cannot provide field data - surface only (also reads .dpv files)


  ///////////////////// WRITE FUNCTIONS ///////////////////////////////////////////
  void save(const string& filename)const;       
  void save_ascii(const string& s)const;
  void save_gifti(const string& s) const;
  void save_vtk(const string& s)const;
  void save_dpv(const string& s)const;
  void save_matrix(const string& s)const;

  // void save_gifti(const string& s,const int& type=GIFTI_ENCODING_B64GZ) const;


  void print(const string& filename);

  //  ostream& operator <<(ostream& flot);
};

/// generic function for projecting points 
 bool SameSide(const NEWMESH::Pt&, const NEWMESH::Pt&, const NEWMESH::Pt&, const NEWMESH::Pt&);  
  bool PointInTriangle(const NEWMESH::Pt& , const NEWMESH::Pt&, const NEWMESH::Pt&, const NEWMESH::Pt&);
 
 NEWMESH::Pt projectPoint(const NEWMESH::Pt &, const NEWMESH::Pt &, const NEWMESH::Pt &, const NEWMESH::Pt &, NEWMESH::Pt &);  
 void projectPoint(const NEWMESH::Pt &, const NEWMESH::Pt &, const NEWMESH::Pt &, NEWMESH::Pt &); 
 void projectVector(const NEWMESH::Pt &, const NEWMESH::Pt &, const NEWMESH::Pt &, const NEWMESH::Pt &, NEWMESH::Pt &); 


}

#endif
