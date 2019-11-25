/*  Relations.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
/*  CLASS FOR STORING MESH VERTEX NEIGHBOURHOOD INFORMATION - distance metric is just angular separation so currently this class assumes that meshes involved are spheres*/
#ifndef RELATIONS_h
#define RELATIONS_h

#define PI    3.14159265359
#define TWOPI 6.28318530718
#define RAD 100.0

#include "newmesh.h"
#include "newmat.h"
#include "newmatio.h"
#include "miscmaths/SpMat.h"
using namespace std;
using namespace NEWMESH;

struct Mem
{
  double rank;
  int ind;
  };

struct SortAsc
{
  bool operator() (const Mem& m1, const Mem& m2)
  { return m1.rank < m2.rank; }
  };

struct SortbyIndDec
{
  bool operator() (const Mem& m1, const Mem& m2)
  { return m1.ind >  m2.ind; }
};


namespace NEWMESH {
 

  class Grid{  /// grid class subdivides all mesh vertices into cells thus speeding up searching of immediate neighbourhood, currently splits data into 30x30 bins unless the target mesh has less vertices than this, probably bin size should be better linked to desired search space
    
    int Nth,Nph;  // number of th and ph bins
    float rth,rph; // angular spacing of bins
    vector< vector<int> > cell2mesh; /// lists all mesh vertices for each 2D bin (indexed by th and ph coords) 
    vector< int >         mesh2cell;  // identifies each mesh vertex with a cell
    vector<boost::shared_ptr<Mpoint> >  m; /// save m in memory (unecessary? should really be a pointer or passed between grid functions)
    bool _debug;
    Matrix vals;
  public:
    
   

    Grid(){_debug=false;};
    
    /////////////////////////// SETTING UP GRID ////////////////////////////
    void Initialize(const newmesh & _mesh,const int& _Nth,const int& _Nph); // creates grid
    void Initialize(const vector<boost::shared_ptr<Mpoint> > &,const int&,const int&); // creates grid
    void create_luts(); // bins all mesh vertices using get_cell  
    int get_cell(const float&,const float&,const float&)const;
    
    /////////////////////////// LOCATING CLOSEST VERTICES IN GRID FOR SOME UNSEEN VERTEX (from another mesh) ////////////////////////////////
    vector<int> get_cells_in_range(const int &, int ,const int &, const int &,bool=false)const;

    vector<int> get_cell_group(int c,float ang)const; /// once grid is created an unseen point can be matched with its nearest neighbours by finding all nearbye cells  
    vector<int> get_cell_group_exact_range(int c,float ang) const;
  
    vector<int> get_points_in_cell(int c) const{return cell2mesh[c];}; 
    vector< pair<float,int> > get_points(const float& x,const float& y,const float& z,const float& ang)const;
    
    vector< vector<int> > return_cell2mesh()const{return cell2mesh;}
    void debug(){_debug=true;}
    void stopdebug(){_debug=false;}

  };
 
  class RELATIONS{  // stores neighbourhood info
  
    vector<vector<int> > mat;   // for every query vertex/triangle this will store the nearest neighbours found using the GRID
   

    Grid _grid;
    float _ang;  /// this determines the search space, all neighours within this angular separation will be saved in mat
    boost::shared_ptr< MISCMATHS::SpMat<int> > ADJ; // sparse adjacency matrix - binary matrix recording which neighbouring pairs are being considered. 
    bool _ASFACES;
    bool _EstAdj;
  public:
    
    
    RELATIONS(bool faces=false){_ang=0.0; _ASFACES=faces; _EstAdj=false;}
    RELATIONS(NEWMESH::newmesh &source,const NEWMESH::newmesh &target, double A,bool faces=false){ _ASFACES=faces; _EstAdj=false; Initialize(source,target,A);}
    RELATIONS(NEWMESH::newmesh &source,const vector<boost::shared_ptr<Mpoint> > target, double A,bool faces=false){ _ASFACES=faces; _EstAdj=false; Initialize(source,target,A);}
    RELATIONS(const string &fname){ load(fname);}
    RELATIONS(const int cols){ vector<int> tmp; for(int i=0;i<cols;i++) mat.push_back(tmp);}
    RELATIONS(const RELATIONS& R){
      mat=R.mat;
      ADJ=R.ADJ;
      _grid=R._grid;
      _ang=R._ang;
     _ASFACES=R._ASFACES;
     _EstAdj=R._EstAdj;
      //*this=R;
    }
    
    ~RELATIONS(){};
    
    inline RELATIONS&  operator=(const RELATIONS& R){ //correct assignment operator
      if (this== &R) return *this;
      else{ mat=R.mat; ADJ=R.ADJ; _grid=R._grid;_ang=R._ang;_EstAdj=R._EstAdj;return *this;}
    };
    
    inline const int &operator() (int i, int j) const{
      if((i <=0 || j<=0 ) || (j > (int) mat.size() || i>(int) mat[j-1].size())){ throw  NEWMESHException(" newmesh::RELATIONS::() Error. Relations matrix dimensions are not compatible with funciton call");}
      return mat[j-1][i-1];
    };

    void estimate_adj(int rows, int cols){ADJ= boost::shared_ptr<MISCMATHS::SpMat<int> >(new MISCMATHS::SpMat<int>(rows,cols)); _EstAdj=true;};
    //////// INITILISATION and FILLING //////////////////////////////////////////
    void Initialize(NEWMESH::newmesh &, const NEWMESH::newmesh &, const double&); // input order: source mesh, target mesh. 
    void Initialize(NEWMESH::newmesh &, const vector<boost::shared_ptr<Mpoint> > &, const double&); // input order: source mesh, target mesh. 

    ///Will set up GRID of target mesh points so that all neighbours of the source mesh vertices can be found
    vector<int> return_cell_group(const Pt &p,const double &range);
    void update_RELATIONS(const NEWMESH::newmesh &); 
    void update_RELATIONSTRI(const NEWMESH::newmesh &,const NEWMESH::newmesh &); 

    void update_RELATIONS_for_ind(int, const NEWMESH::newmesh &, double newang=0);
    void update_RELATIONS_for_tri(int, const NEWMESH::newmesh &, const NEWMESH::newmesh &, double newang=0);

    vector<pair<float,int> > update_RELATIONS_for_ind(const Pt &, double newang=0) const;
    inline void Set(const int& i,const int& j, const int &val){mat[j-1].insert(mat[j-1].begin()+i-1,val);}
    inline void Add(const int& j, const int &val){mat[j-1].push_back(val);}

    int get_closest_point_index(Pt &cr);
    void find_next_closest_neighbour( int  , int &, NEWMESH::Pt &,const NEWMESH::newmesh & );
    void update_w_querypoints_for_ind(const int & index,const vector<int> &neighbours);  // used for gradient based registration 
    //for replacing the nearest neighbours found using GRID with the vertices that make up the faces surrounding a given point 
   

    ////////////////// INVERSION /////////////////////
    RELATIONS  invert_relations(const NEWMESH::newmesh &,const NEWMESH::newmesh &,double ang=0);  /// invert such that you have all source neighbours for each target vertex
    RELATIONS  invert_relationsTR(const NEWMESH::newmesh &,const NEWMESH::newmesh &,double ang=0);  /// invert such that you have all source neighbours for each target vertex

     
  
    /////////////// READ/WRITE/ACCESS FUNCTIONS //////////////////////
    inline bool is_neighbour(int i,int j)const { if(!_EstAdj) throw  NEWMESHException("newmesh::RELATIONS::is_neighbour() Error. adjacency not initialised"); return  (ADJ->Peek(i,j)!=0);};
 
    vector<int> Col(int i)const;  /// returns all neighbours of source point i
    inline const int Ncols() const{  return mat.size();}; // should be same size as the number of vertices of the source mesh 
    inline const int Nrows(const int& i) const { return mat[i-1].size();};  /// there will be a different number of rows for each column as there are different numbers of neighbours within _ang f
    inline const double get_ang() const{return _ang;};
    inline const vector<int> get_cell_members(int c){return _grid.get_points_in_cell(c);}

    void load(const string &fname);
    void Save(const string&)const;
  };
  
  //// HELPER FNS - ALSO USED ELSEWHERE IN CODE ///////
  
  void check_scale(NEWMESH::newmesh& in,const NEWMESH::newmesh& ref);
  bool check_scale(const NEWMESH::newmesh& ,const double &);
  void true_rescale(NEWMESH::newmesh& M, const double &rad);
  Matrix recentre(NEWMESH::newmesh& M);

}


#endif
