/*  meshfns.h

    Emma Robinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
/* MISCELLANEOUS FUNCTIONS FOR MESH PROCESSING USED BY MSM */
#if !defined(__meshfns_h)
#define __meshfns_h

#include "resampler.h"
#include "miscmaths/histogram.h"
#include <boost/shared_ptr.hpp>


using namespace boost;


namespace NEWMESH {
	
  /// TANGENT CLASS, STILL USED FOR CURRENT AFFINE IMPLEMENTATION	
  struct Tangs{
    NEWMESH::Pt e1;
    NEWMESH::Pt e2;
    
    double temp;
    double temp1;
    double temp2;
    double temp3;
    double temp4;
  };
  
  class Tangent {
    
  public:
    
    Tangent(){}
    ~Tangent(){}
    
    Tangs calculate(int,const NEWMESH::newmesh&);
    Tangs calculate_tri(const Pt &);

    Tangs calculate_tri(int,const NEWMESH::newmesh&);

    Tangs calculate2(int,const NEWMESH::newmesh&);

  };


  
  void projectPoint(const NEWMESH::Pt &, const Tangs &, double & , double & ); //checked

  Matrix get_coordinate_transformation(const double &,const double &, ColumnVector &);
  Matrix form_matrix_from_points(Pt, Pt,Pt, bool trans=false); 
  newmesh WLS_gradient(const newmesh &,const RELATIONS &, double percentile=100);
  //newmesh resample_regular_grid(newmesh &,newmesh &,newmesh &, const double &);
  newmesh projectmesh(newmesh, newmesh, newmesh ANAT = newmesh());//,const boost::shared_ptr<RELATIONS> &rel= boost::shared_ptr<RELATIONS>());
  // boost::shared_ptr<RELATIONS > _rel =boost::shared_ptr<RELATIONS >()
  newmesh calculate_strains(int, const vector<int> &, const NEWMESH::newmesh &,const NEWMESH::newmesh &, const boost::shared_ptr<Matrix>&,const boost::shared_ptr<RELATIONS>&);
  newmesh calculate_strains(double, const NEWMESH::newmesh &,const NEWMESH::newmesh &, const boost::shared_ptr<Matrix>&,const boost::shared_ptr<RELATIONS>REL=boost::shared_ptr<RELATIONS >());
  newmesh calculate_triangular_strains(const NEWMESH::newmesh &, const NEWMESH::newmesh &, double MU=0.1, double KAPPA=10.0);

  double calculate_triangular_strain(const Triangle &,const Triangle &,const double &,const double &,bool calc_emery=false, const boost::shared_ptr<ColumnVector> & indexSTRAINS=boost::shared_ptr<ColumnVector> (), bool legacy_strain = false);

  double calculate_triangular_strain(int, const NEWMESH::newmesh &, const NEWMESH::newmesh &,const double &,const double &,bool calc_emery=false, const boost::shared_ptr<ColumnVector> & indexSTRAINS=boost::shared_ptr<ColumnVector> (), bool legacy_strain = false);

  double triangle_strain(const Matrix&, const Matrix &,  const double &,const double &, const boost::shared_ptr<ColumnVector> & strains=boost::shared_ptr<ColumnVector> (), bool emery_strain=false, bool legacy_strain = false);

  void mean_curvature(double, NEWMESH::newmesh &,const boost::shared_ptr<RELATIONS>&);
  ColumnVector eig2(const double &, const double &, const double &);
  /////////////// INTENSITY NORMALISATION /////////////////
  
  Histogram build_histogram(const ColumnVector& M,const string &,const int &);// checked
  Histogram build_histogram(const ColumnVector& ,boost::shared_ptr<NEWMESH::newmesh>,const int &b);

  void histogram_normalization(NEWMESH::newmesh &,const NEWMESH::newmesh &,const string & ,const string &,int=256);// checked
  void multivariate_histogram_normalization(BFMatrix &,BFMatrix &,boost::shared_ptr<NEWMESH::newmesh> ,boost::shared_ptr<NEWMESH::newmesh>,bool=false); // checked
  void get_range(const int &, const BFMatrix &,const ColumnVector &,double & , double& );
  void set_range(const int &, BFMatrix &,const ColumnVector &, double & , double& );
 
  ///////////////////// AFFINE ROTATIONS /////////////////////////////////////////////////

  NEWMAT::ReturnMatrix  rotate_euler(const ColumnVector &,const double &,const double &,const double &);

  //// UTILITIES FOR MATRIX DATA ////
  bool set_data(const  string &,boost::shared_ptr<BFMatrix>&,NEWMESH::newmesh &, bool issparse=false);

  void logtransform(BFMatrix &);// checked
  
  void normalise(BFMatrix &);// checked
   
  void log_transform_and_normalise(BFMatrix &);

  //////////////////// MESH UNFOLD////////////////////////////////////

  bool check_for_intersections(const int , const double &, NEWMESH::newmesh &);
  NEWMESH::Pt spatialgradient(const int&,const NEWMESH::newmesh &);

  void unfold(NEWMESH::newmesh &); 
 
  //////////////////// UTILITIES FOR CALCULATION UNFOLDING GRADIENTS - based on yeo spherical demons implementation //////////////////////////////// 

  void computeNormal2EdgeOfTriangle(const NEWMESH::Pt &,const  NEWMESH::Pt &,const  NEWMESH::Pt &,  NEWMESH::Pt &); // barycentric similarity gradient calculations

  NEWMESH::Pt computeGradientOfBarycentricTriangle(const  NEWMESH::Pt &, const  NEWMESH::Pt &, const  NEWMESH::Pt&);

  ColumnVector barycentricSurfaceGrad(const int &,const NEWMESH::newmesh &); // new approach barycentric gradient  

  double barycentricGradforInd(const NEWMESH::Pt &, const NEWMESH::Pt &,const NEWMESH::Pt &,const NEWMESH::Pt &, const double &,const double &,const double &, ColumnVector &);

  /////////////////MESH CHECKING/////////////////
  ColumnVector getarealdistortion(const NEWMESH::newmesh  &, NEWMESH::newmesh  &);
  ColumnVector getarealseparation(const NEWMESH::newmesh  &, NEWMESH::newmesh  &);
  double getarealseparation(const int &, const NEWMESH::newmesh  &, NEWMESH::newmesh  &);

  ColumnVector getarealdistortionFACES(const NEWMESH::newmesh  &, NEWMESH::newmesh  &);

  /////////////////// ROTATE POINTS ON MESH ///////////////////
  Matrix estimate_rotation_matrix(const Pt &,const Pt& );
  Matrix estimate_axis_of_rotation(const Pt &,const Pt&, int & );
  Matrix rodriguez_rot(const double &,const Matrix &, const int &);
   /////////////////// GENERIC /////////////////////////////////
  bool get_all_neighbours(const int &index, vector<int> &,const NEWMESH::Pt & point, const int&,const NEWMESH::newmesh &,const RELATIONS &, SpMat<int> &); // temporary
  
  void parallel_mesh_copy(NEWMESH::newmesh &,const NEWMESH::newmesh &);

  NEWMESH::newmesh binarise_cfweighting(const NEWMESH::newmesh &);

  /////////////////////////////////////////////////////////////////////////////////////////////
  NEWMESH::newmesh create_exclusion(NEWMESH::newmesh &,const Matrix &,const float &,const float &);
  vector<string>  read_ascii_list(string);
  template <typename Iterator> inline bool next_combination(const Iterator first, Iterator k, const Iterator last) {
      if ((first == last) || (first == k) || (last == k))
         return false;
      Iterator itr1 = first;
      Iterator itr2 = last;
      ++itr1;
      if (last == itr1)
         return false;
      itr1 = last;
      --itr1;
      itr1 = k;
      --itr2;
      while (first != itr1)
      {
         if (*--itr1 < *itr2)
         {
            Iterator j = k;
            while (!(*itr1 < *j)) ++j;
            std::iter_swap(itr1,j);
            ++itr1;
            ++j;
            itr2 = k;
            std::rotate(itr1,j,last);
            while (last != j)
            {
               ++j;
               ++itr2;
            }
            std::rotate(k,itr2,last);
            return true;
         }
      }
      std::rotate(first,k,last);
      return false;
   }
}

#endif
