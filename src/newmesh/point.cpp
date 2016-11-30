/*  Copyright (C) 1999-2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include "point.h"

namespace NEWMESH {

const double operator|(const Pt &v1, const Pt &v2)
{
  return v1.X*v2.X+ v1.Y*v2.Y+ v1.Z*v2.Z;
}

const Pt operator*(const Pt &v1, const Pt &v2)
{
 
  return(Pt(v1.Y * v2.Z - v1.Z * v2.Y, 
	     v2.X * v1.Z - v2.Z * v1.X, 
	    v1.X * v2.Y - v2.X * v1.Y));  ///cross product
}

const Pt operator+(const Pt &v1, const Pt &v2)
{
  return Pt(v1.X+v2.X, v1.Y+v2.Y, v1.Z+v2.Z);
}

const Pt operator-(const Pt &v1, const Pt &v2)
{
  return Pt(v1.X-v2.X, v1.Y-v2.Y, v1.Z-v2.Z);
}

const Pt operator/(const Pt &v, const double &d)
{
  if (d!=0)
    {
      return(Pt(v.X/d, v.Y/d, v.Z/d));
    }
  else {cerr<<"newmesh::point division by zero"<<endl; return v;}
} 

const Pt operator*(const Pt &v, const double &d)
{
  return(Pt(v.X*d, v.Y*d, v.Z*d));
} 

const Pt operator*(const Matrix &M,const Pt &v)
{
  if(M.Ncols()!=3 || M.Nrows()!=3){
    cout <<M.Ncols() <<  " NEWMESH::newmesh Pt matrix multiply error: matrix should be 3x3" <<M.Nrows() << endl; exit(1);
  }

  return(Pt(M(1,1)*v.X+M(1,2)*v.Y+M(1,3)*v.Z,M(2,1)*v.X+M(2,2)*v.Y+M(2,3)*v.Z,M(3,1)*v.X+M(3,2)*v.Y+M(3,3)*v.Z));
} 

const Pt operator*(const Pt &v,const Matrix &M)
{
  if(M.Ncols()!=3 || M.Nrows()!=3){cout << " NEWMESH::newmesh Pt matrix multiply error: matrix should be 3x3" <<endl; exit(1);}

  return(Pt(v.X*M(1,1)+v.Y*M(2,1)+v.Z*M(3,1),v.X*M(1,2)+v.Y*M(2,2)+v.Z*M(3,2),v.X*M(1,3)+v.Y*M(2,3)+v.Z*M(3,3)));
}
}
