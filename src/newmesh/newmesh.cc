/*  newmesh.cc

    Emma Robinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "newmesh.h"
#include <climits>

namespace NEWMESH{
//////// CSVMESH
// 1: ascii, 2:vtk, 3: gii, 4 as .txt (for supplyng data in a text file) -1: unknown
  int  meshFileType(const string& filename){
    if(filename.size()<=5){return -1;}
    string last_3 = filename.substr(filename.size()-3, 3);
    if( last_3 == ".gz" ){
      last_3 = filename.substr(filename.size()-6, 3);
      if( last_3 == "gii" ) {return NEW_GIFTI;}
    }
    
    if( last_3 == "gii" ) {return NEW_GIFTI;}
    else if( last_3 == "txt" ){return NEW_MATRIX;}
    else if ( last_3=="dpv" ){return NEW_DPV;}
    else if ( last_3=="asc" ){return NEW_ASCII;}
    else if ( last_3=="vtk" ){return NEW_VTK;}

    ifstream f(filename.c_str());
    //reading the header
    string header;
    getline(f, header);
    {
      string::size_type pos = header.find("# vtk DataFile Version");
      if (pos!= string::npos) {
	f.close();
	return NEW_VTK;
      }
    }
    {
      string::size_type pos = header.find("#!ascii");
      if (pos != string::npos) {
	f.close();
	return NEW_ASCII;
      }
    }
    return -1;
  }

  bool meshExists(const string& filename){
    int type = meshFileType(filename);
    if(type>0){return true;}
    return false;
  }


  ////// MPOINT FNS ////////////////////

  bool Mpoint::is_neighbour(const int& j){   //  used as part of push_triangle
    bool found=false;
    vector<int>::iterator it;
    it=find (_nID.begin(), _nID.end(), j);   
    
    if (it!=_nID.end()) found=true;
    
    return found;
  }  
  
  bool Mpoint::is_triangle(const int& j){   //  used as part of push_triangle
    bool found=false;
    vector<int>::iterator it;
    it=find (_trID.begin(), _trID.end(), j);   
    
    if (it!=_trID.end()) found=true;
    
    return found;
  }  
  /////         TRIANGLE FNS /////////////////

  Triangle Triangle::copy() const{  
    
    Triangle t;

    for(unsigned int i=0;i<_vertice.size();i++){
      boost::shared_ptr<Mpoint> pt=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(*_vertice[i]));
      t._vertice.push_back(pt);
    }
    t._no=_no;

    return t;
  
  }
  
  //Triangle& Triangle::operator=(const Triangle& t){
  //  _vertice.clear();
    
  //  for(unsigned int i=0;i<_vertice.size();i++){
  //    boost::shared_ptr<Mpoint> pt=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(t.get_vertice(i)));
  //    _vertice.push_back(pt);
  //  }
    
  //  _no=t._no;
  //  return *this;
  // }


  const vector<double> Triangle::get_angles() const{ // get angles in order of vertex: 0,1,2
    Pt v0,v1,v2;
    vector<double> face_angles;

    v0=_vertice[2]->get_coord() -_vertice[0]->get_coord();  /// edge from vertex 0 to 2
    v1=_vertice[1]->get_coord() -_vertice[0]->get_coord();  /// edge from vertex 0 to 1
    v2=_vertice[2]->get_coord() -_vertice[1]->get_coord();  /// edge from vertex 1 to 2
    double dot01=v0|v1; 
    double dot02=v0|v2;
    double dot12=v1|v2;

    double angle0=acos(dot01/(v0.norm()*v1.norm())); face_angles.push_back(angle0);
    double angle1=acos(dot02/(v0.norm()*v2.norm())); face_angles.push_back(angle1);
    double angle2=acos(dot12/(v2.norm()*v1.norm())); face_angles.push_back(angle2);
    return face_angles;
  }

  bool Triangle::isinside(const NEWMESH::Pt& x)const{
    Pt v0,v1,v2;
    v0=_vertice[2]->get_coord() -_vertice[0]->get_coord();
    v1=_vertice[1]->get_coord() -_vertice[0]->get_coord();
    v2=x-_vertice[0]->get_coord();
    double dot00=v0|v0;
    double dot01=v0|v1;
    double dot02=v0|v2;
    double dot11=v1|v1;
    double dot12=v1|v2;
    double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    
    // Check if point is in triangle
    return (u > 0) && (v > 0) && (u + v < 1);
    
  }

  double Triangle::dist_to_point(const NEWMESH::Pt& x0)const{
    double d;
    Pt x1(_vertice[0]->get_coord().X,_vertice[0]->get_coord().Y,_vertice[0]->get_coord().Z);
    Pt x2(_vertice[1]->get_coord().X,_vertice[1]->get_coord().Y,_vertice[1]->get_coord().Z);
    Pt x3(_vertice[2]->get_coord().X,_vertice[2]->get_coord().Y,_vertice[2]->get_coord().Z);
    Pt u;
    double dmin=1000000;
    // test edges
    u=x2-x1;
    if( ((x0-x1)|u)>0 && ((x0-x2)|u)<0 ){
      d=(((x0-x1)*(x0-x2)).norm()/(x2-x1).norm());
      if(d<dmin)dmin=d;
    }
    u=x3-x1;
    if( ((x0-x1)|u)>0 && ((x0-x3)|u)<0 ){
      d=(((x0-x1)*(x0-x3)).norm()/(x3-x1).norm());  
      if(d<dmin)dmin=d;
    }
    u=x3-x2;
    if( ((x0-x2)|u)>0 && ((x0-x3)|u)<0 ){
      d=(((x0-x2)*(x0-x3)).norm()/(x3-x2).norm());
      if(d<dmin)dmin=d;
    }
  
    d=(x0-x1).norm();if(d<dmin)dmin=d;
    d=(x0-x2).norm();if(d<dmin)dmin=d;
    d=(x0-x3).norm();if(d<dmin)dmin=d;
    return dmin;
  }

  // Saad's implementations - not used in MSM instead return_closest_triangle is used to find individual vertices that intersect the plane see return_closest_triangle
  // algorithm from:
  // http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()

  const bool Triangle::intersect(const vector<NEWMESH::Pt> & p) const {
    Pt    u,v,n;   // triangle vectors
    Pt    dir,w0,w; // ray vectors
    double r, a, b;              // params to calc ray-plane intersect

    // check if point is one the vertices
    for(int ii=0;ii<=2;ii++){
      if((*_vertice[ii])==p[0])return true;
      if((*_vertice[ii])==p[1])return true;
    }

    // get triangle edge vectors and plane normal
    u = _vertice[1]->get_coord()-_vertice[0]->get_coord();
    v = _vertice[2]->get_coord()-_vertice[0]->get_coord();
    n = u*v;             // cross product
    if (n.norm()==0) // triangle is degenerate
      return false;                 
    

    dir = p[1]-p[0];             // ray direction vector
    w0 = p[0]-_vertice[0]->get_coord();
    a = -(n|w0)/n.norm()/w0.norm();
    b = (n|dir)/n.norm()/dir.norm();
    if (fabs(b) < 0.001) { // ray is parallel to triangle plane
      if (fabs(a) < 0.001)                 // ray lies in triangle plane
	return true;
      else return false;             // ray disjoint from plane
    }
    
    // get intersect point of ray with triangle plane
    r = a / b;
    if (r < 0.0)                   // ray goes away from triangle
      return false;                  // => no intersect
    if(r > 1.0)
      return false;
    // for a segment, also test if (r > 1.0) => no intersect
    Pt I;
    Pt tmp;
    //tmp=r * dir;
    I = p[0] + (dir*r);           // intersect point of ray and plane
    
    // is I inside T?
    double    uu, uv, vv, wu, wv, D;
    uu = (u|u);
    uv = (u|v);
    vv = (v|v);
    w = I - _vertice[0]->get_coord();
    wu = (w|u);
    wv = (w|v);
    D = uv * uv - uu * vv;
    
    // get and test parametric coords
    double s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0)        // I is outside T
      return false;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  // I is outside T
      return false;
    
    return true;                      // I is in T
    
  }
  
  
  const bool Triangle::intersect(const vector<NEWMESH::Pt> & p,int& ind) const {
    Pt    u,v,n;   // triangle vectors
    Pt    dir,w0,w; // ray vectors
    double r, a, b;              // params to calc ray-plane intersect

    // check if point is one the vertices
    for(int ii=0;ii<=2;ii++){
      if((*_vertice[ii])==p[0]){ind=ii;return true;}
      if((*_vertice[ii])==p[1]){ind=ii;return true;}
    }

    // get triangle edge vectors and plane normal
    u = _vertice[1]->get_coord()-_vertice[0]->get_coord();
    v = _vertice[2]->get_coord()-_vertice[0]->get_coord();
    n = u*v;             // cross product
    if (n.norm()==0) // triangle is degenerate
      return false;                 
    

    dir = p[1]-p[0];             // ray direction vector
    w0 = p[0]-_vertice[0]->get_coord();
    a = -(n|w0);
    b = (n|dir);
    if (fabs(b) < 0.001) { // ray is parallel to triangle plane
      if (fabs(a) < 0.001)                 // ray lies in triangle plane
	{ind=0;return true;}
      else return false;             // ray disjoint from plane
    }
    
    // get intersect point of ray with triangle plane
    r = a / b;
    if (r < 0.0)                   // ray goes away from triangle
      return false;                  // => no intersect
    if(r > 1.0)
      return false;
    // for a segment, also test if (r > 1.0) => no intersect
    Pt I;
    I = p[0] + dir*r;           // intersect point of ray and plane
    
    // is I inside T?
    double    uu, uv, vv, wu, wv, D;
    uu = (u|u);
    uv = (u|v);
    vv = (v|v);
    w = I - _vertice[0]->get_coord();
    wu = (w|u);
    wv = (w|v);
    D = uv * uv - uu * vv;
    
    // get and test parametric coords
    double s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0)        // I is outside T
      return false;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  // I is outside T
      return false;

    // which vertex is closest to where the segment intersects?
    float x=uu-2*wu,y=vv-2*wv;
    if( x<0 ){
      if( x<y ) ind=1;
      else ind=2;
    }
    else{
      if( y<0 ) ind=2;
      else ind=0;
    }
    
    return true;                      // I is in T
    
  }

  
  /// ---operators .-----------

  const bool operator ==(const NEWMESH::Mpoint &p2, const NEWMESH::Mpoint &p1){
   
    return (fabs(p1.get_coord().X- p2.get_coord().X)<1e-4 && fabs(p1.get_coord().Y - p2.get_coord().Y)<1e-4 && fabs(p1.get_coord().Z - p2.get_coord().Z)<1e-4);}


  const bool operator ==(const NEWMESH::Mpoint &p2, const NEWMESH::Pt &p1){ return (fabs(p1.X- p2.get_coord().X)<1e-3 && fabs(p1.Y - p2.get_coord().Y)<1e-3 && fabs(p1.Z - p2.get_coord().Z)<1e-3);}

  const Pt operator -(const NEWMESH::Mpoint &p1, const NEWMESH::Mpoint &p2){
   return Pt (p1.get_coord().X - p2.get_coord().X,p1.get_coord().Y - p2.get_coord().Y,p1.get_coord().Z - p2.get_coord().Z );
 }

  const Pt operator -(const NEWMESH::Pt &p1, const NEWMESH::Mpoint &p2){
   return Pt (p1.X - p2.get_coord().X,p1.Y - p2.get_coord().Y,p1.Z - p2.get_coord().Z );
 }

  const bool operator ==(const NEWMESH::newmesh &M1, const NEWMESH::newmesh &M2){ 
    bool isequal=1; Mpoint p1,p2; 
    if(M1.nvertices()!=M2.nvertices()) return false;
    for (int i=0;i< M1.nvertices();i++){
      p1=M1.get_point(i);p2=M2.get_point(i); 
      if(p1==p2) isequal=1;else{ isequal=0; break;}
    } 
    return isequal;}


  //////////// NEWMESH ///////////////////////////
  newmesh::newmesh(){
    /// SET DEFAULT COORD SYSTEM
    global_defaultcoord.clear();  
    string default_space="NIFTI_XFORM_TALAIRACH";
    std::vector<double> transform(16,0);
    //////// set as identity
    transform[0]=1; transform[5]=1; transform[10]=1; transform[15]=1;
    GIFTIcoordinateSystem default_coord(default_space, default_space, transform); // save out a default coord system
    global_defaultcoord.push_back(default_coord); 
  }


  newmesh::newmesh(const newmesh &m):_normals(m._normals),_pvalues(m._pvalues),_tvalues(m._tvalues),global_metaData(m.global_metaData),global_Attributes(m.global_Attributes),global_GIFTIlabels(m.global_GIFTIlabels),global_defaultcoord(m.global_defaultcoord){ // do neigglobal_metaData(m.global_metaData),global_Attributes(m.global_Attributes),global_GIFTIlabels(m.global_GIFTIlabelhbours get copied across


    _points.clear();
    _triangles.clear(); 
    _normals.clear(); 
    
    for (vector<boost::shared_ptr<Mpoint> >::const_iterator p= m._points.begin(); p!=m._points.end(); p++)
      {
	boost::shared_ptr<Mpoint> pt=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(**p));
	_points.push_back(pt);
	
      }

    for (vector<Triangle>::const_iterator t= m._triangles.begin(); t!=m._triangles.end(); t++)
      {
	int v0 = t->get_vertex_no(0), v1 = t->get_vertex_no(1), v2 = t->get_vertex_no(2);
	Triangle tr(_points[v0], _points[v1], _points[v2],t->get_no());
	_triangles.push_back(tr);
      }

    for(int i=0;i<(int) _points.size();i++)
      _normals.push_back(local_normal(i));
    
  }

  newmesh& newmesh::operator=(const newmesh& m){

  
    _points.clear();
    _triangles.clear();  
    _normals=m._normals;
    _pvalues=m._pvalues;
    _tvalues=m._tvalues;


    global_metaData=m.global_metaData;
    global_Attributes=m.global_Attributes;
    global_GIFTIlabels=m.global_GIFTIlabels;
    global_defaultcoord=m.global_defaultcoord;
   

    for (vector<boost::shared_ptr<Mpoint> >::const_iterator p= m._points.begin(); p!=m._points.end(); p++)
    {
      boost::shared_ptr<Mpoint> pt=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(**p));
      _points.push_back(pt);

    }

    for (vector<Triangle>::const_iterator t= m._triangles.begin(); t!=m._triangles.end(); t++)
      {
	int v0 = t->get_vertex_no(0), v1 = t->get_vertex_no(1), v2 = t->get_vertex_no(2);
	Triangle tr(_points[v0], _points[v1], _points[v2],t->get_no());
	_triangles.push_back(tr);
      }

    return(*this);
  }



  void newmesh::copy_meta(const newmesh& m){

    global_metaData=m.global_metaData;
    global_Attributes=m.global_Attributes;
    global_GIFTIlabels=m.global_GIFTIlabels;
    global_defaultcoord=m.global_defaultcoord;
  }

  void newmesh::make_mesh_from_icosa(int n)   
  {
    clear();
  
    const double tau=0.8506508084;	
    const double one=0.5257311121;
  
    //creates regular icosahedron ////////////  
    boost::shared_ptr<Mpoint> ZA=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint( tau , one,  0,    0));
    boost::shared_ptr<Mpoint> ZB=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint( -tau, one,  0,    1));
    boost::shared_ptr<Mpoint> ZC=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint( -tau, -one, 0,    2));
    boost::shared_ptr<Mpoint> ZD=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint( tau,  -one, 0,    3));
    boost::shared_ptr<Mpoint> YA=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint( one,  0 ,   tau,  4));
    boost::shared_ptr<Mpoint> YB=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint( one,  0,    -tau, 5));
    boost::shared_ptr<Mpoint> YC=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint( -one, 0,    -tau, 6));
    boost::shared_ptr<Mpoint> YD=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint( -one, 0,    tau,  7));
    boost::shared_ptr<Mpoint> XA=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint( 0 ,   tau,  one,  8));
    boost::shared_ptr<Mpoint> XB=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint( 0,    -tau, one,  9));
    boost::shared_ptr<Mpoint> XC=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint( 0,    -tau, -one, 10));
    boost::shared_ptr<Mpoint> XD=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint( 0,    tau,  -one, 11));
    
    
    Triangle t0(YD, XA, YA,0);
    Triangle t1(XB, YD, YA,1);
    Triangle t2(XD, YC, YB,2);
    Triangle t3(YC, XC, YB,3);
    Triangle t4(ZD, YA, ZA,4);
    Triangle t5(YB, ZD, ZA,5);
    Triangle t6(ZB, YD, ZC,6);
    Triangle t7(YC, ZB, ZC,7);
    Triangle t8(XD, ZA, XA,8);
    Triangle t9(ZB, XD, XA,9);
    Triangle t10(ZD, XC, XB,10);
    Triangle t11(XC, ZC, XB,11);
    Triangle t12(ZA, YA, XA,12);
    Triangle t13(YB, ZA, XD,13);
    Triangle t14(ZD, XB, YA,14);
    Triangle t15(XC, ZD, YB,15);
    Triangle t16(ZB, XA, YD,16);
    Triangle t17(XD, ZB, YC,17);
    Triangle t18(XB, ZC, YD,18);
    Triangle t19(ZC, XC, YC,19);
 
    _points.push_back(ZA);
    _points.push_back(ZB);
    _points.push_back(ZC);
    _points.push_back(ZD);
    _points.push_back(YA);
    _points.push_back(YB);
    _points.push_back(YC);
    _points.push_back(YD);
    _points.push_back(XA);
    _points.push_back(XB);
    _points.push_back(XC);
    _points.push_back(XD);
    
    
    push_triangle(t0);
    push_triangle(t1);
    push_triangle(t2);
    push_triangle(t3);
    push_triangle(t4);
    push_triangle(t5);
    push_triangle(t6);
    push_triangle(t7);
    push_triangle(t8);
    push_triangle(t9);
    push_triangle(t10);
    push_triangle(t11);
    push_triangle(t12);
    push_triangle(t13);
    push_triangle(t14);
    push_triangle(t15);
    push_triangle(t16);
    push_triangle(t17);
    push_triangle(t18);
    push_triangle(t19);
    
    for (vector<NEWMESH::Triangle>::iterator i= tbegin(); i!=tend(); i++){
      i->swap();  //changes triangle orientation
    }
    
    
    for (int io = 0; io<n ; io++)
      {
	//re-tesselates (subdivides)
	retessellate();
	
	
      }
    
    
    vector<float> tmp_pvalues(_points.size(),0);
    _pvalues.push_back(tmp_pvalues);


    vector<float> tmp_tvalues(_triangles.size(),0);
    _tvalues.push_back(tmp_tvalues);

    
  }


  int newmesh::get_ico_resolution() const {
    int ico=0;
    int n=_points.size();
   
    if(n==42)
      ico=1;
    if(n==162)
      ico=2;
    if(n==642)
      ico=3;
    if(n==2562)
      ico=4;
    if(n==10242)
      ico=5;
    if(n==40962)
      ico=6;

    if(ico==0)  { throw  NEWMESHException(" NEWMESH::newmesh::get_ico_resolution mesh is not an icosphere or has not been initialised ");	}
    return ico;   
  }
 
  void newmesh::retessellate() {
    
    vector<boost::shared_ptr<Mpoint> > added_points; // need added points to point to same memory as _points
    vector<NEWMESH::Triangle> tr = _triangles;
    
    added_points.clear();
    
    int count=0;
    int tot_triangles=0;
    
    _triangles.clear();
    
    for (unsigned int i=0;i<_points.size();i++)
      _points[i]->_nID.clear();
    
    for (unsigned int i=0;i<_points.size();i++)
      _points[i]->_trID.clear();
  
    
    for (vector<NEWMESH::Triangle>::iterator t=tr.begin(); t!=tr.end(); t++)
      {
	
	boost::shared_ptr<Mpoint> v0=t->_vertice[0];
	boost::shared_ptr<Mpoint> v1=t->_vertice[1];
	boost::shared_ptr<Mpoint> v2=t->_vertice[2];
	NEWMESH::Pt pt0((v1->_coord.X + v2->_coord.X)/2, (v1->_coord.Y + v2->_coord.Y)/2, (v1->_coord.Z + v2->_coord.Z)/2); // adds new points at the center of existing faces
	NEWMESH::Pt pt2((v0->_coord.X + v1->_coord.X)/2, (v0->_coord.Y + v1->_coord.Y)/2, (v0->_coord.Z + v1->_coord.Z)/2);
	NEWMESH::Pt pt1((v0->_coord.X + v2->_coord.X)/2, (v0->_coord.Y + v2->_coord.Y)/2, (v0->_coord.Z + v2->_coord.Z)/2);
	
	
	boost::shared_ptr<Mpoint> p1,p2,p0;
	
	bool b0=true, b1=true, b2=true;
	count=0;
	int index=0;
	for (vector<boost::shared_ptr<Mpoint> >::const_iterator mi=added_points.begin(); mi!=added_points.end(); mi++)
	  {
	    index++;
	    NEWMESH::Pt current = (*mi)->_coord;
	    if (pt0==current) {b0=false; p0=*mi;}
	    if (pt1==current) {b1=false; p1=*mi;}
	    if (pt2==current) {b2=false; p2=*mi;}
	  }
	
	
	if (b0) {p0=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(pt0, nvertices() + count)); count++;};
	if (b1) {p1=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(pt1, nvertices() + count)); count++;};
	if (b2) {p2=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(pt2, nvertices() + count)); count++;}; 
	
	if (b0) {_points.push_back(p0);added_points.push_back(p0);}  // ER moved, neighbours must be added before triangles are pushed back
	if (b1) {_points.push_back(p1);added_points.push_back(p1);}
	if (b2) {_points.push_back(p2);added_points.push_back(p2);}
	
	Triangle t0(p2, p0, p1,tot_triangles);  tot_triangles++;
	Triangle t1(p1, v0, p2,tot_triangles);  tot_triangles++;
	Triangle t2(p0, v2, p1,tot_triangles);  tot_triangles++;
	Triangle t3(p2, v1, p0,tot_triangles);  tot_triangles++;

	//Triangle t0(p1, p0, p2,tot_triangles);  tot_triangles++;
	//Triangle t1(p2, v0, p1,tot_triangles);  tot_triangles++;
	//Triangle t2(p1, v2, p0,tot_triangles);  tot_triangles++;
	//Triangle t3(p0, v1, p2,tot_triangles);  tot_triangles++;

	push_triangle(t0);
	push_triangle(t1);
	push_triangle(t2); 
	push_triangle(t3);
	
	
    }
    
    for (vector<boost::shared_ptr<Mpoint> >::const_iterator i=vbegin(); i!=vend(); i++)
	(*i)->normalize();
    
  }
  
  void newmesh::retessellate(vector<vector<int> > &OLD_tr_neighbours) {
    
    vector<boost::shared_ptr<Mpoint> > added_points; // need added points to point to same memory as _points
    vector<NEWMESH::Triangle> tr = _triangles;
    vector<int> tmp;

    added_points.clear();
    
    int count=0;
    int tot_triangles=0;
    
    _triangles.clear();
    _pvalues.clear(); 

    for (unsigned int i=0;i<_points.size();i++)
      _points[i]->_nID.clear();
    
    for (unsigned int i=0;i<_points.size();i++)
      _points[i]->_trID.clear();
  
    
    for (vector<NEWMESH::Triangle>::iterator t=tr.begin(); t!=tr.end(); t++)
      {
	OLD_tr_neighbours.push_back(tmp);
	boost::shared_ptr<Mpoint> v0=t->_vertice[0];
	boost::shared_ptr<Mpoint> v1=t->_vertice[1];
	boost::shared_ptr<Mpoint> v2=t->_vertice[2];
	NEWMESH::Pt pt0((v1->_coord.X + v2->_coord.X)/2, (v1->_coord.Y + v2->_coord.Y)/2, (v1->_coord.Z + v2->_coord.Z)/2); // adds new points at the center of existing faces
	NEWMESH::Pt pt2((v0->_coord.X + v1->_coord.X)/2, (v0->_coord.Y + v1->_coord.Y)/2, (v0->_coord.Z + v1->_coord.Z)/2);
	NEWMESH::Pt pt1((v0->_coord.X + v2->_coord.X)/2, (v0->_coord.Y + v2->_coord.Y)/2, (v0->_coord.Z + v2->_coord.Z)/2);
	
	
	boost::shared_ptr<Mpoint> p1,p2,p0;
	
	bool b0=true, b1=true, b2=true;
	count=0;
	int index=0;
	for (vector<boost::shared_ptr<Mpoint> >::const_iterator mi=added_points.begin(); mi!=added_points.end(); mi++)
	  {
	    index++;
	    NEWMESH::Pt current = (*mi)->_coord;
	    if (pt0==current) {b0=false; p0=*mi;}
	    if (pt1==current) {b1=false; p1=*mi;}
	    if (pt2==current) {b2=false; p2=*mi;}
	  }
	
	
	if (b0) {p0=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(pt0, nvertices() + count)); count++;};
	if (b1) {p1=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(pt1, nvertices() + count)); count++;};
	if (b2) {p2=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(pt2, nvertices() + count)); count++;}; 
	
	if (b0) {_points.push_back(p0);added_points.push_back(p0);}  // ER moved, neighbours must be added before triangles are pushed back
	if (b1) {_points.push_back(p1);added_points.push_back(p1);}
	if (b2) {_points.push_back(p2);added_points.push_back(p2);}
	
	Triangle t0(p2, p0, p1,tot_triangles); 
	OLD_tr_neighbours[t->get_no()].push_back(tot_triangles); 
	tot_triangles++; 

	Triangle t1(p1, v0, p2,tot_triangles);  
	OLD_tr_neighbours[t->get_no()].push_back(tot_triangles); 
	tot_triangles++;

	Triangle t2(p0, v2, p1,tot_triangles);  
	OLD_tr_neighbours[t->get_no()].push_back(tot_triangles);
	tot_triangles++;

	Triangle t3(p2, v1, p0,tot_triangles);  
	OLD_tr_neighbours[t->get_no()].push_back(tot_triangles);
	tot_triangles++;
	//Triangle t0(p1, p0, p2,tot_triangles);  tot_triangles++;
	//Triangle t1(p2, v0, p1,tot_triangles);  tot_triangles++;
	//Triangle t2(p1, v2, p0,tot_triangles);  tot_triangles++;
	//Triangle t3(p0, v1, p2,tot_triangles);  tot_triangles++;

	push_triangle(t0);
	push_triangle(t1);
	push_triangle(t2); 
	push_triangle(t3);
	

    }
    
    for (vector<boost::shared_ptr<Mpoint> >::const_iterator i=vbegin(); i!=vend(); i++)
	(*i)->normalize();
    
  }
  
  void newmesh::estimate_normals(){
    
    _normals.clear();
    for (int i=0;i<(int) _points.size();i++)
      _normals.push_back(local_normal(i));
    
  }

  Pt newmesh::local_normal(const int& pt)const{     
    Pt v(0,0,0);
    for(int i=0;i<_points[pt]->ntriangles();i++){      
      v+=_triangles[ _points[pt]->get_trID(i) ].normal();
    }
    v.normalize();
    return v;    
  }


  //----------ACCESS --------------------//

  const vector<vector<double> > newmesh::get_face_angles()const {

    vector<vector<double> > all_angles(0,vector<double>());

    for (int i=0;i<(int) _triangles.size();i++){
      all_angles.push_back(_triangles[i].get_angles());

    }
    return all_angles;
  }

  const int newmesh::get_triangleID_from_vertexIDs(const int & n0, const int & n1,const int & n2)const{
    boost::shared_ptr<Mpoint>  mp0,mp1,mp2;
    vector<int> IDs_0, IDs_10;
    int ID=-1;
    
    mp0=_points[n0]; mp1=_points[n1]; mp2=_points[n2];
    for (int i=0;i<mp0->ntriangles();i++){
      IDs_0.push_back(mp0->get_trID(i));
    }

    for (int i=0;i<mp1->ntriangles();i++){
      int tmpID=mp1->get_trID(i);

      for (int j=0;j<(int) IDs_0.size();j++){
	if(tmpID==IDs_0[j]){IDs_10.push_back(tmpID); }
	//if(IDs_10.size()==2) break;
      }
    }

    for (int i=0;i<mp2->ntriangles();i++){

      int tmpID=mp2->get_trID(i);
      for (int j=0;j<(int) IDs_10.size();j++){
	if(tmpID==IDs_10[j]){ ID=tmpID; break;}
      }
    }

    if(ID==-1){cout <<" newmesh::get_triangleID_from_3vertices error: cannot find triangle ID " << endl; exit(1); }
    return ID;
  }

  /// the following 4 functions are used during write
  vector< float > newmesh::getPointsAsVectors()const{   
    vector< float > ret;
    for(int i=0;i<nvertices();i++){
      ret.push_back( get_point(i).get_coord().X );
      ret.push_back( get_point(i).get_coord().Y );
      ret.push_back( get_point(i).get_coord().Z );
    }
    return ret;
  }
  vector< float > newmesh::getValuesAsVectors(const int & d)const{
    vector< float > ret;

    for(unsigned int i=0;i< _pvalues[d].size();i++){
      ret.push_back( get_pvalue(i,d) );
    }
    return ret;
  }

  vector< vector<unsigned int> > newmesh::getTrianglesAsVectors()const{
    
    vector< vector<unsigned int> > ret;
    for(int i=0;i<ntriangles();i++){
      vector<unsigned int> tmp(3);
      tmp[0]=(unsigned int)_triangles[i].get_vertice(0).get_no();
      tmp[1]=(unsigned int)_triangles[i].get_vertice(1).get_no();
      tmp[2]=(unsigned int)_triangles[i].get_vertice(2).get_no();
      ret.push_back(tmp);
    }
    return ret;
  }

  vector<int> newmesh::getTrianglesAsVector()const{
     vector<int> ret;
     for(int i=0;i<ntriangles();i++){
      ret.push_back(_triangles[i].get_vertice(0).get_no());
      ret.push_back(_triangles[i].get_vertice(1).get_no());
      ret.push_back(_triangles[i].get_vertice(2).get_no());
      //  cerr << i << " " << _triangles[i].get_vertice(0).get_no() << " " <<  _triangles[i].get_vertice(1).get_no() << " " << _triangles[i].get_vertice(2).get_no() << endl;
    }
    return ret;
  }


  Matrix  newmesh::get_pvalues() const{

    Matrix M(0,0);
    if(_pvalues.size()>=1){
      M.ReSize(_pvalues.size(),_pvalues[0].size());
      for(unsigned int dim=0;dim<_pvalues.size();dim++){
	if(_pvalues[dim].size()!=_pvalues[0].size()){throw NEWMESHException("get_pvalues: inconsistent dimensions");}
	for(unsigned int i=0;i<_pvalues[dim].size();i++){
	  M(dim+1,i+1)=_pvalues[dim][i];
	}
      }
    }
    return M;
  }



 vector<unsigned int> newmesh::cluster(const double threshold, const int field) {  //clusters all fields ( hopefully! )
  unsigned int currentCluster = 1;//give each cluster a different value, including across fields, starting at 1
  vector<unsigned int> clusterSizes(1,0); //size of cluster "0" is 0;
  int startingField(0), finalField(get_dimension()-1);
  if ( field > -1 )
    startingField=finalField=field;
  for (int field = startingField; field <= finalField; ++field) { //This needs to go across all fields
    vector<int> clusteredField(nvertices(), 0), marked(nvertices(), 0);
    for (int i = 0; i < nvertices(); ++i)
      if (get_pvalue(i,field) > threshold) 
	marked[i] = 1;     
    for (int i = 0; i < nvertices(); ++i) {
      if (marked[i]) {//node is above threshold
	vector<int> connectedNodes(1,i);
	marked[i] = 0;//unmark it when added to list to prevent multiples
	for (int index = 0; index < (int)connectedNodes.size(); ++index) { //NOTE: vector grows inside loop
	  int node = connectedNodes[index];//keep list around so we can put it into the output immediately if it is large enough
	  int numNeigh = get_total_neighbours(node);
	  for (int n = 0; n < numNeigh; ++n) {
	    const int32_t& neighbor =  get_neighbour(node, n);
	    if (marked[neighbor]) {		       
	      connectedNodes.push_back(neighbor);
	      marked[neighbor] = 0;
	    }
	  }
	} //searched nodes, now label

	if (currentCluster > INT_MAX) 
	  throw Exception("too many clusters, unable to mark them uniquely");
	int clusteredNodes((int)connectedNodes.size());
	for (int index = 0; index < clusteredNodes; ++index)
	  clusteredField[connectedNodes[index]] = currentCluster;
	if(clusteredNodes)
	  clusterSizes.push_back(clusteredNodes);
	++currentCluster;
      }
    }
    for (int i = 0; i < nvertices(); ++i)
      set_pvalue(i,clusteredField[i],field);
  } //Clustered field c
  return clusterSizes;
}

 /// find face that vb is in
  ReturnMatrix newmesh::closest_triangle(int ind, const NEWMESH::Pt& vb) const { //checked 
    
    NEWMAT::ColumnVector CV(3);
    CV = -1;
      
     
    for ( vector<int>::const_iterator i=_points[ind]->_trID.begin(); i!=_points[ind]->_trID.end(); i++)
      {
  	
	NEWMESH::Pt v1 =  _triangles[*i].get_vertex_coord(0),  v2 = _triangles[*i].get_vertex_coord(1), v3 = _triangles[*i].get_vertex_coord(2);
		
	NEWMESH::Pt PP;
	
	projectPoint(vb,v1,v2,v3,PP);
	
	if (PointInTriangle(PP,v1,v2,v3)) {
	  	  
	  CV(1) = _triangles[*i].get_vertex_no(0);
	  CV(2) = _triangles[*i].get_vertex_no(1);
	  CV(3) = _triangles[*i].get_vertex_no(2);
	  break;
	}
	
      }
    
    CV.Release();
    return CV;
  }


// find face that vb is in
  int newmesh::closest_triangle_ID(int ind, const NEWMESH::Pt& vb) const { //checked 
    
    int  c=-1;
      
     
    for ( vector<int>::const_iterator i=_points[ind]->_trID.begin(); i!=_points[ind]->_trID.end(); i++)
      {
  	
	NEWMESH::Pt v1 =  _triangles[*i].get_vertex_coord(0),  v2 = _triangles[*i].get_vertex_coord(1), v3 = _triangles[*i].get_vertex_coord(2);
		
	NEWMESH::Pt PP;
	
	projectPoint(vb,v1,v2,v3,PP);
	
	if (PointInTriangle(PP,v1,v2,v3)) {
	  c=*i;	  
	  break;
	}
	
      }
   
    return c;
  }
  /// The below functions are standard formula for estimating whether a point is actually within a triangle (note point is already on same plane due to projection in closest_triangle...)
  bool newmesh::return_closest_points(int closestvertex, NEWMESH::Pt &a1, NEWMESH::Pt &a2, NEWMESH::Pt &a3,const NEWMESH::Pt &p, int &n1, int &n2, int &n3) const {
    
    NEWMAT::ColumnVector tr(3); tr = -1;
    bool found;
    tr = closest_triangle(closestvertex, p);
    
    n1 = (int) tr(1);
    n2 = (int) tr(2);
    n3 = (int) tr(3);
    
    NEWMAT::ColumnVector tr0(3);
    
    tr0 = -1;
    
    //NOW: if tr == tr0 then we need to find another starting point
    
    if (tr == tr0) {
      
     
      found =0;
     

    } else{
     
      a1 = _points[n1]->get_coord();
      a2 = _points[n2]->get_coord();
      a3 = _points[n3]->get_coord();
      found=1;

    }
   
    return found;
  }
  
 

  ///// --------READ/WRITE--------------------///

  void newmesh::load(const string& filename, const bool loadSurfaceData, const bool appendFieldData){
    int type=meshFileType(filename);
    if(type==NEW_ASCII){
      load_ascii(filename,loadSurfaceData,appendFieldData);
    }
    else if(type==NEW_VTK){
      load_vtk(filename);
    }
    else if(type==NEW_GIFTI){
      load_gifti(filename,loadSurfaceData,appendFieldData);

    }
    else if(type==NEW_MATRIX || type==NEW_DPV){
      load_matrix(filename,type);

    }
    else{
    cerr<<"Mesh::load:error reading file: "<<filename<<"  ... Unknown format"<<endl;
    exit(1);
    }
    
  }

  void newmesh::load_gifti(const string& filename, const bool loadSurfaceData, const bool appendFieldData) {
    GIFTIwrapper reader;
    reader.readGIFTI(filename);

    if ( !appendFieldData ) {
      _pvalues.clear();
      _tvalues.clear();
    }

    if ( loadSurfaceData ) {

      _points.clear();
      _triangles.clear();


      vector<GIFTIfield> surfaceData=reader.returnSurfaceFields();
      global_defaultcoord=surfaceData[0].getCoordSystems();  // save out a default coord system

      for (int point = 0; point < surfaceData[0].getDim(0) ; point++ ) {
	vector<float> coords(surfaceData[0].fVector(point));
	boost::shared_ptr<Mpoint> m=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(coords[0],coords[1],coords[2], point));
	_points.push_back(m);
      }


      for ( int triangle = 0; triangle < surfaceData[1].getDim(0) ; triangle++ ) {
	vector<int> points(surfaceData[1].iVector(triangle));
	Triangle t(_points[points[0]],_points[points[1]],_points[points[2]], triangle);
	push_triangle(t);
      }      
    }
   
    vector<GIFTIfield> nonSurfaceData=reader.returnNonSurfaceFields();
    int point;
    ////////////// READ IN FEATURES  ///////////////////////////
    
    for (unsigned int dim=0; dim< nonSurfaceData.size();dim++){      
      /////// NEED TO CHECK FOR tvalues or pvalues ///////////////////
      if(_points.size()==0 || nonSurfaceData[dim].getDim(0)==(int)_points.size()){
	vector<float> tmp_pvalues;
	for (point = 0; point < nonSurfaceData[dim].getDim(0) ; point++ ) {
	  tmp_pvalues.push_back(nonSurfaceData[dim].fScalar(point));
	  
	}
	_pvalues.push_back(tmp_pvalues);
      }
      else if(nonSurfaceData[dim].getDim(0)==(int) _triangles.size()){ /// NOT TESTED  - makes assumption that gifti allows for storing value for faces in the same way as for vertices
	vector<float> tmp_tvalues;
	for (point = 0; point < nonSurfaceData[dim].getDim(0) ; point++ ) {
	  tmp_tvalues.push_back(nonSurfaceData[dim].fScalar(point));
	  
	}
	
      }else {throw NEWMESHException(" mismatch between data and surface dimensions");}
      
    }

    //////////////// SAVE METADATA ///////////////////
    global_metaData=reader.metaData;
    global_Attributes=reader.extraAttributes; //Any other attributes that the tag has
    global_GIFTIlabels=reader.GIFTIlabels;

  }

  void newmesh::load_vtk(const string& filename) {
    _points.clear();
    _triangles.clear();

    ifstream f(filename.c_str());
    if (f.is_open())
      {	
      //reading the header
      string header;
      getline(f, header);
      string::size_type pos = header.find("# vtk DataFile Version");
      if (pos == string::npos) {
	cerr<<"Mesh::load_vtk:error in the header"<<endl;exit(1);
      }
      getline(f,header);
      getline(f,header);
      getline(f,header);
      int NVertices, NFaces;
      f>>header>>NVertices>>header;	  
      //reading the points     
      for (int i=0; i<NVertices; i++)
	{
	  double x, y, z;
	  f>>x>>y>>z;
	  boost::shared_ptr<Mpoint> m=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(x,y,z,i));
	  _points.push_back(m);
	}
      f>>header>>NFaces>>header;
   
      //reading the triangles
      for (int i=0; i<NFaces; i++)
	{
	  int p0, p1, p2;
	  int j;
	  f>>j>>p0>>p1>>p2;
	  Triangle t(_points[p0],_points[p1],_points[p2],i);
	  push_triangle(t);
	 
	}
      f>>header>>header;
      f>>header>>header>>header;
      f>>header>>header;
      //reading the values

      for (int i=0; i<NVertices; i++)
	{
	  int val;
	  f>>val;	      
	  set_pvalue(i,val);
	}      	  
      f.close();
    }
  else {cout<<"Mesh::error opening file: "<<filename<<endl; exit(1);}
  }

  ////// if loadSurfaceData = true then pvalues are ignored (for compatability with gifti, where .surf and .func are loaded separately)
  void newmesh::load_ascii(const string& filename,const bool loadSurfaceData, const bool appendFieldData) { //load a freesurfer ascii mesh
     
    if(loadSurfaceData){_triangles.clear();
      _points.clear();}

    if ( !appendFieldData ) {
      _pvalues.clear();
      _tvalues.clear();
    }
    
    ifstream f(filename.c_str());
    if (f.is_open())
      {	
	//reading the header
	string header;
	getline(f, header);
	string::size_type pos = header.find("#!ascii");
	if (pos == string::npos) {
	  cerr<<"Mesh::load_ascii:error in the header"<<endl;exit(1);
	}
	


	 //reading the size of the mesh
	int NVertices, NFaces;
	f>>NVertices>>NFaces;
	
	if(!loadSurfaceData){
	  vector<float> tmp_pvalues(NVertices,0);
	  _pvalues.push_back(tmp_pvalues);
	}
	//reading the points
	for (int i=0; i<NVertices; i++)
	  {
	    double x, y, z;
	    float val;
	    f>>x>>y>>z>>val;
	    
	    if(loadSurfaceData){boost::shared_ptr<Mpoint> m=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(x,y,z,i));
	      _points.push_back(m);}
	    else _pvalues[_pvalues.size()-1][i]=val;
	  }      
	//reading the triangles
	
	if(!loadSurfaceData){vector<float> tmp_tvalues(NFaces,0);
	  _tvalues.push_back(tmp_tvalues);
	}

	for (int i=0; i<NFaces; i++)
	{
	  int p0, p1, p2;
	  float val;
	  f>>p0>>p1>>p2>>val;

	  if(loadSurfaceData){Triangle t(_points[p0],_points[p1],_points[p2],i); push_triangle(t);}
	  else _tvalues[_tvalues.size()-1][i]=val;
	}
	f.close();
      }
    else {cout<<"Mesh::load_ascii:error opening file: "<<filename<<endl; exit(1);}
  }


 
  void newmesh::load_ascii_file(const string& filename) { //load a freesurfer ascii mesh for pvalues only
    clear();
  
     ifstream f(filename.c_str());
     if (f.is_open())
       {	
	 //reading the header
	 string header;
	 getline(f, header);
	 string::size_type pos = header.find("#!ascii");
	 if (pos == string::npos) {
	   cerr<<"Mesh::load_ascii:error in the header"<<endl;exit(1);
	 }
	 
	 //reading the size of the mesh
	int NVertices, NFaces;
	f>>NVertices>>NFaces;
	
	vector<float> tmp_pvalues;
	
	//reading the points


	for (int i=0; i<NVertices; i++)
	  {
	    double x, y, z;
	    float val;
	    f>>x>>y>>z>>val;
	    boost::shared_ptr<Mpoint> m=boost::shared_ptr<Mpoint>(new NEWMESH::Mpoint(x,y,z,i));
	    _points.push_back(m);
	    tmp_pvalues.push_back(val);
	  }      
	_pvalues.push_back(tmp_pvalues);
	//reading the triangles
	
	vector<float> tmp_tvalues;
	_tvalues.push_back(tmp_tvalues);
	  

	for (int i=0; i<NFaces; i++)
	{
	  int p0, p1, p2;
	  float val;
	  f>>p0>>p1>>p2>>val;
	  Triangle t(_points[p0],_points[p1],_points[p2],i);
	  push_triangle(t);
	  tmp_tvalues.push_back(val);
	}
	_tvalues.push_back(tmp_tvalues);
	f.close();
  
       }else {cout<<"Mesh::load_ascii:error opening file: "<<filename<<endl; exit(1);}
  }
 
  void newmesh::load_matrix(const string& filename, const int& type) {  // for pvalues only - when data is held in a textfile
    Matrix tmp,tmp2;
    tmp=read_ascii_matrix(filename);
    if(type==NEW_DPV){
      tmp2=tmp;
      tmp.ReSize(tmp.Nrows(),1);
    //  if(tmp(1,1)!= 001)
    
     if(tmp2.Ncols()!=5)  {cout<<"Mesh::load_dpv:error opening file (wrong format) : "<<filename<<endl; exit(1);}
      for (int i=1;i<=tmp.Nrows();i++){
	if(tmp2(i,1)!=i-1) {cout<<i << " " <<  tmp2(i,1) << "Mesh::load_dpv:error opening file (wrong format 2) : "<<filename<<endl; exit(1);}
	tmp(i,1)=tmp2(i,5);	
}
    }
    if(tmp.Nrows()==0 && tmp.Ncols()==	0) {cout<<"Mesh::load_txt:error opening file (wrong format). Note matrix file must be delimited with spaces ?? "<<endl; exit(1);}
    if(tmp.Ncols()==5){ int ind=0;
      for (int i=1;i<=tmp.Nrows();i++) 
	if(tmp(i,1)==i-1) ind++;
			
      if(ind==tmp.Nrows()){cout<< "WARNING: this looks like a dpv file but is being read as a text file!! "<<endl;} 
    }
    set_pvalues(tmp);
  }

  void newmesh::save(const string& filename)const{

    int type=meshFileType(filename);
    switch(type)
      {
      case NEW_DPV:
	save_dpv(filename);
	break;
      case NEW_MATRIX:
	save_matrix(filename);
	break;
      case NEW_ASCII:
	save_ascii(filename);
	break;
      case NEW_VTK:
	save_vtk(filename);break;
      default:
	save_gifti(filename);break;
      }
  }
  
  

  
  void newmesh::save_gifti(const string& s)const{
    string filename(s);
    string last_3 = filename.substr(filename.size()-3, 3);
    GIFTIwrapper writer;

    string subtype;

    if( last_3 != "gii" ){
      if(last_3 != ".gz"){
	filename=filename+".gii";
      }
      subtype = filename.substr(filename.size()-9, 5);

    }else  subtype = filename.substr(filename.size()-9, 5);


    //////////////// SAVE METADATA ///////////////////

    writer.metaData=global_metaData;
    writer.extraAttributes=global_Attributes; //Any other attributes that the tag has
    writer.GIFTIlabels=global_GIFTIlabels;
    //writer.report();
    if(subtype == ".surf"){
   
      vector<GIFTIfield> surfaceData;
      vector<int> dims(2,3); //default as 3-vector for surface
      dims[0]=_points.size();

      vector<float> pointData(getPointsAsVectors()); 
      surfaceData.push_back(GIFTIfield(NIFTI_INTENT_POINTSET,NIFTI_TYPE_FLOAT32,2,dims.data(),pointData.data(),GIFTI_IND_ORD_ROW_MAJOR,global_defaultcoord));


      dims[0]=_triangles.size();
      vector<int> triangleData(getTrianglesAsVector());

      surfaceData.push_back(GIFTIfield(NIFTI_INTENT_TRIANGLE,NIFTI_TYPE_INT32,2,dims.data(),triangleData.data(),GIFTI_IND_ORD_ROW_MAJOR,global_defaultcoord));
      writer.allFields=surfaceData; 
      writer.writeGIFTI(filename,GIFTI_ENCODING_B64GZ);
    }
    else if(subtype == ".func" || subtype == "shape"){
      vector<int> scalardims(1);
      vector<GIFTIfield> fieldData;

      for (unsigned int dim=0;dim<_pvalues.size();dim++){
	scalardims[0]=_pvalues[dim].size();
	vector<float> valueData(getValuesAsVectors(dim));
	fieldData.push_back(GIFTIfield(NIFTI_INTENT_NONE,NIFTI_TYPE_FLOAT32,1,scalardims.data(),valueData.data(),GIFTI_IND_ORD_ROW_MAJOR));
      }
    ////// NEED TO IMPLEMENT OPTION FOR TVALUES  ///////////////////
      writer.allFields=fieldData; 
      writer.writeGIFTI(filename,GIFTI_ENCODING_B64GZ);
    }


  }
  


  void newmesh::save_vtk(const string& s)const{
    string filename(s);
    string last_3 = filename.substr(filename.size()-3, 3);
    if( last_3 != "vtk" ){
      filename=filename+".vtk";
    }
    
    ofstream flot(filename.c_str());
    if(flot.is_open()){
      flot<<"# vtk DataFile Version 3.0"<<endl
	  <<"surface file"<<endl
	  <<"ASCII"<<endl
	  <<"DATASET POLYDATA"<<endl
	  <<"POINTS ";
      flot<<_points.size()<<"  float"<<endl;
      
      for (unsigned int i =0; i<_points.size();i++)  { 
	//	flot.precision(6);
	flot<<_points[i]->get_coord().X<<" "
	    <<_points[i]->get_coord().Y<<" "
	    <<_points[i]->get_coord().Z<<endl;
#ifdef PPC64
	if ((n++ % 20) == 0) flot.flush();
#endif
      }
      flot<<"POLYGONS "<<_triangles.size()<<" "<<_triangles.size()*4<<endl;
      for ( unsigned int i=0; i<_triangles.size(); i++) 
	flot<<"3 "
	    <<_triangles[i].get_vertice(0).get_no()<<" "
	    <<_triangles[i].get_vertice(1).get_no()<<" "
	    <<_triangles[i].get_vertice(2).get_no()<<" "<<endl;
#ifdef PPC64
      if ((n++ % 20) == 0) flot.flush();
#endif
    }
    else{
      cerr<<"::save_vtk:error opening file "<<filename<<" for writing"<<endl;
      exit(1);
    }
  }
  
  void newmesh::save_ascii(const string& s)const{
    string filename(s);
    string last_3 = filename.substr(filename.size()-3, 3);
    if( last_3 != "asc" ){filename=filename+".asc";}
   
    ofstream f(filename.c_str());
    stringstream flot;
    if (f.is_open())
      {
	int ptcount(0), tricount(0);
	for(unsigned int i=0;i<_points.size();i++){
	  float val;
	  if(_pvalues.size()==0) val=0;
	  else val=_pvalues[0][i]; // value of first column

	  flot<<_points[i]->get_coord().X<<" "
	      <<_points[i]->get_coord().Y<<" "
	      <<_points[i]->get_coord().Z<<" "
	      <<val<<endl; 	  
	  ptcount++;
	}
	for(unsigned int i=0;i<_triangles.size();i++){
	  float val;
	  if(_tvalues.size()==0) val=0;
	  else val=_tvalues[0][i]; // value of first column

	  flot<<_triangles[i].get_vertice(0).get_no()<<" "
	      <<_triangles[i].get_vertice(1).get_no()<<" "
	      <<_triangles[i].get_vertice(2).get_no()<<" "<<val<<endl;
	  tricount++;
	}
	f<<"#!ascii from Mesh"<<endl;
	f<<ptcount<<" "<<tricount<<endl<<flot.str();
	f.close();
      }
    else cerr<<"Mesh::save_ascii:error opening file for writing: "<<s<<endl;

  }

  void newmesh::save_dpv(const string& s)const{
    string filename(s);
    string last_3 = filename.substr(filename.size()-3, 3);
    if( last_3 != "dpv" ){filename=filename+".dpv";}
   
    ofstream f(filename.c_str());
    if (f.is_open())
      {
	if(_pvalues.size()==0){ throw NEWMESHException("newmesh::save_dpv, cannot write out as dpv as there is no data");}
	else{
	  if(_points.size()!=_pvalues[0].size()){ throw NEWMESHException("newmesh::save_dpv, data and mesh dimensions do not agree");}

	  for(unsigned int i=0;i<_pvalues[0].size();i++){
	    
	    float val=_pvalues[0][i];  /// outputs value for first column only
	    if(i<100)
	      f << setfill('0') << setw(3) << i;
	    else 
	      f << i ;

	    if(_points.size()){
	      f<< " " << _points[i]->get_coord().X
	       << " " <<_points[i]->get_coord().Y
	       << " " <<_points[i]->get_coord().Z
	       << " " <<val<<endl; 	  
	    
	    }
	    else{
	      f<< " 0 0 0 " <<val<<endl; 

	    }
	 
	}
	  f.close();
      }
      }
    else cerr<<"Mesh::save_ascii:error opening file for writing: "<<s<<endl;


  }

    void newmesh::save_matrix(const string& s)const{
    string filename(s);
    string last_3 = filename.substr(filename.size()-3, 3);
    if( last_3 != "txt" ){filename=filename+".txt";}
   
    ofstream f(filename.c_str());
    if (f.is_open())
      {
	if(_pvalues.size()==0){ throw NEWMESHException("newmesh::save_matrix, cannot write out matrix as there is no data");}
	else{
	  if(_points.size()!=_pvalues[0].size()){ throw NEWMESHException("newmesh::save_matrix, data and mesh dimensions do not agree");}

	  for(unsigned int k=0;k<_pvalues.size();k++){
	    for(unsigned int i=0;i<_pvalues[k].size();i++){
	    
	      f<<  _pvalues[k][i] << " " ; 
	    
	    }
	    f << endl;
	  }
	  f.close();
	}
      }
    else cerr<<"Mesh::save_ascii:error opening file for writing: "<<s<<endl;


  }
  ////////////////////////ASSIGNMENT ///////////////////////////////////
 
  void newmesh::set_pvalue(const unsigned int& i,const float& val, const int dim){

    if(_pvalues.size()==0){ cout << "newmesh::set_pvalue, warning, pvalues have not been initialised. Generating pvalues field of length=number of vertices" << endl;  
      initialize_pvalues(dim+1,true);}
    else if((int) _pvalues.size()<dim+1){ cout << "newmesh::set_pvalue, warning, index exceeds known pvalues dimension. Appending pvalues" << endl;  
      initialize_pvalues(dim-_pvalues.size()+1,true); }
    
    if(i>=_pvalues[dim].size()) {cout << i << " dim " << dim << " " << _pvalues[dim].size() << endl;  throw NEWMESHException("newmesh::set_pvalue, index is incompatible with data dimensions");} 
    
    _pvalues[dim][i]=val;} // this will be problematic is i is > the size of the vector
 
  void newmesh::set_tvalue(const unsigned int& i,const float& val, const int dim){
    if(_tvalues.size()==0){ cout << "newmesh::set_tvalue, warning, tvalues have not been initialised. Generating tvalues field of length=number of faces" << endl;  
      initialize_tvalues(dim+1,true);}
    else if((int)_tvalues.size()<dim){ cout << "newmesh::set_tvalue, warning, index exceeds known tvalues dimension. Appending tvalues" << endl;  
      initialize_tvalues(dim-_tvalues.size(),true); }
      
    if(i>=_tvalues[dim].size()){throw NEWMESHException("newmesh::set_tvalue, index is incompatble with data dimensions");} 
    _tvalues[dim][i]=val;}

  void newmesh::set_pvalues(const vector<int>& ids,const float& val,const int dim){
    for(unsigned int i=0;i<ids.size();i++){_pvalues[dim][ids[i]]=val;}
  }
  
  
  void newmesh::set_pvalues(const Matrix& M, const bool appendFieldData){ 
     bool verticesAreColumns=false;

     if( _points.size()>0){ // check that if points have been supplied, that the data has the same length
       if(M.Ncols()==(int) _points.size()) verticesAreColumns=true;
       else if (M.Nrows()==(int) _points.size()) verticesAreColumns=false;
       else{cout << "(M.Ncols()" << M.Ncols() << " (M.Ncols() " <<M.Nrows() << endl;   throw NEWMESHException(" Cannot assign data to mesh. Dimensions do not match that of mesh");}
     }
     else if(_pvalues.size() > 0){
       if(M.Ncols()== (int) _pvalues[0].size()) verticesAreColumns=true;
       else if(M.Nrows()== (int) _pvalues[0].size())  verticesAreColumns=false;
       else{ throw NEWMESHException(" Cannot assign data to mesh. Dimensions do not match that of mesh");}
     }
     else{
       if(M.Ncols()>M.Nrows()) verticesAreColumns=true;
     }

     if(! appendFieldData)
     _pvalues.clear();

     if(verticesAreColumns){   
       for(int i=1;i<=M.Nrows();i++){
	vector<float> tmp_pvalues;
	for(int j=1;j<=M.Ncols();j++){
	  tmp_pvalues.push_back(M(i,j));
	}
	_pvalues.push_back(tmp_pvalues);
      }
     }else{
       for(int i=1;i<=M.Ncols();i++){
	 vector<float> tmp_pvalues;
	 for(int j=1;j<=M.Nrows();j++){
	   tmp_pvalues.push_back(M(j,i));
	 }
	 
	_pvalues.push_back(tmp_pvalues);

      }
    }
  }

  void newmesh::initialize_pvalues(const int dim,const bool appendFieldData)
  {    
    if ( !appendFieldData ) {
      _pvalues.clear();
    }
    
    vector<float> tmp(_points.size(),0.0);
    for (int i=0;i<dim;i++)
      _pvalues.push_back(tmp);
  } 

  void newmesh::initialize_tvalues(const int dim,const bool appendFieldData)
  {   
    if ( !appendFieldData ) {
      _tvalues.clear();
    }
    vector<float> tmp(_triangles.size(),0.0);
    for (int i=0;i<dim;i++)
    _tvalues.push_back(tmp);
  } 

  void newmesh::push_triangle(const NEWMESH::Triangle& t){

    vector<int> _nos;
    
    _triangles.push_back(t);
    
   
    for(int i=0;i<3;i++){
      _points[t.get_vertice(i).get_no()]->push_triangle(t.get_no());
      _nos.push_back(t.get_vertice(i).get_no());
        
    }
    
   
    if(!_points[_nos[0]]->is_neighbour(_nos[1])) _points[_nos[0]]->push_neighbour(_nos[1]); // added by emma to ensure vertex neighbourhoods are saved
    if(!_points[_nos[0]]->is_neighbour(_nos[2])) _points[_nos[0]]->push_neighbour(_nos[2]); 

    if(!_points[_nos[1]]->is_neighbour(_nos[0])) _points[_nos[1]]->push_neighbour(_nos[0]);
    if(!_points[_nos[1]]->is_neighbour(_nos[2])) _points[_nos[1]]->push_neighbour(_nos[2]); 

    if(!_points[_nos[2]]->is_neighbour(_nos[0])) _points[_nos[2]]->push_neighbour(_nos[0]); 
    if(!_points[_nos[2]]->is_neighbour(_nos[1])) _points[_nos[2]]->push_neighbour(_nos[1]);
  
   
}

  Pt newmesh::estimate_origin(){
    vector<Pt> p(4,Pt());
    srand (time(NULL));
    double m11,m12,m13,m14,m15;
   
    for (int i=1;i<=4;i++){  
      p[i-1]=_points[floor(_points.size()/i)-1]->get_coord(); /// might need a better selection criteria?
    }

    Matrix a(4,4);
    Pt c;
    double r;
    //equation of a sphere as determinant of its variables and 4 sampled points
    for (int i=0;i<4;i++){         //find minor 11
      a(i+1,1) = p[i].X;
      a(i+1,2) = p[i].Y;
      a(i+1,3) = p[i].Z;
      a(i+1,4) = 1;
    }
   
    m11= a.Determinant();
    

    for (int i=0;i<4;i++){         //find minor 12
      a(i+1,1) = p[i].X*p[i].X + p[i].Y*p[i].Y + p[i].Z*p[i].Z;
      a(i+1,2) = p[i].Y;
      a(i+1,3) = p[i].Z;
      a(i+1,4) = 1;
    }

    m12= a.Determinant();

    for (int i=0;i<4;i++){         //find minor 13
      a(i+1,1) = p[i].X*p[i].X + p[i].Y*p[i].Y + p[i].Z*p[i].Z;
      a(i+1,2) = p[i].X;
      a(i+1,3) = p[i].Z;
      a(i+1,4) = 1;
    }
    m13 = a.Determinant();

    for (int i=0;i<4;i++){   //find minor 14
      a(i+1,1) = p[i].X*p[i].X + p[i].Y*p[i].Y + p[i].Z*p[i].Z;
      a(i+1,2) = p[i].X;
      a(i+1,3) = p[i].Y;
      a(i+1,4) = 1;
    }
    m14  = a.Determinant();

    for (int i=0;i<4;i++){   // find minor 15
      a(i+1,1) = p[i].X*p[i].X + p[i].Y*p[i].Y + p[i].Z*p[i].Z;
      a(i+1,2) = p[i].X;
      a(i+1,3) = p[i].Y;
      a(i+1,4) = p[i].Y;
    }
    m15 =a.Determinant();

    if (m11 == 0)
      r = 0;
    else{
      c.X =  0.5 * (m12 / m11); //center of sphere
      c.Y = -0.5 * (m13 / m11);
      c.Z =  0.5 * (m14 / m11);
      r   = sqrt( c.X*c.X + c.Y*c.Y + c.Z*c.Z - (m15/m11) );
    }

    return c;
  }

// calculate on what side of a surface a step goes to
// a step here always starts at a vertex (vertind)
// the sign corresponds to the sign of the dot-product with the 
// normal to te closest tile
  int newmesh::step_sign(const int& vertind,const NEWMESH::Pt& step)const{
    int trid;
    float d=0,dmin=0;
    for(int i=0;i<_points[vertind]->ntriangles();i++){
      trid=_points[vertind]->get_trID(i);    
      d=(_triangles[trid].normal()|step);
      if(i==0||(fabs(d)<fabs(dmin))){dmin=d;}
  }
    return (dmin>0?1:-1);
  }



  void newmesh::print(const string& filename){
    ofstream fs(filename.c_str());
    fs<<_points.size()<<" vertices"<<endl;
    for(unsigned int i=0;i<_points.size();i++){
      fs<<_points[i]->get_coord().X<<" "
	  <<_points[i]->get_coord().Y<<" "
	  <<_points[i]->get_coord().Z<<endl;
    }
    fs<<_triangles.size()<<" triangles"<<endl;
    for(unsigned int i=0;i<_triangles.size();i++){
      fs<<_triangles[i].get_vertice(0).get_no()<<" "
	  <<_triangles[i].get_vertice(1).get_no()<<" "
	  <<_triangles[i].get_vertice(2).get_no()<<endl;
    }
    fs<<_pvalues[0].size()<<" scalars"<<endl;
    for(unsigned int i=0;i<_pvalues[0].size();i++){
      if(_pvalues[0][i]!=0){fs<<_pvalues[0][i]<<endl;}
    }
    fs.close();
  }

  bool SameSide(const NEWMESH::Pt& p1, const NEWMESH::Pt& p2, const NEWMESH::Pt& a, const NEWMESH::Pt& b) {//checked 
    
    NEWMESH::Pt cp1 = (b-a)*(p1-a);
    NEWMESH::Pt cp2 = (b-a)*(p2-a);
 
    if ((cp1|cp2) >= -1E-6) {return true;} else {return false;}
    
    
  }

  bool PointInTriangle(const NEWMESH::Pt& p, const NEWMESH::Pt& a, const NEWMESH::Pt& b, const NEWMESH::Pt& c) {

    if (SameSide(p,a,b,c) && SameSide(p,b,c,a) && SameSide(p,c,a,b)) return true;  else return false;
  
  }
  
  
 
 // project onto triangle
  NEWMESH::Pt projectPoint(const NEWMESH::Pt &vb, const NEWMESH::Pt &v1, const NEWMESH::Pt &v2, const NEWMESH::Pt &v3,NEWMESH::Pt &PP) 
  {    
    NEWMESH::Pt s1 = v3 - v1; s1.normalize();
    NEWMESH::Pt s2 = v2 - v1; s2.normalize();
   
    NEWMESH::Pt s3 = s1*s2;/// s3= normal
    s3.normalize(); 
    
    double si = (s3|v1)/(s3|vb); // formula for line plane intersection s3.(v1-sivb)=0 i.e. (v1-sivb) should be perpendicular to plane normal s3
       
     PP = vb*si;  // PP therefore where vb intersects plane
     return s3;
     
  }


// project onto triangle
  void projectPoint(const NEWMESH::Pt &v1, const NEWMESH::Pt &v2, const NEWMESH::Pt &v3,NEWMESH::Pt &PP) 
  {    
    NEWMESH::Pt s1 = v3 - v1; s1.normalize();
    NEWMESH::Pt s2 = v2 - v1; s2.normalize();
   
    NEWMESH::Pt s3 = s1*s2;/// s3= normal
    s3.normalize(); 
    
    double si = (s3|v1)/(s3|PP); // formula for line plane intersection s3.(v1-sivb)=0 i.e. (v1-sivb) should be perpendicular to plane normal s3
    
     PP = PP*si;  // PP therefore where vb intersects plane
    
     
  }


  void projectVector(const NEWMESH::Pt &VEC, const NEWMESH::Pt &v1, const NEWMESH::Pt &v2, const NEWMESH::Pt &v3,NEWMESH::Pt &Pvec) 
  {    
    NEWMESH::Pt s1 = v3 - v1; s1.normalize();
    NEWMESH::Pt s2 = v2 - v1; s2.normalize();
   
    NEWMESH::Pt s3 = s1*s2;/// s3= normal
    s3.normalize(); 

    Pvec=s3*(VEC*s3);
   
     
  }
}

