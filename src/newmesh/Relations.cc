/*  Relations.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */

#include "Relations.h"

namespace NEWMESH {

  /// -----------------GRID FUNCTIONS ---------------------//
  void Grid::Initialize(const newmesh & _mesh,const int& _Nth,const int& _Nph){
    Initialize(_mesh.get_points(),_Nth,_Nph);
  }

  void Grid::Initialize(const vector<boost::shared_ptr<Mpoint> > & mpoints,const int& _Nth,const int& _Nph){
    Nth=_Nth;
    Nph=_Nph;
    m=mpoints;
    mesh2cell.reserve(m.size());
    rth=float(Nth)/PI;rph=float(Nph)/TWOPI;
    cell2mesh.resize(Nth*Nph);
    create_luts();
    
    
  }

  void Grid::create_luts(){    
    int k=0;
    vals.ReSize(1,m.size());
    for (vector<boost::shared_ptr<NEWMESH::Mpoint> >::const_iterator i=m.begin(); i!=m.end(); i++){
      Pt p=(*i)->get_coord();
      int cell=get_cell(p.X,p.Y,p.Z); // find cell for every vertex of m (target mesh)
      mesh2cell.push_back(cell);
      if(cell >= (int) cell2mesh.size() || k>(int) m.size()){ cout << " cell >= cell2mesh.size() " << endl; }
      cell2mesh[cell].push_back(k);k++;
      vals(1,k)=cell;
    }  
   
    
  }
  
  int Grid::get_cell(const float& x,const float& y,const float& z)const{  // divide sphere into grid spaced by angle rth and rph
    float th,ph;ColumnVector v(3);
    int cell;
    v<<x<<y<<z;
    cart2sph(v,th,ph);
   
    //// need typecasting to float to avoid errors when th or ph is max value
    if(std::floor(float(rth*th))==Nth){ 
      if(std::floor(float(rph*(ph+PI)))==Nph) 
	cell=int((std::floor(float(rth*th))-1)*Nph+std::floor(float(rph*(ph+PI))))-1;
      else
	cell=int((std::floor(float(rth*th))-1)*Nph+std::floor(float(rph*(ph+PI))));
      if( cell >= (int) cell2mesh.size()) cout << " cell >= cell2mesh.size() "  << endl; } // enforce round down in event th==PI;
    else if (std::floor(float(rph*(ph+PI)))==Nph){cell=int(std::floor(float(rth*th))*Nph+int(std::floor(float(rph*(ph+PI)))-1)); if(cell >= (int) cell2mesh.size()) cout << " cell >= cell2mesh.size() "  << endl; }  
    else{ 
      if((ph+PI)<0){ /// rounding of PI error
	cell=int(std::floor(float(rth*th))*Nph); 
	if(cell >= (int) cell2mesh.size()) cout << " cell >= cell2mesh.size() "  << endl; 
      }
      else{
	cell=int(std::floor(float(rth*th))*Nph+std::floor(float(rph*(ph+PI))));  
	if(cell >= (int) cell2mesh.size()) cout << " cell >= cell2mesh.size() "  << endl; 
      } 
    }
    
    return cell;  
  }
  
  vector<int> Grid::get_cells_in_range(const int &nth2, int nph2, const int &i, const int &j, bool debug)const {

    int ia,ib,nph;
    vector<int> ret;
    bool theta_out_of_range=false;
    

    for(int a=-nth2;a<=nth2;a++){
    
      ia=i+a;


      if(ia<0){ ia=-ia-1;  theta_out_of_range=true; std::ceil(nph=float(Nph/2)); } 
      else if(ia>=Nth){ ia=2*Nth-ia-1;  theta_out_of_range=true; std::ceil(nph=float(Nph/2));} 
      else nph=nph2;


      for(int b=-nph;b<=nph;b++){
	ib=j+b;

	if(theta_out_of_range==true){  ib=Nph-ib; } 

	if(ib<0){ib+=Nph;}
	else if(ib>=Nph){ib-=Nph;}
		
       	if( ia*Nph+ib  <0 ) { cout << ia*Nph+ib << " ia " << ia << "  " << Nph  << " " << Nth << " ib " <<  ib << " rth " << rth << " i "  << i << " nth2 " << nth2 << " rph " << rph << " nph2 " << nph2 <<" PI " << PI << " debug " << debug <<   endl; exit(1);
	}
	ret.push_back(ia*Nph+ib );
	
	
      }
      theta_out_of_range=false;
    }
    return ret;
  }

  vector<int> Grid::get_cell_group(int c,float ang)const{  /// find all cells around c that are within angle _ang
    int i=std::floor(float(c)/Nph);int j=c-i*Nph; 

   
    int nth2=std::ceil(float(ang*rth));
    double th=i/rth;
    int nph2;


    if(th==0 || abs(th-PI) < 1e-3) std::ceil(nph2=float(Nph/2));
    else nph2=std::ceil(float((1/abs(sin(th)))*ang*rph));

    if(abs(nph2) >= 0.5*Nph )  std::ceil(nph2=float(Nph/2));
  
    vector<int> ret2=get_cells_in_range(nth2,nph2,i,j);
    
    return ret2;
  }
  
  vector<int> Grid::get_cell_group_exact_range(int c,float ang)const{  /// find all cells around c that are within angle _ang
    int i=std::floor(float(c)/Nph);int j=c-i*Nph; 

    double th=i/rth;
    int nth2=MISCMATHS::round(float(ang*rth));
    int nph2;

   
    if(th==0 || abs(th-PI) < 1e-3) MISCMATHS::round(nph2=float(Nph/2));
    else nph2=MISCMATHS::round(float((1/abs(sin(th)))*ang*rph));

    if(abs(nph2) >= 0.5*Nph )  MISCMATHS::round(nph2=float(Nph/2));
    vector<int> ret2=get_cells_in_range(nth2,nph2,i,j,true);


    return ret2;
  }
  
  vector< pair<float,int> > Grid::get_points(const float& x,const float& y,const float& z,const float& ang)const{
    int c=get_cell(x,y,z);

    vector<int> cg=get_cell_group(c,ang);
  
    ColumnVector selected(m.size()); selected=0.0;
    vector< pair<float,int> > pts; pair<float,int> tmp;
    Pt p(x,y,z);
    Pt p2;
    p.normalize();
    // start= clock();

    for(unsigned int i=0;i<cg.size();i++){
      for (unsigned int j=0;j<cell2mesh[cg[i]].size();j++){
	tmp.second=cell2mesh[cg[i]][j];
	p2=m[tmp.second]->get_coord(); 
	p2.normalize();
	if((p2|p) >= cos(ang)){  /// check proximity using angular separation, if within ang keep 
	  if(selected(tmp.second+1)==0){
	    tmp.first=(p2-p).norm();
	    pts.push_back(tmp);
	    selected(tmp.second+1)=1;
	  }
	 
	}
      }
    }
    sort(pts.begin(),pts.end());

    return pts;
  }

  
  ///////    FNS FOR EXTRACTING NEIGHBOURHOODS //////////////////////
  void  RELATIONS::Initialize(NEWMESH::newmesh &source,const NEWMESH::newmesh & target, const double& A)
  { 
    mat.clear();
    int total=0;
    check_scale(source,target);

    if(target.nvertices()/900<10){
      _grid.Initialize(target,floor(sqrt(target.nvertices()/10)),floor(sqrt(target.nvertices()/10)));

    }
    else
    _grid.Initialize(target,30,30);

    vector<int> tmp; 
    _ang=A;

    if(_ASFACES){total=source.ntriangles(); }
    else{ total=source.nvertices(); }

    for (int i=0;i<total;i++){
      mat.push_back(tmp);
    }
  }
 
  void  RELATIONS::Initialize(newmesh & source,const vector<boost::shared_ptr<Mpoint> > &  targetpoints, const double& A)
  { 
    mat.clear();
    int total=0;

    if(targetpoints.size()/900<10){
      _grid.Initialize(targetpoints,floor(sqrt(targetpoints.size()/10)),floor(sqrt(targetpoints.size()/10)));

    }
    else
    _grid.Initialize(targetpoints,30,30);

    vector<int> tmp; 
    _ang=A;

    if(_ASFACES){total=source.ntriangles(); }
    else{ total=source.nvertices(); }

    for (int i=0;i<total;i++){
      mat.push_back(tmp);
    }
  }

  vector<int> RELATIONS::return_cell_group(const Pt &p,const double & range){
    int c=_grid.get_cell(p.X, p.Y, p.Z);
    vector<int> cg=_grid.get_cell_group_exact_range(c,range);
    return cg;
  }

  void RELATIONS::update_RELATIONS(const NEWMESH::newmesh &source){
    //
    int total=0;
    vector<int> tmp; 
 
    if(_ASFACES){ throw  NEWMESHException(" newmesh::RELATIONS::update error. Use update_RELATIONSTRI for face neighbours");}
    mat.clear();  /// resize mat
    
    total=source.nvertices();

    for (int i=0;i<total;i++)
      mat.push_back(tmp);
  
    if(_ASFACES){total=source.ntriangles(); }


    for (int i = 1; i <= total; i++){
      update_RELATIONS_for_ind(i,source);
    }	
    
        
  }

 

  void RELATIONS::update_RELATIONSTRI(const NEWMESH::newmesh &source, const NEWMESH::newmesh &target){
    //
    int total=0;
    vector<int> tmp; 
 
    mat.clear();  /// resize mat
    total=source.ntriangles();
    

    for (int i=0;i<total;i++)
      mat.push_back(tmp);
  

    for (int i = 1; i <= total; i++){
      update_RELATIONS_for_tri(i,source,target);	
     
     
    }
   
   
  }

  void RELATIONS::update_RELATIONS_for_ind(int index,const NEWMESH::newmesh &source,double newang){
    Pt cr;

    if(newang>PI) { cout << " newang " << newang << endl;throw  NEWMESHException(" newmesh::update_RELATIONS_for_ind::() angle error");}
    cr = source.get_coord(index-1);

    if(cr.norm()>0){
      if(newang>0) _ang=newang;
      mat[index-1].clear(); 
      vector< pair<float,int> >  V=update_RELATIONS_for_ind(cr);

      for(unsigned int j=0;j<V.size();j++){
	mat[index-1].push_back(V[j].second+1);  
       	if(_EstAdj) ADJ->Set(V[j].second+1,index,1);
      }
    }
  }

  /// (predominantly FOR STRAINS BASED REGUm_SOURCELARISER) I would prefer not to have to use this as the dot product checks aren't stable
  void RELATIONS::update_RELATIONS_for_tri(int index,const NEWMESH::newmesh &source,const NEWMESH::newmesh &target,double newang){
    Pt cr,s0,s1,s2;
    Pt d0,d1,d2;
    bool found=false;
    ColumnVector NN(source.nvertices()); found=0;
    if(newang>PI) { cout << " newang " << newang << endl;throw  NEWMESHException(" newmesh::update_RELATIONS_for_ind::() angle error");}
      Pt v0,v1,v2;
      v0=source.get_triangle_vertex(index-1,0);
      v1=source.get_triangle_vertex(index-1,1);
      v2=source.get_triangle_vertex(index-1,2);
      cr=(v0+v1+v2)/3;
     
      s1=v1-v0; s1.normalize();
      s2=v2-v0;  s2.normalize();
      s0=v2-v1; s0.normalize();

      
      if(cr.norm()>0){
      if(newang>0) _ang=newang;
      mat[index-1].clear(); 

      vector< pair<float,int> >  V=update_RELATIONS_for_ind(cr);

      for(unsigned int j=0;j<V.size();j++){
	found=false;
	Pt p=target.get_coord(V[j].second);
	  d0=p-v0; d1=p-v1;d2=p-v2;
	  d0.normalize();	  d1.normalize(); d2.normalize();
	 
	  /// check whether lies on face edge or within triangle

 	  if((d2|s2) < -0.98 && (d0|s2) > 0.98){ found=true;} //  0.98 based on comparison of ico1 faces to ico2 faces
	  else if ((d0|s1) > 0.98 && (d1|s1) < -0.98 ){ found=true;}
	  else if ((d2|s0) < -0.98 && (d1|s0) > 0.98) {  found=true;}
	  else if(PointInTriangle(p,v0, v1,v2)) found =true;

	  if(found==true){
	    mat[index-1].push_back(V[j].second+1);  
	    if(_EstAdj) ADJ->Set(V[j].second+1,index,1);
	  }
	  
	}

    }
  
  }
  
  
  vector< pair<float,int> >   RELATIONS::update_RELATIONS_for_ind(const Pt &cr, double  newang) const{
    vector< pair<float,int> >  V; 
    if(newang==0) newang=_ang; 
    if(cr.norm()>0){
      if(_EstAdj && (!ADJ.get())){ throw  NEWMESHException(" RELATIONS::Adjacency matrix has not been initialised! ");}
      
      while(V.size()<7){ // 6 being max number of neighbours of any vertex
	V = _grid.get_points(cr.X, cr.Y, cr.Z, newang); /// I would like a faster function to replace this one 
	if(V.size() < 7){
	  newang=2*newang;
	}
      }
     
      
    }
    return V;
  }
 
  
  int RELATIONS::get_closest_point_index(Pt &cr){
    vector< pair<float,int> >  V; 
    while (V.size()==0){
      V = _grid.get_points(cr.X, cr.Y, cr.Z, _ang);
      if(V.size()==0) _ang=_ang*2;
    }
    return V[0].second+1;
  }
 
  void RELATIONS::find_next_closest_neighbour(int index, int & t_ind, NEWMESH::Pt &point,const NEWMESH::newmesh & REF){
    ////// use neighbours of closest point to identify the next closest (assumes no folding)
    NEWMESH::Pt v0,v1,v2;
    double min_dist;
    ColumnVector dist(3), num(3);
    
    int NN;
    NN=t_ind;
    
    v0= REF.get_coord(t_ind-1);
    min_dist=(point-v0).norm();
    
    for (vector<int>::const_iterator j =REF.tIDbegin(t_ind-1); j!=REF.tIDend(t_ind-1); j++){
      
      NEWMESH::Triangle tri=REF.get_triangle(*j);
      v0 = tri.get_vertex_coord(0);  v1 = tri.get_vertex_coord(1); v2 = tri.get_vertex_coord(2);
      
      
      dist(1)=(point-v0).norm(); dist(2)=(point-v1).norm(); dist(3)=(point-v2).norm();
      num(1)=tri.get_vertex_no(0); num(2)=tri.get_vertex_no(1); num(3)=tri.get_vertex_no(2);  
      
      for (int i=1;i<=3;i++){
	if(dist(i)<min_dist){
	  min_dist=dist(i);
	  NN=(int) num(i)+1;
	}
      }
    }
    
    if(NN==t_ind){cout << index << "closest neighbour was not within nearest triangles, search whole mesh space "  << endl; exit(0);}
    else mat[index-1][0]=NN;
    
    t_ind=NN;
  }
  
  void RELATIONS::update_w_querypoints_for_ind(const int &index,const vector<int> &neighbours){
    
    if(_EstAdj && (!ADJ.get())){ throw  NEWMESHException(" RELATIONS::Adjacency matrix has not been initialised! ");}
 
    ColumnVector newadj; //(ADJ->Nrows()); 
    int t_ind= mat[index-1][0];
    
    if(_EstAdj) {ADJ->Set(t_ind,index,1);newadj.ReSize(ADJ->Nrows());}
    mat[index-1].erase(mat[index-1].begin()+1,mat[index-1].end());
    
    for(unsigned int j=0;j < neighbours.size();j++){
      if(neighbours[j]!= t_ind-1){
	mat[index-1].push_back(neighbours[j]+1);
	  
	if(_EstAdj) ADJ->Set(neighbours[j]+1,index,1);
      }
    }
    
  }
  
 RELATIONS  RELATIONS::invert_relations(const NEWMESH::newmesh &newsource, const NEWMESH::newmesh &newtarget,double ang){
 //////// INVERSE CP neighbourhoods /////////////////
    if(ang==0) ang=_ang;
    newmesh TMP=newsource;

    RELATIONS NEWrel(TMP,newtarget,ang);

    Pt CP,tp;
    double dist;

    if((int) mat.size()!=newtarget.nvertices()) { cout << mat.size() << " " << newtarget.nvertices() << " " << newsource.nvertices() << endl; throw  NEWMESHException(" newmesh::RELATIONS::invert Error. Either the original relations matrix has not been initialised or the supplied target matrix has a diffferent number of vertices to the original");}
    for (int k=1;k<=newtarget.nvertices();k++){
      for(int i=0;i<(int)mat[k-1].size();i++){
	CP=newsource.get_coord(mat[k-1][i]-1);
	dist=asin((CP-newtarget.get_coord(k-1)).norm()/RAD);
	tp=newtarget.get_coord(k-1);
	CP.normalize(); tp.normalize();
	if((CP|tp)>cos(ang)){
	  NEWrel.Add(mat[k-1][i],k);

	}
      }
    }

   
    return NEWrel;
  }

  RELATIONS  RELATIONS::invert_relationsTR(const NEWMESH::newmesh &newsource, const NEWMESH::newmesh &newtarget,double ang){
    //////// INVERSE get neighbours for each FACE of the newsource /////////////////
    newmesh TMP=newsource;
    if(ang==0) ang=_ang;
    RELATIONS NEWrel(TMP,newtarget,ang,true);

    Pt v0,v1,v2,tmp;
    int n0,n1,n2,newID;

    TMP=newtarget;
    for (int i=0;i<TMP.nvertices();i++){
      TMP.set_pvalue(i,0);
    }
    
    //newtarget.save("newtarget.surf");
    //newsource.save("newsource.surf");

    if((int) mat.size()!=newtarget.nvertices()) {cout << mat.size() << " " << newtarget.nvertices() << " " << newsource.nvertices() << endl; throw  NEWMESHException(" newmesh::RELATIONS::invert Error. Either the original relations matrix has not been initialised or the supplied target matrix has a diffferent number of vertices to the original");}
    for (int k=1;k<=newtarget.nvertices();k++){
      tmp=newtarget.get_coord(k-1);
      
      for(int i=0;i<(int) mat[k-1].size();i++){

	if(newsource.return_closest_points(mat[k-1][i]-1,v0,v1,v2,tmp,n0,n1,n2)){
	  newID=newsource.get_triangleID_from_vertexIDs(n0,n1,n2);
	  projectPoint(v0,v1,v2,tmp);
	  NEWrel.Add(newID+1,k);
	  break;
	}
      }
    }

    

    
    return NEWrel;
  }

  vector<int> RELATIONS::Col(int i)const {
    
   vector<int> col; 
    for(unsigned int j=0;j< mat[i-1].size();j++)
	col.push_back(mat[i-1][j]);
       
    return col;
  }

  void RELATIONS::load(const string &fname){
    
    ifstream f(fname.c_str());
    mat.clear();
    
    if (f.is_open())
      {
	string header;
	getline(f, header);
	string::size_type pos = header.find("#!RELATIONS");
	if (pos == string::npos) { throw  NEWMESHException(" MESHREG::RELATIONS::error in the header. ");	}
	
	//reading the size of the mesh
	int Ncols,Nrows,val;
	string row,cell;
	
	vector<int> tmp;
	f>>Ncols>>Nrows>>_ang>>_ASFACES;
	getline(f, row);

	//ADJ= boost::shared_ptr<MISCMATHS::SpMat<int> >(new MISCMATHS::SpMat<int>(Nrows,Ncols)); 
	//	ColumnVector adjcol(Nrows); 
	//reading the points
	for (int i=0; i<Ncols; i++)
	  {
	    mat.push_back(tmp);
	    getline(f, row);
	    stringstream lstr(row);
	    //  adjcol=0;
	    while(getline(lstr,cell,',')){
	      val=atoi(cell.c_str());
	      if(val!=0){
		mat[i].push_back(val);
		//	adjcol(val)=1;
	      }
	    }
	    // ADJ->SetColumn(i+1,adjcol);
	    
	  }
	
	
	
	//reading the triangles
      }
    
  }
  
  void RELATIONS::Save(const string & s)const {
    ofstream out;
    out.open(s.c_str());
    
    out << "#!RELATIONS\n" ;
    out << mat.size() << " "  << _ang <<  " " << _ASFACES <<"\n" ;

    for(unsigned int i=0;i<mat.size();i++){
      for(unsigned int j=0;j<mat[i].size();j++){
	out << mat[i][j] << "," ;
      }
      out << "\n" ;
    }
    out.close();
  }

    ////////////////////// HELPER FUNCTIONS //////////////////////////
  
  void check_scale(NEWMESH::newmesh& in,const  NEWMESH::newmesh& ref){//checked 
    
    Pt cr2,cr3;
    Pt cr = in.get_coord(0);
    
    cr2 = in.get_coord(1);
    cr3 = ref.get_coord(1);
          
    if(abs(cr.norm()-cr2.norm())>1e-3|| abs(cr.norm()-cr3.norm())>1e-3 || abs(cr2.norm()-cr3.norm())>1e-3){
      true_rescale(in,cr3.norm());
    }
  }
  
 
  bool check_scale(const NEWMESH::newmesh& in,const double &scale){//checked 
    
    Pt cr2;
    Pt cr = in.get_coord(0);
    bool actual_scale=true;
    int num;

    if(in.nvertices()> 1000) num=1000;
    else num=in.nvertices();

    for (int i=0;i<num;i++){ // don't need to check every vertex
      cr2 = in.get_coord(0);
      if(cr2.norm()!=scale)
	actual_scale=false;
    }
    return actual_scale;
  }
  /// rescales sphere vertices to have equal radii 
  
  void true_rescale(NEWMESH::newmesh& M, const double &rad){//checked 
    
    for (int i = 0; i < M.nvertices(); i++){
      
      Pt cr = M.get_coord(i);
      cr.normalize();
      cr = cr*rad;
      
      M.set_coord(i,cr);   
    }
    
  }

  Matrix recentre(newmesh &sphere){

    Pt mean;
    Matrix TRANSLATE(4,4) ; TRANSLATE=0; 
    
    mean=sphere.estimate_origin();

    if(mean.norm()>1e-2){

     
      ColumnVector P_in(4);
      TRANSLATE(1,1)=1;  TRANSLATE(2,2)=1;  TRANSLATE(3,3)=1;  TRANSLATE(4,4)=1; 

      TRANSLATE(1,4)=-mean.X;  TRANSLATE(2,4)=-mean.Y;  TRANSLATE(3,4)=-mean.Z;  
     
      for ( vector<boost::shared_ptr<NEWMESH::Mpoint> >::const_iterator i= sphere.vbegin(); i!=sphere.vend(); i++){
       NEWMESH::Pt p = (*i)->get_coord();
       NEWMESH::Pt p2;
       if(p.norm()){
	 P_in(1) = p.X; P_in(2) = p.Y; P_in(3) = p.Z; P_in(4) = 1;  
	 P_in = TRANSLATE * P_in;
	 p2.X=P_in(1); p2.Y=P_in(2);  p2.Z=P_in(3);
	 (*i)->set_coord(p2);
       }
      
      }
    }	 
    return TRANSLATE;
  }

 
}
  
