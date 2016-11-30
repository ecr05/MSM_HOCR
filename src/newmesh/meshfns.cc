/*  meshfns.cc

    Emma Robinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "meshfns.h"

namespace NEWMESH {

  ////TANGENT PLANE ESTIMATIONS

// http://stackoverflow.com/questions/5255806/how-to-calculate-tangent-and-binormal
  Tangs Tangent::calculate(int ind,const NEWMESH::newmesh & SPH_in){
    Tangs T;
    NEWMESH::Pt a = SPH_in.local_normal(ind); 

    Pt tmp=SPH_in.get_coord(ind);
    if((a|tmp)<0) a=a*-1;
  
  /* pick one random vector */
    NEWMESH::Pt b(1.0f,0.0f,0.0f);

    NEWMESH::Pt c = a*b;
    double len = SQR(c.X) + SQR(c.Y) + SQR(c.Z);

    if (len == 0.0){
      /* the vector b was parallel to a */
      b.X=0.0f;
      b.Y=1.0f;
      b.Z=0.0f;
      c = a*b;
      len = SQR(c.X) + SQR(c.Y) + SQR(c.Z);
    }
  /* normalize */
    len=std::sqrt(len);
    
    if (len == 0.0)
      {
	cout<<"Tangs Tangent:: first tangent vector of length zero at vertex: "<<ind<<endl;
	len = 1 ;
      }
    
    T.e1.X = c.X / len ;
    T.e1.Y = c.Y / len ;
    T.e1.Z = c.Z / len ;
    
    b = a*c;
    /* normalize */
    len=std::sqrt(SQR(b.X)+SQR(b.Y)+SQR(b.Z));
    
    if (len == 0)
      {
	cout<<"Tangs Tangent::second tangent vector of length zero at vertex: "<<ind<<endl;
	len = 1 ;
      }
    
    T.e2.X = b.X / len ;
    T.e2.Y = b.Y / len ;
    T.e2.Z = b.Z / len ;
    
    
    return T;
  }

  Tangs Tangent::calculate_tri(const Pt &a){
    Tangs T;
  
  /* pick one random vector */
    NEWMESH::Pt b(1.0f,0.0f,0.0f);

    NEWMESH::Pt c = a*b;
    double len = SQR(c.X) + SQR(c.Y) + SQR(c.Z);
  
    if (len == 0.0){
      /* the vector b was parallel to a */
      b.X=0.0f;
      b.Y=1.0f;
      b.Z=0.0f;
      c = a*b;
      len = SQR(c.X) + SQR(c.Y) + SQR(c.Z);
    }
   
  /* normalize */
    len=std::sqrt(len);
    
    if (len == 0.0)
      {
	cout<<"Tangs Tangent:: first tangent vector of length zero at vertex: "<<endl;
	len = 1 ;
      }
    
    T.e1.X = c.X / len ;
    T.e1.Y = c.Y / len ;
    T.e1.Z = c.Z / len ;
    
    b = a*c;
    /* normalize */
    len=std::sqrt(SQR(b.X)+SQR(b.Y)+SQR(b.Z));
    
    if (len == 0)
      {
	cout<<"Tangs Tangent::second tangent vector of length zero at vertex: "<<endl;
	len = 1 ;
      }
    
    T.e2.X = b.X / len ;
    T.e2.Y = b.Y / len ;
    T.e2.Z = b.Z / len ;
    
    
    return T;
  }
  // http://stackoverflow.com/questions/5255806/how-to-calculate-tangent-and-binormal
  Tangs Tangent::calculate_tri(int ind,const NEWMESH::newmesh & SPH_in){
    Tangs T;
    NEWMESH::Pt a = SPH_in.get_triangle_normal(ind); 

    
  /* pick one random vector */
    NEWMESH::Pt b(1.0f,0.0f,0.0f);

    NEWMESH::Pt c = a*b;
    double len = SQR(c.X) + SQR(c.Y) + SQR(c.Z);
   

    if (len == 0.0){
      /* the vector b was parallel to a */
      b.X=0.0f;
      b.Y=1.0f;
      b.Z=0.0f;
      c = a*b;
      len = SQR(c.X) + SQR(c.Y) + SQR(c.Z);
    }
   
  /* normalize */
    len=std::sqrt(len);
    
    if (len == 0.0)
      {
	cout<<"Tangs Tangent:: first tangent vector of length zero at vertex: "<<ind<<endl;
	len = 1 ;
      }
    
    T.e1.X = c.X / len ;
    T.e1.Y = c.Y / len ;
    T.e1.Z = c.Z / len ;
    
    b = a*c;
    /* normalize */
    len=std::sqrt(SQR(b.X)+SQR(b.Y)+SQR(b.Z));
    
    if (len == 0)
      {
	cout<<"Tangs Tangent::second tangent vector of length zero at vertex: "<<ind<<endl;
	len = 1 ;
      }
    
    T.e2.X = b.X / len ;
    T.e2.Y = b.Y / len ;
    T.e2.Z = b.Z / len ;
    
    
    return T;
  }
  Tangs Tangent::calculate2(int ind,const NEWMESH::newmesh & SPH_in){

    // cout << " in calc2 " << endl;
    Tangs T;
    double mag;
    NEWMESH::Pt a = SPH_in.local_normal(ind); 
    Pt tmp=SPH_in.get_coord(ind);
    if((a|tmp)<0) a=a*-1;
  

    if(abs(a.X)>=abs(a.Y) && abs(a.X) >= abs(a.Z)){
      mag=sqrt(a.Z*a.Z + a.Y*a.Y);
      if(mag==0) {
	T.e1.X = 0 ;
	T.e1.Y = 0;
	T.e1.Z = 1;
      }else{
	T.e1.X = 0;
	T.e1.Y = -a.Z /mag ;
	T.e1.Z = a.Y /mag ;
      }
      //cout << " case 1 a " << a.X << " " << a.Y << " " << a.Z << " T.e1 " <<  T.e1.X << " " << T.e1.Y << " " << T.e1.Z << " mag " << mag << endl;

    }
    else if (abs(a.Y)>=abs(a.X) && abs(a.Y) >=abs(a.Z)){
      mag=sqrt(a.Z*a.Z + a.X*a.X);
      if(mag==0) {
	T.e1.X = 0 ;
	T.e1.Y = 0;
	T.e1.Z = 1;
      }else{
	T.e1.X = -a.Z/mag ;
	T.e1.Y = 0 ;
	T.e1.Z = a.X /mag ;
      }
    
    }
    else{
      mag=sqrt(a.Y*a.Y + a.X*a.X);
      if(mag==0) {
	T.e1.X = 1 ;
	T.e1.Y = 0;
	T.e1.Z = 0;
      }else{
	T.e1.X = -a.Y/mag ;
	T.e1.Y = a.X/mag ;
	T.e1.Z = 0 ;
      }

     
    }

    
    
    T.e2 = a*T.e1;
    T.e2.normalize();
    
    
    return T;
  }
  
  // project onto tangent plane
  void projectPoint(const NEWMESH::Pt &vb,const Tangs &T, double &e1coord, double & e2coord)//checked 
  { 
    e1coord=vb|T.e1;
    e2coord=vb|T.e2;
  }

  Matrix get_coordinate_transformation(const double &dNdT1,const double &dNdT2, ColumnVector &Norm){
    Pt G1(1,0,dNdT1);
    Pt G2(0,1,dNdT2);

    Pt G3=G1*G2;

   
    G3=G3/sqrt((G3|G3));

    Matrix G(3,3);

    G(1,1)=G1.X; G(1,2)=G2.X; G(1,3)=G3.X;
    G(2,1)=G1.Y; G(2,2)=G2.Y; G(2,3)=G3.Y;
    G(3,1)=G1.Z; G(3,2)=G2.Z; G(3,3)=G3.Z;
    Norm(1)=G3.X;  Norm(2)=G3.Y;  Norm(3)=G3.Z; 
    return (G.i()).t();
  }

  Matrix form_matrix_from_points(Pt p1, Pt p2,Pt p3, bool trans) {
    Matrix T(3,3);

    T(1,1)=p1.X; T(1,2)=p2.X; T(1,3)=p3.X;
    T(2,1)=p1.Y; T(2,2)=p2.Y; T(2,3)=p3.Y;
    T(3,1)=p1.Z; T(3,2)=p2.Z; T(3,3)=p3.Z;
  
    if(trans)
      return T.t();
    else 
      return T;
  }
 
  
  //// weighted least squares similarity gradient estimation
  newmesh WLS_gradient(const newmesh &SPHERE,const RELATIONS &REL,  double percentile)
  {
   
    NEWMESH::Pt point, v0, v1,v2,cr,grad;
    NEWMESH::newmesh GRADIENT=SPHERE;
    ColumnVector gradtmp;
    Matrix Disp,Dispinv,Data;
    Matrix weights;
    double x11,x21;
    double y11,y21;
    ColumnVector dist(2);
    Pt origin;
    double min_sigma=4;

    for (int index=0;index<SPHERE.nvertices();index++){
      Tangent Tang;
      Tangs T = Tang.calculate(index,SPHERE);
      gradtmp.ReSize(2); gradtmp=0.0;
    
      Disp.ReSize(REL.Nrows(index+1),2);
      weights.ReSize(REL.Nrows(index+1),REL.Nrows(index+1));
      Data.ReSize(REL.Nrows(index+1),1);
      Disp=0; weights=0;
      point=SPHERE.get_coord(index);


      origin=T.e1*T.e2;
      origin.normalize();
      origin=origin*RAD;

      projectPoint(point-origin,T,y11,y21);  
  
      for (int i=1;i<=REL.Nrows(index+1);i++){
	cr=SPHERE.get_coord(REL(i,index+1)-1);
	projectPoint(cr-origin,T,x11,x21);  
	dist(1)=x11-y11;dist(2)=x21-y21; 
    
        
     
       	
	
	Disp(i,1)=x11; Disp(i,2)=x21; 

	//////////// THE DISTANCES SHOULD BE CALCULATED FOR THE TANGENT PLANE///////////// - we ideally want an analytica expression

	weights(i,i)=exp(-(dist(1)*dist(1)+dist(2)*dist(2))/(2*min_sigma*min_sigma)); //
  
	Data(i,1)=SPHERE.get_pvalue(REL(1,index+1)-1)-SPHERE.get_pvalue(REL(i,index+1)-1); // difference of query point to mean similarity	
      }

      //estimates fits a gradient to the set samples and their similarity values using WLS
      Matrix tmp=((Disp.t()*weights)*Disp);
      if(tmp.Determinant()>1e-5){
	Dispinv=(tmp.i()*Disp.t())*weights; 
	gradtmp=Dispinv*Data; 
      }
    
      
      grad.X = (gradtmp(1)*(T.e1).X + gradtmp(2)*(T.e2).X) ;
      grad.Y = (gradtmp(1)*(T.e1).Y + gradtmp(2)*(T.e2).Y) ;
      grad.Z = (gradtmp(1)*(T.e1).Z + gradtmp(2)*(T.e2).Z) ;
      
    
      GRADIENT.set_pvalue(index,grad.norm());
      
    }

    if(percentile<100){
      GRADIENT.save("FULLGRAD.func");
      ColumnVector vals(SPHERE.nvertices()); vals=0;
      double thresh=0;
      for (int i=0;i<SPHERE.nvertices();i++){
	vals(i+1)=GRADIENT.get_pvalue(i); 
      }
      
      Histogram gradHist(vals,256);
      
      thresh=gradHist.getPercentile(percentile);
      for (int i=0;i<SPHERE.nvertices();i++){
	  if(GRADIENT.get_pvalue(i)<thresh) GRADIENT.set_pvalue(i,0);
      }
    }

    return GRADIENT;
    
  }
 

  newmesh projectmesh(newmesh ORIG, newmesh TARGET, newmesh ANAT){//, const boost::shared_ptr<RELATIONS> &rel){
    newmesh TRANS=ORIG;
    resampler R;
    vector<std::map<int,double> > baryweights=R.get_all_barycentric_weights(ORIG,TARGET);
    
     for (int Node = 0; Node < ORIG.nvertices(); ++Node)//this loop can't be parallelized
       {

	 Pt newcoord;
	 for (map<int, double >::iterator iter = baryweights[Node].begin(); iter != baryweights[Node].end(); ++iter)//convert scattering weights to gathering weights
	   {
	     if(ANAT.nvertices()==TARGET.nvertices())
	       newcoord+=ANAT.get_coord(iter->first)*iter->second;
	     else
	       newcoord+=TARGET.get_coord(iter->first)*iter->second;

	     
	   }
	 TRANS.set_coord(Node,newcoord);
       }
     return TRANS;
  }

  ColumnVector calculate_strains(const int &index,const vector<int> &kept, const  NEWMESH::newmesh & ORIG,const NEWMESH::newmesh & FINAL, const boost::shared_ptr<Matrix> &PrincipalStretches){

    ColumnVector STRAINS(4); STRAINS=0;
    Pt Normal_O, Normal_F;
    Matrix TRANS(3,3),principal_STRAINS(2,ORIG.nvertices());
    double T1_coord,T2_coord;
    Pt orig_tang,final_tang,tmp;
    Tangent Tang;  
    Matrix a,b,c,d; // coefficients
    int maxind, minind;

 
   
    Normal_O=ORIG.get_normal(index-1); Normal_F=FINAL.get_normal(index-1);
   
    Tangs T = Tang.calculate2(index-1,ORIG);
   
	  

    TRANS=form_matrix_from_points(T.e1,T.e2,Normal_O,true);
    
    Matrix pseudo_inv, pseudo_V, pseudo_U;
    DiagonalMatrix pseudo_D;
    Matrix alpha(kept.size(),6),N(kept.size(),1);
    Matrix n(kept.size(),1),t1(kept.size(),1),t2(kept.size(),1);
    
    alpha=0; N=0; n=0; t1=0;t2=0;
    
    for(int j=0;j<(int) kept.size();j++){
      tmp=ORIG.get_coord(kept[j])-ORIG.get_coord(index-1);
      projectPoint(tmp,T,T1_coord,T2_coord);
      
      N(j+1,1)=tmp|Normal_O;
      
      //// fit poly 
      alpha(j+1,2)=T1_coord; alpha(j+1,3)=T2_coord; alpha(j+1,4)=0.5*T1_coord*T1_coord; 
      alpha(j+1,5)=0.5*T2_coord*T2_coord; alpha(j+1,6)=T1_coord*T2_coord; 
      /// get transformed coords 
      tmp=FINAL.get_coord(kept[j])-FINAL.get_coord(index-1);
      projectPoint(tmp,T,T1_coord,T2_coord);
	    
      n(j+1,1)=tmp|Normal_O;
      t1(j+1,1)=T1_coord; t2(j+1,1)=T2_coord;
	
      }

    
      SVD(alpha, pseudo_D, pseudo_U, pseudo_V);  
      
      for(int i=1;i<=pseudo_D.Nrows();i++){
	if(pseudo_D(i)!=0)
	  pseudo_D(i)=1/pseudo_D(i);
      }
      pseudo_inv = pseudo_V*pseudo_D*pseudo_U.t();
       
      /// get coeffieicnts for different fits
      a = pseudo_inv *N ;
    
      
      b = pseudo_inv * t1;
      
      c = pseudo_inv * t2;
       
      d = pseudo_inv * n;
       
      double dNdT1, dNdT2, dt1dT1,dt1dT2,dt2dT1,dt2dT2,dndT1,dndT2;
      dNdT1=a(2,1); dNdT2=a(3,1);
      dt1dT1=b(2,1); dt1dT2=b(3,1);
      dt2dT1=c(2,1); dt2dT2=c(3,1);
      dndT1=d(2,1); dndT2=d(3,1);
      
      ColumnVector G3(3);
      Matrix G_cont=get_coordinate_transformation(dNdT1,dNdT2,G3);
      
        
      Pt g1(dt1dT1,dt2dT1,dndT1);
      Pt g2(dt1dT2,dt2dT2,dndT2);
      Pt g3=g1*g2; g3=g3/sqrt((g3|g3));
      
      Matrix g(3,3);
      g(1,1)=g1.X; g(1,2)=g2.X; g(1,3)=g3.X;
      g(2,1)=g1.Y; g(2,2)=g2.Y; g(2,3)=g3.Y;
      g(3,1)=g1.Z; g(3,2)=g2.Z; g(3,3)=g3.Z;
      
     
      Matrix F=g*G_cont.t();
      Matrix C= F.t()*F;
	
     
      DiagonalMatrix Omega;
      Matrix U,Umax,Umin;
      SVD(C,Omega,U); 
      
     
      Matrix mm=G3.t()*U;
    
      if(abs(mm(1,1))>=abs(mm(1,2)) && abs(mm(1,1))>=abs(mm(1,3))){
	if(sqrt(Omega(2)) > sqrt(Omega(3))){
	  maxind=2; minind=3;
	}
	else{   maxind=3; minind=2; }
      }
      else if(abs(mm(1,2))>=abs(mm(1,1)) && abs(mm(1,2))>=abs(mm(1,3))){
	if(sqrt(Omega(1)) > sqrt(Omega(3))){
	  maxind=1; minind=3;
	}
	else{  maxind=3; minind=1;}
      }
      else {
	if(sqrt(Omega(1)) > sqrt(Omega(2))){
	  maxind=1; minind=2;
	}
	else{ maxind=2; minind=1;}
      }
      
      STRAINS(1)=sqrt(Omega(maxind)); STRAINS(2)= sqrt(Omega(minind)); 
      Umax=TRANS.i()*U.SubMatrix(1,3,maxind,maxind);
      Umin=TRANS.i()*U.SubMatrix(1,3,minind,minind);
      if(PrincipalStretches.get()){
	
	(*PrincipalStretches)(index,1)=Umax(1,1); (*PrincipalStretches)(index,2)=Umax(2,1); (*PrincipalStretches)(index,3)=Umax(3,1);
	(*PrincipalStretches)(index,4)=Umin(1,1); (*PrincipalStretches)(index,5)=Umin(2,1); (*PrincipalStretches)(index,6)=Umin(3,1);  
      }
      STRAINS(3)=0.5*( STRAINS(1)*STRAINS(1) - 1); STRAINS(4)=0.5*( STRAINS(2)*STRAINS(2) - 1);
      
      
     
   return STRAINS;
 }
    

  newmesh calculate_strains(double fit_radius, const NEWMESH::newmesh & ORIG, const NEWMESH::newmesh & FINAL, const boost::shared_ptr<Matrix> &PrincipalStretches,const boost::shared_ptr<RELATIONS> REL){

    newmesh STRAIN=FINAL;
    Matrix STRAINS(4,ORIG.nvertices());
    ColumnVector nodestrain;
    double dir_chk1,dir_chk2,dist;

    vector<int> kept;
   
    int neighbour, numneighbours=ORIG.nvertices();
 
    if(PrincipalStretches.get()) PrincipalStretches->ReSize(ORIG.nvertices(),6) ; 
    //    cout << " fit radius " << fit_radius << endl;
    double fit_temp;
    kept.clear();
    for (int index=1; index <= ORIG.nvertices(); index++){
      if(REL.get()) numneighbours=REL->Nrows(index);
      fit_temp=fit_radius;
      while (kept.size()<=8){
		  kept.clear();
		  for (int j=1; j <=numneighbours; j++){

			if(REL.get()) neighbour=(*REL)(j,index); 
			else neighbour=j;
      
			dist=(ORIG.get_coord(index-1)-ORIG.get_coord(neighbour-1)).norm();
     
			/// reject points with normals in opposite directions
			dir_chk1=ORIG.get_normal(neighbour-1)|ORIG.get_normal(index-1);
			dir_chk2=1; //Normals_F[j-1]|Normals_F[index-1]; ////SET to 1?

			if(dist<=fit_temp && dir_chk1>=0 && dir_chk2>=0){
			kept.push_back(neighbour-1);
	 
			}
			}
		
    

		if(kept.size()>8){
			nodestrain=calculate_strains(index,kept,ORIG,FINAL,PrincipalStretches);
   
			STRAINS(1,index)=nodestrain(1);STRAINS(2,index)=nodestrain(2);STRAINS(3,index)=nodestrain(3);STRAINS(4,index)=nodestrain(4);
		}else{ fit_temp+=0.5;  }

     }
    kept.clear();
    }	
    STRAIN.set_pvalues(STRAINS);
    

    return STRAIN;

  }

   
  newmesh calculate_triangular_strains(const NEWMESH::newmesh & ORIG, const NEWMESH::newmesh & FINAL, double MU, double KAPPA){

    newmesh STRAIN=FINAL;
    Matrix STRAINS(3,ORIG.ntriangles());
   
    //// get normals //////////////////
    for (int index=0; index < ORIG.ntriangles(); index++){
      boost::shared_ptr<ColumnVector> strainstmp=boost::shared_ptr<ColumnVector>(new ColumnVector (2));
      STRAINS(3,index+1)=calculate_triangular_strain(index,ORIG,FINAL,MU,KAPPA,strainstmp);
      STRAINS(1,index+1)=(*strainstmp)(1);       STRAINS(2,index+1)=(*strainstmp)(2); 


      
    }

    for (int index=0; index < ORIG.nvertices(); index++){
      for (int j=0;j<3;j++){
	double SUM=0.0;
	for ( vector<int>::const_iterator i=ORIG.tIDbegin(index); i!=ORIG.tIDend(index); i++){
	     SUM+=STRAINS(j+1,*i+1);


	}
	SUM/= ORIG.get_total_triangles(index);
	STRAIN.set_pvalue(index,SUM, j);
      }
    }
    // 

    return STRAIN;
  }

  double calculate_triangular_strain(const Triangle & ORIG_tr, const Triangle & FINAL_tr, const double & mu, const double &kappa, bool calc_emery,  const boost::shared_ptr<ColumnVector> & indexSTRAINS){

    Pt Normal_O,Normal_F;
    Tangent Tang;  
    Matrix ORIG3D(3,3), FINAL3D(3,3);
    Matrix ORIG2D(3,3), FINAL2D(3,3);
    Matrix TRANS, TRANS2;

    Normal_O=ORIG_tr.normal(); Normal_F=FINAL_tr.normal() ;

  
    Tangs T = Tang.calculate_tri(Normal_O);
    Tangs T_trans = Tang.calculate_tri(Normal_F);
    double W;
    Matrix TMP;
    DiagonalMatrix D,D2;
   
    TRANS=form_matrix_from_points(T.e1,T.e2,Normal_O);
    if(TRANS.Determinant()<0){
      TMP=TRANS;
      TMP(1,1)=TRANS(1,2); TMP(2,1)=TRANS(2,2); TMP(3,1)=TRANS(3,2);
      TMP(1,2)=TRANS(1,1); TMP(2,2)=TRANS(2,1); TMP(3,2)=TRANS(3,1);
      TRANS=TMP;
    }
    TRANS2=form_matrix_from_points(T_trans.e1,T_trans.e2,Normal_F);
    if(TRANS.Determinant()<0){
      TMP=TRANS2;
      TMP(1,1)=TRANS2(1,2); TMP(2,1)=TRANS2(2,2); TMP(3,1)=TRANS2(3,2);
      TMP(1,2)=TRANS2(1,1); TMP(2,2)=TRANS2(2,1); TMP(3,2)=TRANS2(3,1);
      TRANS2=TMP;
    }
    
    for(int i=0;i<3;i++){
      Pt vertex=ORIG_tr.get_vertex_coord(i); 
      ORIG3D(i+1,1)=vertex.X;	ORIG3D(i+1,2)=vertex.Y;	ORIG3D(i+1,3)=vertex.Z;
  
      vertex=FINAL_tr.get_vertex_coord(i); //get_triangle_vertex(index,i);
      FINAL3D(i+1,1)=vertex.X;	FINAL3D(i+1,2)=vertex.Y;	FINAL3D(i+1,3)=vertex.Z;

    }

    
  
    ORIG2D=ORIG3D*TRANS;
    FINAL2D=FINAL3D*TRANS2;

  
   
    W=triangle_strain(ORIG2D,FINAL2D,mu,kappa,indexSTRAINS,calc_emery);
 
  
    return W;
  }

  double calculate_triangular_strain(int index, const NEWMESH::newmesh & ORIG, const NEWMESH::newmesh & FINAL, const double & mu, const double &kappa, bool calc_emery,  const boost::shared_ptr<ColumnVector> & indexSTRAINS){

    Pt Normal_O,Normal_F;
    Tangent Tang;  
    Matrix ORIG3D(3,3), FINAL3D(3,3);
    Matrix ORIG2D(3,3), FINAL2D(3,3);
    Matrix TRANS, TRANS2;
    Normal_O=ORIG.get_triangle_normal(index); Normal_F=FINAL.get_triangle_normal(index);
  
    Tangs T = Tang.calculate_tri(index,ORIG);
    Tangs T_trans = Tang.calculate_tri(index,FINAL);
    double W;
    Matrix TMP;
    DiagonalMatrix D,D2;
  

    TRANS=form_matrix_from_points(T.e1,T.e2,Normal_O);
    if(TRANS.Determinant()<0){
      TMP=TRANS;
      TMP(1,1)=TRANS(1,2); TMP(2,1)=TRANS(2,2); TMP(3,1)=TRANS(3,2);
      TMP(1,2)=TRANS(1,1); TMP(2,2)=TRANS(2,1); TMP(3,2)=TRANS(3,1);
      TRANS=TMP;
    }
    TRANS2=form_matrix_from_points(T_trans.e1,T_trans.e2,Normal_F);
    if(TRANS.Determinant()<0){
      TMP=TRANS2;
      TMP(1,1)=TRANS2(1,2); TMP(2,1)=TRANS2(2,2); TMP(3,1)=TRANS2(3,2);
      TMP(1,2)=TRANS2(1,1); TMP(2,2)=TRANS2(2,1); TMP(3,2)=TRANS2(3,1);
      TRANS2=TMP;
    }
  
    for(int i=0;i<3;i++){
      Pt vertex=ORIG.get_triangle_vertex(index,i);
      ORIG3D(i+1,1)=vertex.X;	ORIG3D(i+1,2)=vertex.Y;	ORIG3D(i+1,3)=vertex.Z;
  
      vertex=FINAL.get_triangle_vertex(index,i);
      FINAL3D(i+1,1)=vertex.X;	FINAL3D(i+1,2)=vertex.Y;	FINAL3D(i+1,3)=vertex.Z;

    }

 
    ORIG2D=ORIG3D*TRANS;
    FINAL2D=FINAL3D*TRANS2;

   
    W=triangle_strain(ORIG2D,FINAL2D,mu,kappa,indexSTRAINS,calc_emery);
 
    
    return W;
  }

  /// taken from matlab code Alina Oltean, Oct. 2014, alinabme[at]gmail.com
  // http://web.stanford.edu/class/cs205b/lectures/lecture7.pdf
  //  for  principal strains: http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/plane_strain_principal.cfm - 2d eigenvalue solution to matrix of normal and shear strains



  double triangle_strain(const Matrix& AA, const Matrix & BB, const double &MU,const double &KAPPA, const boost::shared_ptr<ColumnVector> & strains, bool emery_strain){
    double c0,c1,c2,c3,c4,c5;
    double c0c,c1c,c2c,c3c,c4c,c5c;
    double I1,I3,I1st,J,W;
    Matrix Edges(2,2),edges(2,2);
    Matrix F,F3D(3,3),F3D_2;

    F3D(3,1)=0; F3D(3,2)=0; F3D(3,3)=1;

    // all t0 distances
    c0=AA(2,1)-AA(1,1); //x2-x1 t0
    c1=AA(2,2)-AA(1,2); //y2-y1
    c2=AA(3,1)-AA(2,1); //x3-x2 
    c3=AA(3,2)-AA(2,2); // y3-y2
    c4=AA(3,1)-AA(1,1); // x3-x1
    c5=AA(3,2)-AA(1,2); // y3-y1 
    
    //all tf distances
    c0c=BB(2,1)-BB(1,1); //x2-x1 tf
    c1c=BB(2,2)-BB(1,2); // y2-y1
    c2c=BB(3,1)-BB(2,1); // x3-x2 
    c3c=BB(3,2)-BB(2,2); // y3-y2
    c4c=BB(3,1)-BB(1,1);  // x3-x1
    c5c=BB(3,2)-BB(1,2); // y3-y1 
    
    
    Edges(1,1)=c0; Edges(1,2)=c4; Edges(2,1)=c1; Edges(2,2)=c5;
    edges(1,1)=c0c; edges(1,2)=c4c; edges(2,1)=c1c; edges(2,2)=c5c;

    
    F = edges*Edges.i();

  
    F3D(1,1)=F(1,1); F3D(1,2)=F(1,2); F3D(1,3)=0; 
    F3D(2,1)=F(2,1); F3D(2,2)=F(2,2); F3D(2,3)=0; 
  
    F3D_2=F3D.t()*F3D;
    I1=F3D_2.Trace();
    I3=F3D_2.Determinant();
    I1st = I1*std::pow(I3,(-1.0/3.0));
    J = sqrt(I3);
    
    if(emery_strain==true){
      W=0.5*MU*(I1-3-2*log(J))+0.5*KAPPA*log(J)*log(J);
    }
    else
      W = 0.5*(MU*(I1st - 3) + KAPPA*(J - 1)*(J - 1));

   
    //// calculating prinicipal strains (for testing)
    if(strains.get()){
      double a11,a21,a31,a12,a22,a32,a13,a23,a33;
      double A, A11,A21,A31,A12,A22,A32,A13,A23,A33;
      double B1,B2,B3;
      double e11,e22,e12,X,Y;
      //2*dx^2, then 2*dy^2, then 2*dx*dy
      a11=2.*c0*c0 ;
      a21=2.*c2*c2 ;
      a31=2.*c4*c4;
      a12=2.*c1*c1 ;
      a22=2.*c3*c3 ;
      a32=2.*c5*c5;
      a13=4.*c0*c1 ;
      a23=4.*c2*c3;
      a33=4.*c4*c5;

      
      A=a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a23*a32*a11-a33*a21*a12;

      A11=(a22*a33-a32*a23)/A;
      A12=(a13*a32-a12*a33)/A;
      A13=(a12*a23-a13*a22)/A;
      A21=(a23*a31-a21*a33)/A;
      A22=(a11*a33-a13*a31)/A;
      A23=(a13*a21-a11*a23)/A;
      A31=(a21*a32-a22*a31)/A;
      A32=(a12*a31-a11*a32)/A;
      A33=(a11*a22-a12*a21)/A;


      // deformed distances between points
      B1=c0c*c0c+c1c*c1c-c0*c0-c1*c1;
      B2=c2c*c2c+c3c*c3c-c2*c2-c3*c3;
      B3=c4c*c4c+c5c*c5c-c4*c4-c5*c5;

      //     strains wrt x,y coords

      e11=B1*A11+B2*A12+B3*A13;
      e22=B1*A21+B2*A22+B3*A23;
      e12=B1*A31+B2*A32+B3*A33;

    
      X=e11+e22 ;
      Y=e11-e22;

      (*strains)(1)=X/2+sqrt((Y/2)*(Y/2)+(e12)*(e12)) ;																					
      (*strains)(2)=X/2-sqrt((Y/2)*(Y/2)+(e12)*(e12));  

     
    }

    return W;
  }

  void mean_curvature(double fit_radius, NEWMESH::newmesh & MESH, const boost::shared_ptr<RELATIONS> &REL){

    double dir_chk1,dist;
    double T1_coord,T2_coord;
    double Dxx,Dyy,Dxy,Cmean;
    Pt orig_tang,final_tang,tmp;
    vector<int> kept;
    vector<Pt> Normals_O;
    Tangent Tang;  
    ColumnVector Lambdas; 
    Matrix a;
    int neighbour, numneighbours=MESH.nvertices();

    //// get normals //////////////////
    for (int index=0; index < MESH.nvertices(); index++){
      Pt N= MESH.local_normal(index);
         Normals_O.push_back(N);
   
    }

   
   
    for (int index=1; index <= MESH.nvertices(); index++){
      Tangs T = Tang.calculate2(index-1,MESH);
     
    
      kept.clear();
    
      if(REL.get()) numneighbours=REL->Nrows(index);
      
      for (int j=1; j <=numneighbours; j++){

	if(REL.get()) neighbour=(*REL)(j,index); 
	else neighbour=j;

	dist=(MESH.get_coord(index-1)-MESH.get_coord(neighbour-1)).norm();

	/// reject points with normals in opposite directions
	dir_chk1=Normals_O[neighbour-1]|Normals_O[index-1];

	if(dist<=fit_radius && dir_chk1>=0 ){
	  kept.push_back(neighbour-1);
	  
	  
	}
      }
     

      if(kept.size()>=8){
	
	Matrix pseudo_inv, pseudo_V, pseudo_U;
	DiagonalMatrix pseudo_D;
	Matrix alpha(kept.size(),6),N(kept.size(),1);
	Matrix n(kept.size(),1),t1(kept.size(),1),t2(kept.size(),1);
	
	alpha=0; N=0; n=0; t1=0;t2=0;
      
	for(int j=0;j<(int) kept.size();j++){
	  tmp=MESH.get_coord(kept[j])-MESH.get_coord(index-1);
	  projectPoint(tmp,T,T1_coord,T2_coord);
	  
	  N(j+1,1)=tmp|Normals_O[index-1];
	
	  //// fit poly 
	  alpha(j+1,1)=T1_coord*T1_coord; alpha(j+1,2)=T2_coord*T2_coord; alpha(j+1,3)=T1_coord*T2_coord; 
	  alpha(j+1,4)=T1_coord; alpha(j+1,5)=T2_coord; alpha(j+1,6)=1; 
	
	/// get transformed coords 
	 	  
	}


	SVD(alpha, pseudo_D, pseudo_U, pseudo_V);  


	for(int i=1;i<=pseudo_D.Nrows();i++){
	  if(pseudo_D(i)!=0)
	    pseudo_D(i)=1/pseudo_D(i);
	}
	 
	pseudo_inv = pseudo_V*pseudo_D*pseudo_U.t();

	/// get coeffieicnts for different fits
	a = pseudo_inv *N ;

	// Make hessian
	Dxx=2*a(1,1); Dxy=a(3,1); Dyy=2*a(2,1);
	Lambdas=eig2(Dxx,Dyy,Dxy);
	Cmean=(Lambdas(1)+Lambdas(2))/2;
	if(fabs(Cmean)>=10) {cout << index << " Cmean out of range kept.size() " << kept.size() << " Cmean " << Cmean << endl;
	  cout << a << endl;
	}
	MESH.set_pvalue(index-1,Cmean);
      }else cout << "kept.size <8 " << endl;

    }
  }

  /// taken from matlab patchcurvature - only need eigenvalues
  ColumnVector eig2(const double & Dxx,const double & Dyy,const double & Dxy){
    ColumnVector Lambda(2); 
    double tmp,mu1,mu2;
   
    tmp = sqrt((Dxx - Dyy)*(Dxx - Dyy) + 4*Dxy*Dxy);
   

    // Compute the eigenvalues
    mu1 = (0.5*(Dxx + Dyy + tmp));
    mu2 = (0.5*(Dxx + Dyy - tmp));

    //Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
    if(fabs(mu1)<fabs(mu2)){
      Lambda(1)=mu1;
      Lambda(2)=mu2;
      //  I2=[v1x v1y]; // directions not needed for mean curvature?
      //I1=[v2x v2y];
    }
    else{
      Lambda(1)=mu2;
      Lambda(2)=mu1;
      // I2=[v2x v2y];
      //I1=[v1x v1y];
    }
    return Lambda;
  }
  //####################################################################
  Histogram build_histogram(const ColumnVector& M,const string &excl,const int &bins){
    
    NEWMESH::newmesh EXC;
    ColumnVector excludedvalues(M.Nrows());
    excludedvalues=1;
    
    if(excl!=""){ EXC.load(excl);} 
    for (int i=1;i<=M.Nrows();i++){
      if(excl!="" ){
	excludedvalues(i)=EXC.get_pvalue(i-1);
      }
      
    }
    
    Histogram Hist(M,bins);
    Hist.generate(excludedvalues);
    Hist.setexclusion(excludedvalues);
    return Hist;
  }
  
  //####################################################################
  Histogram build_histogram(const ColumnVector& M,boost::shared_ptr<NEWMESH::newmesh> EXCL,const int &bins){
    
    
    ColumnVector excludedvalues(M.Nrows());
    excludedvalues=1;
    
    for (int i=1;i<=M.Nrows();i++){
      if(EXCL.get()){
	excludedvalues(i)=EXCL->get_pvalue(i-1);
      }
    
      
    }
    
    Histogram Hist(M,bins);
    Hist.generate(excludedvalues);
    Hist.setexclusion(excludedvalues);
    return Hist;
  }

  void histogram_normalization(NEWMESH::newmesh  & IN, const NEWMESH::newmesh & REF,const string &excl_in,const string &excl_ref, int numbins){
    
    ColumnVector CDF_in,CDF_ref;
    Histogram Hist_in,Hist_ref; 
    ColumnVector datain(IN.nvertices()),dataref(REF.nvertices());
    
    for (int i=0;i<IN.nvertices();i++)
      datain(i+1)=IN.get_pvalue(i);

    for (int i=0;i<REF.nvertices();i++)
      dataref(i+1)=REF.get_pvalue(i);

    Hist_in=build_histogram(datain,excl_in,numbins);
    Hist_ref=build_histogram(dataref,excl_ref,numbins);
  
    
    Hist_in.generateCDF(); Hist_ref.generateCDF();
    
    Hist_in.match(Hist_ref);   /// performs histogram normalisation by finding intensity bin for reference that will generate same CDF as for current input bin and then mapping reference intensity across histogram http://en.wikipedia.org/wiki/Histogram_matching

    
    datain =Hist_in.getsourceData();
    Hist_in.generateCDF();
    CDF_in=Hist_in.getCDF();

    for (int i=1;i<=IN.nvertices();i++)
      IN.set_pvalue(i-1,datain(i));
    
  }
  
  /// multivariate does the same as above for each feature separately
  void multivariate_histogram_normalization(BFMatrix &IN , BFMatrix & REF,boost::shared_ptr<NEWMESH::newmesh> EXCL_IN,boost::shared_ptr<NEWMESH::newmesh> EXCL_REF, bool rescale){
    
    ColumnVector CDF_in,CDF_ref;
    ColumnVector datain(IN.Ncols()),dataref(REF.Ncols());
    double max,min;
    ColumnVector excluded_in(IN.Ncols()); excluded_in=1;
    ColumnVector excluded_ref(REF.Ncols());excluded_ref=1;
    int numbins=256;
   
    for(int d=1; d<=(int) IN.Nrows();d++){      

      for (unsigned int i=1;i<=IN.Ncols();i++){
	datain(i)=IN.Peek(d,i);

	if(EXCL_IN.get()){ // if using an exclusion mask these values will be eliminated from the histogram matching
	  if(EXCL_IN->get_dimension()>=d)
	    excluded_in(i)=EXCL_IN->get_pvalue(i-1,d-1);
	  else
	    excluded_in(i)=EXCL_IN->get_pvalue(i-1);

	}
      }
    

      for (unsigned int i=1;i<=REF.Ncols();i++){	

	dataref(i)=REF.Peek(d,i);

	if(EXCL_REF.get()){
	  if(EXCL_REF->get_dimension()>=d)	    
	    excluded_ref(i)=EXCL_REF->get_pvalue(i-1,d-1);
	  else
	    excluded_ref(i)=EXCL_REF->get_pvalue(i-1);
	}
		
      }

      Histogram Hist_in(datain,numbins),Hist_ref(dataref,numbins); 
      Hist_in.generate(excluded_in);
      Hist_in.setexclusion(excluded_in);
      Hist_ref.generate(excluded_ref);
      Hist_ref.setexclusion(excluded_ref);
      Hist_in.generateCDF(); Hist_ref.generateCDF();
      Hist_in.match(Hist_ref);

      datain =Hist_in.getsourceData();
      // Hist_in.generateCDF();
      CDF_in=Hist_in.getCDF();
      for (unsigned int i=1;i<=IN.Ncols();i++){
	IN.Set(d,i,datain(i));

      }

      if(rescale){
	if(d>1){
	  set_range(d,IN,excluded_in,max,min);
	  set_range(d,REF,excluded_ref,max,min);
	  
	}else{	
	  get_range(d,IN,excluded_in,max,min);
	}
      }
    }
  }

  void get_range(const int &dim, const BFMatrix &M, const ColumnVector & excluded, double & min, double& max){
    max=0;
    int n=0;
    for (unsigned int i=1;i<=M.Ncols();i++){
      if(M.Peek(dim,i)>max && excluded(i)) max=M.Peek(dim,i);
      if(!excluded(i)) n++;
    }
    min=max;
    for (unsigned int i=1;i<=M.Ncols();i++){
      if(M.Peek(dim,i)<min && excluded(i)) min=M.Peek(dim,i);
    }
  }

  void set_range(const int &dim, BFMatrix &M, const  ColumnVector & excluded, double & min, double& max){

    double nmin,nmax=0;
    double nrange,range;
    double nval, val;
    int n=0;
    range=max-min;
    for (unsigned int i=1;i<=M.Ncols();i++){
      if(M.Peek(dim,i)>nmax && excluded(i)) nmax=M.Peek(dim,i);
    }
    nmin=nmax;
    for (unsigned int i=1;i<=M.Ncols();i++){
      if(M.Peek(dim,i)<nmin && excluded(i)) nmin=M.Peek(dim,i);
      if(!excluded(i)) n++;
    }
    nrange=nmax-nmin;

    for (unsigned int i=1;i<=M.Ncols();i++){
      if(excluded(i)){
	nval=(M.Peek(dim,i)-nmin)/nrange;
	val=nmin+nval*range;
	M.Set(dim,i,val);
      }
    }
  }

  ///// get rotation using euler angles, used for current affine implementation
  NEWMAT::ReturnMatrix  rotate_euler(const ColumnVector &V,const double &w1,const double &w2,const double &w3){
    
    NEWMAT::Matrix Rot(3,3);
    Rot << cos(w2)*cos(w3) << -cos(w1)*sin(w3) + sin(w1)*sin(w2)*cos(w3) << sin(w1)*sin(w3) + cos(w1)*sin(w2)*cos(w3)
	<< cos(w2)*sin(w3) << cos(w1)*cos(w3) + sin(w1)*sin(w2)*sin(w3) << -sin(w1)*cos(w3) + cos(w1)*sin(w2)*sin(w3)
	<< -sin(w2) << sin(w1)*cos(w2) << cos(w1)*cos(w2);
    Rot = Rot.t();
    
    NEWMAT::ColumnVector Vrot(3);
    
    Vrot = Rot*V;
    Vrot.Release();
    return Vrot;
    
  }


  //////////////////////// BASIC DATA PROCESSING /////////////////////////
  bool set_data(const  string &dataname,boost::shared_ptr<BFMatrix> &BF,NEWMESH::newmesh &M, bool issparse){

    
    bool isfunc=false;

    int filetype=meshFileType(dataname); // will detect file type from ending, data matrices need to be ended .txt at the moment

   
    if(issparse){
      SpMat<double> sparse_mat(dataname);
      BF=boost::shared_ptr<BFMatrix > (new SparseBFMatrix<double>(sparse_mat)); // may not be desirable to do it this way if dimensions are v high

    }
    else{
      M.load(dataname,false,false);
      Matrix tmp;
      double sum=0;
      tmp=M.get_pvalues();
      sum=0;
      for(int i=0;i<M.npvalues();i++){
	sum+=M.get_pvalue(i);
      }
      if(sum==0)throw  NEWMESHException("No data has been suppled");
      isfunc=true;

      BF=boost::shared_ptr<BFMatrix > (new FullBFMatrix(tmp)); 
      
    }
   
    if((int) BF->Ncols()!=M.npvalues()){ 
      if((int) BF->Nrows()!=M.npvalues())throw  NEWMESHException("data does not match mesh dimensions");
      else
	BF=BF->Transpose();
    }
    return isfunc;
  }

  void logtransform(BFMatrix &features){
    
    
    //#pragma omp parallel for 
    for(int i=1;i<=(int) features.Ncols();i++){

      for (BFMatrixColumnIterator it = features.begin(i); it != features.end(i); ++it) {
	features.Set(it.Row(),i,log10(*it+1));
	    
      }

      
    }
    
  }
  
  void normalise(BFMatrix &features){
    
    double size=0, mean;
    for(int i=1;i<= (int) features.Ncols();i++){
	size=0;
	mean=0;
	for (BFMatrixColumnIterator it = features.begin(i); it != features.end(i); ++it) 
	  mean+=(*it);

	mean/=features.Nrows();

	for (BFMatrixColumnIterator it = features.begin(i); it != features.end(i); ++it) 
	  size+=(*it-mean)*(*it-mean);
		
	size=sqrt(size);
	
	for (BFMatrixColumnIterator it = features.begin(i); it != features.end(i); ++it) {
	  if(*it>0) features.Set(it.Row(),i,*it/size);
	}
	
      }
          
  }
  

  void  log_transform_and_normalise( BFMatrix  &data){
    logtransform(data);
    normalise(data);    
  }

  
 ///////////////////////////////   UNFOLDING METHOD //////////////
   bool check_for_intersections(const int ind, const double &eps, NEWMESH::newmesh &IN){  // not needed if we have negdet?
    
    NEWMESH::Pt c = IN.get_coord(ind);
    NEWMESH::Pt V=c;
    V.normalize(); 
    int a = 0;
         

    NEWMESH::Triangle tr=IN.get_triangle_from_vertex(ind,0);
    c = IN.get_coord(ind);

    NEWMESH::Pt N = tr.normal();
    for (vector<int>::const_iterator j=IN.tIDbegin(ind) ; j!=IN.tIDend(ind); j++){
      NEWMESH::Triangle tr2=IN.get_triangle(*j);
      NEWMESH::Pt N2 = tr2.normal();

      a = a || ((N|N2) <= 0.5);

      if(a==1){
	break;
      }
    }


    return a;
  }

  /////////// spatial gradient for vertex trianglular mesh
  NEWMESH::Pt spatialgradient(const int &index,const  NEWMESH::newmesh &_SOURCE){
    
    NEWMESH::Pt v0,v1,v2;
    NEWMESH::Pt dA,grad;
    NEWMESH::Pt ci  = _SOURCE.get_coord(index);
    
    for ( vector<int>::const_iterator j=_SOURCE.tIDbegin(index); j!=_SOURCE.tIDend(index); j++){
      
      v0 = _SOURCE.get_triangle_vertex(*j,0);  v1 =_SOURCE.get_triangle_vertex(*j,1); v2 =_SOURCE.get_triangle_vertex(*j,2);

      
      if((ci-v0).norm()==0){
	dA=computeGradientOfBarycentricTriangle(v1, v2, v0);
      }
      else if((ci-v1).norm()==0){
	dA=computeGradientOfBarycentricTriangle(v2, v0, v1);
      }
      else{
	dA=computeGradientOfBarycentricTriangle(v0, v1, v2);
      }
      
      grad=grad+dA;
    }
    
   
    return grad;
    
  }

 
  void unfold(NEWMESH::newmesh &_SOURCE){

    /// check for intersections
    vector<int> foldedvertices;
    vector<NEWMESH::Pt> foldinggradients; 
    NEWMESH::Pt grad,ppp,ci;
    double current_stepsize;
    NEWMESH::newmesh tmp=_SOURCE;
    NEWMESH::Pt pp;
    bool folded=true;
    int it=0;
  
    while (folded){

      foldedvertices.clear();
      foldinggradients.clear();
   
      for (int i=0;i<_SOURCE.nvertices();i++){
       if(check_for_intersections(i,0,_SOURCE)){  /// identfy points that are flipped.
        foldedvertices.push_back(i);
        }
      }
    
      if(foldedvertices.size()==0) folded=false;

      else{ if(it%100==0) cout << " mesh is folded, total folded vertices =" << foldedvertices.size() << " unfold " <<  endl;}

      for (unsigned int i=0;i<foldedvertices.size();i++){
        grad=spatialgradient(foldedvertices[i],_SOURCE);  /// estimates gradient of triangle area for all vertices that are detected as folded
        foldinggradients.push_back(grad);
      }
  
  
      for (unsigned int i=0;i<foldedvertices.size();i++){

	current_stepsize=1;
	ci=_SOURCE.get_coord(foldedvertices[i]);
	  grad=foldinggradients[i];
        do{
            pp =ci - grad*current_stepsize; // gradient points in direction of increasing area therefore opposite to desired direction
            pp.normalize();
            _SOURCE.set_coord(foldedvertices[i],pp*RAD);
            current_stepsize/=2.0;
           
         
        }while(check_for_intersections(foldedvertices[i],0,_SOURCE)==1 && current_stepsize > 1e-3);
        
        _SOURCE.set_coord(foldedvertices[i],pp*RAD);
    
      }
      
     it++;

     if(it==1000) break;
     
    }
  }

  // inspired by SPHERICAL DEMONS MATLAB CODE YEO MARS_linearInterp.h 
  // Find the unit normal to the vertex v0 then cross with edge v1v0 to find vector in the plane of triangle perpendicular to edge v1v0
  // This is the gradient of triangle area for the triangle v2v1v0 where v2 is the moving point.
   //
  void computeNormal2EdgeOfTriangle(const NEWMESH::Pt &v0,const NEWMESH::Pt &v1,const NEWMESH::Pt &v2, NEWMESH::Pt &norm2edge)
  {
    
    NEWMESH::Pt s1 = v2 - v0;
    NEWMESH::Pt s2 = v1 - v0;
    
    if(s1.norm()>1e-10)  s1.normalize(); else s1=NEWMESH::Pt(0,0,0);
    if(s2.norm()>1e-10)  s2.normalize(); else s2=NEWMESH::Pt(0,0,0); /// if these are very very small (rounding errors) all this goes horribly wrong
    
    NEWMESH::Pt norm2triangle = s1*s2;/// s3= normal to triangle at v0
    
    
    if(norm2triangle.norm()>1e-10) norm2triangle.normalize();else norm2triangle=NEWMESH::Pt(0,0,0);
    
    
    norm2edge=s2*norm2triangle; //Then take cross product of norm and the edge v1-v0 to get vector perpendicular to v1-v0
    
    if((s1|norm2edge) < 0)  /*first edge is in wrong direction, flip it.*/
      {
	norm2edge= norm2edge*-1;
      }    
  }
  
  
  // Find the gradient of triangle area associated with edge v1v0 (end moving point v2) points away from v1v0 with magnitude equal to half the length of v1v0
  // used for unfolding
  NEWMESH::Pt computeGradientOfBarycentricTriangle(const NEWMESH::Pt &v0, const NEWMESH::Pt &v1, const NEWMESH::Pt &v2)
  {
    NEWMESH::Pt v0v1, norm2edge;
    double base;
    
    computeNormal2EdgeOfTriangle(v0, v1, v2, norm2edge);
    v0v1 = v1 - v0;
    
    base = v0v1.norm();
    
    return norm2edge*0.5*base;
 
  } 
  
  ColumnVector barycentricSurfaceGrad(const int &index,const NEWMESH::newmesh &IN)	     
  {
    int n0,n1,n2;
    NEWMESH::Pt point, v0, v1,v2,cr;
    ColumnVector temp_grad,grad;
    double d0,d1,d2;
    
    double area, total_area=0;
    
    n0=0;n1=0;n2=0;
    
    point=IN.get_coord(index);
    grad.ReSize(3); grad=0;

    
    /// problem if closest point coincides with an Input vertex
    // if so this version assume points match exactly
    for (vector<int>::const_iterator j=IN.tIDbegin(index);j!=IN.tIDend(index); j++){
      
      
      v0 = IN.get_triangle_vertex(*j,0);  v1 = IN.get_triangle_vertex(*j,1); v2 = IN.get_triangle_vertex(*j,2);
      
      n0 = IN.get_triangle_vertexID(*j,0);
      n1 = IN.get_triangle_vertexID(*j,1);
      n2=  IN.get_triangle_vertexID(*j,2);
    
      
      d0 = IN.get_pvalue(n0); // gets the similarity terms for these
      d1 = IN.get_pvalue(n1); 
      d2 = IN.get_pvalue(n2); 
      
      point=(v0+v1+v2)/3 ;//// find the barycentric graient for each of the centre triangles and then take a weighted (by area) average
      area=barycentricGradforInd(point,v0,v1,v2,d0,d1,d2,temp_grad);
      total_area+=area;
      
      
      grad+=area*temp_grad; 
      
     
    }
    grad/=total_area;
    
    
    
  return grad;
  }
  
  //// gets barycentric gradient for a given point (at the barycentre of a triangle) or equivalently for one vertex of a triangle
  double barycentricGradforInd(const NEWMESH::Pt &point, const NEWMESH::Pt &v0,const NEWMESH::Pt &v1,const NEWMESH::Pt &v2, const double &d0,const double &d1,const double &d2, ColumnVector &GRAD)	    
  {
    
    double A0, A1, A2, totalA;
    NEWMESH::Pt projected_point, norm, dA0, dA1, dA2;
    ColumnVector temp_grad(3);
    Matrix temp_mat(3,3); 
    Matrix I(3,3),tmp(3,3);
    Matrix test;
    
    norm=projectPoint(point, v0, v1, v2, projected_point);
    
    dA2=computeGradientOfBarycentricTriangle(v0, v1, projected_point); // should be the gradient of triangle pp0p1 according to text in yeo
    
    dA1=computeGradientOfBarycentricTriangle(v2, v0,projected_point);
    
    dA0=computeGradientOfBarycentricTriangle(v1, v2, projected_point);
    
    A0 = computeArea(projected_point, v1, v2);
    A1 = computeArea(projected_point, v0, v2);
    A2 = computeArea(projected_point, v0, v1);    
    
    totalA = A0 + A1 + A2;
    
    /// Yeo had something here for multidimensional data ......
    
    temp_grad(1) = (dA0.X*d0 + dA1.X*d1 + dA2.X*d2)/totalA;
    temp_grad(2) = (dA0.Y*d0 + dA1.Y*d1 + dA2.Y*d2)/totalA;
    temp_grad(3) = (dA0.Z*d0 + dA1.Z*d1 + dA2.Z*d2)/totalA;
    

    float C1, C2;
    C1 = (v0|norm)/(point|norm);
    C2 = C1/(point|norm);
    
    I << C1 << 0 << 0
      << 0 << C1 << 0
      << 0 << 0 << C1;
   
    tmp << point.X*norm.X << point.X*norm.Y << point.X*norm.Z
	<< point.Y*norm.X << point.Y*norm.Y << point.Y*norm.Z
	<< point.Z*norm.X << point.Z*norm.Y << point.Z*norm.Z;
    
    
    temp_mat=I-C2*tmp;
    
    
    GRAD=(temp_grad.t()*temp_mat).t();

    
    return totalA;
  }
  
  /////////////////////////////////// MESH CHECKING /////////////////////////////
  double getarealseparation(const int &index, const NEWMESH::newmesh &ORIG,  NEWMESH::newmesh &OUT){
    
    double distortion=0;
    Pt v0,v1,v2;
    Pt a_v0,a_v1,a_v2;

    double area1=0, area2=0;
    double meandistortion=0.0;

    if(ORIG.nvertices()!=OUT.nvertices()){ throw  NEWMESHException("NEWMESH ERROR in meshfns:: getarealdistortion initial and distorted surfaces should have them same number of vertices");}
    
    for ( vector<int>::const_iterator j=ORIG.tIDbegin(index); j!=ORIG.tIDend(index); j++)
      {
	v0 = OUT.get_triangle_vertex(*j,0);  v1 =OUT.get_triangle_vertex(*j,1); v2 =OUT.get_triangle_vertex(*j,2);
	a_v0 = ORIG.get_triangle_vertex(*j,0);  a_v1 =ORIG.get_triangle_vertex(*j,1); a_v2 =ORIG.get_triangle_vertex(*j,2);

	area2=computeArea(v0, v1, v2);
	area1=computeArea(a_v0, a_v1, a_v2);
	distortion=((area2-area1)/area1)*((area2-area1)/area1); // areal distortion measure used in caret
	meandistortion+=distortion; 

      }
    
    meandistortion/=OUT.get_total_triangles(index);
  
    return meandistortion;
  }

  ColumnVector getarealdistortion(const NEWMESH::newmesh &ORIG,  NEWMESH::newmesh &OUT){
    
    cout << " here " << endl;
    double distortion=0;
    double area1=0, area2=0;
    ColumnVector meandistortion(OUT.nvertices());
    ColumnVector weights(OUT.nvertices());
    int faces=0;
    meandistortion=0;weights=0;
    if(ORIG.nvertices()!=OUT.nvertices()){ throw  NEWMESHException("NEWMESH ERROR in meshfns:: getarealdistortion initial and distorted surfaces should have them same number of vertices");}
    int ind=0;
    for ( vector<NEWMESH::Triangle>::const_iterator i=ORIG.tbegin() ; i!=ORIG.tend(); i++)
      {
	NEWMESH::Pt v1 = (*i).get_vertex_coord(0),  v2 = (*i).get_vertex_coord(1), v3 = (*i).get_vertex_coord(2);
	  ind++;
	area1=computeArea(v1, v2, v3);
	if(area1==0) cout <<ind << "area1 v1-v2 " << (v1-v2).norm() << " (v3-v2).norm() " << (v3-v2).norm() << "v1-v3 " << (v1-v3).norm() << " aera1 " << area1 <<" area2 " << area2 <<  endl;
	
	v1=OUT.get_coord((*i).get_vertex_no(0)); v2=OUT.get_coord((*i).get_vertex_no(1)); v3=OUT.get_coord((*i).get_vertex_no(2)); 
	if(area2==0) cout << ind << "area2 v1-v2 " << (v1-v2).norm() << " (v3-v2).norm() " << (v3-v2).norm() << "v1-v3 " << (v1-v3).norm() << " aera1 " << area1 << " area2 " << area2 << endl;

	area2=computeArea(v1, v2, v3);
	
	distortion=log2(area2/area1); // areal distortion measure used in caret
	//if(distortion!=distortion) cout << area1 << " inf " << area2 << " distortion " << distortion << endl;
	meandistortion((*i).get_vertex_no(0)+1)+=0.3*area1*distortion; 
	meandistortion((*i).get_vertex_no(1)+1)+=0.3*area1*distortion; 
	meandistortion((*i).get_vertex_no(2)+1)+=0.3*area1*distortion; 
	weights((*i).get_vertex_no(0)+1)+=0.3*area1;
	weights((*i).get_vertex_no(1)+1)+=0.3*area1;
	weights((*i).get_vertex_no(2)+1)+=0.3*area1;
	//cout << area1 << " " << area2 << " distortion " << distortion <<   endl;
	faces++;
    }
    
    for (int i=1;i<=OUT.nvertices();i++){
		if(weights(i)>0)
			meandistortion(i)=meandistortion(i)/weights(i);
    }
   
    return meandistortion;
  }

  ColumnVector getarealseparation(const NEWMESH::newmesh &ORIG,  NEWMESH::newmesh &OUT){
    
    double distortion=0;
    double area1=0, area2=0;
    ColumnVector meandistortion(OUT.nvertices());
    int faces=0;
    meandistortion=0;
    if(ORIG.nvertices()!=OUT.nvertices()){ throw  NEWMESHException("NEWMESH ERROR in meshfns:: getarealdistortion initial and distorted surfaces should have them same number of vertices");}
    
    for ( vector<NEWMESH::Triangle>::const_iterator i=ORIG.tbegin() ; i!=ORIG.tend(); i++)
      {
	NEWMESH::Pt v1 = (*i).get_vertex_coord(0),  v2 = (*i).get_vertex_coord(1), v3 = (*i).get_vertex_coord(2);
	  
	area1=computeArea(v1, v2, v3);
	v1=OUT.get_coord((*i).get_vertex_no(0)); v2=OUT.get_coord((*i).get_vertex_no(1)); v3=OUT.get_coord((*i).get_vertex_no(2)); 
	area2=computeArea(v1, v2, v3);
	distortion=((area2-area1)/area1)*((area2-area1)/area1); // areal distortion measure used in caret
	
	meandistortion((*i).get_vertex_no(0)+1)+=distortion; 
	meandistortion((*i).get_vertex_no(1)+1)+=distortion; 
	meandistortion((*i).get_vertex_no(2)+1)+=distortion; 
	
	faces++;
      }
    
    for (int i=1;i<=OUT.nvertices();i++){
      meandistortion(i)=meandistortion(i)/OUT.get_total_triangles(i-1);
    }
   
    return meandistortion;
  }

  ColumnVector getarealdistortionFACES(const NEWMESH::newmesh &ORIG,  NEWMESH::newmesh &OUT){
    
    double area1=0, area2=0;
    ColumnVector distortion(OUT.ntriangles());
  
    if(ORIG.nvertices()!=OUT.nvertices()){ throw  NEWMESHException("NEWMESH ERROR in meshfns:: getarealdistortion initial and distorted surfaces should have them same number of vertices");}
    
    for ( vector<NEWMESH::Triangle>::const_iterator i=ORIG.tbegin() ; i!=ORIG.tend(); i++)
      {
	NEWMESH::Pt v1 = (*i).get_vertex_coord(0),  v2 = (*i).get_vertex_coord(1), v3 = (*i).get_vertex_coord(2);
	  
	area1=computeArea(v1, v2, v3);
	v1=OUT.get_coord((*i).get_vertex_no(0)); v2=OUT.get_coord((*i).get_vertex_no(1)); v3=OUT.get_coord((*i).get_vertex_no(2)); 
	area2=computeArea(v1, v2, v3);
	distortion(i->get_no()+1)=log2(area2/area1); // areal distortion measure used in caret
	
      }
   
    return distortion;
  }

  Matrix estimate_rotation_matrix(const Pt &p1,const Pt &p2) // ci is start point index is end point
  {
    Matrix I(3,3),u(3,3),R(3,3);
    double theta;
    Pt ci=p1;
    Pt index=p2;
    ci.normalize(); index.normalize();
    
    theta=acos(ci|index);

  

    if(theta > PI ) {
      cout << " rotation angle is greater than 90 degrees " << endl;
      exit(1);
    }
    I << 1 << 0 << 0
      << 0 <<  1 << 0
      << 0 <<  0 << 1;
    
    //calculate axis of rotation from cross product between centre point and neighbour - could also use barycentres
    Pt cross=ci*index; cross.normalize();
  
  
    if(fabs(1-(ci|index))<EPSILON){
      R=I;
      
    }
    else if(cross.norm()==0) {
      R<< -1 << 0 << 0
       << 0 <<  -1 << 0
       << 0 <<  0 << -1;

    }
    else{
      
      u << 0 << -cross.Z << cross.Y
	<< cross.Z <<  0 << -cross.X
	<< -cross.Y <<  cross.X << 0;

      //cout << u << endl;
      //// now use rodriguez formula to rotate all data points that contribute the similarity of this control point
      if(fabs(-1-(ci|index))< EPSILON){// numerical problems at theta cose to PI
	Matrix outerprod(3,3); 

	outerprod << cross.X*cross.X << cross.X*cross.Y << cross.X*cross.Z
	  << cross.Y*cross.X <<  cross.Y*cross.Y << cross.Y*cross.Z
	  << cross.Z*cross.X <<  cross.Z*cross.Y << cross.Z*cross.Z;

	R=2*outerprod-I;
	
	if(R.Trace()!=R.Trace()){ cout << " R=2*outerprod-I " << " (ci|index) " << (ci|index) << endl; cout << outerprod << endl;}
      }else{
      

	R=I+u*sin(theta)+(1-cos(theta))*(u*u);
      }

      if(R.Trace()!=R.Trace()){ cout << " theta " << theta << " (ci|index) " << (ci|index) << " ci.X " << ci.X << " ci.Y " << ci.Y << " ci.Z " << ci.Z << " index.X " << index.X << " index.Y " << index.Y << " index.Z " << index.Z << endl; }

    }
  
    
    return R;
  }

  Matrix rodriguez_rot(const double &theta,const Matrix &AXIS, const int & index){
    Matrix R(3,3),I(3,3);

    I << 1 << 0 << 0
      << 0 <<  1 << 0
      << 0 <<  0 << 1;
    
    switch (index) {
    case 1:
      R=I;
      break;
    case 2:
      R << -1 << 0 << 0
       << 0 <<  -1 << 0
       << 0 <<  0 << -1;
      break;
    case 3:
	R=I+AXIS*sin(theta)+(1-cos(theta))*(AXIS*AXIS);
      break;
    case 4:

	R=2*AXIS-I;
	break;
    default:
      cout << " rodriguez_rot:: case is not found " << endl;
      exit(1);
    }
	
    if(R.Trace()!=R.Trace()){ cout << " R.Trace()!=R.Trace())" << endl; exit(1);}

    return R;
  }

  Matrix estimate_axis_of_rotation(const Pt &p1,const Pt &p2, int & result) // ci is start point index is end point
  {
    Matrix I(3,3),u(3,3); 
    Pt ci=p1;
    Pt index=p2;
    ci.normalize(); index.normalize();
    result=1;
    //calculate axis of rotation from cross product between centre point and neighbour - could also use barycentres
    Pt cross=ci*index; cross.normalize();


    if((fabs(1-(ci|index)) < EPSILON)){
      result=1;
    }else if(cross.norm()==0){
      result=2;
    }else{
      
      if(fabs(-1-(ci|index))< EPSILON){
	u << cross.X*cross.X << cross.X*cross.Y << cross.X*cross.Z
	  << cross.Y*cross.X <<  cross.Y*cross.Y << cross.Y*cross.Z
	  << cross.Z*cross.X <<  cross.Z*cross.Y << cross.Z*cross.Z;
	result=4;
	cout << " angle approx PI" << endl; exit(1); 
      }else{ u << 0 << -cross.Z << cross.Y
	<< cross.Z <<  0 << -cross.X
	<< -cross.Y <<  cross.X << 0;
	result=3;
      }
      //// now use rodriguez formula to rotate all data points that contribute the similarity of this control point
    }
  
    
    return u;
  }


  bool get_all_neighbours(const int &index, vector<int> & N, const  NEWMESH::Pt & point, const int &n, const NEWMESH::newmesh &REF, const RELATIONS& _rel, SpMat<int> & found) {
  
    int n0,n1,n2;
    bool update=false;
    
    
    for (vector<int>::const_iterator j=REF.tIDbegin(n); j!=REF.tIDend(n); j++){
      n0 = REF.get_triangle_vertexID(*j,0); n1 = REF.get_triangle_vertexID(*j,1); n2=  REF.get_triangle_vertexID(*j,2);
      
      if(_rel.is_neighbour(n0+1,index+1)==0 || _rel.is_neighbour(n1+1,index+1)==0 || _rel.is_neighbour(n2+1,index+1)==0){update=true;}
      
      if(found.Peek(n0+1,1)==0){ N.push_back(n0); found.Set(n0+1,1,1); }
      
      if(found.Peek(n1+1,1)==0){ N.push_back(n1); found.Set(n1+1,1,1);}
      
      if(found.Peek(n2+1,1)==0){ N.push_back(n2); found.Set(n2+1,1,1);} 
    }
   	   
    return update;
      
  }

  NEWMESH::newmesh  binarise_cfweighting(const NEWMESH::newmesh &CFWEIGHTING ){
    Matrix DATA=CFWEIGHTING.get_pvalues();
    newmesh exclusion=CFWEIGHTING;
    Matrix newdata(1,CFWEIGHTING.nvertices()); newdata=0;
    for(int i=1;i<=CFWEIGHTING.nvertices();i++){
      for(int d=1;d<=DATA.Nrows();d++){
	if(DATA(d,i)!=0)
	  newdata(1,i)=1;
      }
    }
     exclusion.set_pvalues(newdata);
     return exclusion;
  }

 
  NEWMESH::newmesh create_exclusion(NEWMESH::newmesh &IN,const Matrix & DATA,const float &thrl,const float &thru){
    Matrix cfweighting(1,IN.npvalues());
    NEWMESH::newmesh EXCL=IN;
    int flag;
    cfweighting=1;
    
    for (int i=1;i<=IN.npvalues();i++){
      flag=0;
      for (int j=1;j<=DATA.Nrows();j++){
	if(DATA(j,i) >= thrl-EPSILON && DATA(j,i) <= thru+EPSILON){ /// only exclude if cfweighting for all feat dimensions is zero
	  flag=1;
	 
	}else
	  flag=0;
	   
      }
      if(flag==1) cfweighting(1,i)=0;
    }

    EXCL.set_pvalues(cfweighting);

    return EXCL;
  }

  vector<string>  read_ascii_list(string filename)
  {
    vector<string > list;
    ifstream fs(filename.c_str());
    string tmp;
    if(fs){
      fs>>tmp;
      do{
	list.push_back(tmp);

	fs>>tmp;
      }while(!fs.eof());
    }
    else{
      cerr<<filename<<" does not exist"<<endl;
      exit(0);
    }
  
    return list;
  }

 
  
}



  
