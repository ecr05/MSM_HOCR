/*  ContinuosOpt.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  

    Some sections of code inspired by Alek Petrovic.*/
/*  CCOPYRIGHT  */
#include "ContinuosOpt.h"

namespace MESHREG {
  
  void MeshCF::Initialize(){

    //////// temporary REGISTRATION PARAMETERS ///////////////  
    current_sim.ReSize(_SOURCE.nvertices());
   
  
    MVD=Calculate_MVD(_SOURCE);
   
    min_sigma=MVD;
    _rel=boost::shared_ptr<RELATIONS >(new RELATIONS (_SOURCE,_TARGET,2*asin(4*MVD/(2*RAD))));  
    _rel->estimate_adj(_TARGET.nvertices(),_SOURCE.nvertices());
    sim.Resize(_TARGET.nvertices(),_SOURCE.nvertices());  /// creates sparse sim kernel to fill with similarities of nearest neighbours only
   	
    sim.set_reference(FEAT->get_reference_data()); sim.set_input(FEAT->get_input_data()); 
    sim.Initialize(_simmeasure);
    sim.set_relations(_rel);
    update_similarity();
    
  }
   
  //// weighted least squares similarity gradient estimation
  ColumnVector MeshCF::WLS_simgradient(const Tangs  &T, int index, const vector<int> &querypoints)
  {
    
    NEWMESH::Pt point, v0, v1,v2,cr;
    ColumnVector grad,gradtmp;
    SpMat<int> found(_TARGET.nvertices(),1);
    Matrix Disp,Dispinv,Data;
    Matrix weights;
    double x11,x21;
    double y11,y21;
    double SUM=0.0;
    double JPsim=0.0; 
    ColumnVector dist(2);
    Pt origin;

   
    gradtmp.ReSize(2); gradtmp=0.0;
    
    Disp.ReSize(querypoints.size(),2);
    weights.ReSize(querypoints.size(),querypoints.size());
    Data.ReSize(querypoints.size(),1);
    Disp=0; weights=0;
    point=_SOURCE.get_coord(index);


    origin=T.e1*T.e2;
    origin.normalize();
    origin=origin*RAD;

    projectPoint(point-origin,T,y11,y21);  
    //#pragma omp critical
  
    for (unsigned int i=0;i<querypoints.size();i++){
      cr=_TARGET.get_coord(querypoints[i]);
      projectPoint(cr-origin,T,x11,x21);  
      dist(1)=x11-y11;dist(2)=x21-y21; 
    
        
      if((dist(1)*dist(1)+dist(2)*dist(2))>0){
       	
	
	//	double si=2*RAD*asin(dist.norm()/(2*RAD))/dist.norm();  // geodesic distance - use this instead?
	Disp(i+1,1)=x11; Disp(i+1,2)=x21; 

	//////////// THE DISTANCES SHOULD BE CALCULATED FOR THE TANGENT PLANE///////////// - we ideally want an analytica expression

	weights(i+1,i+1)=_refweight(1,querypoints[i]+1)*exp(-(dist(1)*dist(1)+dist(2)*dist(2))/(2*min_sigma*min_sigma)); //
  
	SUM+=weights(i+1,i+1); 
	JPsim+=sim.Peek(querypoints[i]+1,index+1)*weights(i+1,i+1);

      }

    }
 

    if(SUM>0)  
      JPsim/=SUM;
    
    
    for (unsigned int i=0;i<querypoints.size();i++){
      Data(i+1,1)=sim.Peek(querypoints[i]+1,index+1)-JPsim; // difference of query point to mean similarity
   
    }

    //estimates fits a gradient to the set samples and their similarity values using WLS
    Matrix tmp=((Disp.t()*weights)*Disp);
    if(tmp.Determinant()>1e-5){
      Dispinv=(tmp.i()*Disp.t())*weights; 
      gradtmp=Dispinv*Data; 
    }
    
    
    grad.ReSize(3); grad=0.0;
    grad(1) += (gradtmp(1)*(T.e1).X + gradtmp(2)*(T.e2).X) ;
    grad(2) += (gradtmp(1)*(T.e1).Y + gradtmp(2)*(T.e2).Y) ;
    grad(3) += (gradtmp(1)*(T.e1).Z + gradtmp(2)*(T.e2).Z) ;
      
    
    current_sim(index+1)=JPsim;

    return grad;
    
    }
  
 
  ////////////////////////////////// UPDATES ////////////////////////////////////////
    
  void MeshCF::update_similarity()
  {
      
    

    for (int i=1;i<=_SOURCE.nvertices();i++)
      _rel->update_RELATIONS_for_ind(i,_SOURCE); 
        
    sim.set_input(FEAT->get_input_data()); 
    sim.set_reference(FEAT->get_reference_data());
    

    for(int i=1;i<=_SOURCE.nvertices();i++){
      sim.calculate_sim_column(i); 
    }
    
  }
    
  void MeshCF::update_similarity(const int &i,vector<int> & points)
  {
    
    _rel->update_w_querypoints_for_ind(i+1,points);   
    sim.calculate_sim_column(i+1); 
      
  }
  
  ColumnVector MeshCF::Evaluate_SIMGradient(int i,const Tangs &T){
    
    SpMat<int> found(_TARGET.nvertices(),1);
    int n0,n1,n2,t_ind;
    NEWMESH::Pt point, v0, v1,v2;
    vector<int> querypoints;
    bool update=false;
    ColumnVector vecnew;
    
    point=_SOURCE.get_coord(i);
    
    if(_rel->Nrows(i+1)){
      t_ind=(int) (*_rel)(1,i+1);
    
      vecnew.ReSize(3); vecnew=0;
      do{ 
	querypoints.clear();
	if(_TARGET.return_closest_points(t_ind-1,v0,v1,v2,point,n0,n1,n2)){ // finds nearest point in target mesh
	  /// obtains all the vertices connected to this point
	  if((_refweight(1,n0+1) > 0) && (_refweight(1,n1+1) > 0) && ( _refweight(1,n2+1) > 0) ){
	    if(get_all_neighbours(i,querypoints, point,n0,_TARGET,*_rel,found))  update=true;
	    if(get_all_neighbours(i,querypoints, point,n1,_TARGET,*_rel,found) || update==true ) update=true;
	    if(get_all_neighbours(i,querypoints, point,n2,_TARGET,*_rel,found) || update==true ) update=true;
	    
	    if(update){
	      update_similarity(i,querypoints) ;
	    }
	    /// calculates similarity gradient
	    vecnew=WLS_simgradient(T,i,querypoints);
	    vecnew=_inweight(1,i+1)*_refweight(1, t_ind)*vecnew;
	    update=false;
	    
	  }else{ update=false; current_sim(i+1)=0;}
	}else{update=true; _rel->find_next_closest_neighbour(i+1,t_ind,point,_TARGET);}
      }while(update==true);
    } 
   
    return vecnew;
  }


   /////////////// AFFINE REGISTRATION FUNCTIONS - from Alek Petrovic //  
  void affineMeshCF::Rotate_IN_mesh(const double &a1, const double &a2,const double & a3){
    
   
    for (int index=0; index < _SOURCE.nvertices(); index++){
      
      NEWMESH::Pt cii = _SOURCE.get_coord(index);
      NEWMAT::ColumnVector V(3), VR(3);
      V(1) = cii.X; V(2) = cii.Y; V(3) = cii.Z;
      VR = rotate_euler(V,a1,a2,a3); // in meshregfns
      
      NEWMESH::Pt ppc(VR(1), VR(2), VR(3));
      _SOURCE.set_coord(index,ppc);
    }
    
  }
  
  double affineMeshCF::Affine_cost_mesh(const double & dw1, const double & dw2,const double &dw3){
    
    double SUM = 0;
    NEWMESH::newmesh tmp=_SOURCE;
    Tangent Tang;  // holds tangent plane basis at each point

    Rotate_IN_mesh(dw1,dw2,dw3);
    
    for (int index=0; index < _SOURCE.nvertices(); index++){
	if(_rel->Nrows(index+1)){
	  Tangs T = Tang.calculate(index,_SOURCE);
	  Evaluate_SIMGradient(index,T); /// only need mean similarity from this most of this code is now redundant!
	  SUM+=current_sim(index+1);
	}
    }
    
     _SOURCE=tmp;
    return SUM;
  }
  
  NEWMESH::newmesh affineMeshCF::run(){
    if(_verbosity)
    cout<<"Affine registration started"<<endl;
    double Euler1,Euler2,Euler3;
    double mingrad_zero;
    NEWMAT::Matrix Rot(3,3), Rot_f(3,3);
    Rot_f = 0; Rot_f(1,1) = 1; Rot_f(2,2) = 1; Rot_f(3,3) = 1;
    Euler1 = 0; Euler2 = 0; Euler3 = 0;
    double per ; //1000 - for the artificial ex.
    ColumnVector E(3);
    double RECinit, RECfinal;
    RECfinal=0;
    double grad_zero = Affine_cost_mesh(Euler1, Euler2, Euler3); /// tries different small rotations and gets mesh similarity for each of these points
    int min_iter=0;
   

    mingrad_zero=grad_zero; RECinit=grad_zero;
    int loop=0;
    newmesh tmp;
   while (_spacing > 0.05){ 
     double  step=_stepsize;

     per =_spacing; 
    for (int _iter=1; _iter <= _iters; _iter++){
      //Initially Euler angles are zero
      Euler1 = 0; Euler2 = 0; Euler3 = 0;
           

      NEWMESH::Pt grad;
      
      grad.X = (Affine_cost_mesh(Euler1+per, Euler2, Euler3)-grad_zero)/per;
      
      grad.Y = (Affine_cost_mesh(Euler1, Euler2+per, Euler3)-grad_zero)/per;
      
      grad.Z = (Affine_cost_mesh(Euler1, Euler2, Euler3+per)-grad_zero)/per;
	
      grad.normalize();
      Euler1 = Euler1 + step*grad.X;
      Euler2 = Euler2 + step*grad.Y;
      Euler3 = Euler3 + step*grad.Z;
      if(_verbosity){  
	cout<<"******:   "<< _iter << " per " << per << " loop" << loop << "  (loop*_iters)+ _iter " << (loop*_iters)+ _iter <<  endl;
	cout<<"grad_zero:   "<<grad_zero<<endl;
	cout<<"mingrad_zero:   "<<mingrad_zero<< "min_iter " <<  min_iter << endl;
	cout<<"stepsize:   "<<step<<endl;

      } 
      tmp=_SOURCE;
      Rotate_IN_mesh(Euler1,Euler2,Euler3);  
      grad_zero = Affine_cost_mesh(Euler1, Euler2, Euler3);
          
      if(grad_zero > mingrad_zero){ mingrad_zero = grad_zero; min_iter=(loop*_iters)+_iter; RECfinal=mingrad_zero;}

      if((loop*_iters)+_iter -min_iter>0){ step*=0.5; _SOURCE=tmp;}
	  if(step<1e-3) break;
	 
    }
    loop++;
    _spacing*=0.5;
    
}

   
    E(1) = Euler1;
    E(2) = Euler2;
    E(3) = Euler3;

    if(_verbosity && (RECfinal >0) )
      cout<<"Affine improvement: "<< abs(((RECfinal-RECinit))/RECfinal*100)<<"%"<<  endl;
    
   
    return _SOURCE;
   }


}

