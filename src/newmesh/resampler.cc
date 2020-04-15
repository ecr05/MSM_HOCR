/*  resampler.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford

    Some sections of code inspired by A. Petrovic.
*/
/*  CCOPYRIGHT  */
/*  Adaptive Barycentric code used with the permission of Tim Coulson under the below licence:
 *  Copyright (C) 2014  Washington University School of Medicine
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included
 *  in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 *  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 *  CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 *  TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.*/

#include "resampler.h"

namespace NEWMESH {

  enum StringValue { evNotDefined,
		     evStringValue1,
		     evStringValue2,
		     evStringValue3,
		     evStringValue4,
		     evEnd };
  // Map to associate the strings with the enum values
  static std::map<std::string, StringValue> s_mapStringValues;


  //thin-plate spline
  double phi(double r){

    double u = 0;
    if (r != 0) u = r*r*log(r);
    return u;
  }


  //gaussian
  double phi_g(double r, double c){

    double u = 0;
    if (r != 0) u = exp(-c*r*r); else u = 1;
    return u;
  }



  void RBF::set_point(double X, double Y, double Z, double fn){

    int NR = cent.Nrows();

    cent.Release();
    Matrix cent_temp = cent;
    cent.ReSize(NR+1,3);

    f.Release();
    ColumnVector f_temp = f;
    f.ReSize(NR+1);

    if (NR > 0)
      cent.Rows(1,NR) << cent_temp;

    cent(NR+1,1) = X;
    cent(NR+1,2) = Y;
    cent(NR+1,3) = Z;

    if (NR > 0)
      f.Rows(1,NR) = f_temp;
    f(NR+1) = fn;

  }


  ReturnMatrix RBF::direct(Matrix& cent){

    // Input
    //         cent  n*dim matrix  coordinates of centers
    //               of reals
    // Output  AI     (n+dim+1)*   Symmetric matrix of the linear
    //               (n+dim+1)    system obtained if we solve
    //               array of     for the radial basis function
    //                            interpolant directly.
    //
    // Write the matrix A in the form
    //             A    P
    //     AI  =
    //             P^t  O
    // where P is the polynomial bit.
    int n = cent.Nrows();
    Matrix A(n,n);
    A = 0;

    for (int i=1; i<=n; i++){
      for (int j=1; j<=i; j++){

	NEWMESH::Pt Vi(cent(i,1),cent(i,2),cent(i,3));
	NEWMESH::Pt Vj(cent(j,1),cent(j,2),cent(j,3));
	double r = (Vi - Vj).norm();
	double temp = phi(r);

	A(i,j) = temp;
	A(j,i) = temp;
      }
    }
    Matrix P(n, 4);
    P = 1;
    P.Columns(2,4) = cent;

    Matrix AI(n+4,n+4);
    AI = 0;
    AI.SubMatrix(1,n,1,n) = A;
    AI.SubMatrix(n+1,n+4,1,n) = P.t();
    AI.SubMatrix(1,n,n+1,n+4) = P;

    AI.Release();
    return AI;
  }

  ReturnMatrix RBF::fitit(Matrix& cent, ColumnVector& f){

    // Input
    //    n       is the number of points participating in the interpolation
    //    cent    n by dim array of centres
    //    f       n by 1 vector of values at centres
    //
    // Output
    //    coeff   (n + 3) by 1 vector of coefficients with the
    //            coefficients of 1, x and y (2D) or 1, x, y
    //            and z (3D) last.
    int n = cent.Nrows();
    Matrix A(n+4,n+4);

    A = direct(cent);

    f.Release();
    ColumnVector X=f;
    f.ReSize(n+4);
    f = 0;
    f.Rows(1,n)=X.Rows(1,n);

    ColumnVector coeff;

    if (A.Determinant() < 1e-5) {
      cout<<"Matrix is singular 2 "<< A.Determinant() << " " << A.Nrows() << endl;
      exit(1);
    }

    coeff = A.i()*f;

    coeff.Release();
    return coeff;
  }

  double RBF::eval_direct(Matrix& cent, ColumnVector coeff, double uX, double uY, double uZ){

    // Inputs
    //    cent   n by dim matrix of coordinates for the centres
    //    coeff  n+dim +1 vector of coefficients with the linear
    //           polynomial part last.
    //    u      row vector    point at which to evaluate
    //
    // Output
    //    s    Value of the RBF at position u

    int n = cent.Nrows();
    double s = 0;
    NEWMESH::Pt U(uX, uY, uZ);

    double MAX_C = 1;
    MAX_C = 1;

    for (int i=1; i<=n; i++){

      NEWMESH::Pt Vi(cent(i,1),cent(i,2),cent(i,3));
      s += coeff(i)*phi((U - Vi).norm());
     }

    s += coeff(n+1);
    s += coeff(n+2)*uX + coeff(n+3)*uY + coeff(n+4)*uZ;
    return s;
  }


  double RBF::interpolate(double X, double Y, double Z, bool estimate_coeff){

    if(estimate_coeff) coeff = fitit(cent, f);
    return eval_direct(cent, coeff, X, Y, Z);

  }

  ////////////////////////////////// RESAMPLER METHODS /////////////////////////////////////
  void resampler:: set_method(const string &M){
    s_mapStringValues["NN"] = evStringValue1;
    s_mapStringValues["LINEAR"] = evStringValue2;
    s_mapStringValues["GAUSSIAN"] = evStringValue3;
    s_mapStringValues["ADAP_BARY"] = evStringValue4;
    s_mapStringValues["BARYCENTRIC"] = evEnd;
    _method=M;
  }


  RELATIONS resampler::initialize_relations(const double& ang, NEWMESH::newmesh & in,const  NEWMESH::newmesh & ref){
    RELATIONS REL(in,ref,ang);
    REL.update_RELATIONS(in);
    return REL;
  }

  double resampler::calc_weight(const double &dist , const double &kernel_width){

    double weight=1.0;
    switch(s_mapStringValues[_method])
      {
      case evStringValue2:
	if(dist<=kernel_width)
	  weight=1-dist/kernel_width;
        break;
      case evStringValue3:
	weight=(1/sqrt(2*PI*kernel_width*kernel_width))*exp(-(dist*dist)/(2*kernel_width*kernel_width));  ///
        break;
      default:
        cout << "'" <<_method
	     << "' is an invalid string. s_mapStringValues now contains "
             << s_mapStringValues.size()
             << " entries." << endl;
        break;
      }
    return weight;
  }

  double resampler::guess_angular_spacing(const int &minvertexnum){
    double MVD;
    NEWMESH::newmesh icotmp;
    int ico=1;
    if(minvertexnum>= 42)
      ico=1;
    if(minvertexnum>= 162)
      ico=2;
    if(minvertexnum>= 642)
      ico=3;
    if(minvertexnum>=2562)
      ico=4;
    if(minvertexnum>= 10242)
      ico=5;
    if(minvertexnum>=40962)
      ico=6;

    icotmp.make_mesh_from_icosa(ico); true_rescale(icotmp,RAD);
    MVD=Calculate_MVD(icotmp);

    return 2*asin(MVD/(RAD));  // make angular spacing twice that of a regular grid of similar size
  }

  ///// generic resample function: note EXCL is an exclusion masking or confidence weighting function (though is only used to mask here), in is the to be resampled mesh, ref is the reference surface, data should have as many columns as input vertices (check here??) and sigma is the size of the smoothing kernel
  void resampler::resampledata(NEWMESH::newmesh& in,const NEWMESH::newmesh& ref, boost::shared_ptr<NEWMESH::newmesh> &EXCL, boost::shared_ptr<BFMatrix> &data, const double &sigma, boost::shared_ptr<RELATIONS > _rel){
    if(EXCL.get()!=0 && EXCL->npvalues()!=in.nvertices()){ throw  NEWMESHException("RESAMPLER:: exclusion mask does not have the same number of vertices as the surface mesh.");}

    if(_method=="BARYCENTRIC" || _method=="ADAP_BARY"){
      barycentric_interpolation(in,ref,data,EXCL,_rel);
    }
    else{
      if(((_rel.get()) &&_rel->Ncols()==in.nvertices()) || ((!_rel.get()) && (in.nvertices()>ref.nvertices()))){
      	downsample_w_interpolation(in,ref,sigma,data,EXCL,_rel);
      }
      else{
	upsample_w_interpolation(in,ref,sigma,data,EXCL,_rel);
      }
    }
  }

  void resampler::resampledata_to_coords(NEWMESH::newmesh& in,boost::shared_ptr<BFMatrix> &targetcoords, boost::shared_ptr<BFMatrix> &data, boost::shared_ptr<RELATIONS > _rel){

    //cout << "resampledata_to_coords - barycentric only" << endl;
    barycentric_interpolation_to_coords(in,targetcoords,data,_rel);



  }

  void resampler::resampledata(NEWMESH::newmesh& in,const NEWMESH::newmesh& ref, boost::shared_ptr<BFMatrix> &data, const double &sigma, boost::shared_ptr<RELATIONS > _rel){
    boost::shared_ptr<NEWMESH::newmesh> EXCL;
    if(_method=="BARYCENTRIC" || _method=="ADAP_BARY"){
      barycentric_interpolation(in,ref,data,EXCL,_rel);
    }
    else{
      if(((_rel.get()) &&_rel->Ncols()==in.nvertices()) || ((!_rel.get()) && (in.nvertices()>ref.nvertices())))
	downsample_w_interpolation(in,ref,sigma,data,EXCL,_rel);
      else
	upsample_w_interpolation(in,ref,sigma,data,EXCL,_rel);
    }
  }



  void resampler::nearest_neighbour_interpolation(NEWMESH::newmesh &in,const  NEWMESH::newmesh &SPH, boost::shared_ptr<BFMatrix> &data,  boost::shared_ptr<NEWMESH::newmesh> &EXCL, boost::shared_ptr<RELATIONS > _rel){
    Matrix newdata(data->Nrows(),SPH.nvertices()); newdata=0;
    newmesh exclusion=SPH;
    check_scale(in,SPH);
    if(!_rel.get()||_rel->Ncols()!=SPH.nvertices()){
    double ang=2*asin(Calculate_MaxVD(SPH)/(2*RAD));
    _rel = boost::shared_ptr<RELATIONS > ( new RELATIONS(exclusion,in,ang));
    _rel->update_RELATIONS(exclusion);
    }


    for (int i = 1; i <=  SPH.nvertices(); i++){
      exclusion.set_pvalue(i-1,0);
      if(_rel->Nrows(i)>0){

	if(EXCL.get()==0 || EXCL->get_pvalue((*_rel)(1,i)-1)){

	    if(EXCL.get()) exclusion.set_pvalue(i-1,EXCL->get_pvalue((*_rel)(1,i)-1));
	    for (BFMatrixColumnIterator it = data->begin((*_rel)(1,i)); it != data->end((*_rel)(1,i)); ++it) {
	      newdata(it.Row() ,i)=*it;
	  }

	}
      }
    }
    if(EXCL.get())*EXCL=exclusion;
    boost::shared_ptr<FullBFMatrix > pin =boost::dynamic_pointer_cast<FullBFMatrix>(data);
    if(pin.get()){      data = boost::shared_ptr<BFMatrix >(new FullBFMatrix (newdata));}
    else {data = boost::shared_ptr<BFMatrix >(new  SparseBFMatrix<double>  (newdata));}

  }

  /*  ADAPTIVE BARYCENTRIC CODE SUPPLIED BY TIM COULSON AT WASHU Copyright (C) 2014  Washington University School of Medicine */
  /* ORIGINAL CODE: https://github.com/Washington-University/workbench/blob/master/src/Files/SurfaceResamplingHelper.cxx */

  void resampler::barycentric_interpolation(NEWMESH::newmesh &IN, const NEWMESH::newmesh  &ref, boost::shared_ptr<BFMatrix> &data, boost::shared_ptr<NEWMESH::newmesh> &EXCL, boost::shared_ptr<RELATIONS > _rel)
  {
    double val,excl_val;
    NEWMESH::Pt ci;
    Matrix newdata(data->Nrows(),ref.nvertices());
    double excl_sum;
    NEWMESH::newmesh exclusion=ref;
    vector<map<int,double> > adapt=get_all_adap_barycentric_weights(IN,exclusion,EXCL,_rel);

      //// take barycentric weights and use these to resample whilst excluding 0 mask values

    for (int k=0;k<(int) ref.nvertices();k++){
      exclusion.set_pvalue(k,0);
      for (int i=1;i<=(int) data->Nrows();i++){
      	val=0; excl_val=0; excl_sum=0;
      	for (std::map<int,double>::iterator it=adapt[k].begin(); it!=adapt[k].end(); ++it){
      	   if(EXCL.get()==0 || EXCL->get_pvalue(it->first)){
      	      val+=it->second*data->Peek(i,it->first+1);


	         if(i==1 && EXCL.get()) {
	            excl_val+=it->second*EXCL->get_pvalue(it->first);
	         }else{val+=0;
	            if(i==1) excl_val+=0;
	    	    }
	         }
	        }

	      newdata(i,k+1)=val;
	      if(EXCL.get()&& i==1){

	    exclusion.set_pvalue(k,excl_val);

	}
	else if(i==1)
	  exclusion.set_pvalue(k,1);

	if(exclusion.get_pvalue(k)!=exclusion.get_pvalue(k))
	  cout <<" exclusion.get_pvalue(k)!=exclusion.get_pvalue(k) " << endl;
      }
    }

    if(EXCL.get())*EXCL=exclusion;
    boost::shared_ptr<FullBFMatrix > pin =boost::dynamic_pointer_cast<FullBFMatrix>(data);

    if(pin.get()){data = boost::shared_ptr<BFMatrix >(new FullBFMatrix (newdata));}
    else {data = boost::shared_ptr<BFMatrix >(new  SparseBFMatrix<double>  (newdata));}

  }

  void resampler::barycentric_interpolation_to_coords(NEWMESH::newmesh &IN, const boost::shared_ptr<BFMatrix>  &targetcoords, boost::shared_ptr<BFMatrix> &data,  boost::shared_ptr<RELATIONS > _rel)
  {
    double val;
    Matrix newdata(targetcoords->Nrows(),data->Nrows());
    vector<map<int,double> >  bary=get_all_barycentric_weights(IN,targetcoords,_rel);
    //cout << " barycentric_interpolation_to_coords " << targetcoords->Nrows() << " " << targetcoords->Ncols() << " " << newdata.Nrows() << newdata.Ncols () <<  endl;
    //// take barycentric weights and use these to resample whilst excluding 0 mask values

    for (int k=0;k<(int) targetcoords->Nrows();k++){

      for (int i=1;i<=(int) data->Nrows();i++){
        val=0;
        for (std::map<int,double>::iterator it=bary[k].begin(); it!=bary[k].end(); ++it)
            val+=it->second*data->Peek(i,it->first+1);

        newdata(k+1,i)=val;
      }
    }
    boost::shared_ptr<FullBFMatrix > pin =boost::dynamic_pointer_cast<FullBFMatrix>(data);

    if(pin.get()){data = boost::shared_ptr<BFMatrix >(new FullBFMatrix (newdata));}
    else {data = boost::shared_ptr<BFMatrix >(new  SparseBFMatrix<double>  (newdata));}

  }

  //// describes all vertices in mesh A in terms of barycentric coordinates for surrounding triangles defined in MESHB
  vector<std::map<int,double> > resampler::get_all_barycentric_weights(NEWMESH::newmesh &MESHA, NEWMESH::newmesh &MESHB, boost::shared_ptr<RELATIONS > _rel){
      vector<map<int,double> > weights;
      map<int,double> indexweight;
      double ang=0;

      if(!_rel.get() ||_rel->Ncols()!=MESHA.nvertices())
        {
  	       if(MESHB.nvertices() > MESHA.nvertices()){ ang=guess_angular_spacing(MESHA.nvertices());}
  	        else{ ang=2*asin(Calculate_MaxVD(MESHB)/(RAD));}

  	        _rel = boost::shared_ptr<RELATIONS > ( new RELATIONS(MESHA,MESHB,ang));

  	        _rel->update_RELATIONS(MESHA);
        }

      for (int k=0;k<(int) MESHA.nvertices();k++){
          indexweight=get_barycentric_weight_for_ind(k,ang,MESHA,MESHB,*_rel);
          weights.push_back(indexweight);
        }
    return weights;
  }

  //// describes all vertices in mesh A in terms of barycentric coordinates for surrounding triangles defined in MESHB
  vector<std::map<int,double> > resampler::get_all_barycentric_weights(NEWMESH::newmesh &in, const boost::shared_ptr<BFMatrix>  &targetcoords, boost::shared_ptr<RELATIONS > _rel){
      vector<map<int,double> > weights;
      map<int,double> indexweight;
      double ang=0;
      //cout << " get_all_barycentric_weights " << endl;

      if(!_rel.get() ||_rel->Ncols()!=in.nvertices())
        {
           ang=guess_angular_spacing(in.nvertices());
           _rel = boost::shared_ptr<RELATIONS > ( new RELATIONS(in,in,ang));
           //cout << " initialise rel " << _rel->Ncols() << " " << in.nvertices() << " " << targetcoords->Nrows() << endl;
        }

      for (int k=1;k<=(int) targetcoords->Nrows();k++){
          Pt ci(targetcoords->Peek(k,1),targetcoords->Peek(k,2),targetcoords->Peek(k,3));
          //cout << k << " " << ci.X << " " << ci.Y << " " <<  ci.Z << " " <<  endl;
          indexweight=get_barycentric_weight_for_ind(k,ang,ci,in,*_rel);
          weights.push_back(indexweight);
        }
    return weights;
  }

  void resampler::reverse_barycentric_weights(vector<std::map<int,double> > & weights, const int &origsize, const int &newsize){
    vector<map<int,double> > reordered_weights;
    reordered_weights.resize(newsize);
    for (int oldNode = 0; oldNode < origsize; ++oldNode)//this loop can't be parallelized
      {

	       for (map<int, double >::iterator iter = weights[oldNode].begin(); iter != weights[oldNode].end(); ++iter)//convert scattering weights to gathering weights
	        {
	           reordered_weights[iter->first][oldNode] = iter->second;

	          }
          }

    weights=reordered_weights;
  }

  vector<std::map<int,double> > resampler::get_all_adap_barycentric_weights(NEWMESH::newmesh &IN, NEWMESH::newmesh &ref,  boost::shared_ptr<NEWMESH::newmesh> &EXCL, boost::shared_ptr<RELATIONS > _rel){
    vector<map<int,double> > forward, reverse,reverse_reorder,adapt;
    vector<vector<int> > testinvert;
    double  ang;
    boost::shared_ptr<RELATIONS> REL,_reverseREL;

    if(ref.nvertices() > IN.nvertices()){ang=guess_angular_spacing(IN.nvertices()); }
    else{ang=guess_angular_spacing(ref.nvertices()); }
    if(!_rel.get() || (_rel->Ncols()!=IN.nvertices() && _rel->Ncols()!=ref.nvertices())){


      REL = boost::shared_ptr<RELATIONS > ( new RELATIONS(ref,IN,ang));
      REL->update_RELATIONS(ref);

      if(_method=="ADAP_BARY") {
	       _reverseREL=boost::shared_ptr<RELATIONS > ( new RELATIONS (REL->invert_relations(IN,ref,ang)));
      }
    }
    else{
      if(_rel->Ncols()==ref.nvertices()){ REL=_rel;
 	      if(_method=="ADAP_BARY") _reverseREL=boost::shared_ptr<RELATIONS > ( new RELATIONS (_rel->invert_relations(IN,ref,ang)));
      }
      else{
	       _reverseREL=_rel;
	       REL=boost::shared_ptr<RELATIONS > ( new RELATIONS (_rel->invert_relations(ref,IN,ang)));
      }
    }
    check_scale(IN,ref);

    forward=get_all_barycentric_weights(ref,IN,REL);  // get forward weights i.e. find for each reference vertex which source mesh vertices surround it and what are their barycentric weights
    if(_method=="BARYCENTRIC"){
      adapt=forward;

      }else{ // use adaptive barycentric weighting

      //	_reverseREL=boost::shared_ptr<RELATIONS > ( new RELATIONS (_rel->invert_relations(IN,exclusion,ang)));
      reverse=get_all_barycentric_weights(IN,ref,_reverseREL); // calculate reverse weights i.e. weights for each source mesh vertex

      int numOldNodes=IN.nvertices();
      int numNewNodes=ref.nvertices();
      vector<double> newAreas(numNewNodes,0);
      vector<double> oldAreas(numOldNodes,0);
      vector<double> correction(numOldNodes,0);

      //// now reorder reverse weights so they can be used to sample forwards
      reverse_reorder.resize((int) forward.size());
      adapt.resize((int) forward.size());
      for (int oldNode = 0; oldNode < numOldNodes; ++oldNode)//this loop can't be parallelized
	{
	  oldAreas[oldNode]=computevertexArea(oldNode,IN);
	  for (map<int, double >::iterator iter = reverse[oldNode].begin(); iter != reverse[oldNode].end(); ++iter)//convert scattering weights to gathering weights
	    {
	      reverse_reorder[iter->first][oldNode] = iter->second;
	    }
	}
      for (int newNode = 0; newNode < numNewNodes; ++newNode)
	{
	  if(EXCL.get()==0 || EXCL->get_pvalue((*REL)(1,newNode+1)-1)){
	    newAreas[newNode]=computevertexArea(newNode,ref);

	    if ( reverse_reorder[newNode].size() <=forward[newNode].size() )  // choose whichever has the most samples as this will be the higher res mesh (this differs from source code)
	      adapt[newNode] = forward[newNode];
	    else
	      adapt[newNode] = reverse_reorder[newNode];


	    for (map<int, double>::iterator iter = adapt[newNode].begin(); iter != adapt[newNode].end(); ++iter)//begin area correction by multiplying by target node area
	      {
	       	iter->second *= newAreas[newNode];//begin the process of area correction by multiplying by gathering node areas

		correction[iter->first] += iter->second;//now, sum the scattering weights to prepare for first normalization
	      }
	    }
	}

      for (int newNode = 0; newNode < numNewNodes; ++newNode)
	{
	  double weightsum = 0.0f;
	  if(EXCL.get()==0 || EXCL->get_pvalue((*REL)(1,newNode+1)-1)){

	    for (map<int, double>::iterator iter = adapt[newNode].begin(); iter != adapt[newNode].end(); ++iter)//begin area correction by multiplying by target node area
	      {
		iter->second *= oldAreas[iter->first] / correction[iter->first];//divide the weights by their scatter sum, then multiply by current areas
                weightsum += iter->second;//and compute the sum

	      }
	      if (weightsum != 0.0)//this shouldn't happen unless no nodes remain due to roi, or node areas can be zero
		{
		  for (map<int, double>::iterator iter = adapt[newNode].begin(); iter != adapt[newNode].end(); ++iter)
		    {
		      iter->second /= weightsum;//and normalize to a sum of 1
		    }
		}

	  }
	}
    }
    return adapt;
  }

  map<int,double> resampler::get_barycentric_weight_for_ind(int k,double ang, NEWMESH::newmesh &MESHA, NEWMESH::newmesh &MESHB, RELATIONS &_rel){

    map<int,double> vertex_weights;
    bool closestfound=false, notfound=false;

    Pt ci,v0,v1,v2;
    int n0,n1,n2,ind,t_ind;
    double dist=0,maxsepA=0;

    ci =MESHA.get_coord(k);
    ind=1;

    while(closestfound==false){
      if(ind>_rel.Nrows(k+1)) {
	       if(notfound){ang=2*ang; } // break;
	        else{
	           for (vector<int>::const_iterator it=MESHA.nbegin(k); it!=MESHA.nend(k); it++){
	              dist=(ci-MESHA.get_coord(*it)).norm();
	              if(dist>maxsepA)
	               maxsepA=dist;
	            }
	            ang=2*asin(maxsepA/(RAD));
	         }
	          _rel.update_RELATIONS_for_ind(k+1,MESHA,ang);
	           maxsepA=0; ind=1;
	           notfound=true;

        }

      t_ind=_rel(ind,k+1)-1;
      if(MESHB.return_closest_points(t_ind,v0,v1,v2,ci,n0,n1,n2)){
	vertex_weights=get_barycentric_weights(v0,v1,v2,ci,n0,n1,n2);  /// barycentric weights stored as a map with first term = indices of neighbouring points and second = barycentric weights

	closestfound=true;
      }
      ind++;

    }
    return vertex_weights;
  }

  map<int,double> resampler::get_barycentric_weight_for_ind(const int &k,const double &ang,const Pt &ci,const NEWMESH::newmesh &MESHB,const  RELATIONS &_rel){

    map<int,double> vertex_weights;
    bool closestfound=false;
    double newang=ang;


    vector< pair<float,int> > neighbours=_rel.update_RELATIONS_for_ind(ci,ang);
    //cout << " get_barycentric_weight_for_ind( " << ci.X << " " << ci.Y << " mesh vertices " << MESHB.nvertices() <<" neighbours size " << neighbours.size() <<  " " << ang << endl;


    Pt v0,v1,v2;
    int n0,n1,n2,ind,t_ind;

    ind=1;

    while(closestfound==false){
      if(ind>(int) neighbours.size()) {
	       newang=2*newang;
	        neighbours=_rel.update_RELATIONS_for_ind(ci,newang);
	         ind=1;
           //cout << " ind>(int) neighbours.size() " << ang << endl;
      }

      t_ind=neighbours[ind-1].second;
      //cout << t_ind << "  neighbour " << ind << endl;
      if(MESHB.return_closest_points(t_ind,v0,v1,v2,ci,n0,n1,n2)){
         //cout << " (MESHB.return_closest_points(t_ind,v0,v1,v2,ci,n0,n1,n2) = TRUE " << n0 << " " << n1 << " " << n2 << endl;
	       vertex_weights=get_barycentric_weights(v0,v1,v2,ci,n0,n1,n2);  /// barycentric weights stored as a map with first term = indices of neighbouring points and second = barycentric weights
         closestfound=true;
      }else
      ind++;

    }

    return vertex_weights;
  }


    ///// if upsampling then target mesh will be higher resolution it is quicker to calculate the neighbourhood mask  from the high res mesh to the low res mesh so upsampling is performed by finding all neighbours of the target and summing their contribution and downsampling is found by instead finding all neighbours of the source

  void resampler::upsample_w_interpolation(NEWMESH::newmesh &in,const NEWMESH::newmesh &SPH,  const double &maxdist, boost::shared_ptr<BFMatrix> &data, boost::shared_ptr<NEWMESH::newmesh> &EXCL, boost::shared_ptr<RELATIONS > _rel){

    if(_method=="NN"){
      nearest_neighbour_interpolation(in,SPH, data, EXCL);
    }
    else{
      check_scale(in,SPH);

      Matrix newdata(data->Nrows(),SPH.nvertices()); newdata=0;

      NEWMESH::newmesh exclusion=SPH;

      if(!_rel.get()||_rel->Ncols()!=SPH.nvertices()){ double ang=4*asin(maxdist/(2*RAD));  _rel = boost::shared_ptr<RELATIONS > ( new RELATIONS(exclusion,in,ang));
      _rel->update_RELATIONS(exclusion);}

      ColumnVector V;

      NEWMESH::Pt cr, ci;
      double norm,weight;
      double SUM,excl_sum;
      for (int i = 1; i <=  SPH.nvertices(); i++){
	exclusion.set_pvalue(i-1,0);

	ci = SPH.get_coord(i-1);

	SUM=0;excl_sum=0.0;
	if(EXCL.get()==0 || (_rel->Nrows(i) >0 && EXCL->get_pvalue((*_rel)(1,i)-1) > 0 )){
	  for (int j = 1; j<= _rel->Nrows(i); j++){


	    cr=in.get_coord((*_rel)(j,i)-1);

	    norm = (ci - cr).norm();
	    norm=2*RAD*asin(norm/(2*RAD)); // geodesic dist

	    weight=calc_weight(norm,maxdist);
	    excl_sum+=weight;

	    if(EXCL.get()) weight=EXCL->get_pvalue((*_rel)(j,i)-1)*weight;

	    SUM+=weight;

	    for (BFMatrixColumnIterator it = data->begin((*_rel)(j,i)); it != data->end((*_rel)(j,i)); ++it) {
	      newdata(it.Row(),i)=newdata(it.Row(),i)+*it*weight;
	    }

	  }
	  if(excl_sum)
	    if(EXCL.get()){ exclusion.set_pvalue(i-1,SUM/excl_sum);}

	  for (int it=1; it <= newdata.Nrows();it++){
	    if(SUM!=0)
	      newdata(it,i)=newdata(it,i)/SUM;

	  }


	}else  exclusion.set_pvalue(i-1,0);



      }

      if(EXCL.get()) *EXCL=exclusion;

      boost::shared_ptr<FullBFMatrix > pin =boost::dynamic_pointer_cast<FullBFMatrix>(data);
      if(pin.get()){data = boost::shared_ptr<BFMatrix >(new FullBFMatrix (newdata));}
      else {data = boost::shared_ptr<BFMatrix >(new  SparseBFMatrix<double>  (newdata)); }


    }


  }

  void resampler::downsample_w_interpolation(NEWMESH::newmesh &in,const NEWMESH::newmesh &SPH, const double &maxdist,  boost::shared_ptr<BFMatrix> &data, boost::shared_ptr<NEWMESH::newmesh> &EXCL, boost::shared_ptr<RELATIONS > _rel){


    if(_method=="NN" ){
      nearest_neighbour_interpolation(in,SPH,data,EXCL);
    }
    else{

      check_scale(in,SPH);
      Matrix newdata(data->Nrows(),SPH.nvertices()); newdata=0;
      ColumnVector  SUM(SPH.nvertices()),excl_sum(SPH.nvertices()); // as neighbourhoods are built for the source the SUM of all source points contributing to each target is a vector and is gradually filled as we iterate through the source data
      ColumnVector exclude(SPH.nvertices());
      ColumnVector check(SPH.nvertices()); check=0;

      SUM=0 ;excl_sum=0; exclude=0;
       NEWMESH::newmesh exclusion=SPH;
      if(!_rel.get() ||_rel->Ncols()!=in.nvertices() ){ double ang=4*asin(maxdist/(2*RAD));  _rel = boost::shared_ptr<RELATIONS > ( new RELATIONS(in,SPH,ang));
      _rel->update_RELATIONS(in);}


      ColumnVector V;

      NEWMESH::Pt cr, ci;
      double norm,weight,val;

      for (int i = 0; i <  in.nvertices(); i++){

	ci = in.get_coord(i);

	for (int j = 1; j <= _rel->Nrows(i+1); j++){
	  cr=SPH.get_coord((*_rel)(j,i+1)-1);
	  norm = (ci - cr).norm();
	  norm=2*RAD*asin(norm/(2*RAD)); // geodesic dist
	  check((*_rel)(j,i+1))=1;

	  if(!EXCL.get() || EXCL->get_pvalue(i) > 0 ){


	    weight=calc_weight(norm,maxdist);
	    if(EXCL.get()) exclude((*_rel)(j,i+1))+=weight*EXCL->get_pvalue(i);

	         excl_sum((int) (*_rel)(j,i+1))+=weight;

	    if(EXCL.get()) weight=EXCL->get_pvalue(i)*weight;

	    SUM((*_rel)(j,i+1))+=weight;

	    for (BFMatrixColumnIterator it = data->begin(i+1); it != data->end(i+1); ++it) {
	      val=newdata(it.Row(),(*_rel)(j,i+1))+*it*weight;
	      newdata(it.Row(),(*_rel)(j,i+1))=val;
	    }

	  }
	}

      }




      for (int i = 1; i <=  SPH.nvertices(); i++){
	exclusion.set_pvalue(i-1,0);

	if(SUM(i)>0) exclude(i)/=SUM(i);
	if(!EXCL.get() || exclude(i)>0.5){
	  if(excl_sum(i)!=0)
	    exclusion.set_pvalue(i-1,SUM(i)/excl_sum(i));
	  for (int it=1;it<=newdata.Nrows();it++){ //BFMatrixColumnIterator it = data->begin(i); it != data->end(i); ++it){
	    if(SUM(i)!=0)
	      newdata(it,i)=newdata(it,i)/SUM(i);
	    if(newdata(it,i)!=newdata(it,i)){ cout << "Nan " << endl;exit(1) ; }
	  }

	}else {
	  exclusion.set_pvalue(i-1,0);
	  for (int it=1;it<=newdata.Nrows();it++)
	    newdata(it,i)=0;

	}
      }
       if(EXCL.get()) *EXCL=exclusion;
      boost::shared_ptr<FullBFMatrix > pin =boost::dynamic_pointer_cast<FullBFMatrix>(data);
      if(pin.get()){data = boost::shared_ptr<BFMatrix >(new FullBFMatrix (newdata));}
      else {data = boost::shared_ptr<BFMatrix >(new  SparseBFMatrix<double>  (newdata)); }

    }


  }

  void resampler::smooth_data(const double & sigma,  boost::shared_ptr<BFMatrix>  &data, NEWMESH::newmesh  &IN)
  {
    boost::shared_ptr<NEWMESH::newmesh> EXCL;

    upsample_w_interpolation(IN,IN,sigma,data,EXCL);
  }

  void resampler::resampledata(NEWMESH::newmesh& in,const NEWMESH::newmesh& ref,boost::shared_ptr<NEWMESH::newmesh> &EXCL, Matrix &data, const double &sigma, boost::shared_ptr<RELATIONS > _rel){
    boost::shared_ptr<BFMatrix> BFDATA;
    bool transpose=false;

    BFDATA=boost::shared_ptr<BFMatrix> (new  FullBFMatrix (data));
    if((int) BFDATA->Ncols()!=in.nvertices()){
      if((int) BFDATA->Nrows()==in.nvertices()){transpose=true; BFDATA=BFDATA->Transpose();}
      else{
	cout << "resampler::error data has incompatible dimensions " << BFDATA->Ncols() << " " <<  BFDATA->Nrows() << " in.nvertices " << in.nvertices() <<  endl;
	exit(2);
      }
    }
    resampledata(in,ref,EXCL,BFDATA,sigma,_rel)  ;
    if(transpose) BFDATA=BFDATA->Transpose();
    data=BFDATA->AsMatrix();
  }

   void resampler::resampledata(NEWMESH::newmesh& in,const NEWMESH::newmesh& ref, Matrix &data, const double &sigma, boost::shared_ptr<RELATIONS > _rel){
    boost::shared_ptr<NEWMESH::newmesh> EXCL;

    resampledata(in,ref,EXCL,data,sigma,_rel)  ;

  }


  void resampler::resampledata(NEWMESH::newmesh& in,const NEWMESH::newmesh& ref, boost::shared_ptr<NEWMESH::newmesh> &EXCL, SpMat<double> &data, const double &sigma, boost::shared_ptr<RELATIONS > _rel){
    boost::shared_ptr<BFMatrix> BFDATA;

    BFDATA=boost::shared_ptr<BFMatrix> (new  SparseBFMatrix<double> (data));
    resampledata(in,ref,EXCL,BFDATA,sigma,_rel)  ;
    SpMat<double> tmp(BFDATA->Nrows(),BFDATA->Ncols());
    data=tmp;
    for (unsigned int i=1;i<=BFDATA->Nrows();i++){
      for (unsigned int j=1;j<=BFDATA->Ncols();j++)
	data.Set(i,j,BFDATA->Peek(i,j));
    }

  }

  // scalar data
  void resampler::resample_scalar(NEWMESH::newmesh &in,const NEWMESH::newmesh &SPH,const double & maxdist,boost::shared_ptr<NEWMESH::newmesh> EXCL){


    boost::shared_ptr<BFMatrix> scalardata;
    scalardata =boost::shared_ptr<BFMatrix> (new FullBFMatrix (1,in.nvertices()));
    NEWMESH::newmesh output=SPH;

    for (int i=0;i<in.nvertices();i++)
      scalardata->Set(1,i+1,in.get_pvalue(i));

    resampledata(in,SPH,EXCL,scalardata,maxdist)  ;

    for (int i=0;i<SPH.nvertices();i++){
      output.set_pvalue(i,scalardata->Peek(1,i+1));
    }

    in=output;

  }


  ////  ANATOMICAL MESH DOWNSAMPLING /////////////////

  newmesh mesh_resample(const NEWMESH::newmesh &ANAT, NEWMESH::newmesh SPH_ORIG, NEWMESH::newmesh & SPH_LOW, bool _adapbary, boost::shared_ptr<RELATIONS > _rel){
    vector<map<int,double> > baryweights;
    Pt newPt, newPtsphere, ci;
    resampler R;
    newmesh ANAT_LOW=SPH_LOW;
    boost::shared_ptr<NEWMESH::newmesh> EXCL;
    vector<int> neighbours;


    baryweights=R.get_all_barycentric_weights(SPH_LOW,SPH_ORIG);

    for(int i=0;i<SPH_LOW.nvertices();i++){
      newPt*=0.0;newPtsphere*=0.0;
      ci=SPH_LOW.get_coord(i); neighbours.clear();
      for (std::map<int,double>::iterator it=baryweights[i].begin(); it!=baryweights[i].end(); ++it){
	newPt+=ANAT.get_coord(it->first)*it->second;
	Pt P=ANAT.get_coord(it->first);
      }


      ANAT_LOW.set_coord(i,newPt);
    }

    return ANAT_LOW;
  }
  ////////////////////////////////////  MESH INTERPOLATION FUNCTIONS /////////////////////////////////

  //// barycentric mesh interpolation is like surface-project-reproject in Caret
   void barycentric_mesh_interpolation(NEWMESH::newmesh &SPH_up, NEWMESH::newmesh &SPH_low_init,NEWMESH::newmesh &SPH_low_final,boost::shared_ptr<RELATIONS > _rel){
     vector<map<int,double> >  baryweights;
     map<int,double> weight;
     Pt newPt,ci;
     resampler R;
     double ang;
     if(!_rel.get() ||_rel->Ncols()!=SPH_up.nvertices())
      {
	if(SPH_low_init.nvertices() > SPH_up.nvertices()){ ang=R.guess_angular_spacing(SPH_up.nvertices());}
	else{ ang=2*asin(Calculate_MaxVD(SPH_low_init)/(RAD));}
	_rel = boost::shared_ptr<RELATIONS > ( new RELATIONS(SPH_up,SPH_low_init,ang));
	_rel->update_RELATIONS(SPH_up);
      }else ang=_rel->get_ang();

    for(int i=0;i<SPH_up.nvertices();i++){
      weight=R.get_barycentric_weight_for_ind(i,_rel->get_ang(),SPH_up,SPH_low_init,*_rel);
      newPt*=0.0;
      ci=SPH_up.get_coord(i);
      for (std::map<int,double>::iterator it=weight.begin(); it!=weight.end(); ++it){

	newPt+=SPH_low_final.get_coord(it->first)*it->second;

      }


      newPt.normalize();
      newPt=newPt*100;
      SPH_up.set_coord(i,newPt);
    }


   }


  // Here reference is the untransformed regular mesh and IN is the transformed mesh
  void upsample_transform_RBF(NEWMESH::newmesh &SPH_up, NEWMESH::newmesh &SPH_ref,NEWMESH::newmesh &SPH_in, const double & ang){


    RELATIONS _rel(SPH_up,SPH_ref,ang);
    _rel.update_RELATIONS(SPH_up);
    upsample_transform_RBF(SPH_up,SPH_ref,SPH_in,_rel);

  }

  ////// THIS IS THIN PLATE SPLINE INTERPOLATION, following Alex Petrovic implementation following Bookstein
  void upsample_transform_RBF(NEWMESH::newmesh &SPH_up, NEWMESH::newmesh &SPH_ref,NEWMESH::newmesh &SPH_in, const RELATIONS & _rel){
    if(!check_scale(SPH_ref,RAD)) true_rescale(SPH_ref,RAD);
    if(!check_scale(SPH_ref,RAD)) true_rescale(SPH_in,RAD);
    if(!check_scale(SPH_ref,RAD)) true_rescale(SPH_up,RAD);



    RBF CRx;
    RBF CRy;
    RBF CRz;
    for (int i = 0; i < SPH_up.nvertices(); i++){
      NEWMESH::Pt cr = SPH_up.get_coord(i);


      if(cr.norm()){

	for (int j = 1; j <= _rel.Nrows(i+1); j++){
	  NEWMESH::Pt rm = SPH_ref.get_coord(_rel(j,i+1)-1);
	  NEWMESH::Pt im = SPH_in.get_coord(_rel(j,i+1)-1);
	  CRx.set_point(rm.X,rm.Y,rm.Z,im.X);
	  CRy.set_point(rm.X,rm.Y,rm.Z,im.Y);
	  CRz.set_point(rm.X,rm.Y,rm.Z,im.Z);

	}
	double x = CRx.interpolate(cr.X, cr.Y, cr.Z);
	double y = CRy.interpolate(cr.X, cr.Y, cr.Z);
	double z = CRz.interpolate(cr.X, cr.Y, cr.Z);

	NEWMESH::Pt pp(x,y,z);
	pp.normalize();
	pp = pp*RAD;


	SPH_up.set_coord(i,pp);

	CRx.clean();
	CRy.clean();
	CRz.clean();

      }
    }


  }

  /////////////////////// THE FOLLOWING FUNCTIONS ARE FOR PIECEWISE AFFINE MESH INTERPOLATIONS /////////////////////

  void project(NEWMESH::newmesh & mesh_new, NEWMESH::newmesh & mesh_A,const NEWMESH::newmesh & mesh_B,const double &ang){

    RELATIONS _rel(mesh_new,mesh_A,ang); _rel.update_RELATIONS(mesh_new);

    for ( vector<boost::shared_ptr<NEWMESH::Mpoint> >::const_iterator i= mesh_new.vbegin(); i!=mesh_new.vend(); i++){
      NEWMESH::Pt p = (*i)->get_coord();
      NEWMESH::Pt p2;
      if(p.norm()){

	p2=  project(_rel.Col((*i)->get_no()+1),p,mesh_A, mesh_B);

	(*i)->set_coord(p2);
      }

    }


  }

  ////// this function is used to approximate the affine transform between two meshes by using linear regression (used to initialise --trans option in  msm)
  Matrix affine_transform(NEWMESH::newmesh & mesh_new, NEWMESH::newmesh mesh_A, NEWMESH::newmesh mesh_B){

     vector<int> list;
     ColumnVector P_in(4);
     Matrix AFFINETRANS;


     for (int i=1;i<=mesh_A.nvertices();i++)
       list.push_back(i);

     AFFINETRANS=return_affine_transform(list,mesh_A,mesh_B);


     for ( vector<boost::shared_ptr<NEWMESH::Mpoint> >::const_iterator i= mesh_new.vbegin(); i!=mesh_new.vend(); i++){
       NEWMESH::Pt p = (*i)->get_coord();
       NEWMESH::Pt p2;
       if(p.norm()){
	 P_in(1) = p.X; P_in(2) = p.Y; P_in(3) = p.Z; P_in(4) = 1;
	 P_in = AFFINETRANS * P_in;
	 p2.X=P_in(1); p2.Y=P_in(2);  p2.Z=P_in(3);
	 (*i)->set_coord(p2);
       }

     }

   return AFFINETRANS;
  }

  void affine_transform(NEWMESH::newmesh & mesh_new, Matrix AFFINETRANS ){

     vector<int> list;
     ColumnVector P_in(4);


     for ( vector<boost::shared_ptr<NEWMESH::Mpoint> >::const_iterator i= mesh_new.vbegin(); i!=mesh_new.vend(); i++){
       NEWMESH::Pt p = (*i)->get_coord();
       NEWMESH::Pt p2;
       if(p.norm()){
		P_in(1) = p.X; P_in(2) = p.Y; P_in(3) = p.Z; P_in(4) = 1;
		P_in = AFFINETRANS * P_in;
		p2.X=P_in(1); p2.Y=P_in(2);  p2.Z=P_in(3);
		(*i)->set_coord(p2);
       }

     }


  }


  Pt project(const vector<int>& R,const NEWMESH::Pt& cr,const NEWMESH::newmesh & mesh_A,const NEWMESH::newmesh & mesh_B){//  int ind_IN, int neigh){

    Matrix T=return_affine_transform(R,mesh_A,mesh_B);

    ColumnVector P_in(4), P_out(4);
    Pt pnew;
    P_in(1) = cr.X; P_in(2) = cr.Y; P_in(3) = cr.Z; P_in(4) = 1;
    P_out = T * P_in;
    pnew.X=P_out(1); pnew.Y=P_out(2);  pnew.Z=P_out(3);

    P_out.Release();

    //OUT(P_out);
    return pnew;
  }

  //// mesh A and mesh B represent the same mesh before and after some warp R lists the vertices from which the affine transformation will be approximate, this might be just a subsection of the mesh (for piecewise non linear transformation) or the whole mesh (global affine)
  Matrix return_affine_transform(const vector<int>& R,const NEWMESH::newmesh & mesh_A,const NEWMESH::newmesh & mesh_B){

    ColumnVector tilda(R.size()*3); tilda = 0;
    ColumnVector a(12); a = 0;
    Matrix T(4,4); T = 0;
    Matrix pseudo_inv(12,R.size()*3); pseudo_inv = 0;
    Matrix A(R.size()*3,12); A = 0;
    for (unsigned int i=1; i <= R.size(); i++){
      Pt p = mesh_B.get_coord(R[i-1]-1);
      tilda(3*(i-1)+1) = p.X;
      tilda(3*(i-1)+2) = p.Y;
      tilda(3*(i-1)+3) = p.Z;
    }
    for (unsigned int i=1; i <= R.size(); i++){

      Pt p = mesh_A.get_coord(R[i-1]-1);

      A(3*(i-1)+1,1) = p.X;
      A(3*(i-1)+1,2) = p.Y;
      A(3*(i-1)+1,3) = p.Z;
      A(3*(i-1)+1,4) = 1;
      A(3*(i-1)+1,5) = 0;
      A(3*(i-1)+1,6) = 0;
      A(3*(i-1)+1,7) = 0;
      A(3*(i-1)+1,8) = 0;
      A(3*(i-1)+1,9) = 0;
      A(3*(i-1)+1,10) = 0;
      A(3*(i-1)+1,11) = 0;
      A(3*(i-1)+1,12) = 0;

      A(3*(i-1)+2,1) = 0;
      A(3*(i-1)+2,2) = 0;
      A(3*(i-1)+2,3) = 0;
      A(3*(i-1)+2,4) = 0;
      A(3*(i-1)+2,5) = p.X;
      A(3*(i-1)+2,6) = p.Y;
      A(3*(i-1)+2,7) = p.Z;
      A(3*(i-1)+2,8) = 1;
      A(3*(i-1)+2,9) = 0;
      A(3*(i-1)+2,10) = 0;
      A(3*(i-1)+2,11) = 0;
      A(3*(i-1)+2,12) = 0;

      A(3*(i-1)+3,1) = 0;
      A(3*(i-1)+3,2) = 0;
      A(3*(i-1)+3,3) = 0;
      A(3*(i-1)+3,4) = 0;
      A(3*(i-1)+3,5) = 0;
      A(3*(i-1)+3,6) = 0;
      A(3*(i-1)+3,7) = 0;
      A(3*(i-1)+3,8) = 0;
      A(3*(i-1)+3,9) = p.X;
      A(3*(i-1)+3,10) = p.Y;
      A(3*(i-1)+3,11) = p.Z;
      A(3*(i-1)+3,12) = 1;

    }
    pseudo_inv = ((A.t()*A).i())*A.t();
    a = pseudo_inv * tilda;
    T(1,1) = a(1); T(1,2) = a(2); T(1,3) = a(3); T(1,4) = a(4);
    T(2,1) = a(5); T(2,2) = a(6); T(2,3) = a(7); T(2,4) = a(8);
    T(3,1) = a(9); T(3,2) = a(10); T(3,3) = a(11); T(3,4) = a(12);
    T(4,1) = 0; T(4,2) = 0; T(4,3) = 4; T(4,4) = 1;

    return T;
  }


 //// Calculate affine transformation between pointsets START and FINISH
  Matrix return_affine_transform(const vector<Pt>& START,const vector<Pt>& FINISH){

    ColumnVector tilda(START.size()*3); tilda = 0;
    ColumnVector a(12); a = 0;
    Matrix T(4,4); T = 0;
    Matrix pseudo_inv(12,START.size()*3); pseudo_inv = 0;
    Matrix A(START.size()*3,12); A = 0;
    if(FINISH.size()!=START.size()){ throw  NEWMESHException("RESAMPLER:: return_affine_transform, points sets are not the same size");}

    for (unsigned int i=0; i < FINISH.size(); i++){
      Pt p = FINISH[i];
      tilda(3*i+1) = p.X;
      tilda(3*i+2) = p.Y;
      tilda(3*i+3) = p.Z;
    }

    for (unsigned int i=0; i < START.size(); i++){

      Pt p = START[i];

      A(3*i+1,1) = p.X;
      A(3*i+1,2) = p.Y;
      A(3*i+1,3) = p.Z;
      A(3*i+1,4) = 1;
      A(3*i+1,5) = 0;
      A(3*i+1,6) = 0;
      A(3*i+1,7) = 0;
      A(3*i+1,8) = 0;
      A(3*i+1,9) = 0;
      A(3*i+1,10) = 0;
      A(3*i+1,11) = 0;
      A(3*i+1,12) = 0;

      A(3*i+2,1) = 0;
      A(3*i+2,2) = 0;
      A(3*i+2,3) = 0;
      A(3*i+2,4) = 0;
      A(3*i+2,5) = p.X;
      A(3*i+2,6) = p.Y;
      A(3*i+2,7) = p.Z;
      A(3*i+2,8) = 1;
      A(3*i+2,9) = 0;
      A(3*i+2,10) = 0;
      A(3*i+2,11) = 0;
      A(3*i+2,12) = 0;

      A(3*i+3,1) = 0;
      A(3*i+3,2) = 0;
      A(3*i+3,3) = 0;
      A(3*i+3,4) = 0;
      A(3*i+3,5) = 0;
      A(3*i+3,6) = 0;
      A(3*i+3,7) = 0;
      A(3*i+3,8) = 0;
      A(3*i+3,9) = p.X;
      A(3*i+3,10) = p.Y;
      A(3*i+3,11) = p.Z;
      A(3*i+3,12) = 1;

    }
    pseudo_inv = ((A.t()*A).i())*A.t();
    a = pseudo_inv * tilda;
     T(1,1) = a(1); T(1,2) = a(2); T(1,3) = a(3); T(1,4) = a(4);
    T(2,1) = a(5); T(2,2) = a(6); T(2,3) = a(7); T(2,4) = a(8);
    T(3,1) = a(9); T(3,2) = a(10); T(3,3) = a(11); T(3,4) = a(12);
    T(4,1) = 0; T(4,2) = 0; T(4,3) = 4; T(4,4) = 1;
    return T;
  }


  // this should be suitable for regular mesh only
  double Calculate_MVD(const NEWMESH::newmesh &SPH_in){

    //----mean_inter_vertex_distance----

    /// calculate from first vertex only so not true unless regularised!!!!
    double k = 0, kr = 0;
    double MVD;

    for (int i=0;i<SPH_in.nvertices();i++){
      NEWMESH::Pt Sp = SPH_in.get_coord(i);

      for (vector<int>::const_iterator j=SPH_in.nbegin(i); j!=SPH_in.nend(i); j++){
	k++;

	kr += (SPH_in.get_coord(*j) - Sp).norm();
      }

    }
    MVD = kr/k;

    return MVD;
  }

  double Calculate_MaxVD(const NEWMESH::newmesh &SPH_in){

      double dist=0,MAX=0;
      Pt CP;
      int tot=0;

      for (int k=0;k<SPH_in.nvertices();k++){

	CP=SPH_in.get_coord(k);
	for (vector<int>::const_iterator it=SPH_in.nbegin(k);it !=SPH_in.nend(k);it++){

	  dist= 2*RAD*asin((CP-SPH_in.get_coord(*it)).norm()/(2*RAD));
	  tot++;
	  if(dist > MAX)
	    MAX=dist;
	}
      }
    return MAX;
  }

  double Calculate_MaxVD(const NEWMESH::newmesh &SPH_in, ColumnVector &SEP){

    double dist=0, MAX=0;
    Pt CP;
    int tot=0;
    for (int k=0;k<SPH_in.nvertices();k++){

      CP=SPH_in.get_coord(k);
      for (vector<int>::const_iterator it=SPH_in.nbegin(k);it !=SPH_in.nend(k);it++){

	dist= 2*RAD*asin((CP-SPH_in.get_coord(*it)).norm()/(2*RAD));
	tot++;
	if(dist > SEP(k+1))
	 SEP(k+1)=dist;
	if(SEP(k+1)>MAX) MAX=SEP(k+1);
      }
    }
    return MAX;
  }


  double computevertexArea(const int & ind,const NEWMESH::newmesh &MESH)
  {
    double sum=0;
    for ( vector<int>::const_iterator i=MESH.tIDbegin(ind); i!=MESH.tIDend(ind); i++){
	NEWMESH::Pt v1 = MESH.get_triangle_vertex(*i,0),  v2 = MESH.get_triangle_vertex(*i,1), v3 = MESH.get_triangle_vertex(*i,2);
	sum+=computeArea(v1,v2,v3);
    }

    return sum/MESH.get_total_triangles(ind);
  }

  double computeArea(const NEWMESH::Pt &v0, const NEWMESH::Pt &v1, const NEWMESH::Pt &v2)
  {
    NEWMESH::Pt res, v0v1, v0v2;

    v0v1 = v1 - v0;
    v0v2 = v2 - v0;
    res=v0v1*v0v2;
    return 0.5*res.norm();
  }

  // for barycentric interpolation
  double barycentric(const NEWMESH::Pt& v1, const NEWMESH::Pt& v2, const NEWMESH::Pt& v3, const NEWMESH::Pt& vref, const double &va1, const double &va2,const double &va3){


    double A, Aa,Ab,Ac;

    Aa=computeArea(vref, v2,v3);
    Ab=computeArea(vref, v1,v3);
    Ac=computeArea(vref, v1,v2);

    A=Aa+Ab+Ac;
    Aa=Aa/A;
    Ab=Ab/A;
    Ac=Ac/A;
    return  Aa*va1+Ab*va2+Ac*va3;
  }

 // for barycentric interpolation
  Pt barycentric(const NEWMESH::Pt& v1, const NEWMESH::Pt& v2, const NEWMESH::Pt& v3, const NEWMESH::Pt& vref, const NEWMESH::Pt& va1, const NEWMESH::Pt& va2, const NEWMESH::Pt& va3){


    double A, Aa,Ab,Ac;


    Aa=computeArea(vref, v2,v3);
    Ab=computeArea(vref, v1,v3);
    Ac=computeArea(vref, v1,v2);


    A=Aa+Ab+Ac;
    Aa=Aa/A;
    Ab=Ab/A;
    Ac=Ac/A;

    return  va1*Aa+va2*Ab+va3*Ac;
  }
  // for barycentric interpolation // this should be bf matrices and columnvectors!
  SpMat<double> barycentric(const NEWMESH::Pt& v1, const NEWMESH::Pt& v2, const NEWMESH::Pt& v3, const NEWMESH::Pt& vref, const int &n1,const int &n2,const int &n3,const SpMat<double> &DATA){

    map<int,double> weights=get_barycentric_weights(v1,v2,v3,vref,n1,n2,n3);
    SpMat<double> NEW(DATA.Nrows(),DATA.Ncols());
    for(unsigned int i=1;i<=DATA.Nrows();i++){
      NEW.Set(i,n1,weights[n1]*DATA.Peek(i,n1)+weights[n2]*DATA.Peek(i,n2)+weights[n3]*DATA.Peek(i,n3));
    }
    return  NEW;
  }


  void barycentric2(const NEWMESH::Pt& v1, const NEWMESH::Pt& v2, const NEWMESH::Pt& v3, const NEWMESH::Pt& vref, const int &n1,const int &n2,const int &n3,const BFMatrix &DATA,vector<double> &DATA2){

    map<int,double> weights=get_barycentric_weights(v1,v2,v3,vref,n1,n2,n3);
    for(unsigned int i=1;i<=DATA.Nrows();i++){
      DATA2.push_back(weights[n1]*DATA.Peek(i,n1)+weights[n2]*DATA.Peek(i,n2)+weights[n3]*DATA.Peek(i,n3));
    }

  }

  std::map<int,double> get_barycentric_weights(const NEWMESH::Pt& v1, const NEWMESH::Pt& v2, const NEWMESH::Pt& v3, const NEWMESH::Pt& vref, const int &n1,const int &n2,const int &n3){

    float A, Aa,Ab,Ac;
    Pt PP;

    map<int,double> weights;

    projectPoint(vref,v1,v2,v3,PP);
    Aa=computeArea(PP, v2,v3);
    Ab=computeArea(PP, v1,v3);
    Ac=computeArea(PP, v1,v2);

    A=Aa+Ab+Ac;

    weights[n1]=Aa/A;
    weights[n2]=Ab/A;
    weights[n3]=Ac/A;


    return weights;
  }



}
