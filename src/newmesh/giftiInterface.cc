/*  giftiInterface.cc
    Matthew Webster (FMRIB) 
    Copyright (C) 2013 University of Oxford  */
/*  CCOPYRIGHT  */

#include "giftiInterface.h"
#include <algorithm>
#include <cstring>
#include "boost/lexical_cast.hpp"
using namespace std;

vector<GIFTIfield> GIFTIwrapper::returnSurfaceFields() { //returns the first NIFTI_INTENT_POINTSET and NIFTI_INTENT_TRIANGLE ( in that order )
  
  bool foundPointSet(false);
  bool foundTriangle(false);
  vector<GIFTIfield> surfaceFields;
  for ( unsigned int field = 0 ; field < allFields.size() && !foundPointSet; field++ ) {
    if ( allFields[field].intent == NIFTI_INTENT_POINTSET ) { 
      surfaceFields.push_back(allFields[field]);
      foundPointSet=true;
    }
  }
  for ( unsigned int field = 0 ; field < allFields.size()  && !foundTriangle; field++ ) {
    if ( allFields[field].intent == NIFTI_INTENT_TRIANGLE ) {
      surfaceFields.push_back(allFields[field]);
      foundTriangle=true;
    }
  }
  if ( !foundPointSet || !foundTriangle )
    throw GiftiException("Unable to locate surface information in file");
  return surfaceFields;
}

vector<GIFTIfield> GIFTIwrapper::returnNonSurfaceFields() { //removes a field array minus the first NIFTI_INTENT_POINTSET and NIFTI_INTENT_TRIANGLE fields ( if extant )
  
  bool removedPointset(false);
  bool removedTriangle(false);
  vector<GIFTIfield> nonSurfaceFields(allFields);
  for ( unsigned int field = 0; !removedPointset && field < nonSurfaceFields.size() ; field++ ) { //only removes first pointset - this code is not suitable for multiple surfaces in a file
    if ( nonSurfaceFields[field].intent == NIFTI_INTENT_POINTSET ) {
      nonSurfaceFields.erase(nonSurfaceFields.begin()+field);
      removedPointset=true;
    }
  }
  for ( unsigned int field = 0; !removedTriangle && field < nonSurfaceFields.size() ; field++ ) { //only removes first triangle field - this code is not suitable for multiple surfaces in a file
    if ( nonSurfaceFields[field].intent == NIFTI_INTENT_TRIANGLE ) {
      nonSurfaceFields.erase(nonSurfaceFields.begin()+field);
      removedTriangle=true;
    }
  }
  return nonSurfaceFields;
}

int GIFTIwrapper::readGIFTI(const string filename)
{
  gifti_image* rawSurface = gifti_read_image(filename.c_str(), 1);
  //Label data
  int Nentries = rawSurface->labeltable.length;
  int* key = rawSurface->labeltable.key;
  float* rgba = rawSurface->labeltable.rgba;
  for ( int i = 0 ; i < Nentries ; ++i,++key)
    {
      if ( rawSurface->labeltable.rgba != NULL ) {
	GIFTIlabels[*key].RGBA=vector<float>(rgba,rgba+4);
	rgba+=4;
      }
      GIFTIlabels[*key].name=string(rawSurface->labeltable.label[i]);
    }
  //Read data

  int NumberOfDataArrays = rawSurface->numDA;
  for ( int array = 0 ; array < NumberOfDataArrays ; array++ )
      allFields.push_back(GIFTIfield( rawSurface->darray[array]));

  for ( int tag=0;tag<rawSurface->meta.length;tag++) 
    metaData.push_back(GIFTImeta(rawSurface->meta.name[tag],rawSurface->meta.value[tag]));
  for ( int tag=0;tag<rawSurface->ex_atrs.length;tag++) 
    extraAttributes.push_back(GIFTImeta(rawSurface->ex_atrs.name[tag],rawSurface->ex_atrs.value[tag]));

  gifti_free_image(rawSurface);
  return 0;
}

int GIFTIwrapper::writeGIFTI(const string filename, int encoding) const
{
  gifti_set_verb(0);
  gifti_image *rawGIFTI =gifti_create_image( -1, -1, -1, -1, NULL, -1 );
  gifti_add_empty_darray(rawGIFTI,  allFields.size());
  if ( rawGIFTI==NULL )
    throw GiftiException("Error in Gii create!!!");

  for ( unsigned int tag=0;tag<metaData.size();tag++) {
      gifti_add_to_meta(&rawGIFTI->meta,metaData[tag].name.c_str(),metaData[tag].value.c_str(), 1);
      //cout << metaData[tag].name << " " << metaData[tag].value << endl;

  }
  for ( unsigned int tag=0;tag<extraAttributes.size();tag++)
    gifti_add_to_nvpairs(&rawGIFTI->ex_atrs,extraAttributes[tag].name.c_str(),extraAttributes[tag].value.c_str());

  rawGIFTI->labeltable.length=GIFTIlabels.size();
  rawGIFTI->labeltable.key = (int *)malloc(rawGIFTI->labeltable.length*sizeof(int));
  rawGIFTI->labeltable.label = (char **)malloc(rawGIFTI->labeltable.length * sizeof(char *)); 
  if ( GIFTIlabels.begin()->second.RGBA.size()!=0)
    rawGIFTI->labeltable.rgba = (float *)malloc(rawGIFTI->labeltable.length*4*sizeof(float));
  int count=0;
  for (map<int,GIFTIlabel>::const_iterator it=GIFTIlabels.begin(); it!=GIFTIlabels.end(); ++it,count++) {
    rawGIFTI->labeltable.key[count]=it->first;
    for ( int j =0; j < 4 && it->second.RGBA.size()!=0; j++ ) 
      rawGIFTI->labeltable.rgba[j+count*4]=it->second.RGBA[j];
    rawGIFTI->labeltable.label[count] = gifti_strdup(it->second.name.c_str());
  }

  for ( unsigned int field = 0 ; field < allFields.size() ; field++ )
  {
    for ( unsigned int tag=0;tag<allFields[field].metaData.size();tag++) 
      gifti_add_to_meta(&rawGIFTI->darray[field]->meta,allFields[field].metaData[tag].name.c_str(),allFields[field].metaData[tag].value.c_str(), 1);
    for ( unsigned int tag=0;tag<allFields[field].extraAttributes.size();tag++)
      gifti_add_to_nvpairs(&rawGIFTI->darray[field]->ex_atrs,allFields[field].extraAttributes[tag].name.c_str(),allFields[field].extraAttributes[tag].value.c_str());

    rawGIFTI->darray[field]->numCS=0;	 
    if ( allFields[field].intent == NIFTI_INTENT_POINTSET )  { //write out any coordinate systems
      if ( allFields[field].coordSystems.size() == 0 ) 
	throw GiftiException("Field "+boost::lexical_cast<string>(field)+" has intent NIFTI_INTENT_POINTSET but no coordSystem.");
      for ( unsigned int system=0;system<allFields[field].coordSystems.size();system++) {
	gifti_add_empty_CS(rawGIFTI->darray[field]);
	rawGIFTI->darray[field]->coordsys[system]->dataspace = gifti_strdup(allFields[field].coordSystems[system].dataSpace.c_str());
	rawGIFTI->darray[field]->coordsys[system]->xformspace = gifti_strdup(allFields[field].coordSystems[system].transformSpace.c_str());
	for(int row=0;row<4;row++) {
	  for(int column=0;column<4;column++)
	    rawGIFTI->darray[field]->coordsys[system]->xform[row][column] = allFields[field].coordSystems[system].transform[column+row*4];
	}
      }
    }

    rawGIFTI->darray[field]->intent = allFields[field].intent;
    rawGIFTI->darray[field]->encoding = encoding;	
    rawGIFTI->darray[field]->num_dim=allFields[field].dims.size();
    long long nvals = rawGIFTI->darray[field]->dims[0] = allFields[field].dims[0];
    for(unsigned char dim=1;dim<allFields[field].dims.size();dim++) {
      nvals*=allFields[field].dims[dim];
      rawGIFTI->darray[field]->dims[dim] = allFields[field].dims[dim];
    }
    rawGIFTI->darray[field]->nvals = nvals;
    rawGIFTI->darray[field]->datatype = allFields[field].dataType;

    if ( allFields[field].dataType == NIFTI_TYPE_FLOAT32 ) {
      rawGIFTI->darray[field]->nbyper = sizeof(float);
      rawGIFTI->darray[field]->data=malloc(sizeof(float)*nvals);
      memcpy(rawGIFTI->darray[field]->data,allFields[field].fdata.data(), sizeof(float)*nvals);
    }
    if ( allFields[field].dataType ==NIFTI_TYPE_INT32 ) {
      rawGIFTI->darray[field]->nbyper = sizeof(int); 
      rawGIFTI->darray[field]->data=malloc(sizeof(int)*nvals);
      memcpy(rawGIFTI->darray[field]->data,allFields[field].idata.data(), sizeof(int)*nvals);
    }
    if ( allFields[field].dataType == NIFTI_TYPE_UINT8 ) {
      rawGIFTI->darray[field]->nbyper = sizeof(char);
      rawGIFTI->darray[field]->data=malloc(sizeof(char)*nvals);
      memcpy(rawGIFTI->darray[field]->data,allFields[field].bdata.data(), sizeof(char)*nvals);

    }
    rawGIFTI->darray[field]->ind_ord = allFields[field].ordering;
    rawGIFTI->darray[field]->ext_offset=allFields[field].externalOffset;
    rawGIFTI->darray[field]->ext_fname = gifti_strdup(allFields[field].externalFilename.c_str());
 }
 
  gifti_write_image( rawGIFTI , filename.c_str(), 1) ;
  gifti_free_image(rawGIFTI);		

  return 0;
}

void GIFTIwrapper::report() const
{
  cout << "Number of master metadata pairs: " << metaData.size() << endl;
  for ( unsigned int tag=0;tag<metaData.size();tag++)
    cout << metaData[tag].name << " " << metaData[tag].value << endl;
  cout << "Number of master extra attributes: " << extraAttributes.size() << endl;
  for ( unsigned int tag=0;tag<extraAttributes.size();tag++)
    cout << extraAttributes[tag].name << " " << extraAttributes[tag].value << endl;
  cout << "Number of labels: " << GIFTIlabels.size() << endl;
  for (map<int,GIFTIlabel>::const_iterator it=GIFTIlabels.begin(); it!=GIFTIlabels.end(); ++it) {
    cout << it->first << " ";
    for ( int j =0; j < 4 && it->second.RGBA.size()!=0; j++ ) 
      cout <<  it->second.RGBA[j] << " ";
    cout << it->second.name << endl;
  }
  for ( unsigned int field = 0 ; field < allFields.size() ; field++ ) {
    cout << "DATA FIELD: " << field << endl;
    allFields[field].report();
    cout << " after field " << endl;
  }
  cout << "END" << endl;
}

void GIFTIfield::setIScalar(const size_t location, const int value) 
{
  if ( dims.size() != 1 )
    throw GiftiException("Tried to set non-scalar GIFTI field to scalar");
  if ( dataType != NIFTI_TYPE_INT32 )
    throw GiftiException("Tried to set non-int GIFTI field with int data");
  idata[location]=value;
}

void GIFTIfield::setFScalar(const size_t location, const float value) 
{
  if ( dims.size() != 1 )
    throw GiftiException("Tried to set non-scalar GIFTI field to scalar");
  if ( dataType != NIFTI_TYPE_FLOAT32 )
    throw GiftiException("Tried to set non-float GIFTI field with float data");
  fdata[location]=value;
}

void GIFTIfield::setBScalar(const size_t location, const char value) 
{
  if ( dims.size() != 1 )
    throw GiftiException("Tried to set non-scalar GIFTI field to scalar");
  if ( dataType != NIFTI_TYPE_UINT8 )
    throw GiftiException("Tried to set non-char GIFTI field with char data");
  bdata[location]=value;
}

int GIFTIfield::asIscalar(const size_t location) const
{
  if ( dims.size() != 1 )
    throw GiftiException("Requested scalar data from non-scalar GIFTI field");
  if ( dataType == NIFTI_TYPE_FLOAT32 )
    return (int)fdata[location];
  if ( dataType == NIFTI_TYPE_UINT8 )
    return (int)bdata[location];
  return idata[location];
}

float GIFTIfield::asFscalar(const size_t location) const
{
  if ( dims.size() != 1 )
    throw GiftiException("Requested scalar data from non-scalar GIFTI field");
  if ( dataType ==NIFTI_TYPE_INT32 )
    return (int)idata[location];
  if ( dataType == NIFTI_TYPE_UINT8 )
    return (int)bdata[location];
  return fdata[location];
}

char GIFTIfield::asBscalar(const size_t location) const
{
  if ( dims.size() != 1 )
    throw GiftiException("Requested scalar data from non-scalar GIFTI field");
  if ( dataType == NIFTI_TYPE_FLOAT32 )
    return (int)fdata[location];
  if ( dataType ==NIFTI_TYPE_INT32 )
    return (int)idata[location];
  return bdata[location];
}


int GIFTIfield::iScalar(const size_t location) const
{
  if ( dims.size() != 1 )
    throw GiftiException("Requested scalar data from non-scalar GIFTI field");
  if ( dataType != NIFTI_TYPE_INT32 )
    throw GiftiException("Requested non-int data from int GIFTI field");
  return idata[location];
}

float GIFTIfield::fScalar(const size_t location) const
{
  if ( dims.size() != 1 )
    throw GiftiException("Requested scalar data from non-scalar GIFTI field");
  if ( dataType != NIFTI_TYPE_FLOAT32 )
    throw GiftiException("Requested non-float data from float GIFTI field");
  return fdata[location];
}

char GIFTIfield::bScalar(const size_t location) const
{
  if ( dims.size() != 1 )
    throw GiftiException("Requested scalar data from non-scalar GIFTI field");
  if ( dataType != NIFTI_TYPE_UINT8 )
    throw GiftiException("Requested non-char data from char GIFTI field");
  return bdata[location];
}

vector<int> GIFTIfield::iVector(const size_t location) const
{
  if ( dims.size() != 2 )
    throw GiftiException("Requested vector data from non-vector GIFTI field");
  if ( dataType != NIFTI_TYPE_INT32 )
    throw GiftiException("Requested non-int data from int GIFTI field");
  vector<int> test;
  for( int element = 0; element < dims[1]; element++) { 
    //cerr << idata[dims[1]*location+element] << " ";
    test.push_back(idata[dims[1]*location+element]);
  }
  //cerr << endl;
  return test;
}

vector<float> GIFTIfield::fVector(const size_t location) const
{
  if ( dims.size() != 2 )
    throw GiftiException("Requested vector data from non-vector GIFTI field");
  if ( dataType != NIFTI_TYPE_FLOAT32 )
    throw GiftiException("Requested float data from non-float GIFTI field");
  vector<float> test;
  for( int element = 0; element < dims[1]; element++) { 
    //  cerr << fdata[dims[1]*location+element] << " ";
    test.push_back(fdata[dims[1]*location+element]);
  }
  // cerr << endl;
  return test;
}

vector<char> GIFTIfield::bVector(const size_t location) const
{
  if ( dims.size() != 2 )
    throw GiftiException("Requested vector data from non-vector GIFTI field");
  if ( dataType != NIFTI_TYPE_UINT8 )
    throw GiftiException("Requested char data from non-cgar GIFTI field");
  throw GiftiException("bVector not implemented yet!!");
  vector<char> test;
  return test;
}

int GIFTIfield::getDim(const unsigned char dim) const
{
  if ( dim >= dims.size() )
    throw GiftiException("Requested invalid dimension size.");
  return dims[dim];
}

int GIFTIfield::getDataType(void) const
{
  return dataType;
}

int GIFTIfield::getIntent(void) const
{
  return intent;
}


void GIFTIfield::swapOrdering()
{
  vector<int> loopdims(6,1);
  for ( unsigned char dim=0; dim<dims.size(); dim++ )
    loopdims[dim]=dims[dim];
  vector<int> coords(6,0);

  for (coords[0]=0;coords[0]<loopdims[0];coords[0]++)
    for (coords[1]=min(coords[0],loopdims[1]-1);coords[1]<loopdims[1];coords[1]++)
      for (coords[2]=min(coords[1],loopdims[2]-1);coords[2]<loopdims[2];coords[2]++)  
	for (coords[3]=min(coords[2],loopdims[3]-1);coords[3]<loopdims[3];coords[3]++)
	  for (coords[4]=min(coords[3],loopdims[4]-1);coords[4]<loopdims[4];coords[4]++)
	    for (coords[5]=min(coords[4],loopdims[5]-1);coords[5]<loopdims[5];coords[5]++) {
	      vector<int> swapCoords(coords);
	      reverse(swapCoords.begin(), swapCoords.begin()+ dims.size());
	      size_t offset1=(((((coords[5]*loopdims[4]+coords[4])*loopdims[3]+coords[3])*loopdims[2]+coords[2])*loopdims[1]+coords[1])*loopdims[0]+coords[0]);
	      size_t offset2=(((((swapCoords[5]*loopdims[4]+swapCoords[4])*loopdims[3]+swapCoords[3])*loopdims[2]+swapCoords[2])*loopdims[1]+swapCoords[1])*loopdims[0]+swapCoords[0]);
	      for ( unsigned char dim=0; dim<dims.size(); dim++ )
                if ( dataType == NIFTI_TYPE_FLOAT32 )
		  swap(fdata[offset1],fdata[offset2]);
		if ( dataType ==NIFTI_TYPE_INT32 ) {
		  swap(idata[offset1],idata[offset2]);
		}
		if ( dataType == NIFTI_TYPE_UINT8 )
		  swap(bdata[offset1],bdata[offset2]);

	    }
  if ( ordering == GIFTI_IND_ORD_ROW_MAJOR )
    ordering = GIFTI_IND_ORD_COL_MAJOR;
  if ( ordering == GIFTI_IND_ORD_COL_MAJOR )
    ordering = GIFTI_IND_ORD_ROW_MAJOR;
}

void GIFTIfield::report() const
{
  cout << "intent : " << gifti_intent_to_string(intent) << endl;
  cout << "datatype : " << gifti_datatype2str(dataType) << endl;
  cout << "dims: " << dims.size() << " : ";
  for ( unsigned char dim=0; dim<dims.size(); dim++ )
    cout << dims[dim] << " ";
  cout << endl;
  if ( dataType == NIFTI_TYPE_FLOAT32 )
    cout << "Nvals : " << fdata.size() << endl;
  if ( dataType ==NIFTI_TYPE_INT32 )
    cout << "Nvals : " << idata.size() << endl;
  if ( dataType == NIFTI_TYPE_UINT8 )
    cout << "Nvals : " << bdata.size() << endl;
  if ( ordering == GIFTI_IND_ORD_UNDEF )
    cout << "Ordering: Unknown " << endl;
  if ( ordering == GIFTI_IND_ORD_ROW_MAJOR )
    cout << "Ordering: Row (C++)" << endl;
  if ( ordering == GIFTI_IND_ORD_COL_MAJOR )
    cout << "Ordering: Column (Matlab)" << endl;


  if ( intent == NIFTI_INTENT_POINTSET ) { //Print out any coordinate systems
    for ( unsigned int system=0;system<coordSystems.size();system++) {
      cout << coordSystems[system].dataSpace << endl;
      cout << coordSystems[system].transformSpace << endl;
      for(int row=0;row<4;row++) {
	for(int column=0;column<4;column++)
	  cout << coordSystems[system].transform[column+row*4] << " ";
	cout << endl;
      }
    }
  }

  cout << "Number of metadata pairs: " << metaData.size() << endl;
  for ( unsigned int tag=0;tag<metaData.size();tag++)
    cout << metaData[tag].name << " " << metaData[tag].value << endl;
  cout << "Number of extra attributes: " << extraAttributes.size() << endl;
  for ( unsigned int tag=0;tag<extraAttributes.size();tag++)
    cout << extraAttributes[tag].name << " " << extraAttributes[tag].value << endl;
}

void GIFTIfield::printAsSixTensor() const
{
  vector<int> loopdims(6,1);
  for ( unsigned char dim=0; dim<dims.size(); dim++ )
    loopdims[dim]=dims[dim];
  for (int dim0=0;dim0<loopdims[0];dim0++)
    for (int dim1=0;dim1<loopdims[1];dim1++)
      for (int dim2=0;dim2<loopdims[2];dim2++)
	for (int dim3=0;dim3<loopdims[3];dim3++)
	  for (int dim4=0;dim4<loopdims[4];dim4++)
	    for (int dim5=0;dim5<loopdims[5];dim5++)
	    {
	      cerr << dim0 << " "  << dim1 << " " << dim2 << " " << dim3 << " " << dim4 << " " << dim5 << " ";
	      if ( dataType == NIFTI_TYPE_FLOAT32 )
		cerr << fdata[(((((dim5*loopdims[4]+dim4)*loopdims[3]+dim3)*loopdims[2]+dim2)*loopdims[1]+dim1)*loopdims[0]+dim0)] << endl;
	      if ( dataType ==NIFTI_TYPE_INT32 )
		cerr << idata[(((((dim5*loopdims[4]+dim4)*loopdims[3]+dim3)*loopdims[2]+dim2)*loopdims[1]+dim1)*loopdims[0]+dim0)] << endl;
	      if ( dataType == NIFTI_TYPE_UINT8 )
		cerr << bdata[(((((dim5*loopdims[4]+dim4)*loopdims[3]+dim3)*loopdims[2]+dim2)*loopdims[1]+dim1)*loopdims[0]+dim0)] << endl;
	    }
}

GIFTIfield::GIFTIfield(const int inputIntent, const int inputDatatype, const int inputNDim, const int* inputDims, const void* inputData, const int inputOrdering, const vector<GIFTIcoordinateSystem>& inputCoordSystems,const vector<GIFTImeta>& inputMeta,const vector<GIFTImeta>& inputExtraAttributes, const long long inputExternalOffset, const string& inputExternalFilename) : intent(inputIntent), dataType(inputDatatype), coordSystems(inputCoordSystems), metaData(inputMeta), extraAttributes(inputExtraAttributes), externalOffset(inputExternalOffset),  externalFilename(inputExternalFilename)
{
  bdata.clear();
  idata.clear();
  fdata.clear();

  dims.resize(inputNDim);
  dims=vector<int>(inputDims,inputDims+inputNDim);
  long long nvals(dims[0]);
  for(unsigned char dim=1;dim<dims.size();dim++)
    nvals*=dims[dim];
  if ( dataType == NIFTI_TYPE_FLOAT32 )
    fdata=vector<float>((float*)inputData,((float*)inputData)+nvals);
  if ( dataType ==NIFTI_TYPE_INT32 )
    idata=vector<int>((int*)inputData,((int*)inputData)+nvals);
  if ( dataType == NIFTI_TYPE_UINT8 )
    bdata=vector<char>((char*)inputData,((char*)inputData)+nvals);
  ordering=inputOrdering;



}

GIFTIfield::GIFTIfield(const giiDataArray* fieldPointer ) : externalOffset(0),  externalFilename("")
{
  bdata.clear();
  idata.clear();
  fdata.clear();
  metaData.clear();
  coordSystems.clear();
  extraAttributes.clear();

  intent=fieldPointer->intent;
  dataType=fieldPointer->datatype;
  dims.resize(fieldPointer->num_dim);
  dims=vector<int>(fieldPointer->dims,fieldPointer->dims+fieldPointer->num_dim);
  long long nvals(dims[0]);
  for(unsigned char dim=1;dim<dims.size();dim++)
    nvals*=dims[dim];
  if ( dataType == NIFTI_TYPE_FLOAT32 )
    fdata=vector<float>((float*)fieldPointer->data,((float*)fieldPointer->data)+nvals);
  if ( dataType ==NIFTI_TYPE_INT32 )
    idata=vector<int>((int*)fieldPointer->data,((int*)fieldPointer->data)+nvals);
  if ( dataType == NIFTI_TYPE_UINT8 )
    bdata=vector<char>((char*)fieldPointer->data,((char*)fieldPointer->data)+nvals);
  ordering=fieldPointer->ind_ord;

  if ( intent == NIFTI_INTENT_POINTSET ) { //Search for coordinate systems
    int nCS(fieldPointer->numCS);
    for ( int system=0;system<nCS;system++) {
      vector<double> matrix;
      for(int row=0;row<4;row++)
	for(int column=0;column<4;column++) 
	  matrix.push_back(fieldPointer->coordsys[system]->xform[row][column]);
      coordSystems.push_back( GIFTIcoordinateSystem(string(fieldPointer->coordsys[system]->dataspace),string(fieldPointer->coordsys[system]->xformspace),matrix));
    }
  }

  for ( int tag=0;tag<fieldPointer->meta.length;tag++) 
    metaData.push_back(GIFTImeta(fieldPointer->meta.name[tag],fieldPointer->meta.value[tag]));
  for ( int tag=0;tag<fieldPointer->ex_atrs.length;tag++) 
    extraAttributes.push_back(GIFTImeta(fieldPointer->ex_atrs.name[tag],fieldPointer->ex_atrs.value[tag]));

  externalOffset=fieldPointer->ext_offset;
  externalFilename=string( fieldPointer->ext_fname ? fieldPointer->ext_fname : "" );

  if ( ordering == GIFTI_IND_ORD_COL_MAJOR ) {
    //  cout << "Swapping ordering to row" << endl;
    swapOrdering();
  }
}


GIFTImeta::GIFTImeta(const char* inputName, const char* inputValue) : name(inputName ? inputName : ""), value(inputValue ? inputValue : "") {
}

GIFTImeta::GIFTImeta(const string inputName, const string inputValue) : name(inputName), value(inputValue) {
}
