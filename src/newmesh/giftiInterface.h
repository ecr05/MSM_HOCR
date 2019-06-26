/*  giftiInterface.h
    Matthew Webster (FMRIB)
    Copyright (C) 2013 University of Oxford  */
/*  CCOPYRIGHT  */

#include <stdlib.h>
#include <map>
#include <vector>
#include <iostream>
extern "C" {
#include <giftiio/gifti_io.h>
}

 struct GiftiException : public std::exception
{
  std::string errorMessage;
  GiftiException(const std::string& error) : errorMessage(error) {}
  ~GiftiException() throw() {}
  const char* what() const throw() { return errorMessage.c_str(); }
};

struct GIFTImeta {
  GIFTImeta(const char* inputName, const char* inputValue);
  GIFTImeta(const std::string inputName, const std::string inputValue);

  std::string name;
  std::string value;
};

struct GIFTIcoordinateSystem {
  GIFTIcoordinateSystem(const std::string& inputDataSpace, const std::string& inputTransformSpace, const std::vector<double>& inputTransform) : dataSpace(inputDataSpace), transformSpace(inputTransformSpace), transform(inputTransform){};
  std::string dataSpace;
  std::string transformSpace;
  std::vector<double> transform;
};

class GIFTIfield { //NOT templated by field type as I want to store all fields in single container
  friend class GIFTIwrapper;
  int intent;
  std::vector<int> dims; // size of this is number of dims
  int dataType;
  int ordering;
  std::vector<char> bdata; //Stores data for byte fields
  std::vector<int>  idata; //Stores data for int fields
  std::vector<float> fdata; //Stores data for float fields;
  void swapOrdering();

  std::vector<GIFTIcoordinateSystem> coordSystems;
  std::vector<GIFTImeta> metaData;
  std::vector<GIFTImeta> extraAttributes; //Any other attributes that the tag has - note these are _not_ saved by the GIFTI library despite being in the struct!!!
  long long externalOffset;
  std::string externalFilename;

public:

  GIFTIfield(const int intent, const int datatype, const int nDim,const int* dims,const void* data, const int ordering, const std::vector<GIFTIcoordinateSystem>& inputCoordSystems=std::vector<GIFTIcoordinateSystem>(), const std::vector<GIFTImeta>& inputMeta=std::vector<GIFTImeta>(), const std::vector<GIFTImeta>& inputExtraAttribures=std::vector<GIFTImeta>(), const long long inputExternalOffset=0, const std::string& inputExternalFilename=std::string());
  GIFTIfield(const giiDataArray* fieldPointer );
  void report() const;
  void printAsSixTensor() const;
  int getDim(const unsigned char dim) const;
  int getDataType(void) const;
  int getIntent(void) const;

  std::vector<GIFTIcoordinateSystem> getCoordSystems() const {return coordSystems;};
  std::vector<GIFTImeta> getMetaData() const {return metaData;};
  std::vector<GIFTImeta> getExtraAttributes() const {return extraAttributes;};
  long long getExternalOffset() const {return externalOffset;};
  std::string getExternalFilename() const {return externalFilename;};


  float fScalar(const size_t location) const; //Returns a uni-dimensional float fields value at location, throws exception for other types
  int iScalar(const size_t location) const; //Returns a uni-dimensional  int fields value at location, throws exception for other types
  char bScalar(const size_t location) const; //Returns a uni-dimensional byte fields value at location, throws exception for other types
  float asFscalar(const size_t location) const; //Returns a uni-dimensional float fields value at location, throws exception for other types
  int asIscalar(const size_t location) const; //Returns a uni-dimensional  int fields value at location, throws exception for other types
  char asBscalar(const size_t location) const; //Returns a uni-dimensional byte fields value at location, throws exception for other types
  void setFScalar(const size_t location, const float value);
  void setIScalar(const size_t location, const int value);
  void setBScalar(const size_t location, const char value );
  std::vector<float> fVector(const size_t location) const;  //Returns 2-dimensional vector field at location, throws exception for other types
  std::vector<int> iVector(const size_t location) const;  //Returns 2-dimensional vector field at location, throws exception for other types
  std::vector<char> bVector(const size_t location) const;  //Returns 2-dimensional vector field at location, throws exception for other types
};

class GIFTIlabel {
public:
  std::vector<float> RGBA;
  std::string name;
};

class GIFTIwrapper {
public:
  std::vector<GIFTIfield> allFields;
  std::vector<GIFTImeta> metaData;
  std::vector<GIFTImeta> extraAttributes; //Any other attributes that the tag has
  std::map<int,GIFTIlabel> GIFTIlabels;
  int readGIFTI(const std::string filename);
  int writeGIFTI(const std::string filename, int encoding=GIFTI_ENCODING_ASCII) const;
  std::vector<GIFTIfield> returnSurfaceFields();
  std::vector<GIFTIfield> returnNonSurfaceFields();
  void report() const;
};
