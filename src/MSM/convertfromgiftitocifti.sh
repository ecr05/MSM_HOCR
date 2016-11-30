#!/bin/bash
#   Emma Robinson, FMRIB Image Analysis Group
#   Copyright (C) 2012 University of Oxford
#
#   SHCOPYRIGHT

Usage() {
  echo "convertfromgiftitocifti <Subject> <dirpath>  <FeatType> <OutputRegName> <palettemax> <palettemin> "
  exit 1 
}

[ "$6" = "" ] && Usage

Subject=$1
dirname=$2
FeatType=$3
OutputRegName=$4
palettemax=$5
palettemin=$6

 wb_command -cifti-create-dense-timeseries ${dirname}/${Subject}.${FeatType}.${OutputRegName}.dtseries.nii -left-metric ${dirname}/${Subject}.L.${FeatType}.${OutputRegName}.func.gii -right-metric ${dirname}/${Subject}.R.${FeatType}.${OutputRegName}.func.gii
 
 wb_command -cifti-convert-to-scalar ${dirname}/${Subject}.${FeatType}.${OutputRegName}.dtseries.nii ROW ${dirname}/${Subject}.${FeatType}.${OutputRegName}.dscalar.nii 
 wb_command -set-map-name  ${dirname}/${Subject}.${FeatType}.${OutputRegName}.dscalar.nii 1 /${Subject}.${FeatType}.${OutputRegName}

 wb_command -cifti-palette  ${dirname}/${Subject}.${FeatType}.${OutputRegName}.dscalar.nii  MODE_USER_SCALE  ${dirname}/${Subject}.${FeatType}.${OutputRegName}.dscalar.nii -pos-user 0 $palettemax -neg-user 0 $palettemin -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false
