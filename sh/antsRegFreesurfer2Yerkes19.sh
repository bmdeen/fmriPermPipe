#!/bin/bash
#

if [ $# = 0 ] 
then
  echo " "
  echo "======================================================= "
  echo " This script performs nonlinear volumetric registration"
  echo " of macaque anatomical images to the Yerkes19 template"
  echo " space ,using the ANTs toolkit."
  echo " "
  echo " ***WARNING*** The input image MUST BE in Freesurfer-"
  echo " conformed space (fake 1mm resolution, 256x256x256, LIA"
  echo " orientation), and the script assumes that the original"
  echo " image was .5mm resolution. Input should be a whole-head"
  echo " image, not bias-corrected(e.g. mri/rawavg.mgz)."
  echo " "
  echo " The brain mask should be tight (ideally manually edited)"
  echo " as this will determine registration quality at the edges."
  echo " "
  echo " If functional template path and transform files are "
  echo " specified, registrations will additionally be computed"
  echo " between functional space and Yerkes19 space, and"
  echo " parcels will be registered to functional space."
  echo " "
  echo "Usage: antsRegFreesurfer2Yerkes19.sh inputPath maskPath outputDir"
  echo " 	inputPath: path to whole-head image to transform"
  echo "	maskPath: path to brain mask"
  echo " 	outputDir (optional): directory for output files"
  echo "		(Default: INPUT_DIR/freesurfer2yerkes19)"
  echo " "
  echo "Optional arguments:"
  echo " 	-f, --func_template : path to functional template space"
  echo "    -fx, --func_template_xfm: FSL-format (.mat) transform"
  echo "		from functional template to freesurfer space."
  echo " "
  echo "======================================================= "
  echo " "
  return 0
fi

origResolution=.5			# Original anatomical resolution, in mm
origResolutionInverse=`echo "scale=4 ; 1 / $origResolution" | bc`
yerkesDir=${BASH_SOURCE%/*}/Yerkes19

# Check for ANTs, c3d, and Freesurfer paths
if [ -z `which antsRegistration` ]; then
	echo "ERROR: ANTs directory not included in PATH."
	return
fi
if [ -z `which c3d_affine_tool` ]; then
	echo "ERROR: c3d directory not included in PATH."
	return
fi
if [ -z `which tkregister2` ]; then
	echo "ERROR: Freesurfer must be sourced."
	return
fi

baseArguments=()
funcTemplatePath=""
funcTemplateXfmPath=""

# Loop through arguments
for arg in "$@"
do
    case $arg in
        -f|--func_template)
        funcTemplatePath="$2"
        shift # Remove argument name from processing
        ;;
        -fx|--func_template_xfm)
        funcTemplateXfmPath="$2"
        shift # Remove argument name from processing
        ;;
        *)
        baseArguments+=("$1")
        shift # Remove generic arguments from processing
        ;;
    esac
done

# Check if required input argument wasn't specified
if [ ${#baseArguments[@]} -eq 0 ]; then
    echo "ERROR: input image path argument is required."
    return
fi
inputPath=${baseArguments[0]}

# Check if required mask argument wasn't specified
if [ ${#baseArguments[@]} -eq 1 ]; then
    echo "ERROR: mask image path argument is required."
    return
fi
maskPath=${baseArguments[1]}

# Define output directory, from argument if specified
if [ -z ${baseArguments[2]} ]; then
	outputDir=$(dirname "${inputPath}")/reg_freesurfer2yerkes19
else
	outputDir=${baseArguments[2]}
fi

# If functional template is specified, transform to freesurfer space must be specified, and vice versa
if [[ ! -z $funcTemplatePath && -z $funcTemplateXfmPath ]]; then
	echo "ERROR: if func_template is specified, func_template_xfm must be defined as well."
	return
fi
if [[ -z $funcTemplatePath && ! -z $funcTemplateXfmPath ]]; then
	echo "ERROR: if func_template_xfm is specified, func_template must be defined as well."
	return
fi

# Check if input files exist, return if not
if [ ! -f $inputPath ]; then
	echo "ERROR: input data file does not exist."
    return
fi
if [ ! -f $maskPath ]; then
	echo "ERROR: mask file does not exist."
    return
fi
if [ ! -z $funcTemplatePath ] && [ ! -f $funcTemplatePath ]; then
	echo "ERROR: functional template does not exist."
    return
fi
if [ ! -z $funcTemplateXfmPath ] && [ ! -f $funcTemplateXfmPath ]; then
	echo "ERROR: functional template xfm does not exist."
    return
fi

# Create output directory, if it doesn't exist
if [ ! -d $outputDir ]; then
	mkdir $outputDir
fi

# Copy input images
cp ${yerkesDir}/Yerkes19_v1.2_T1_brain.nii.gz \
	${outputDir}/yerkes19.nii.gz
# cp ${yerkesDir}/Yerkes19_v1.2_T1.nii.gz \
# 	${outputDir}/yerkes19.nii.gz
mri_convert $inputPath ${outputDir}/freesurfer.nii.gz
mri_convert $maskPath ${outputDir}/freesurferMask.nii.gz

# N4 bias field correction
N4BiasFieldCorrection -d 3 -i ${outputDir}/freesurfer.nii.gz -o \
	${outputDir}/freesurfer.nii.gz

# Apply brain mask
fslmaths ${outputDir}/freesurfer.nii.gz -mul ${outputDir}/freesurferMask.nii.gz ${outputDir}/freesurfer.nii.gz

# Convert to LAS orientation
mri_convert ${outputDir}/freesurfer.nii.gz ${outputDir}/freesurferLAS.nii.gz --out_orientation las

# Convert input image header to original resolution, from fake 1mm Freesurfer-conformed input
hdrInfo=(`fslorient -getqform ${outputDir}/freesurferLAS.nii.gz`)
fslcreatehd 256 256 256 1 ${origResolution} ${origResolution} ${origResolution} 0 ${hdrInfo[3]} \
	${hdrInfo[7]} ${hdrInfo[11]} 16 ${outputDir}/emptyHeader.nii.gz
fslmaths ${outputDir}/emptyHeader.nii.gz -add ${outputDir}/freesurferLAS.nii.gz \
	${outputDir}/freesurferLAS${origResolution}mm.nii.gz
rm -rf ${outputDir}/emptyHeader.nii.gz

# Define registrations between original, LAS, and LAS${origResolution} images
if [ -f ${outputDir}/freesurferLAS.5mm2freesurferLAS_fslformat.mat ]; then
	rm -rf ${outputDir}/freesurferLAS.5mm2freesurferLAS_fslformat.mat
fi
xfmRows=("${origResolutionInverse} 0 0 0" "0 ${origResolutionInverse} 0 0" "0 0 ${origResolutionInverse} 0" "0 0 0 1")
for r in "${xfmRows[@]}"; do
	echo $r >> ${outputDir}/freesurferLAS.5mm2freesurferLAS_fslformat.mat;
done

if [ -f ${outputDir}/freesurferLAS2freesurfer_fslformat.mat ]; then
	rm -rf ${outputDir}/freesurferLAS2freesurfer_fslformat.mat
fi
xfmRows=("1 0 0 0" "0 0 -1 256" "0 1 0 0" "0 0 0 1")
for r in "${xfmRows[@]}"; do
	echo $r >> ${outputDir}/freesurferLAS2freesurfer_fslformat.mat;
done

# DEPRECATED
# Compute registrations between original, LAS, and .5mm images using FLIRT with constraints
# echo "Computing registrations between freesurfer, freesurferLAS, and freesurferLAS${origResolution}mm spaces"
# flirt -in ${outputDir}/freesurferLAS${origResolution}mm.nii.gz -ref ${outputDir}/freesurferLAS.nii.gz -omat \
# 	${outputDir}/freesurferLAS${origResolution}mm2freesurferLAS_fslformat.mat -dof 7 \
# 	-searchrx -2 2 -searchry -2 2 -searchrz -2 2
# flirt -in ${outputDir}/freesurferLAS.nii.gz -ref ${outputDir}/freesurfer.nii.gz -omat \
# 	${outputDir}/freesurferLAS2freesurfer_fslformat.mat -dof 6 \
# 	-searchrx -180 180 -searchry -180 180 -searchrz -180 180

convert_xfm -omat ${outputDir}/freesurferLAS${origResolution}mm2freesurfer_fslformat.mat -concat \
	${outputDir}/freesurferLAS2freesurfer_fslformat.mat \
	${outputDir}/freesurferLAS${origResolution}mm2freesurferLAS_fslformat.mat
convert_xfm -omat ${outputDir}/freesurfer2freesurferLAS${origResolution}mm_fslformat.mat -inverse \
	${outputDir}/freesurferLAS${origResolution}mm2freesurfer_fslformat.mat

# Test transformation
# flirt -in ${outputDir}/freesurfer.nii.gz -ref ${outputDir}/freesurferLAS${origResolution}mm.nii.gz \
# 	-out ${outputDir}/freesurfer2freesurferLAS${origResolution}mm.nii.gz -applyxfm -init \
# 	${outputDir}/freesurfer2freesurferLAS${origResolution}mm_fslformat.mat

# Run ANTs registration
echo "Running ANTs nonlinear registration to template"
. ${BASH_SOURCE%/*}/antsRegWrapper.sh \
	${outputDir}/freesurferLAS${origResolution}mm.nii.gz ${outputDir}/yerkes19.nii.gz \
	${outputDir}/freesurfer2yerkes19_
# Check if output doesn't exist, return if so
if [ ! -f ${outputDir}/freesurfer2yerkes19_Warped.nii.gz ]; then
    return
fi
# Rename ANTs outputs
mv ${outputDir}/freesurfer2yerkes19_affine.mat ${outputDir}/freesurferLAS${origResolution}mm2yerkes19_antsformat.mat
mv ${outputDir}/freesurfer2yerkes19_affine_fslformat.mat \
	${outputDir}/freesurferLAS${origResolution}mm2yerkes19_fslformat.mat
rm -rf ${outputDir}/freesurfer2yerkes19_affine_freesurferformat.dat

# Combine ANTs affine transform with FLIRT-based freesurfer2freesurferLAS${origResolution}mm transform
convert_xfm -omat ${outputDir}/freesurfer2yerkes19_fslformat.mat -concat \
	${outputDir}/freesurferLAS${origResolution}mm2yerkes19_fslformat.mat \
	${outputDir}/freesurfer2freesurferLAS${origResolution}mm_fslformat.mat
convert_xfm -omat ${outputDir}/yerkes19_2freesurfer_fslformat.mat -inverse \
	${outputDir}/freesurfer2yerkes19_fslformat.mat
# Convert to ANTs and Freesurfer formats
c3d_affine_tool -ref ${outputDir}/yerkes19.nii.gz -src ${outputDir}/freesurfer.nii.gz \
	${outputDir}/freesurfer2yerkes19_fslformat.mat -fsl2ras -oitk \
	${outputDir}/freesurfer2yerkes19_antsformat.mat
c3d_affine_tool -src ${outputDir}/yerkes19.nii.gz -ref ${outputDir}/freesurfer.nii.gz \
	${outputDir}/yerkes19_2freesurfer_fslformat.mat -fsl2ras -oitk \
	${outputDir}/yerkes19_2freesurfer_antsformat.mat
tkregister2 --mov ${outputDir}/freesurfer.nii.gz --targ ${outputDir}/yerkes19.nii.gz --fsl \
	${outputDir}/freesurfer2yerkes19_fslformat.mat --reg ${outputDir}/freesurfer2yerkes19_freesurferformat.dat --noedit
tkregister2 --targ ${outputDir}/freesurfer.nii.gz --mov ${outputDir}/yerkes19.nii.gz --fsl \
	${outputDir}/yerkes19_2freesurfer_fslformat.mat --reg ${outputDir}/yerkes19_2freesurfer_freesurferformat.dat --noedit
# If functional template argument is specified, generate transforms to/from this space
if [ ! -z $funcTemplatePath ]; then
	convert_xfm -omat ${outputDir}/funcTemplate2yerkes19_fslformat.mat -concat \
		${outputDir}/freesurfer2yerkes19_fslformat.mat $funcTemplateXfmPath
	convert_xfm -omat ${outputDir}/yerkes19_2funcTemplate_fslformat.mat -inverse \
		${outputDir}/funcTemplate2yerkes19_fslformat.mat
	c3d_affine_tool -ref ${outputDir}/yerkes19.nii.gz -src $funcTemplatePath \
		${outputDir}/funcTemplate2yerkes19_fslformat.mat -fsl2ras -oitk \
		${outputDir}/funcTemplate2yerkes19_antsformat.mat
	c3d_affine_tool -src ${outputDir}/yerkes19.nii.gz -ref $funcTemplatePath \
		${outputDir}/yerkes19_2funcTemplate_fslformat.mat -fsl2ras -oitk \
		${outputDir}/yerkes19_2funcTemplate_antsformat.mat
	tkregister2 --mov $funcTemplatePath --targ ${outputDir}/yerkes19.nii.gz --fsl \
		${outputDir}/funcTemplate2yerkes19_fslformat.mat --reg ${outputDir}/funcTemplate2yerkes19_freesurferformat.dat --noedit
	tkregister2 --targ $funcTemplatePath --mov ${outputDir}/yerkes19.nii.gz --fsl \
		${outputDir}/yerkes19_2funcTemplate_fslformat.mat --reg ${outputDir}/yerkes19_2funcTemplate_freesurferformat.dat --noedit
fi
# Test transformations
# antsApplyTransforms -d 3 -i ${outputDir}/freesurfer.nii.gz -r ${outputDir}/yerkes19.nii.gz -n linear \
# 	-t ${outputDir}/freesurfer2yerkes19_warp.nii.gz -t ${outputDir}/freesurfer2yerkes19_antsformat.mat \
# 	-o ${outputDir}/freesurfer2yerkes19_Warped2.nii.gz
# antsApplyTransforms -d 3 -i ${outputDir}/yerkes19.nii.gz -r ${outputDir}/freesurfer.nii.gz -n linear \
#  	-t ${outputDir}/yerkes19_2freesurfer_antsformat.mat -t ${outputDir}/freesurfer2yerkes19_invwarp.nii.gz \
#  	-o ${outputDir}/yerkes19_2freesurfer_Warped.nii.gz
# antsApplyTransforms -d 3 -i $funcTemplatePath -r ${outputDir}/yerkes19.nii.gz -n linear \
# 	-t ${outputDir}/freesurfer2yerkes19_warp.nii.gz -t ${outputDir}/funcTemplate2yerkes19_antsformat.mat \
# 	-o ${outputDir}/funcTemplate2yerkes19_Warped.nii.gz

# Remove internal registration files
rm -rf ${outputDir}/freesurferLAS.nii.gz ${outputDir}/freesurferLAS${origResolution}mm.nii.gz \
 	${outputDir}/freesurferLAS2freesurfer_fslformat.mat \
 	${outputDir}/freesurferLAS${origResolution}mm2freesurferLAS_fslformat.mat \
 	${outputDir}/freesurferLAS${origResolution}mm2freesurfer_fslformat.mat \
 	${outputDir}/freesurfer2freesurferLAS${origResolution}mm_fslformat.mat \
 	${outputDir}/freesurferLAS${origResolution}mm2yerkes19_antsformat.mat \
 	${outputDir}/freesurferLAS${origResolution}mm2yerkes19_fslformat.mat
 
# Convert parcellations to Freesurfer and funcTemplate spaces
parcPaths=(${yerkesDir}/Yerkes19_Parcellations_v2.nii.gz ${yerkesDir}/D99_atlas_v1.2b.nii.gz)
parcDir=${outputDir}/parcellations_freesurfer
parcDir2=${outputDir}/parcellations_funcTemplate
mkdir $parcDir
if [ ! -z $funcTemplatePath ]; then
	mkdir $parcDir2
fi
for parcPath in "${parcPaths[@]}"
do
	parcName=`basename $parcPath`
	parcPathNew=${outputDir}/${parcName}
	cp $parcPath $parcPathNew
	fslsplit ${parcPathNew} ${parcPathNew//".nii.gz"/}_tmp_
	
	# Transform parcellations to freesurfer space
	for parcVolPath in ${parcPathNew//".nii.gz"/}_tmp_*
	do
		antsApplyTransforms -d 3 --i ${parcVolPath} -r ${outputDir}/freesurfer.nii.gz -n NearestNeighbor \
			-t ${outputDir}/yerkes19_2freesurfer_antsformat.mat -t ${outputDir}/freesurfer2yerkes19_invwarp.nii.gz \
			-o ${parcDir}/$(basename $parcVolPath)
	done
	fslmerge -t ${parcDir}/$parcName ${parcDir}/${parcName//".nii.gz"}_tmp_*
	rm -rf ${parcDir}/${parcName//".nii.gz"}_tmp_*
	
	# Transform parcellations to funcTemplate space
	if [ ! -z $funcTemplatePath ]; then
		for parcVolPath in ${parcPathNew//".nii.gz"/}_tmp_*
		do
			antsApplyTransforms -d 3 --i ${parcVolPath} -r $funcTemplatePath -n NearestNeighbor \
				-t ${outputDir}/yerkes19_2funcTemplate_antsformat.mat -t ${outputDir}/freesurfer2yerkes19_invwarp.nii.gz \
				-o ${parcDir2}/$(basename $parcVolPath)
		done
		fslmerge -t ${parcDir2}/$parcName ${parcDir2}/${parcName//".nii.gz"}_tmp_*
		rm -rf ${parcDir2}/${parcName//".nii.gz"}_tmp_*
	fi
	
	rm -rf $parcPathNew ${parcPathNew//".nii.gz"}_tmp_*
done
