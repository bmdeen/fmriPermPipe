#!/bin/bash
#

if [ $# = 0 ] 
then
  echo " "
  echo "======================================================= "
  echo " This script provides a wrapper for antsRegistration,"
  echo " using default numerical parameters, and outputs FSL- "
  echo " and Freesurfer-style .mat and .dat files."
  echo " "
  echo "Usage: antsRegWrapper.sh inputPath referencePath outputStem"
  echo " 	inputPath: path to image to transform"
  echo "	referencePath: path to target (e.g. template) image "
  echo " 	outputStem: filestem for output files"
  echo " "
  echo "Optional arguments:"
  echo " 	-nhm : no histogram matching (use this flag for multimodal registration)"
  echo " 	-m, --masks : mask string, structured as \"[inputMaskPath,referenceMaskPath]\""
  echo "    -r, --initial-moving-transform : initial transform .mat file"
  echo " "
  echo " Notes: for best results, use masked, bias-corrected"
  echo " images, with the same rough orientation (e.g. LAS)."
  echo " "
  echo "======================================================= "
  echo " "
  return 0
fi

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

useHistogramMatching=1
maskString=""
initXfmString=""
baseArguments=()

# Loop through arguments
for arg in "$@"
do
    case $arg in
        -nhm)
        useHistogramMatching=0
        shift # Remove -nhm from processing
        ;;
        -m|--masks)
        maskString="-x $2"
        shift # Remove argument name from processing
        ;;
        -r|--initial-moving-transform)
        initXfmString="--initial-moving-transform $2"
        shift # Remove argument name from processing
        ;;
        *)
        baseArguments+=("$1")
        shift # Remove generic arguments from processing
        ;;
    esac
done

# Check if required input arguments weren't specified
if [ ${#baseArguments[@]} -lt 3 ]; then
    echo "ERROR: antsRegWrapper has three required arguments: input image path, " \
    "reference image path, and output stem."
    return
fi

inputPath=${baseArguments[0]}
referencePath=${baseArguments[1]}
outputStem=${baseArguments[2]}
outputDir=$(dirname "${outputStem}")

# TEST
#initXfmString="-r [${inputPath},${referencePath},1]"

# Check if input files exist, return if not
if [[ ! -f $inputPath || ! -f $referencePath ]]; then
	echo "ERROR: input data files do not exist."
    return
fi

# Create output directory, if it doesn't exist
if [ ! -d $outputDir ]; then
	mkdir $outputDir
fi

# If initial registration wasn't specified, use image center of gravity matching
if [ -z $initXfmString ]; then
	initXfmString="--initial-moving-transform [${referencePath},${inputPath},1]"
fi

# antsRegistration: robust, 3-stage (rigid, affine, nonlinear) image registration.
echo " "
echo " "
echo "antsRegistration --dimensionality 3 --float 0 --output \
[${outputStem},${outputStem}Warped.nii.gz] --interpolation \
Linear --winsorize-image-intensities [0.005,0.995] --use-histogram-matching $useHistogramMatching \
${initXfmString} --transform Rigid[0.1] \
--metric MI[${referencePath},${inputPath},1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
--shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[0.1] \
--metric MI[${referencePath},${inputPath},1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
--shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform SyN[0.1,3,0] \
--metric CC[${referencePath},${inputPath},1,4] --convergence [100x70x50x20,1e-6,10] \
--shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox ${maskString} ${initXfmString}"
echo " "
echo " "
antsRegistration --dimensionality 3 --float 0 --output \
	[${outputStem},${outputStem}Warped.nii.gz] --interpolation \
	Linear --winsorize-image-intensities [0.005,0.995] --use-histogram-matching $useHistogramMatching \
	${initXfmString} --transform Rigid[0.1] \
	--metric MI[${referencePath},${inputPath},1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
	--shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[0.1] \
	--metric MI[${referencePath},${inputPath},1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
	--shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform SyN[0.1,3,0] \
	--metric CC[${referencePath},${inputPath},1,4] --convergence [100x70x50x20,1e-6,10] \
	--shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox ${maskString}

# Check if output doesn't exist, return if so
if [ ! -f ${outputStem}Warped.nii.gz ]; then
	echo "ERROR: antsRegistration failed."
    return
fi

# Rename antsRegistration outputs
mv ${outputStem}1Warp.nii.gz ${outputStem}warp.nii.gz
mv ${outputStem}1InverseWarp.nii.gz ${outputStem}invwarp.nii.gz
mv ${outputStem}0GenericAffine.mat ${outputStem}affine.mat

# Convert ANTs/ITK affine xfm file to FSL-style xfm
c3d_affine_tool -ref ${referencePath} -src ${inputPath} -itk ${outputStem}affine.mat \
-ras2fsl -o ${outputStem}affine_fslformat.mat
# Convert FSL-style affine xfm to Freesurfer-style xfm
tkregister2 --mov ${inputPath} --targ ${referencePath} --fsl ${outputStem}affine_fslformat.mat \
--reg ${outputStem}affine_freesurferformat.dat --noedit
