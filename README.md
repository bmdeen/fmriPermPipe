
# fMRIPermPipe

fMRIPermPipe (FPP, version 2.0.0) is a MATLAB-based pipeline for fMRI data analysis, optimized for individual subject analyses and multi-echo data, incorporating multiple neuroimaging tools: FSL, Freesurfer, AFNI, Connectome Workbench, tedana, MSM, and dcm2niix.

Steps include dicom conversion, preprocessing of anatomical and fMRI data, and statistical modeling and analysis of fMRI data. Outputs conform to the [Brain Imaging Data Structure (1.4.1)](https://bids.neuroimaging.io/specification.html) specification. Anatomical preprocessing uses a [Human Connectome Project](https://www.humanconnectome.org/)-like pipeline, yielding an accurate cortical surface reconstruction and surface-based registration to the fsLR atlas. Functional preprocessing uses a simple but powerful approach, including motion parameter estimation, despiking, slice-timing correction, single-shot motion and distortion correction and linear registration to a subject-specific template, and multi-echo ICA denoising. Statistical modeling includes both nonparametric (permutation-based) and parametric (FSL's Improved Linear Model-based) methods.

**NOTE: These scripts are currently under active development, and will be changing frequently. A more stable version is expected by 2022. Currently in progress: wrapper scripts using BIDS naming conventions; json metadata for statistical modeling outputs; Docker container to facilitate installation**

Contents:
1. [Installation](https://github.com/bmdeen/fmriPermPipe/tree/fppipe-v2#installation)
2. [Data requirements](https://github.com/bmdeen/fmriPermPipe/tree/fppipe-v2#data-requirements)
3. [Usage](https://github.com/bmdeen/fmriPermPipe/tree/fppipe-v2#usage)
4. [Processing details](https://github.com/bmdeen/fmriPermPipe/tree/fppipe-v2#processing-details)
5. [References](https://github.com/bmdeen/fmriPermPipe/tree/fppipe-v2#installation)
6. [Licenses](https://github.com/bmdeen/fmriPermPipe/tree/fppipe-v2#licenses)
<br />



## Installation

FPP runs on Mac or Linux (not Windows), on MATLAB 2016b+. MATLAB should be opened through a Bash terminal, to ensure that Bash environment variables are defined. FPP incorporates a number of neuroimaging software packages as dependencies, which must be installed and sourced in the local Bash or MATLAB environment:

Bash:
* [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) v6.0
* [Freesurfer](https://surfer.nmr.mgh.harvard.edu/) v7.1.1
* [AFNI](https://afni.nimh.nih.gov/) v20+
* [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) v1.5+
* [tedana](https://github.com/ME-ICA/tedana) v0.0.9a1
* [MSM_HOCR](https://github.com/ecr05/MSM_HOCR)
* [ASHS](https://www.nitrc.org/projects/ashs/) v2.0.0 (if using MTL segmentation)
* dcm2niix, through [MRICroGL](https://www.nitrc.org/projects/mricrogl/) v1.2+

MATLAB:
* [bids-matlab](https://github.com/bids-standard/bids-matlab)
* [cifti-matlab](https://github.com/Washington-University/cifti-matlab)
* [Freesurfer](https://surfer.nmr.mgh.harvard.edu/) MATLAB tools



## Data requirements

The minimum required data for this pipeline include:
* One or more T1-weighted anatomical images, at .7-1mm resolution, preferably MPRAGE
* T2*-weighted functional images, preferably multi-echo

Recommended data for optimal preprocessing additionally includes:
* One or more T2-weighted anatomical images, at .7-1mm resolution
* Spin-echo "field map" images with opposing phase-encode directions (e.g. AP and PA), matched to functional data in matrix size, number of slices, voxel size, echo spacing, iPAT acceleration properties, and phase-encode direction.

For anatomical images, we recommend collecting 2-3 images at relatively high resolution (.7-.8mm), both T1- and T2-weighted. For functional data, we recommend a multi-echo sequence with 2-2.5mm resolution.

Note: if distortion correction based on spin-echo acquisitions is used, the pipeline assumes that all functional data have a common phase-encode direction.



## Usage

FPP may be used as a full pipeline, or in segments, including dicom conversion, preprocessing, and analysis. For partial usage, raw or preprocessed input data must be formatted based on the BIDS specification. The package also contains a number of standalone functions, including wrappers for FSL, Freesurfer, and Connectome Workbench tools that generate BIDS-formatted JSON metadata for outputs, and other image processing utilities.

The functions are currently written to act on specific input files, such as raw or preprocessed image files, and require a wrapper to be applied automatically across subjects, tasks, and runs. Wrapper scripts based on the BIDS filesystem structure are in progress.

For detailed instructions on function usage, check the help information for each function.


### Dicom conversion

This segment converts dicom images to NIFTI files with JSON metadata in BIDS naming convention, using a wrapper for dcm2niix.

Steps:
`fpp.util.convertDCM2BIDS`


### Preprocessing

This segment preprocesses field map, anatomical, and functional data.

Steps:
1. `fpp.anat.preproc` - Preprocess anatomical images
2. Run Freesurfer's `recon-all` (details below)
3. `fpp.anat.postproc` - Postprocess anatomical images
4. `fpp.fmap.preproc` - Preprocess spin-echo "field map" images (optional)
5. `fpp.func.defineTemplate` - Define functional template image
6. `fpp.func.register` - Register functional template to anatomical
7. `fpp.func.preproc` - Preprocess functional data
8. `fpp.func.removeNuisance` - Remove additional nuisance signals (optional)
9. `fpp.func.surfaceResample` - Resample data to cortical surface (optional)


### Running recon-all

At step 2, Freesurfer's recon-all should be run via command line, using preprocessed T1 and T2 images generated by fpp.anat.preproc as inputs. For .7-.8mm resolution data, the -hires flag should be included, along with an export options file specifying "mris_inflate -n 30". Example call:

`recon-all -s mystudy01 -i ${bidsDir}/derivatives/fpp/sub-mystudy01/anat/sub-mystudy01_space-individual_res-p8_desc-preproc_T1w.nii.gz -all -deface [-T2 ${bidsDir}/derivatives/fpp/sub-mystudy01/anat/sub-mystudy01_space-individual_res-p8_desc-preproc_T1w.nii.gz -T2pial] [-hires -expert ${fppDir}/data/recon-all.opts]`


### Analysis

This segment performs whole-brain and region-of-interest-based analysis of fMRI data.

Permutation-based modeling:
1. `fpp.func.modelPerm` - First-level (within-run) permutation-based analysis
2. `fpp.func.model2ndPerm` - Second-level (cross-run, within subject) permutation-based analysis

Parametric (FILM-based) modeling:
1. `fpp.func.modelFeat` - First-level (within-run) FILM-based analysis
2. `fpp.func.model2ndFeat` - Second-level (cross-run, within subject) FILM-based analysis

Region-of-interest analysis:
`fpp.func.roiExtract` - Extract responses from functionally defined ROIs


### Wrappers

FPP includes wrappers for FSL, Freesurfer, and Connectome Workbench tools, which produce BIDS-formatted JSON metadata for the output image.

FSL:
* `fpp.fsl.aff2Rigid` - wrapper for aff2rigid
* `fpp.fsl.concatWarp` - wrapper for convertwarp
* `fpp.fsl.concatXfm` - wrapper for convert\_xfm -concat
* `fpp.fsl.flirt` - wrapper for flirt
* `fpp.fsl.fnirt` - wrapper for fnirt
* `fpp.fsl.invertWarp` - wrapper for invwarp
* `fpp.fsl.invertXfm` - wrapper for convert\_xfm -inverse
* `fpp.fsl.maths` - wrapper for fslmaths
* `fpp.fsl.moveImage` - wrapper for flirt -applyxfm and applywarp
* `fpp.fsl.robustFOV` - wrapper for robustfov

Freesurfer:
* `fpp.fs.mriConvert` - wrapper for mri\_convert
* `fpp.fs.mrisConvert` - wrapper for mris\_convert

Connectome Workbench:
* `fpp.wb.command` - wrapper for wb\_command


### Utilities

FPP also includes a variety of image processing, statistics, and visualization utilities.

Data I/O:
* `fpp.util.readDataMatrix` - read NIFTI/CIFTI data as coordinate by D matrix (wrapper for MRIRead and cifti\_read)
* `fpp.util.writeDataMatrix` - write coordinate by D matrix to NIFTI/CIFTI file  (wrapper for MRIWrite and cifti\_write)
* `fpp.util.fileParts` - wrapper for MATLAB's fileparts, accommodates NIFTI/GIFTI/CIFTI file extensions

Image processing:
* `fpp.util.checkMRIProperty` - determine properties of MRI images from data or metadata (TR, TE, # vols, dimensions, voxel size, phase-encode dir, slice-encode dir, slice timing)
* `fpp.util.getImageOrientation` - determine MRI image orientation (wrapper for fslval)
* `fpp.util.smoothInMask` - smooth a 3D/4D image within a binary mask
* `fpp.util.firFilter` - apply a zero-phase Hamming-window based FIR filter to a matrix of time series
* `fpp.util.mriFilter` - apply firFilter to a NIFTI or CIFTI dataset
* `fpp.util.reorientToStd` - reorient volumetric MRI data to RAS/LAS via rigid rotation, including orientation-related JSON metadata fields (wrapper for fslreorient2std)
* `fpp.util.tsnrMap` - compute tSNR of fMRI dataset (NIFTI)
* `fpp.util.label2ROI` - convert discrete segmentation image to ROI comprised of multiple subregions (NIFTI/CIFTI)
* `fpp.util.fpp.util.defineROI` - generate ROI from search space and statistical map, based on percentage and/or threshold criteria (NIFTI/CIFTI)
* `fpp.util.sphericalROI` - generate spherical volumetric ROI
* `fpp.util.convertToGii` - convert Freesurfer format surface file to GIFTI format with JSON metadata

Statistics:
* `fpp.util.convertFtoT` - convert F-statistic map from 1D contrast to t-statistic map (NIFTI)
* `fpp.util.convertTtoZ` - convert T-statistic map to Z-statistics (NIFTI/CIFTI)

Visualization:
* `fpp.util.carpetPlot` - generate carpet plot of fMRI data and nuisance metrics
* `fpp.util.barColor` - generate bar plot with colored bars and overlaid data points



## Processing details

Detailed descriptions of processing steps in the pipeline are provided below.


### Anatomical pre/postprocessing

The anatomical pipeline uses a Human Connectome Project (HCP)-like approach to generate a highly accurate cortical surface reconstruction from high-resolution T1- and T2-weighted anatomical images. The specific processing steps are adapted from the HCP's PreFreesurer and PostFreesurfer pipelines.

In preprocessing (`fpp.anat.preproc`), anatomical images are registered to one another, averaged, rigidly aligned with MNI152NLin6Asym space, and bias-corrected. This generates a high-resolution anatomical template space for each individual subject.

Optionally, if high-resolution coronal images of the medial temporal lobe were acquired, `fpp.anat.preprocCoronal` can be used to preprocess these images, and perform anatomical segmentation of hippocampal/parahippocampal structures using Automatic Segmentation of Hippocampal Subfields (ASHS). To perform ASHS, the [UPenn PMC Atlas data](https://www.nitrc.org/frs/?group_id=370) (file 	
ashs\_atlas\_upennpmc\_20170810.tar) must be added. The downloaded folder 	
ashs\_atlas\_upennpmc\_20170810 should be placed in the data folder within the FPP script directory. The [ASHS software](https://www.nitrc.org/projects/ashs/) (v2.0.0) must also be downloaded and sourced.

After preprocessing, a cortical reconstruction is generated using Freesurfer's `recon-all`, with the `-hires` flag to accomodate submillimeter resolution when needed.

In postprocessing (`fpp.anat.postproc`), the following steps are performed:
1. Conversion of Freesurfer volumetric files to anatomical template space
2. Mask definition: brain, gm, wm, csf
3. Conversion of Freesurfer surface, metric, and label files to GIFTI/CIFTI format
4. Generation of midthickness, inflated, and very inflated surfaces
5. Surface-based registration to fsLR space, using Multimodal Surface Matching
6. Generation of .spec files for wb\_view
7. Resampling of HCP1200 atlas files and parcellations to Freesurfer native resolution (optional)

In order to perform step 7, the [HCP multimodal parcellation data](https://balsa.wustl.edu/study/show/RVVG) must be added. The downloaded folder should be renamed `HCP_S1200_Atlas` and placed in the data folder within the FPP script directory.

### Functional preprocessing

The fMRI preprocessing pipeline is optimized to boost temporal signal-to-noise ratio and diminish the impact of non-neural signals, while minimizing steps involving spatial smoothing or interpolation, to maintain high spatial resolution. The pipeline produces preprocessed volumetric data, in a subject-specific functional or anatomical template space. Data can be resampled to the cortical surface using `fpp.func.surfaceResample`.

Preprocessing steps:
1. Motion parameter estimation (FSL MCFLIRT)
2. Despiking (AFNI 3dDespike)
3. Slice timing correction (FSL slicetimer)
4. One-shot motion/distortion correction and template registration (FSL MCFLIRT, TOPUP, FLIRT)
5. Multi-echo ICA (tedana) and brain masking
6. Temporal filtering (optional, not recommended if using tedana)
7. Intensity normalization
8. Nuisance signal extraction (wm/csf/global mean, aCompCorr, DVARS)
9. Volumetric spatial smoothing (optional, not recommended)

The order of steps 1-4 is intended to robustly handle interactions between head motion and slice timing. Slice timing correction (STC) involves temporal interpolation to correct for differences in data acquisition time across slices (e.g. Parker & Razlighi, 2019). Logically, this should be applied before head motion correction (MC), because head motion can change the slice that a certain part of the brain is covered by, such that the appropriate time shift for a given part of the brain changes. However, transient head movements can lead to artefactual spikes in MR signal intensity, and applying STC before MC incurs a risk of spreading this artefact to adjacent time points via temporal interpolation. To avoid this risk, and to generally diminish the impact of signal spikes, we implement despiking prior to STC. Despiking (using AFNI's 3dDespike) removes voxelwise signal outliers by replacing them with values from a smooth function fit to the full time series. Applying despiking and/or STC prior to MC can yield inaccurate head motion estimates (Power et al., 2017). Therefore, we first estimate head motion parameters, then apply despiking and STC, and then apply head motion correction. Note that this procedure still leaves artefactual signal associated with head movements in the data. Such signals will be further mitigated by 1) applying multi-echo ICA; 2) using scrubbing (removal of high-motion time points) during data analysis. Additionally, global signal removal may be used when appropriate with `fpp.func.removeNuisance`, to diminish global effects of respiratory variation (Power et al. 2018).

Steps involving spatial interpolation (motion correction, distortion correction, and template registration) are combined into a single resampling step, such that only one spatial interpolation is performed, to minimize spatial blurring. Distortion correction is performed using blip-up blip-down spin echo acquisitions, using FSL's topup tool. One of two subject-specific template spaces may be chosen: 1) a functional template, defined by a native functional image; or 2) an anatomical template, defined by the subject's anatomical image, rotated and resampled to align with 2mm MNI152NLin6Asym space. Data are not volumetrically resampled to stereotactic space. As an alternative method for transforming data to a common space across participants, we recommend using fpp.func.surfaceResample to convert data to a surface-based atlas space, such as fsLR (Coalson et al. 2018).

Multi-echo ICA is implemented by default for multi-echo data. ME-ICA takes advantage of multi-echo data to specifically remove non-BOLD components from the data. This technique provides a powerful method to boost tSNR, without removing signal of interest. When using ME-ICA, it's important to check the component report, to ensure that no BOLD-like components have been removed. See more details at [the tedana documentation](https://tedana.readthedocs.io/en/stable/).

Temporal filtering is not used by default, because it is largely redundant with multi-echo ICA. If ME-ICA is not used, we recommend high-pass filtering with a cutoff of .01 Hz to mitigate the influence of temporal drift.

At step 8, several nuisance signals are computed from the data, but not removed. These include mean signal from white matter, CSF, and the whole brain; anatomical CompCorr regressors, defined as the top 5 principal components of time series from eroded white matter and CSF segments (Behzadi et al. 2007); and motion-related measures including framewise displacement, DVARS, and standardized DVARS (Power et al. 2012, Nichols 2017). If desired, these signals may be removed from the data via linear regression using `fpp.func.removeNuisance`, or included as nuisance covariates at the statistical modeling step in `fpp.func.modelPerm` or `fpp.func.modelFeat`. Nuisance signals are written to a BIDS-style confounds.tsv file. (Note: we only recommend using aCompCorr signals as nuisance covariates if ME-ICA was not used. In this case, these regressors may also supercede the need for temporal filtering.)

Additionally, pairs of time points with a large head movement between them are flagged as outliers, specified in an outliers.tsv file. Motion cutoffs can be based on framewise displacement (FD), or total translation or rotation. The default behavior is to remove pairs of time points with >.5mm FD between them. We recommend using a more stringent criteria (e.g. >.25mm) with data used for functional connectivity analyses. By default, time points defined in an outliers.tsv file will be removed at the analysis stage in fpp.func.modelPerm or fpp.func.modelFeat.

Volumetric spatial smoothing is not used by default, and not recommended. If spatial smoothing is desired, we recommend using surface-based smoothing functionality in `fpp.func.surfaceResample`, to ensure that smoothing respects the geometry of the cortical surface.


### Functional analysis

Two methods are provided for computing whole-brain General Linear Model-based statistics: nonparametric permutation-based analysis, and parametric analysis using FSL's Improved Linear Model.

In computing individual-subject fMRI time-series statistics, the issue of temporal autocorrelation must be addressed. Autocorrelation implies that individual fMRI time points are not independent measurements, violating assumptions of ordinary least squares. This can be addressed in one of two ways: 1) using nonparametric statistics that avoid any assumptions about the structure of temporal autocorrelation; 2) parametrically modeling temporal autocorrelation structure, and using this model to prewhiten linear model residuals before computing statistics. When using a parametric approach, the specific method used to model autocorrelation is also important to consider. While AR(1)-based approaches to modeling temporal autocorrelation have been found to underestimate parameter variance (Eklund et al. 2012), more robust estimation procedures used by FSL's improved linear model have been argued to yield largely unbiased variance estimates (Woolrich et al. 2001).

FPP offers implementations of both types of statistical model: a nonparametric option (`fpp.func.modelPerm` and `fpp.func.model2ndPerm`), and a parametric option, using FSL's Improved Linear Model (`fpp.func.modelFeat` and `fpp.func.model2ndFeat`). When using the parametric option, we recommend validating the statistical modeling by checking null distributions of test statistics by evaluating models in resting-state data, with no task-driven responses.

After running second-level modeling (within-subject, cross-run), correction for multiple comparisons should be implemented before evaluating whole-brain statistical maps. Multiple comparison correction is not currently built into 2nd-level modeling scripts. However, false discovery rate-based correction can be implemented using `fpp.func.analysis.fdrCorrect`.

**Details on permutation-based statistical inference.** `fpp.func.modelPerm` computes within-subject statistics using a permutation test. Specifically, on each iteration (of 5,000 total, by default) the order of blocks within an experiment is randomly permuted, and voxelwise contrast values are computed, in order to build a null distribution of contrast values at each voxel. As with any permutation test, this analysis relies on the assumption of exchangeability (in this case, of block orders): that under the null hypothesis of no condition differences, the distribution of the test statistic is identical regardless of block order. The best way to ensure that this criterion is satisfied is to build it into the design of the experiment by using a randomization scheme to determine block order. For more details of the logic of permutation-based fMRI analysis, see Nichols and Holmes (2001).



## References

FPP uses tools from numerous neuroimaging software packages, choosing the best tool available for each processing step. The below papers should be cited when the corresponding processing steps have been used.

_Dicom conversion_
* Li X, Morgan PS, Ashburner J, Smith J, Rorden C (2016) The first step for neuroimaging data analysis: DICOM to NIfTI conversion. J Neurosci Methods 264: 47-56. doi: 10.1016/j.jneumeth.2016.03.001

_Anatomical pipeline_
* Glasser MF, Sotiropoulos SN, Wilson JA, Coalson TS, Fischl B, Andersson JSL, Xu J, Jbabdi S, Webster M, Polimeni JR, Van Essen DC, Jenkinson M, Wu-Minn HCP Consortium (2013) The minimal preprocessing pipelines for the Human Connectome Project. NeuroImage 80: 105-124. doi: 10.1016/j.neuroimage.2013.04.127
* Robinson EC, Jbabdi S, Glasser MF, Andersson J, Burgess GC, Harms MP, Smith SM, Van Essen DC, Jenkinson M (2014) MSM: a new flexible framework for Multimodal Surface Matching. NeuroImage 100: 414-26. doi: 10.1016/j.neuroimage.2014.05.06
* Robinson EC, Garcia K, Glasser MF, Chen Z, Coalson TS, Makropoulos A, Bozek J, Wright R, Schuh A, Webster M, Hutter J (2018) Multimodal surface matching with higher-order smoothness constraints. NeuroImage 167: 453-65. doi: 10.1016/j.neuroimage.2017.10.037

_Freesurfer recon-all_
* Dale AM, Fischl B, Sereno MI (1999) Cortical surface-based analysis. I. Segmentation and surface reconstruction. NeuroImage 9: 179-194. doi: 10.1006/nimg.1998.0395

_Functional registration: boundary-based registration_
* Greve DN, Fischl B (2009) Accurate and robust brain image alignment using boundary-based registration. NeuroImage 48(1): 63-72. doi: 10.1016/j.neuroimage.2009.06.060

_Functional preprocessing: motion correction, linear registration_
* Jenkinson M, Bannister P, Brady JM, Smith, SM (2002) Improved optimisation for the robust and accurate linear registration and motion correction of brain images. NeuroImage 17(2): 825-841. doi: 10.1016/s1053-8119(02)91132-8

_Functional preprocessing: distortion correction_
* Andersson JLR, Skare S, Ashburner J (2003) How to correct susceptibility distortions in spin-echo echo-planar images: application to diffusion tensor imaging. NeuroImage 20(2): 870-888. doi: 10.1016/S1053-8119(03)00336-7

_Functional preprocessing: multi-echo ICA_
* DuPre E, Salo T, Markello R, Kundu P, Whitaker K, & Handwerker D (2020) ME-ICA/tedana: 0.0.9a (Version 0.0.9a). Zenodo. doi: 10.5281/zenodo.3786890
* Kundu P, Inati SJ, Evans JW, Luh WM, Bandettini PA (2011) Differentiating BOLD and non-BOLD signals in fMRI time series using multi-echo EPI. NeuroImage, 60: 1759-1770. doi: 10.1016/j.neuroimage.2011.12.028

_Functional preprocessing: aCompCor nuisance signals_
* Behzadi Y, Restom K, Liau J, Liu TT (2007) A component based noise correction method (CompCor) for BOLD and perfusion based fMRI. NeuroImage 37(1):9 0-101. doi: 10.1016/j.neuroimage.2007.04.042

_Functional preprocessing: standardized DVARS nuisance signal_
* Nichols TE (2017) Notes on creating a standardized version of DVARS. arXiv: 1704.01469.

_Functional preprocessing: scrubbing artifact time points_
* Power JD, Barnes KA, Snyder AZ, Schlaggar BL, Petersen SE (2012) Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. NeuroImage 59(3): 2142-54. doi: 10.1016/j.neuroimage.2011.10.018
* Power JD, Mitra A, Laumann TO, Snyder AZ, Schlaggar BL, Petersen SE (2014) Methods to detect, characterize, and remove motion artifact in resting state fMRI. NeuroImage 84: 320-41. doi: 10.1016/j.neuroimage.2013.08.048

_Functional preprocessing: carpet plot_
* Power JD (2017) A simple but useful way to assess fMRI scan qualities. NeuroImage 154: 150-8. doi: 10.1016/j.neuroimage.2016.08.009

_Functional analysis: FILM modeling_
* Woolrich MW, Ripley BD, Brady M, Smith SM (2001) Temporal autocorrelation in univariate linear modeling of fMRI data. NeuroImage 14(6): 1370–1386. doi: 10.1006/nimg.2001.0931

_Brain imaging data structure_
* Gorgolewski KJ, Auer T, Calhoun VD, Craddock RC, Das S, Duff EP, Flandin G, Ghosh SS, Glatard T, Halchenko YO, Handwerker DA (2016) The brain imaging data structure, a format for organizing and describing outputs of neuroimaging experiments. Scientific Data 3(1):1-9. doi: 10.1038/sdata.2016.44

_Additional References_
* Coalson TS, Van Essen DC, Glasser MF (2018) The impact of traditional neuroimaging methods on the spatial localization of cortical areas. PNAS 115(27): E6356-65. doi: 10.1073/pnas.1801582115
* Eklund A, Andersson M, Josephson C, Johannesson M, Knutsson H (2012) Does parametric fMRI analysis with SPM yield valid results?—An empirical study of 1484 rest datasets. NeuroImage 61(3): 565-78. doi: 10.1016/j.neuroimage.2012.03.093
* Eklund A, Nichols TE, Knutsson H (2016) Cluster failure: Why fMRI inferences for spatial extent have inflated false-positive rates. PNAS 13(28):7900-5. doi: 10.1073/pnas.1602413113
* Nichols TE, Holmes AP (2002) Nonparametric permutation tests for functional neuroimaging: a primer with examples. Human Brain Mapping 15(1): 1-25. doi: 10.1002/hbm.1058
* Parker DB, Razlighi QR (2019) The benefit of slice timing correction in common fMRI preprocessing pipelines. Frontiers in Neuroscience 13:821. doi: 10.3389/fnins.2019.00821
* Power JD, Plitt M, Kundu P, Bandettini PA, Martin A (2017) Temporal interpolation alters motion in fMRI scans: Magnitudes and consequences for artifact detection. PLOS One 12(9): e0182939. doi: 10.1371/journal.pone.0182939
* Power JD, Plitt M, Gotts SJ, Kundu P, Voon V, Bandettini PA, Martin A (2018) Ridding fMRI data of motion-related influences: Removal of signals with distinct spatial and physical bases in multiecho data. PNAS 115(9): E2105-14. doi: 10.1073/pnas.1720985115



## Version History

Version 1 (2012-2016): early version of pipeline, with more basic preprocessing and analysis functionality.

Version 2 (2019-current): updated version, incorporating optimized fMRI preprocessing, HCP-like anatomical pipeline, and BIDS compatibility.



## Licenses

fMRIPermPipe relies on numerous software packages, each of which has its own license. Additionally, the data folder contains files specifying template images from the Human Connectome Project (`*_space-fsLR_*` and `*_space-fsaverage_*`), and the Montreal Neurological Institute (`*_space-MNI152NLin6Asym_*`). These files have their own licenses, which are reproduced below.


### fMRIPermPipe License

Copyright (c) 2012-2021 Benjamin Deen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


### Human Connectome Project License

Copyright (c) 2011-2019 [The Human Connectome Project][HCP] and [The Connectome Coordination Facility][CCF]

Redistribution and use in source and binary forms, with or without modification,
is permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, 
  this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions, and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* The names of Washington University in St. Louis, the University of Minnesota,
  Oxford University, the Human Connectome Project, or any contributors
  to this software may *not* be used to endorse or promote products derived
  from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


### MNI152 License

Copyright (C) 1993–2009 Louis Collins, McConnell Brain Imaging Centre, Montreal Neurological Institute, McGill University

Permission to use, copy, modify, and distribute this software and its documentation for any purpose and without fee is hereby granted, provided that the above copyright notice appear in all copies. The authors and McGill University make no representations about the suitability of this software for any purpose. It is provided “as is” without express or implied warranty. The authors are not responsible for any data loss, equipment damage, property loss, or injury to subjects or patients resulting from the use or misuse of this software package.


<!-- References -->

[HCP]: https://www.humanconnectome.org
[CCF]: https://www.humanconnectome.org
