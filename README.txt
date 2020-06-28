
[[ fMRIPermPipe ]]

fMRIPermPipe is a MATLAB-based pipeline for fMRI data analysis, optimized for individual subject analyses.  Within-subject statistics are computed using a permutation test, addressing concerns about the validity of statistics based on parametric approaches to autocorrelation correction (Eklund et al., 2012) and cluster correction (Eklund et al., 2016).

**** NOTE: These scripts are currently undergoing a major overhaul, incorporating BIDS naming conventions, improved preprocessing and registration functionality, and expanded multi-echo denoising functionality...stay tuned ****


The scripts provide a full analysis pipeline, from dicom conversion to statistical modeling, with an emphasis on aggressive denoising and permutation-based statistics at the individual subject level. They may be used as a stand-alone pipeline, or pieces of the code may be useful for writing your own analysis scripts. Denoising strategies include early scrubbing of artifact time points (at the preprocessing phase), and removal of PCA-based noise covariates (anatomical CompCorr, Behzadi et al., 2007). The use of permutation-based statistics makes this pipeline particularly well-suited to studies in which establishing a significant response in individual participants is key, such as dense sampling approaches. These scripts rely heavily on the fMRIB Software Library (FSL) for preprocessing and registration steps, and aspects of the preprocessing stream are modeled after FSL's FEAT.


Highlights
- Permutation-based statistics at the individual level
- Integrates Freesurfer for surface-based registration and parcellation
- Aggressive denoising: built-in artifact detection and nuisance regression
- Intuitive, automated filesystem management
- Incorporates multi-echo data



[DEPENDENCIES]
- MATLAB (tested with 2016b+)
- FSL (requires 5+)
- Freesurfer (tested with 4+)
- BXH/XCEDE tools



[INSTALLATION/SETUP]
1) Ensure that the above dependencies are installed and sourced in your bash environment.
2) Organize your input data (dicoms and experimental timing information), as described below in "Filesystem Setup/Structure."
3) Edit the line of code beginning with "dcm = dir" in convertDCM, to fit your system's dicom naming convention (see comments in the script for more information).
4) Open MATLAB and run the scripts.



[BASIC USAGE]

The scripts should be run in the following order:
1) convertDCM: convert DCM images to NIFTI, and reorient to LAS
2) recon-all (freesurfer): generate surface reconstruction
3) preprocMRI: fMRI preprocessing, including artifact detection/scrubbing, motion correction, skull-stripping, spatial smoothing, intensity normalization, high-pass temporal filtering (optional), and registration to a functional template
4) registerMRI: generates registrations between functional template and anatomical image, as well as Freesurfer space and MNI space
5) modelPerm: permutation-based 1st-level (within-run) modeling, including aCompCorr-based nuisance covariates
6) model2ndPerm: permutation-based 2nd-level (cross-run) modeling.
7) randomize (FSL): permutation-based group-level modeling

All MATLAB scripts are run with the following generic format: function(requiredArg1,requiredArg2,...,'variableArgName1',variableArgValue1,'variableArgName2',variableArgValue2,...). For detailed information on running MATLAB scripts, run "help scriptName". For group-level analysis, FSL's randomize can be used to provide permutation-based statistics. As inputs, use MNI-registered copeN_standard.nii.gz images located in 2nd-level analysis (.gperm) directories. For surface reconstruction, run the following bash command for each subject:

export SUBJECTS_DIR=$STUDYDIR/freesurfer
recon-all -s $SUBJECT -i $STUDY_DIR/$SUBJECT/anat/$ANAT_NAME.nii.gz -all

where $STUDYDIR specifies the study directory ($FMRI_BASE_DIR/$STUDY), $SUBJECT specifies the name of the subject, and $ANAT_NAME specifies the name of the anatomical image.



[FILESYSTEM SETUP/STRUCTURE]

The pipeline implements a rigid, intuitive filesystem structure, allowing simple function calls that hide filesystem details. Prior to running the pipeline, a few steps must be taken to properly organize the input data (dicoms, scan information, and experiment timing information):

1) For a given fMRI study, create the study directory STUDY_DIR. This will contain all of the study's data and analyses.

2) Create the freesurfer directory $STUDY_DIR/freesurfer, which will contain surface reconstructions.

3) For each subject in the study, create the subject directory, SUBJ_DIR=$STUDY_DIR/$SUBJECT, where $SUBJECT is a subject identifier of your choice. This will contain all of that subject's data and analyses.

4) For each scan for a given subject, create a dicom directory, labeled $SUBJ_DIR/dicom* (e.g., dicom, dicom2, dicom_scanA, etc). Place all of the dicom images from that scan directly into the dicom directory (no subdirectories).

5) For each scan for a given subject, create a scanlog file, labeled $SUBJ_DIR/scanlog* (e.g., scanlog, scanlog2, scanlog_scanA, etc). These should be named with the same suffices as the corresponding dicom directory for that scan (they must be consistently ordered by MATLAB's dir function). This file has a three column format, with one row for each run to be analyzed. The first column specifies the acquisition number: the order of this run in the sequence of images produced by the scanner, which should be in the DCM filenames. The second column specifies the experiment name: a label chosen by you that will be used to name files and directories produced by the scripts. Anatomical images should contain the string "anat" within their experiment name; DTI data should contain "dti"; and resting-state data should contain "rest" (all case-insensitive). The third column specifies the run number for a given experiment, which can accumulate across scanlog files (i.e., across scans). A sample scanlog file is provided below:

3		anat			1
4		FaceLoc			1
6		FaceLoc			2
8		FaceLoc			3
10		MainExpt		1
12		MainExpt		2
14		MainExpt		3
16		MainExpt		4
18		Rest2mm			1
20		Rest2mm			2
27		dti_60dir		1

7) For each subject in the study, create the regressor directory, REGR_DIR=$SUBJ_DIR/regressors, which contains files specifying event timing for each experiment, run, and condition. These files must be labeled $SUBJECT.$EXPT.$RUN_NUM.$CONDITION_NUM (e.g.: SUB01.FaceLocalizer.1.1), with no file extension. They are specified in FSL-style 3-column format, with a row for each event, the first column specifying event onset in seconds, second column specifying duration in seconds, and third column specifying boxcar regressor height. A sample regressor file is provided below:

18		18		1
54		18		1
90		18		1
162		18		1
198		18		1
234		18		1

----------

After running the scripts, each $SUBJ_DIR will have the following contents:

- dicom*: raw DCM images
- regressors: text files specifying regressor timing in 3-column format
- func: raw functional data in NIFTI format
- anat: raw anatomical data, and mask/parcellation files in anatomical space, in NIFTI format
- dti: raw DTI data in NIFTI format
- reg: registration files between anatomical/Freesurfer/MNI space
-- target_func (or custom target image name): registration files between target functional image and other spaces
- roi
-- target_func (or custom target image name): regions-of-interest in target functional space
- analysis
-- *.prep: output of preprocMRI
-- *.perm: output of modelPerm
-- *.gperm: output of model2ndPerm



[NOTES ON THEORY AND PRACTICE]

[Permutation-based statistical inference.] These scripts compute within-subject statistics using a permutation test. Specifically, on each iteration (of 5,000 total, by default) the order of blocks within an experiment is randomly permuted, and voxelwise contrast values are computed, in order to build a null distribution of contrast values at each voxel. For multiple comparisons correction across voxels, a cluster extent threshold is used, also determined using a block-order permutation test. As with any permutation test, this analysis relies on the assumption of exchangeability (in this case, of block orders): that under the null hypothesis of no condition differences, the distribution of the test statistic is identical regardless of block order. The best way to ensure that this criterion is satisfied is to build it into the design of the experiment by using a randomization scheme to determine block order. For more details of the logic of permutation-based fMRI analysis, see Nichols and Holmes (2001). Another practical consideration is whether to permute rest blocks along with task blocks (controlled by the permuteRest variable argument in modelPerm, which defaults to true). Strictly speaking, this should only be done if rest blocks were treated identically to task blocks in the randomization scheme determining block order; however, in my experience with blocked designs, running these scripts with rest blocks permuted or not typically gives very similar results for balanced contrasts even when this criterion is not met. Note that in order to compute statistics for unbalanced contrasts (e.g., comparing task versus rest), it is necessary to permute rest blocks.

[Preprocessing options.] By default, preprocMRI runs the following preprocessing steps, in order: artifact detection/scrubbing, motion correction, skull-stripping, spatial smoothing, intensity normalization, and registration to a functional template.  Spatial smoothing can be removed by setting fwhm to 0. Slice-timing correction (STC) is not currently incorporated (will be added soon). High-pass temporal filtering is not included by default, but can be added using the tempFilt option, and will occur after intensity normalization. It is not recommended to use temporal filtering when linear and aCompCorr nuisance regressors are included at the modeling phase (as they are by default), because these procedures are often redundant, and nuisance regression typically outperforms temporal filtering. If temporal filtering is used, the tempFilt option should also be set to 1 when running modelPerm, so that regressors are also temporally filtered.

[Artifact detection/scrubbing.] Head movement can induce highly deleterious artifacts in fMRI data.  Because the strongest effects of head movement on fMRI data are local in time, signal artifacts can be diminished by simply removing (or "scrubbing") these time points from the analysis. preprocMRI implements artifact detection and scrubbing by default, removing pairs of time points that are separated by more than .5mm of total translation or .5 degrees of total rotation, as well as time points for which mean signal across the brain is farther than 3.5 standard deviations from the mean. Time point scrubbing can be removed by setting these thresholds to Inf. Relevant variable arguments: transCutoff, rotCutoff, transSingleAxisCutoff, rotSingleAxisCutoff, stdCutoff, tptsAfter, disdaqs.

[Multi-echo data.] Multi-echo pulse sequences, in which data are acquired at muliple TE values at each time point, can boost SNR and diminish signal dropout by weighting signal from different TEs based on local T2* values, which vary across the brain. preprocMRI incorporates multi-echo data by computing an optimal weighting data from different TEs based on locally estimated T2* values. Relevant variable arguments: multiEcho, teVals, echoesToUse.

[Anatomical CompCorr.] The modeling phase employs nuisance covariates computed using a procedure called anatomical CompCorr or aCompCorr (Behzadi et al., 2007), intended to capture noise variance in fMRI datasets unrelated to neural signal. In this procedure, anatomical masks of white matter (WM) and cerebrospinal fluid (CSF) are defined principal components analysis is performed on time series from these masks; and the top N components (5 by default) are included in time series regressions as nuisance covariates. WM and CSF masks are eroded by one voxel, and PCs are computed using unsmoothed data, in order to ensure that gray matter signal is not included as a result of partial voluming or smoothing.

[Registration and data spaces.] This pipeline computes registration files between the middle functional image of each run (example_func), a user-defined target functional-resolution space used to combine data across runs and potentially across experiments (target_func), a high-resolution anatomical image (highres), that same image in Freesurfer's default orientation and sampling (orig), and MNI152 standard stereotaxic space (standard). Target functional to anatomical resolution is computed using Freesurfer's surface-based registration tool (bbregister), and anatomical to MNI space registration is computed using FSL's nonlinear registration tool (FNIRT). A single anatomical image should be used for both surface reconstruction and registration. The motivation to use a "target functional" is to perform analyses in a space that is similar to native functional space (and not so high resolution as to make computations burdensome), but consistent across runs to enable cross-run analysis. By default, preprocMRI will search for an image called target_func.nii.gz in the $SUBJ_DIR/func directory and use this as the target functional space; if this image doesn't exist, the middle functional image of the first run analyzed will be copied to this path. You may wish to specify your own target functional image. One such scenario would be if you wanted to use different target functional images for different experiments, e.g. if they differ in resolution (although, note that it will require extra effort to compare within-subject results across experiments). To use a custom target functional image, specify the targetName variable argument in preprocMRI, registerMRI, modelPerm, and model2ndPerm (registerMRI can be run multiple times for different target functional images).



[REFERENCES]
- Behzadi Y et al. (2007), "A component based noise correction method (CompCorr) for BOLD and perfusion based fMRI." NeuroImage 37(1), 90-101.
- Eklund A et al. (2016), "Cluster failure: Why fMRI inferences for spatial extent have inflated false-positive rates." PNAS 113(28), 7900-7905.
- Eklund A et al. (2012), "Does parametric fMRI analysis with SPM yield valid results? An empirical study of 1484 rest datasets." NeuroImage 61(3): 565-578.
- Nichols and Holmes (2001), "Nonparametric Permutation Tests For Functional Neuroimaging: A Primer with Examples." Human Brain Mapping 15(1): 1-25.
