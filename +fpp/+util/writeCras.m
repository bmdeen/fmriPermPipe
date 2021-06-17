
% fpp.util.writeCras(inputPath,crasPath)
%
% Function to write c_ras matrix from NIFTI header. Used to add intercept
% to Freesurfer surface coordinates when converting to GIFTI.

function writeCras(inputPath,crasPath)

% Extract intercept (c_ras) info
[~,crasX] = fpp.util.system(['MatrixXYZ=`mri_info --cras ' inputPath '`;'...
	' echo ${MatrixXYZ} | awk ''{print $1;}''']);
crasX = str2num(strtrim(crasX));
[~,crasY] = fpp.util.system(['MatrixXYZ=`mri_info --cras ' inputPath '`;'...
	' echo ${MatrixXYZ} | awk ''{print $2;}''']);
crasY = str2num(strtrim(crasY));
[~,crasZ] = fpp.util.system(['MatrixXYZ=`mri_info --cras ' inputPath '`;'...
	' echo ${MatrixXYZ} | awk ''{print $3;}''']);
crasZ = str2num(strtrim(crasZ));

% Define matrix to add intercept
crasMat = eye(4);
crasMat(1:3,4) = [crasX crasY crasZ]';

% Write output file
fid = fopen(crasPath,'wt');
fprintf(fid,'%f %f %f %f\n',crasMat');
fclose(fid);

end