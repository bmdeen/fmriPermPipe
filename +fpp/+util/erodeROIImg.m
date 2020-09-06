
% erodeROIImg
%
% Function to remove voxels on the boundary of an ROI.
%
% Arguments
% img: 3D binary image matrix.

function newimg = erodeROIImg(img)

if ~(isnumeric(img) || islogical(img)) || numel(size(img))~=3
    disp('Input must be a 3D img matrix.');
    return;
end

if numel(img)~=numel(img(ismember(img,[0 1])))
    disp('Input img must be a binary mask.');
    return;
end

dims = size(img);
newimg = img;

for i = 2:(dims(1)-1)
    for j = 2:(dims(2)-1)
        for k = 2:(dims(3)-1)
            if sum(sum(sum(img((i-1):(i+1),(j-1):(j+1),(k-1):(k+1)))))~=27
                newimg(i,j,k)=0;
            end
        end
    end
end


end