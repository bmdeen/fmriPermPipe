
% d = fpp.util.regExpDir(dirInput,regExp)
%
% Function to return contents of a directory that contain a certain MATLAB
% regular expression (case-insensitive).

function d = regExpDir(dirInput,regExp)

d = dir(dirInput);

i = 1;
while i <= length(d)
    if isempty(regexpi(d(i).name,regExp)), d(i)=[];
    else, i = i+1; end
end

end
    
    