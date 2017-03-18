
% Function to determine TR to use for an image.

function tr = determineTR(funcPath,forceTR,defaultTR,fslPrefix)

if forceTR
    tr = defaultTR;
    fprintf('%s\n',['Forcing default TR of ' num2str(tr) 's.']);
else
    [~, tr] = system([fslPrefix ' fslval ' funcPath ' pixdim4']);
    [~, tu] = system([fslPrefix ' fslval ' funcPath ' time_units']);
    tu = strtrim(tu);
    goodTU = 1;
    if strcmp(tu,'ms')
        tr = str2num(tr)/1000;
    elseif strcmp(tu,'s')
        tr = str2num(tr);
    else
        goodTU = 0;
    end;
    if goodTU && tr>0
        fprintf('%s\n',['Setting TR to ' num2str(tr) 's.']);
    else
        fprintf('%s\n',['Functional data has no TR value or unknown time_units '...
            'field in header. Using default TR of ' num2str(defaultTR) '.']);
        tr = defaultTR;
    end;
end;

end