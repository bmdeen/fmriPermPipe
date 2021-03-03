
% Saves structure st into filename
%
% fpp.bids.tsvWrite(filename,st)
%
% st is a structure created with st=tsvread('file.tsv');
%
% st is a structure with several fields.  Each field is a vector of
% numbers, a cell array of strings, or a matrix of char. Field names are 
% used as headers of each column.
%
% Warning: General format %.16g is used for numerical values. It works fine most
% of the time. Some applications may need to change the output format.
%
% Rafael Palacios, Oct 2009
% Modified by Ben Deen, Jul 2020
%
function tsvWrite(filename,st)

%%Error checking
error(nargchk(2, 2, nargin));  %2 arguments required, 2 maximum
if (~ischar(filename))
    error('First argument must be the name of the file');
end
if (~isstruct(st))
    error('Second argument must be a structure');
end
%Field names
names=fieldnames(st);
% Convert row to column vectors
for i=1:length(names)
    thisField = getfield(st,names{i});
    if size(thisField,1)<size(thisField,2)
        st.(names{i}) = st.(names{i})';
    end
end
rows=size(getfield(st,names{1}),1);
for j=2:length(names)
    if (rows~=size(getfield(st,names{j}),1))
        error('Field $s has a different length than first field (%s)',names{j},names{1});
    end
end

%%
[fp,message]=fopen(filename,'w');
if (fp==-1)
    error('Error opening file: %s',message);
end
%header
fprintf(fp,'%s',names{1});
fprintf(fp,'\t%s',names{2:end});
fprintf(fp,'\n');
%values
for i=1:rows
    for j=1:length(names)
        if (j~=1)
            fprintf(fp,'\t');
        end
        v=getfield(st,names{j});
        if iscell(v)
            if ischar(v{i})
                fprintf(fp,'%s',v{i});
            else
                fprintf(fp,'%.16g',v{i});  %general format
            end
        elseif ischar(v(1,1))
            fprintf(fp,'%s',v(i,:));
        else
            fprintf(fp,'%.16g',v(i));  %general format
        end
    end
    if i~=rows, fprintf(fp,'\n'); end
end
fclose(fp);
