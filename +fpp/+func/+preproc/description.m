
% Function to output a description string for preprocessed fMRI data, given
% an intro sentence and a set of steps.
%
% output = fpp.func.preproc.description(intro,steps)
%
% Arguments:
% - intro (string): Intro sentence string
% - steps (cell array of strings): processing steps
%

function output = description(intro,steps)

while strcmp(intro(end),' '), intro = intro(1:end-1); end
if ~strcmp(intro(end),'.'), intro = [intro '.']; end
output = intro;

if ~isempty(steps)
    output = [output ' Steps include: '];
    for s=1:length(steps)
        output = [output steps{s} '; '];
        if s==length(steps)-1
            output = [output 'and '];
        end
    end
    output = [output(1:end-2) '.'];
end

end