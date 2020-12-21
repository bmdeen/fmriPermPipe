
% fpp.util.system(cmd,verbose)
%
% Wrapper for MATLAB's system function, which throws a MATLAB error if the
% system command errors.

function [status,cmdOut] = system(cmd,verbose)

if ~exist('verbose','var') || isempty(verbose)
    verbose = 0;
end

if verbose
    [status,cmdOut] = system(cmd,'-echo');
else
    [status,cmdOut] = system(cmd);
end

if status~=0 || contains(cmdOut,'Image Exception : #')
    error('%s\n\n%s\n\n','Error executing system command:',cmdOut);
end

end