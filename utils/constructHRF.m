
% Function to construct hemodynamic response function, to use as a kernel
% for convolution with boxcar regressors.
%
% tr: Temporal sampling rate
%
% hrfType: 1 or 2
% 1) Double gamma, with parameters matched to FSL/SPM defaults
% 2) Single gamma, with parameters matched to FsFast defaults

function [y,t] = constructHRF(tr,hrfType)

d = 0;                  % Temporal offset (s)

switch hrfType
    case 1
        
        sigma1=2.449;   % First gamma parameters
        mu1=6;
        sigma2=4;       % Second gamma parameters
        mu2=16;
        ratio=6;        % hrf = gammapdf1 - gammapdf2/ratio;
        
        alpha1 = mu1^2/sigma1^2; tau1 = sigma1^2/mu1;
        alpha2 = mu2^2/sigma2^2; tau2 = sigma2^2/mu2;
        
        t = 0:tr:500;
        y = zeros(length(t),1);
        tShift = t-d;
        y(tShift>0) = tShift(tShift>0).^(alpha1-1).*exp(-tShift(tShift>0)/tau1)/(tau1^alpha1*gamma(alpha1)) ...
            - tShift(tShift>0).^(alpha2-1).*exp(-tShift(tShift>0)/tau2)/(tau2^alpha2*gamma(alpha2))/ratio;
        y = y/max(y);
        
    case 2
        
        tau = 1.25;     % Gamma parameters
        alpha = 2;
        
        t = 0:tr:500;
        y = zeros(length(t),1);
        tShift = t-d;
        y(tShift>0) = ((tShift(tShift>0)/tau).^alpha.*exp(-tShift(tShift>0)/tau));
        y = y/max(y);
        
end;

end