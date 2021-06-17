
% outStr = fpp.util.numPad(num,len)
%
% Function to convert integer to string padded with zeros

function outStr = numPad(num,len)

outStr = int2str(num);

while length(outStr)<len
    outStr = ['0' outStr];
end

end