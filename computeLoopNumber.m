%==========================================================================
%  Computer Loop Number for RANSAC/MSAC Algorithm
%==========================================================================
function N = computeLoopNumber(sampleSize, confidence, pointNum, inlierNum)

% Copyright 2016-2020 The MathWorks, Inc.

%#codegen
pointNum   = cast(pointNum, 'like', inlierNum);
sampleSize = cast(sampleSize, 'like', inlierNum);
inlierProbability = (inlierNum/pointNum)^sampleSize;

if inlierProbability < eps(class(inlierNum))
    N = intmax('int32');
else
    conf = cast(0.01 * confidence, 'like', inlierNum);
    one  = ones(1,    'like', inlierNum);
    num  = log10(one - conf);
    den  = log10(one - inlierProbability);
    N    = int32(ceil(num/den));
end 