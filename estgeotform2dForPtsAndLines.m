function [tform, inlierIndex, status] = estgeotform2dForPtsAndLines(matchedPoints1, matchedPoints2,...
    matchedLines1,matchedLines2)
% Brief: Estimate 2-D geometric rigid transformation from matching point pairs and 
% matching line pairs.
% Details:
%  Similar to build-in `estgeotform2d` function,only for 2-D rigid transformations 
% (i.e., only rotations and translations),but with matching lines added for 
% simultaneous estimation, and other parameters consistent with the built-in function.
% 
% Syntax:  
%     [tform, inlierIndex, status] = estgeotform2dForPtsAndLines(matchedPoints1,...
% matchedPoints2,matchedLines1,matchedLines2)
% 
% Inputs:
%    matchedPoints1 - [m,2] size,[double] type,Matched points from the 
% first image, specified as an m-by-2 matrix of m number of [x y] coordinates.
%
%    matchedPoints2 - [m,2] size,[double] type,Matched points from the
% second image, specified as an m-by-2 matrix of m number of [x y] coordinates.
%
%    matchedLines1 - [n,4] size,[double] type,Matched lines from the 
% first image, specified as an n-by-4 matrix of n number of [x1,y1,x2,y2] 
% coordinates,represents the (x1,y1),(x2,y2) coordinates of the two
% endpoints of the line segment.
%
%    matchedLines2 - [n,4] size,[double] type,Matched lines from the 
% second image, specified as an n-by-4 matrix of n number of [x1,y1,x2,y2] 
% coordinates,represents the (x1,y1),(x2,y2) coordinates of the two
% endpoints of the line segment.
%
%    options - [1,1] size,[name-value] type,See follow arguments input block, 
% which current R2023b does not support codegen for entry-point function name-value
% argument validation.
% 
% Outputs:
%    tform - [2,3] size,[double] type,like
%    [cos(theta),-sin(theta),tx;sin(theta),cos(theta),ty] format.
%
%    inlierIndex - [m+n,1] size,[logical] type,whether it is an inlier point,
% 1 is an inlier point, otherwise it is an outlier point.
%
%    status - [1,1] size,[double] type,0,1,2,Consistent with the MATLAB built-in
% function `estgeotform2d` return value status.
% 
% Notes:
% 1.Matching lines only represent matching consistency in direction only,
% it does not mean that the two endpoints of the line also overlap to match!
%
% 2.Each line of matchedLines1 corresponds to each line of matchedLines2,
% and the order of the two defined points should match in the direction of 
% each other, and there should not be any 180-degree matching in the opposite
% direction.
%
% 3.if no matched points or lines,please pass the parameters as [].
%
% codegen command Example:
%    cfg = coder.config('lib','ecoder',true);
%    cfg.FilePartitionMethod = "SingleFile";
%    cfg.Toolchain  = "CMake";% since R2023a，Generate generic CMakeLists.txt file
%    cfg.GenerateComments = false;
%    cfg.HardwareImplementation.ProdHWDeviceType = "Texas Instruments->C6000";
%    cfg.GenCodeOnly = true;
%    cfg.CppNamespace = "estgeotform2dForPtsAndLines";
%
%    pts1 = coder.typeof(double(0),[inf,2],[1,0]);
%    pts2 = coder.typeof(double(0),[inf,2],[1,0]);
%    matchedLines1 = coder.typeof(0,[inf,4],[1,0]);
%    matchedLines2 = coder.typeof(0,[inf,4],[1,0]);
%
%    codegen -config cfg estgeotform2dForPtsAndLines -args {pts1,pts2,matchedLines1,matchedLines2}  -lang:c++ -report
%
% Example: 
%   theta = 30;
%   t = [100,200];
% 
%   R_theory = [cosd(theta),-sind(theta);
%     sind(theta),cosd(theta)]
%   t_theory = t
%   numPts = 1000;
%   matchedPts1 = 100*rand(numPts,2);
%   matchedPts2 = R_theory*matchedPts1'+t_theory'; 
%   matchedPts2 = matchedPts2';
% 
%   matchedL1 = [0,0,100,0;
%     0,0,100,100*tand(30)];
%   matchedL2 = [0,0,100,100*tand(30);
%     0,100,50,100+50*tand(60)];
%
%   [tform,inds,status] = estgeotform2dForPtsAndLines(matchedPts1,...
%   matchedPts2,matchedL1,matchedL2);
% 
% See also: 
%    estgeotform2d

% Author:                          cuixingxing
% Email:                           cuixingxing150@gmail.com
% Created:                         14-Sep-2023 16:08:40
% Version history revision notes:
%                                  2023.10.25 fix some issues.
% Implementation In Matlab R2023a
% Copyright © 2023 TheMatrix.All Rights Reserved.
%
%#codegen
arguments
    matchedPoints1 (:,2) double
    matchedPoints2 (:,2) double
    matchedLines1 (:,4) double
    matchedLines2 (:,4) double
    % options.MaxNumTrials (1,1) {mustBePositive}=1000;
    % options.Confidence (1,1) {mustBePositive}=99; % Confidence of finding the maximum number of inliers, specified as a positive numeric scalar in the range (0, 100). 
    % options.MaxDistance (1,1) =1.5;
end

m = size(matchedPoints1,1);
n = size(matchedLines1,1);
assert(m==size(matchedPoints2,1)&&n==size(matchedLines2,1));

% least elements for model
cond1 = m>=2;
cond2 = (m==1)&&(n>=1);
cond3_1 = n>=2;
cond3 = false;
if cond3_1
    cond3 = judgeLinesValid(matchedLines1)||judgeLinesValid(matchedLines2);
end

% initialize
tform = [1,0,0;
    0,1,0];
inlierIndex = false(m+n,1);
status = 0;
if ~(cond1||cond2||cond3)
    status = 1; % Lack of sufficient matching points or lines
    return;
end

% MSAC algorithm
maxNumTrials = 1000;
maxSkipTrials = maxNumTrials * 10;
MaxDistance = 1.5;
Confidence = 99;

idxTrial = int32(1);
numTrials = int32(maxNumTrials);
skipTrials = 0;
sampleSize = 2;
maxDis    = cast(MaxDistance*(m+n), "double");
bestDis   = maxDis;

bestInliers = false(m+n, 1);
while idxTrial <= numTrials && skipTrials < maxSkipTrials
    % Random selection without replacement
    indices = randperm(m+n, sampleSize);
    
    % belong which model
    if all(indices<=m)
        matchedPts1 = matchedPoints1(indices,:);
        matchedPts2 = matchedPoints2(indices,:);
        modelParams = computeRigidTformFromPts(matchedPts1,matchedPts2);
    elseif all(indices>m)
        matchedL1 = matchedLines1(indices-m,:);
        matchedL2 = matchedLines2(indices-m,:);
        modelParams = computeRigidTformFromLine(matchedL1,matchedL2);
    else
        indices = sort(indices);
        matchedPt1 = matchedPoints1(indices(1),:);
        matchedPt2 = matchedPoints2(indices(1),:);
        matchedL1 = matchedLines1(indices(2)-m,:);
        matchedL2 = matchedLines2(indices(2)-m,:);
        modelParams = computeRigidTformFromPtAndLine(matchedPt1,matchedPt2,matchedL1,matchedL2);
    end
    
    % Validate the model
    isValidModel = checkTForm(modelParams);
    
    if isValidModel
        % Evaluate model with truncated loss
        [dis, accDis] = evaluateModel(modelParams,...
            matchedPoints1,matchedPoints2,...
            matchedLines1,matchedLines2,...
            MaxDistance);
        
        % Update the best model found so far
        if accDis < bestDis
            bestDis = accDis;
            bestInliers = dis < MaxDistance;
            tform = modelParams;
            inlierIndex = bestInliers;

            inlierNum = cast(sum(inlierIndex),"double");
            num = computeLoopNumber(sampleSize, ...
                Confidence, m+n, inlierNum);
            numTrials = min(numTrials, num);
        end
        
        idxTrial = idxTrial + 1;
    else
        skipTrials = skipTrials + 1;
    end
end

isFound = sum(bestInliers(:)) >= sampleSize;
if isFound
    inliersPts1 = matchedPoints1(bestInliers(1:m),:);
    inliersPts2 = matchedPoints2(bestInliers(1:m),:);
    inliersL1 = matchedLines1(bestInliers(m+1:end),:);
    inliersL2 = matchedLines2(bestInliers(m+1:end),:);
    modelParams = svdFitForInliersPointsAndLines(inliersPts1,inliersPts2,inliersL1,inliersL2);

    dis = evaluateModel(modelParams,...
            matchedPoints1,matchedPoints2,...
            matchedLines1,matchedLines2,...
            MaxDistance);

    isValidModel = checkTForm(modelParams);
    inlierIndex = (dis < MaxDistance);
    if ~isValidModel || ~any(inlierIndex)
        inlierIndex = false(m+n, 1);
    else
        tform = modelParams;
    end

    if numTrials >= int32(maxNumTrials)
        coder.internal.warning('vision:ransac:maxTrialsReached');
    end
else
    status = 2;
    inlierIndex = false(m+n, 1);
end

end

function flag = judgeLinesValid(lines)
num = size(lines,1);
flag = false;
for i = 1:num-1
    for j = i+1
        line1 = lines(i,:);
        line2 = lines(j,:);
        a = [line1(3)-line1(1),line1(4)-line1(2)];
        b = [line2(3)-line2(1),line2(4)-line2(2)];
        distA = sqrt(a(1).^2+a(2).^2);
        distB = sqrt(b(1).^2+b(2).^2);
        val = max(min(dot(a,b)./(distA.*distB),1),-1);
        theta = real(acosd(val));
        if theta>3 % 若满足，即可以看作平行线
            flag = true;
        end
    end
end
end

function tform = svdFitForInliersPointsAndLines(inliersPts1,inliersPts2,...
    inliersL1,inliersL2)
% points
theta1 = 0;
t1 = zeros(2,1);
theta1_flag = false;
t1_flag = false;
if size(inliersPts1,1)>=2
    [R1,t1] = computeRigidTransform(inliersPts1, inliersPts2);
    theta1 = rotm2eul(blkdiag(R1,1));
    theta1 = theta1(1)*180/pi;

    theta1_flag = true;
    t1_flag = true;
elseif size(inliersPts1,1)==1
    t1 = inliersPts2'-inliersPts1';

    t1_flag = true;
end

% lines
theta2 = 0;
t2 = zeros(2,1);
theta2_flag = false;
t2_flag = false;
if size(inliersL1,1)>=2
    a = [inliersL1(:,3)-inliersL1(:,1),inliersL1(:,4)-inliersL1(:,2)];
    b = [inliersL2(:,3)-inliersL2(:,1),inliersL2(:,4)-inliersL2(:,2)];
    distA = sqrt(a(:,1).^2+a(:,2).^2);
    distB = sqrt(b(:,1).^2+b(:,2).^2);
    innerDot = a(:,1).*b(:,1)+a(:,2).*b(:,2);
    val = max(min(innerDot./(distA.*distB),1),-1);
    val = rmoutliers(acosd(val));
    theta2 = mean(val);

    indexs = nchoosek(1:size(inliersL1,1),2);
    pnums = size(indexs,1);
    pnums = min(pnums,50);% avoid to multiply times
    indexsSub = randperm(size(indexs,1),pnums);
    indexs = indexs(indexsSub,:);

    intersectPts1 = zeros(0,2);
    intersectPts2 = zeros(0,2);
    coder.varsize("intersectPts1","intersectPts2",[inf,2],[1,0]);
    for idx = 1:pnums
        idx1 = indexs(idx,1);
        idx2 = indexs(idx,2);
        intesectPt1 = computeIntesect(inliersL1(idx1,:),inliersL1(idx2,:));
        intesectPt2 = computeIntesect(inliersL2(idx1,:),inliersL2(idx2,:));
        if any(isnan(intesectPt1))||any(isnan(intesectPt2))
            continue;
        end
        intersectPts1 = [intersectPts1;intesectPt1];
        intersectPts2 = [intersectPts2;intesectPt2];
    end

    intersectPts = rmoutliers([intersectPts1,intersectPts2]);
    intersectPts1 = intersectPts(:,1:2);
    intersectPts2 = intersectPts(:,3:4);
    if ~isempty(intersectPts1)&&~isempty(intersectPts2)
        t2 = intersectPts2'-[cosd(theta2),-sind(theta2);sind(theta2),cosd(theta2)]*intersectPts1';
        t2 = mean(t2,2);
        t2_flag = true;
    else
        t2_flag = false;
    end
    theta2_flag = true;
elseif size(inliersL1,1)==1
    a = [inliersL1(1,3)-inliersL1(1,1),inliersL1(1,4)-inliersL1(1,2)];
    b = [inliersL2(1,3)-inliersL2(1,1),inliersL2(1,4)-inliersL2(1,2)];
    distA = sqrt(a(1).^2+a(2).^2);
    distB = sqrt(b(1).^2+b(2).^2);
    val = max(min(dot(a,b)./(distA.*distB),1),-1);
    theta2 = real(acosd(val));

    theta2_flag = true;
end

if theta1_flag&&theta2_flag
    theta = (theta1+theta2)/2;
elseif theta1_flag&&(~theta2_flag)
    theta = theta1;
elseif (~theta1_flag)&&theta2_flag
    theta = theta2;
else
    theta = 0;
end

if t1_flag&&t2_flag
    t = (t1+t2)/2;
elseif t1_flag&&(~t2_flag)
    t = t1;
elseif (~t1_flag)&&t2_flag
    t = t2;
else
    t = [0;0];
end

tform = [cosd(theta),-sind(theta),t(1);
    sind(theta),cosd(theta),t(2)];
end


function tform = computeRigidTformFromPts(pts1,pts2)
a = [pts1(2,1)-pts1(1,1),pts1(2,2)-pts1(1,2)];
b = [pts2(2,1)-pts2(1,1),pts2(2,2)-pts2(1,2)];
distA = sqrt(a(1).^2+a(2).^2);
distB = sqrt(b(1).^2+b(2).^2);
val = max(min(dot(a,b)./(distA.*distB),1),-1);
rotateT = acosd(val);
R3 = [cosd(rotateT),-sind(rotateT);
    sind(rotateT),cosd(rotateT)];
P = mean(pts1)';
Q = mean(pts2)';
t3 = Q- R3*P;
tform = [R3,t3];
end

function tform = computeRigidTformFromPtAndLine(pt1,pt2,line1,line2)
a = [line1(1,3)-line1(1,1),line1(1,4)-line1(1,2)];
b = [line2(1,3)-line2(1,1),line2(1,4)-line2(1,2)];
distA = sqrt(a(1).^2+a(2).^2);
distB = sqrt(b(1).^2+b(2).^2);
val = max(min(dot(a,b)./(distA.*distB),1),-1);
rotateT = acosd(val);
R = [cosd(rotateT),-sind(rotateT);
    sind(rotateT),cosd(rotateT)];
t = pt2'-R*pt1';
tform = [R,t];
end

function intesectPt = computeIntesect(line1,line2)
% 过两点(x1,y1),(x2,y2)的直线line1与过两点(m1,n1),(m2,n2)的直线line2求交点(x,y),记为intesectPt
% syms x y x1 y1 x2 y2 m1 n1 m2 n2
% equ1 = (x-x1)*(y2-y1)==(y-y1)*(x2-x1);
% equ2 = (x-m1)*(n2-n1)==(y-n1)*(m2-m1);
% [x,y] = solve([equ1,equ2],x,y)
%
arguments
    line1 (1,4) double % 形如[x1,y1,x2,y2]两个端点的直线
    line2 (1,4) double  % 形如[m1,n1,m2,n2]两个端点的直线
end
x1 = line1(1,1);
y1 = line1(1,2);
x2 = line1(1,3);
y2 = line1(1,4);

m1 = line2(1,1);
n1 = line2(1,2);
m2 = line2(1,3);
n2 = line2(1,4);

d = m1*y1 - n1*x1 - m1*y2 - m2*y1 + n1*x2 + n2*x1 + m2*y2 - n2*x2;
intesectPt = nan(1,2);
if abs(d)>1e-5
    intesectPt(1) = (m1*n2*x1 - m2*n1*x1 - m1*n2*x2 + m2*n1*x2 - m1*x1*y2 + m1*x2*y1 + m2*x1*y2 - m2*x2*y1)/d;
    intesectPt(2) = (m1*n2*y1 - m2*n1*y1 - m1*n2*y2 + m2*n1*y2 - n1*x1*y2 + n1*x2*y1 + n2*x1*y2 - n2*x2*y1)/d;
end
end

function tform = computeRigidTformFromLine(line1,line2)
arguments
    line1 (2,4) double % 每行代表一条直线的两个端点，共两条直线，与line2匹配
    line2 (2,4) double % 每行代表一条直线的两个端点，共两条直线，与line1匹配
end
intesectPt1 = computeIntesect(line1(1,:),line1(2,:));
intesect_x1 = intesectPt1(1);
intesect_y1 = intesectPt1(2);

intesectPt2 = computeIntesect(line2(1,:),line2(2,:));
intesect_x2 = intesectPt2(1);
intesect_y2 = intesectPt2(2);

if any(isnan(intesectPt1))||any(isnan(intesectPt2))
    tform = [1,0,0;
        0,1,0];
else
    % rotate
    a = [line1(1,3)-line1(1,1),line1(1,4)-line1(1,2)];
    b = [line2(1,3)-line2(1,1),line2(1,4)-line2(1,2)];
    distA = sqrt(a(1).^2+a(2).^2);
    distB = sqrt(b(1).^2+b(2).^2);
    val = max(min(dot(a,b)./(distA.*distB),1),-1);
    rotateT1 = acosd(val);

    a = [line1(2,3)-line1(2,1),line1(2,4)-line1(2,2)];
    b = [line2(2,3)-line2(2,1),line2(2,4)-line2(2,2)];
    distA = sqrt(a(1).^2+a(2).^2);
    distB = sqrt(b(1).^2+b(2).^2);
    val = max(min(dot(a,b)./(distA.*distB),1),-1);
    rotateT2 = acosd(val);

    rotateT = (rotateT1+rotateT2)/2;
    R = [cosd(rotateT),-sind(rotateT);
        sind(rotateT),cosd(rotateT)];
    t = [intesect_x2;intesect_y2]-R*[intesect_x1;intesect_y1];
    tform = [R,t];
end

end

function tf = checkTForm(tform)
tf = all(isfinite(tform(:)));
end

function [distances, sumDistances] = evaluateModel(modelIn,matchedPoints1,matchedPoints2,...
    matchedLines1,matchedLines2, threshold) %#codegen
m = size(matchedPoints1,1);
n = size(matchedLines1,1);

dist1 = zeros(m,1);
if m~=0
    pt1 = [matchedPoints1';
        ones(1,m)];
    calPt2 = modelIn*pt1;
    delta = calPt2'-matchedPoints2;
    dist1 = hypot(delta(:,1),delta(:,2));
    dist1(dist1>threshold)=threshold;
end

dist2 = zeros(n,1);
if n~=0
    pt1 = [reshape(matchedLines1',2,2*n);
        ones(1,2*n)];
    pt2 = modelIn*pt1; 
    calPt2 = reshape(pt2,4,n);
    calPt2 = calPt2';

    % 求点(m1,n1),(m2,n2)分别到直线(y2-y1)*x+(x1-x2)*y+x2*y1-x1*y2==0的距离
    x1 = matchedLines2(:,1);
    y1 = matchedLines2(:,2);
    x2 = matchedLines2(:,3);
    y2 = matchedLines2(:,4);
    D = sqrt((y2-y1).^2+(x2-x1).^2);
    m1 = calPt2(:,1);
    n1 = calPt2(:,2);
    m2 = calPt2(:,3);
    n2 = calPt2(:,4);
    distance1 = abs((y2-y1).*m1+(x1-x2).*n1+x2.*y1-x1.*y2)./D;% 第一个端点到n条线的距离
    distance2 = abs((y2-y1).*m2+(x1-x2).*n2+x2.*y1-x1.*y2)./D;% 第二个端点到n条线的距离
    dist2 = (distance1+distance2)/2;
    dist2 = dist2*0.3; % 线段人为考虑比点阈值小30%
    dist2(dist2>threshold)=threshold;
end

distances = [dist1;dist2];
sumDistances = sum(distances);
end
