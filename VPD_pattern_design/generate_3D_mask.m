% =========================================================================
% Random Variable-Density Sampling Mask Generation
%
% This script generates multiple 2D k-space sampling masks using a
% Gaussian-weighted random sampling strategy with a fully sampled center.
% The masks are designed for accelerated MRI acquisitions.
%
% Author:
%   Ruiyu Cao, Fudan University
%   h_li@fudan.edu.cn
%   07/2024
% =========================================================================

clc; clear; close all;
%% ------------------------------------------------------------------------
% Basic mask parameters
% -------------------------------------------------------------------------

maskSize = [132,36];            % [kx, ky] size of sampling mask

USF = 5;                        % Desired undersampling factor

patternNum = 28;                % Number of masks to generate

% Number of samples per mask
totalPoints = ceil(maskSize(1)*maskSize(2)/USF);

%% ------------------------------------------------------------------------
% Gaussian sampling distribution parameters
% -------------------------------------------------------------------------

% Standard deviation controlling sampling density
sigmaX = maskSize(1)/6;         % Gaussian width in kx
sigmaY = maskSize(2)/3;         % Gaussian width in ky

% Fully sampled central ellipse
centerRadiusX = ceil(0.46*sigmaX);
centerRadiusY = ceil(0.46*sigmaY);

%% ------------------------------------------------------------------------
% Sampling control parameters
% -------------------------------------------------------------------------

% Minimum distance scaling factor between sampled points
% Larger value → more uniform sampling
minDistFactor = 2.1;

% Probability decay controlling repetition of samples across masks
probDecay = 0.5;

%% ------------------------------------------------------------------------
% Generate sampling masks
% -------------------------------------------------------------------------

masks = fun_maskGen3D_v2( ...
            maskSize, ...
            centerRadiusX, centerRadiusY, ...
            totalPoints, ...
            patternNum, ...
            sigmaX, sigmaY, ...
            minDistFactor, ...
            probDecay);

%% Display all generated masks

figure
imshow3D(permute(masks,[2,1,3]))
title('All Masks')

%% ------------------------------------------------------------------------
% Evaluate distribution using randomly selected masks
% -------------------------------------------------------------------------

combinedNum = 6;                         % Number of masks to combine
selectedIndices = randperm(patternNum,combinedNum);

combinedMask = sum(masks(:,:,selectedIndices),3);

%% Combined mask with part of samples
figure
imshow(permute(combinedMask>0,[2,1]),[])
title(sprintf('Combined Mask (%d patterns)',combinedNum))

%% ------------------------------------------------------------------------
% Visualize sampling density distribution
% -------------------------------------------------------------------------

[sortedValues,~,idx] = unique(combinedMask);

newValues = 1:length(sortedValues);

combinedMask = reshape(newValues(idx),size(combinedMask));

figure
imshow(permute(combinedMask,[2,1]),[])
colormap jet

cmap = colormap;
cmap(1,:) = [0 0 0];
colormap(cmap)

%% ------------------------------------------------------------------------
% Combine all masks to analyze global sampling distribution
% -------------------------------------------------------------------------

combinedMask = sum(masks,3);

figure
imshow(permute(combinedMask>0,[2,1]),[])

[sortedValues,~,idx] = unique(combinedMask);

newValues = 1:length(sortedValues);

combinedMask = reshape(newValues(idx),size(combinedMask));

figure
imshow(permute(combinedMask,[2,1]),[])
colormap jet

cmap = colormap;
cmap(1,:) = [0 0 0];
colormap(cmap)

title('All Sampling Points Distribution')

%% ------------------------------------------------------------------------
% Calculate effective undersampling factor
% -------------------------------------------------------------------------

underFactor = (maskSize(1)*maskSize(2)*patternNum) / sum(masks(:)>0);

disp("Undersampled Factor: " + underFactor)

%% ------------------------------------------------------------------------
% Store mask in standard format
% -------------------------------------------------------------------------

pattern.mask = permute(masks,[2,1,3]);

%% ------------------------------------------------------------------------
% Visualize masks in temporal order
% -------------------------------------------------------------------------

all_phs = [];

for i = 1:size(masks,3)

    all_phs = cat(2,all_phs,permute(masks(:,:,i),[2,1]));

end

% figure, imshow(all_phs,[])

%% ------------------------------------------------------------------------
% Split mask into three random subsets
% -------------------------------------------------------------------------

mask_3phs = splitTo3Parts(pattern.mask);


%% ========================================================================
% Function: split mask into 3 random subsets
% ========================================================================

function patternsCell = splitTo3Parts(pattern)

% Extract linear indices of sampled points
nonZeroIdx = find(pattern);

% Total number of sampled points
numOnes = numel(nonZeroIdx);

% Size of each subset
oneThird = floor(numOnes/3);

% Random permutation
randomOrder = randperm(numOnes);

% Divide indices
group1Idx = randomOrder(1:oneThird);
group2Idx = randomOrder(oneThird+1:2*oneThird);
group3Idx = randomOrder(2*oneThird+1:end);

% Initialize output masks
patternsCell = cell(1,3);

for i = 1:3
    patternsCell{i} = false(size(pattern));
end

% Assign sampled locations
patternsCell{1}(nonZeroIdx(group1Idx)) = true;
patternsCell{2}(nonZeroIdx(group2Idx)) = true;
patternsCell{3}(nonZeroIdx(group3Idx)) = true;

% Verify reconstruction
assert(isequal(pattern,patternsCell{1} | patternsCell{2} | patternsCell{3}), ...
       'Output patterns do not match the input pattern.');

end
