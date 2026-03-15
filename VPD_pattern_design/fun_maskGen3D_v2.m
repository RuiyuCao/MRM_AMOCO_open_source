function masks = fun_maskGen3D_v2(maskSize, centerRadiusX, centerRadiusY, totalPoints, patternNum, sigmaX, sigmaY, minDistFactor, decayRate)
% fun_maskGen3D_v2
% =========================================================================
% Generate multiple 2D sampling masks using Gaussian-weighted random
% sampling with adaptive minimum-distance constraints.
%
% The generated masks follow a variable-density sampling strategy:
%   • High sampling density near the k-space center
%   • Lower sampling density toward the periphery
%
% A fully sampled elliptical central region is enforced to preserve
% low-frequency information. Remaining samples are generated using a
% Gaussian probability distribution and a minimum-distance constraint
% to ensure spatial uniformity.
%
% Additionally, the sampling probability of previously selected points
% gradually decays to reduce repeated sampling across different masks.
%
% INPUT
%   maskSize      : [width height] size of sampling mask
%   centerRadiusX : semi-major axis of fully sampled ellipse (X direction)
%   centerRadiusY : semi-minor axis of fully sampled ellipse (Y direction)
%   totalPoints   : total number of sampling points per mask
%   patternNum    : number of masks to generate
%   sigmaX        : Gaussian distribution width along X direction
%   sigmaY        : Gaussian distribution width along Y direction
%   minDistFactor : scaling factor controlling minimum distance between
%                   sampled points
%   decayRate     : probability decay factor used to suppress repeated
%                   sampling across masks
% OUTPUT
%   masks : sampling masks
%           size = [height width patternNum]
% Author
%   Ruiyu Cao
%   2024.09.02
% =========================================================================


width = maskSize(1);
height = maskSize(2);

masks = zeros(height, width, patternNum);

% Create initial weight matrix
initialWeight = create_gaussian_weight_matrix(maskSize, sigmaX, sigmaY);

% Calculate center coordinates
centerX = width / 2;
centerY = height / 2;

weight = initialWeight;

% Elliptical region constraint (covers the entire mask size)
[X, Y] = meshgrid(1:width, 1:height);
ellipseConstraint = ((X - centerX) / centerX).^2 + ((Y - centerY) / centerY).^2 <= 1;

for p = 1:patternNum

    mask = zeros(height, width);

    % Ensure the center region is fully sampled
    centerEllipse = ((X - centerX) / centerRadiusX).^2 + ((Y - centerY) / centerRadiusY).^2 <= 1;
    mask(centerEllipse) = 1;

    % Random sampling with weighted probability within the elliptical constraint
    points = random_sampling(maskSize, totalPoints - sum(centerEllipse(:)), weight, minDistFactor, initialWeight, ellipseConstraint);

    % Set the sampled points in the mask
    for i = 1:size(points, 1)
        x = points(i, 1);
        y = points(i, 2);
        if x > 0 && x <= width && y > 0 && y <= height && mask(y, x) == 0
            mask(y, x) = 1;

            % Calculate distance from center
            distFromCenter = sqrt((x - centerX)^2 + (y - centerY)^2) / sqrt(centerX^2 + centerY^2);

            % Determine decay factor based on distance
            if centerEllipse(y, x)
                localDecayFactor = 1;  % Inside the central ellipse, no decay
            else
                localDecayFactor = exp(-decayRate * distFromCenter);
            end

            % Update weight with the calculated local decay factor
            weight(y, x) = weight(y, x) * localDecayFactor;
        end
    end

    % Ensure the mask has the required number of points
    currentPoints = sum(mask(:));
    while currentPoints < totalPoints
        additionalPoints = random_sampling(maskSize, totalPoints - currentPoints, weight, minDistFactor, initialWeight, ellipseConstraint);
        for i = 1:size(additionalPoints, 1)
            x = additionalPoints(i, 1);
            y = additionalPoints(i, 2);
            if x > 0 && x <= width && y > 0 && y <= height && mask(y, x) == 0
                mask(y, x) = 1;

                % Calculate distance from center
                distFromCenter = sqrt((x - centerX)^2 + (y - centerY)^2) / sqrt(centerX^2 + centerY^2);

                % Determine decay factor based on distance
                if centerEllipse(y, x)
                    localDecayFactor = 1;  % Inside the central ellipse, no decay
                else
                    localDecayFactor = exp(-decayRate * distFromCenter);
                end

                % Update weight with the calculated local decay factor
                weight(y, x) = weight(y, x) * localDecayFactor;
            end
        end
        currentPoints = sum(mask(:));
    end

    % Store the mask
    masks(:, :, p) = mask;
end
end

function weight = create_gaussian_weight_matrix(maskSize, sigmaX, sigmaY)
% Create a Gaussian weight matrix where the weight decreases with distance from center
width = maskSize(1);
height = maskSize(2);
centerX = width / 2;
centerY = height / 2;

% Create a grid of distances from the center
[X, Y] = meshgrid(1:width, 1:height);
distancesX = (X - centerX).^2 / (2 * sigmaX^2);
distancesY = (Y - centerY).^2 / (2 * sigmaY^2);

% Create Gaussian weight matrix
weight = exp(-(distancesX + distancesY));
end

function points = random_sampling(maskSize, totalPoints, weight, minDistFactor, initialWeight, ellipseConstraint)
% Random sampling algorithm with weighted probability and dynamically changing minimum distance constraint
width = maskSize(1);
height = maskSize(2);

% Normalize weight to sum to 1
weight = weight / sum(weight(:));

% Create cumulative distribution function (CDF) for weighted sampling
cdf = cumsum(weight(:));

% Generate random points based on the CDF
points = zeros(totalPoints, 2);
count = 0;
while count < totalPoints
    r = rand();
    index = find(cdf >= r, 1, 'first');
    [y, x] = ind2sub([height, width], index);
    newPoint = [x, y];

    % Check if the new point is within the elliptical constraint
    if ellipseConstraint(y, x)
        % Calculate minDist based on initialWeight
        gaussianWeight = initialWeight(y, x);
        factor = (1 - gaussianWeight) / 2 + 0.5;
        minDist = minDistFactor * factor;

        % Check if the new point is at least minDist away from existing points
        if count == 0 || all(sqrt(sum((points(1:count, :) - newPoint).^2, 2)) >= minDist)
            count = count + 1;
            points(count, :) = newPoint;
        end
    end
end
end
