function out = gradEntropy(in,gw)
% GRADENTROPY calculates the gradient entropy of an image
%    GRADENTROPY(IN,GW) returns the gradient entropy of the
%    image IN. Gradient entropy is a metric used to quantify
%    MR motion corruption in the image space [1]. Parameter 
%    GW is a vector with the same length as the number of
%    dimensions of IN. The vector GW determines the order of
%    the separable hanning filter in each direction.
%
%    [1] McGee et al. "Image metric-based correction (autocorrection) 
%         of motion effects: analysis of image metrics. JMRI 2000; 
%          11:174-181.
%
%
% (c) Joseph Y Cheng (jycheng@mrsrl.stanford.edu) 2012
    
% SVN info: 
%   Date:     $Date: 2012-02-11 14:56:23 -0800 (Sat, 11 Feb 2012) $
%   Revision: $Revision: 1071 $
%   Author:   $Author: jycheng $
%   Id:       $Id: gradEntropy.m 1071 2012-02-11 22:56:23Z jycheng $
    
    DEBUG = 0;
    
    if DEBUG
        tic;
    end
    
    %% Dimensions of input
    dims  = size(in);
    ndims = length(dims);
    
    %% Calculate image gradients (not MR imaging gradients!)
    [gs(1:ndims).FIELD] = gradient(abs(in));
    g = zeros(dims);
    for D=1:ndims
        g = g + abs(gs(D).FIELD);
    end
    
    %% First step of gradient entropy calculation
    glog2g = g.*log2(g+eps);
    
    %% Construct separable hanning filters
    hfilts = {};
    perm_i = [1:ndims];
    perm_p = [ndims 1:ndims-1];
    for D=1:ndims
        hfilts{D} = permute(hanning(gw(D)),perm_i);
        perm_i = perm_i(perm_p);
    end
    
    %% Post-convolution weights (fix edges)
    w = ones(size(in));
    
    %% Conovlution: image-based
    glog2gf = glog2g;
    gf      = g;    
    for D=1:ndims
        w       = convn(w,      hfilts{D},'same');
        glog2gf = convn(glog2gf,hfilts{D},'same');
        gf      = convn(gf,     hfilts{D},'same');
    end
    glog2gf = glog2gf./w;
    gf = gf./w;
        
    %% Final calculation and write to output
    out = -(glog2gf./gf-log2(gf+eps));
    
    if DEBUG
        dbdisp(sprintf('Done calculation: %0.4f s', toc));
    end
    