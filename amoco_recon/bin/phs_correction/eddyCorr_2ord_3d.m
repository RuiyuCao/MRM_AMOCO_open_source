function [vDir_eddyCorr, vDir_Eddy] = eddyCorr_2ord_3d(vDir, perc)
% eddyCorr_2ord_3d
% -------------------------------------------------------------------------
% Eddy current correction using 3D second‑order polynomial fitting.
%
% A polynomial surface is fitted to static tissue regions and removed
% from the velocity phase image to eliminate background phase offsets.
%
% INPUT
%   vDir : velocity phase image
%          size = [Nx Ny Nz Nt]
%
%   perc : threshold used to identify static tissue
%          voxels with temporal std < perc * max(std) are considered static
%
% OUTPUT
%   vDir_eddyCorr : corrected velocity phase image
%   vDir_Eddy     : estimated background phase component
%
% METHOD
%   Fit a 3D second‑order polynomial:
%
%   f(x,y,z) =
%       a0 + a1*x + a2*y + a3*z
%       + a4*x^2 + a5*y^2 + a6*z^2
%       + a7*x*y + a8*x*z + a9*y*z
%
% -------------------------------------------------------------------------

[Nx, Ny, Nz, Nt] = size(vDir);

vDir_sd = std(vDir,1,4);
thresh = perc * max(vDir_sd(:));

[xBig,yBig,zBig] = ndgrid(1:Nx,1:Ny,1:Nz);

mask_static = vDir_sd < thresh;

x = xBig(mask_static);
y = yBig(mask_static);
z = zBig(mask_static);

A = [ones(size(x)), ...
     x, y, z, ...
     x.^2, y.^2, z.^2, ...
     x.*y, x.*z, y.*z];

vDir_Eddy = zeros(Nx,Ny,Nz,Nt);
vDir_eddyCorr = zeros(Nx,Ny,Nz,Nt);

for t = 1:Nt

    v_t = vDir(:,:,:,t);
    v_static = v_t(mask_static);

    coeff = pinv(A) * v_static(:);

    E = coeff(1) ...
      + coeff(2)*xBig + coeff(3)*yBig + coeff(4)*zBig ...
      + coeff(5)*xBig.^2 + coeff(6)*yBig.^2 + coeff(7)*zBig.^2 ...
      + coeff(8)*xBig.*yBig + coeff(9)*xBig.*zBig + coeff(10)*yBig.*zBig;

    vDir_Eddy(:,:,:,t) = E;

    vDir_eddyCorr(:,:,:,t) = angle(exp(1i*v_t) ./ exp(1i*E));

end

end