function [data_corr, MRA] = fun_phaseCorr_2ord_3d(data,perc)
% fun_phaseCorr_2ord_3d
% -------------------------------------------------------------------------
% Phase correction for 4D Flow MRI using 3D second‑order polynomial fitting.
%
% This function estimates and removes background phase offsets caused by
% eddy currents and noise. The correction is applied to the velocity
% encoding directions using a polynomial model fitted in static tissue.
%
% INPUT
%   data : complex 5D dataset
%          size = [Nx Ny Nz Nt 4]
%          (reference magnitude + 3 velocity encodings)
%
%   perc : two-element vector controlling threshold selection
%          perc(1) = threshold for static tissue detection
%                    smaller value → more background removed
%                    default = 0.02
%
%          perc(2) = threshold for air/noise removal
%                    smaller value → more aggressive noise removal
%                    default = 0.2
%
% OUTPUT
%   data_corr : phase-corrected dataset
%               size = [Nx Ny Nz Nt 4]
%
%   MRA : magnitude angiography image generated from corrected velocities
%
% NOTE
%   Velocity directions correspond to:
%       AP : Anterior‑Posterior
%       FH : Foot‑Head
%       RL : Right‑Left
%
% -------------------------------------------------------------------------

if nargin < 2
    perc = [0.02,0.2];
end

%% Phase difference calculation (4‑point encoding)

diff_data(:,:,:,:,1) = data(:,:,:,:,2)./data(:,:,:,:,1);
diff_data(:,:,:,:,2) = data(:,:,:,:,3)./data(:,:,:,:,1);
diff_data(:,:,:,:,3) = data(:,:,:,:,4)./data(:,:,:,:,1);

phs_data = angle(diff_data);

%% Phase correction

recon_corr_data = phaCorr_4dflow_2ord_3d( ...
                    abs(data(:,:,:,:,1)), ...
                    phs_data(:,:,:,:,1), ...
                    phs_data(:,:,:,:,2), ...
                    phs_data(:,:,:,:,3), ...
                    perc);

%% Visualization of correction result

map_show_be = cat(2,phs_data(:,:,:,:,1),phs_data(:,:,:,:,2),phs_data(:,:,:,:,3));
map_show_af = cat(2,recon_corr_data{2},recon_corr_data{3},recon_corr_data{4});
map_show_diff = abs(map_show_be-map_show_af);

map_show = cat(1,map_show_be,map_show_af,map_show_diff);

figure
imshow3D(squeeze(map_show(:,:,:,floor(size(map_show,4)/2))))

figure
imshow3D(squeeze(map_show(:,:,floor(size(map_show,3)/2),:)))

%% Generate MRA image

MRA = abs(recon_corr_data{1} .* ...
    sqrt(recon_corr_data{2}.^2 + recon_corr_data{3}.^2 + recon_corr_data{4}.^2));

figure
imshow3D(squeeze(max(MRA,[],3)),[])
drawnow

%% Output corrected data

data_corr(:,:,:,:,1) = recon_corr_data{1}; % magnitude
data_corr(:,:,:,:,2) = recon_corr_data{2}; % FH
data_corr(:,:,:,:,3) = recon_corr_data{3}; % RL
data_corr(:,:,:,:,4) = recon_corr_data{4}; % AP

end
