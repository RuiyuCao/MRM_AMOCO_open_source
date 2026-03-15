function recon_corr = phaCorr_4dflow_2ord_3d(recon_ref, ph_AP_recon, ph_FH_recon, ph_RL_recon,perc)
% phaCorr_4dflow_2ord_3d
% -------------------------------------------------------------------------
% Perform phase correction for 4D Flow MRI using 3D polynomial fitting.
%
% The correction includes two steps:
%   1) Eddy current correction using second‑order polynomial fitting
%   2) Background noise removal using temporal standard deviation
%
% INPUT
%   recon_ref   : magnitude reference image [Nx Ny Nz Nt]
%   ph_AP_recon : AP velocity phase image
%   ph_FH_recon : FH velocity phase image
%   ph_RL_recon : RL velocity phase image
%   perc        : threshold parameters
%
% OUTPUT
%   recon_corr  : cell array containing corrected data
%                 {1} magnitude
%                 {2} AP velocity
%                 {3} FH velocity
%                 {4} RL velocity
% -------------------------------------------------------------------------

recon_corr{1} = recon_ref;
recon_corr{2} = ph_AP_recon;
recon_corr{3} = ph_FH_recon;
recon_corr{4} = ph_RL_recon;

[Nx, Ny, Nz, Nt] = size(recon_ref);

%% Step 1: Eddy current correction

v_ap = ph_AP_recon;
v_fh = ph_FH_recon;
v_rl = ph_RL_recon;

[v_ap_corr, ~] = eddyCorr_2ord_3d(v_ap, perc(1));
[v_fh_corr, ~] = eddyCorr_2ord_3d(v_fh, perc(1));
[v_rl_corr, ~] = eddyCorr_2ord_3d(v_rl, perc(1));

%% Step 2: Background noise removal

std_ap = std(v_ap_corr, 1, 4);
std_fh = std(v_fh_corr, 1, 4);
std_rl = std(v_rl_corr, 1, 4);

perc_noise = perc(2);

mask_ap = std_ap >= perc_noise * max(std_ap(:));
mask_fh = std_fh >= perc_noise * max(std_fh(:));
mask_rl = std_rl >= perc_noise * max(std_rl(:));

mask_combined = mask_ap + mask_fh + mask_rl;
mask_combined(mask_combined < 3) = 0;
mask_combined(mask_combined == 3) = 1;

mask_all = logical(repmat(mask_combined,[1,1,1,Nt]));

v_ap_corr(mask_all) = 0;
v_fh_corr(mask_all) = 0;
v_rl_corr(mask_all) = 0;

recon_corr{2} = v_ap_corr;
recon_corr{3} = v_fh_corr;
recon_corr{4} = v_rl_corr;

end