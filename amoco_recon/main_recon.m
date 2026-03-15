%{
Github Repository: 

Reference: 
The algorithms implemented below corresponds to research
article published in Magnetic Resonance in Medicine (MRM) journal 
(doi:).

Publication Details:
Title: "Adaptive Respiratory Motion Correction Framework for Ungated Hepatic 4D Flow MRI"
Authors: Ruiyu Ca, Shaohang Li, Aiqi Sun, Bingyi Wang, Yuhui Nie, Hongwei Li, Junpu Hu, Ne Yang,
Jiayu Zhu, Ying-Hua Chu, Ning Jin, Xingfeng Shao, Ruokun Li, He Wang, Hao Li

The article provides an accessible overview of the research 
and may be cited for more detailed information pertaining to the 
algorithmic approach used herein.

How to cite:

% ========================================================================%

Written by:
ruiyu cao (rycao22@m.fudan.edu.cn)
%}

close all; 
clear;clc
addpath(genpath('bin'))

%% ------------------------- load data -------------------------
data_file = 'data/example_data_5T.mat';
% data_file = 'data/example_data_3T.mat';
load(data_file)

%===================================================================================================
% Data input should contain data like example data:
%   example_data:
%   kdata - kspace imaging data [kx, ky, kz, coil, time/phases, vel encoding]
%   maps  - coil sensitivy maps [kx,ky,kz,ciol]
%   respSplit_amoco.mask - samping mask for data for motion correction (with motion bins) 
%                          cell: respSplit_amoco.mask{n} [ky,kz,phase,vel encoding] for motion state n
%   respSplit_amoco.wn_eachbin - soft-gating weights for data for motion correction (with motion bins) 
%                          cell: respSplit_amoco.wn_eachbin{n} [ky,kz,phase,vel encoding] for motion state n
%   tranpNtop1_amoco - motion fields for motion state n to motion state 1 
%                      cell: tranpNtop1_amoco{1} motion fields for motion state 2 to motion state 1
%                            [nx,ny,nz,orientation(3 FH RL AP),direction(2 forward or backward)]
%                      cell: tranpNtop1_amoco{2} motion fields for motion state 3 to motion state 1
%                            [nx,ny,nz,orientation(3 FH RL AP),direction(2 forward or backward)]
%===================================================================================================
%% ------------------------- parameter setting -------------------------
Info.recon_phs_num = length(respSplit_amoco.mask); % number of motion states to recon for amoco
Info.d4itermax = 20; % max number of iteration
Info.d4B = 5; % block size for locally low rank                                     
Info.d4tv = 0.02; % lambda of tv        
Info.d4llr = 0.12; % lambda of llr   
Info.gpu = 1; % use gpu 1, cpu 0
%% ------------------------- prepare data and mask -------------------------
respSplit_uncorrect.mask = splitToNParts(abs(squeeze(kdata(floor(size(kdata,1)/2),:,:,1,:,:)))>0,Info.recon_phs_num);
tranpNtop1_uncorrect = 1;
[nx,ny,nz,ncha,nphs,nset] = size(kdata);
maps = reshape(maps,[nx,ny,nz,1,ncha]);
kdata_nphs = single(zeros([nx,ny,nz,ncha,nphs,nset,Info.recon_phs_num]));
kdata_nphs_amoco = single(zeros([nx,ny,nz,ncha,nphs,nset,Info.recon_phs_num]));
soft_gating_mask = single(zeros([ny,nz,nphs,nset,Info.recon_phs_num]));

for i = 1:Info.recon_phs_num
    kdata_nphs(:,:,:,:,:,:,i) = kdata.*reshape(respSplit_uncorrect.mask{i},[1,ny,nz,1,nphs,nset]);
    kdata_nphs_amoco(:,:,:,:,:,:,i) = kdata.*reshape(respSplit_amoco.wn_eachbin{i},[1,ny,nz,1,nphs,nset]);
    kdata_nphs_amoco(:,:,:,:,:,:,i) = kdata.*reshape(respSplit_amoco.mask{i},[1,ny,nz,1,nphs,nset]);
    soft_gating_mask(:,:,:,:,i) = respSplit_amoco.mask{i};
end
soft_gating_mask = reshape(permute(soft_gating_mask,[1,2,3,5,4]),[1,ny,nz,nphs,1,Info.recon_phs_num,nset]);

% undersampling factor for unccrect data
usf = fun_cal_usfactor(sum(kdata_nphs,7));
disp(['Undersampled Factor of uncrrected data:',num2str(usf)])
% undersampling factor for amoco data
usf_amoco = fun_cal_usfactor(sum(kdata_nphs_amoco,7));
disp(['Undersampled Factor of amoco data:',num2str(usf_amoco)])

%% ------------------------- Reconstruction -------------------------
kdata_1dfft_4d = permute(ifftc(kdata_nphs,1),[1,2,3,5,4,7,6]);
kdata_1dfft_4d_amoco = permute(ifftc(kdata_nphs_amoco,1),[1,2,3,5,4,7,6]);
imgR_uncorrect = single(zeros([nx,ny,nz,nphs,nset]));
imgR_amoco = single(zeros([nx,ny,nz,nphs,nset]));
for i = 1:size(kdata,6)
    % reconstruction for uncorrect image
    [img0,mat,opt] = pre_data(kdata_1dfft_4d(:,:,:,:,:,:,i),Info);
    [imgR_uncorrect(:,:,:,:,i),opt_uncorrect] = ADMM_LLR_TV_uncorrect(mat,maps,tranpNtop1_uncorrect,opt,img0);
    imgR_uncorrect(:,:,:,:,i) = gather(imgR_uncorrect(:,:,:,:,i)*opt_uncorrect.scale);
    clear img0 mat opt opt_uncorrect

    % reconstruction for amoco image
    [img0_amoco,mat_amoco,opt] = pre_data(kdata_1dfft_4d_amoco(:,:,:,:,:,:,i),Info);
    opt.Up = soft_gating_mask(:,:,:,:,:,:,i);
    [imgR_amoco(:,:,:,:,i),opt_amoco] = ADMM_LLR_TV_amoco(mat_amoco,maps,tranpNtop1_amoco,opt,img0_amoco);
    imgR_amoco(:,:,:,:,i) = gather(imgR_amoco(:,:,:,:,i)*opt_amoco.scale);
    clear img0_amoco mat_amoco opt opt_amoco
end
if contains(data_file, '3T')
    idx_cut = (round(size(imgR_uncorrect,1)/4)+1):round(3*size(imgR_uncorrect,1)/4);
    imgR_uncorrect = imgR_uncorrect(idx_cut, :, :, :, :)*1e8;
    imgR_amoco = imgR_amoco(idx_cut, :, :, :, :)*1e8;
    imgR_uncorrect = imgR_uncorrect(end:-1:1, end:-1:1, :, :, :);
    imgR_amoco = imgR_amoco(end:-1:1, end:-1:1, :, :, :);
end
%% ------------------------- coil uniform -------------------------
save_dir = "data/results_"+datestr(now, 'yyyy-mm-dd_HH-MM-SS');
img_show = cat(2,imgR_uncorrect,imgR_amoco);
figure,imshow3D(abs(img_show(:,:,:,floor(nphs/2),1)))

% % coil uniform (need python)
% mkdir(fullfile(save_dir,"coil_uniform"))
% h5_filename_uncorrect = fullfile(save_dir,"coil_uniform","uncorrect.h5");
% imgR_uncorrect = fun_coil_uniform(imgR_uncorrect,h5_filename_uncorrect);
% h5_filename_amoco = fullfile(save_dir,"coil_uniform","amoco.h5");
% imgR_amoco = fun_coil_uniform(imgR_amoco,h5_filename_amoco);
% img_show = cat(2,imgR_uncorrect,imgR_amoco);
% figure,imshow3D(abs(img_show(:,:,:,floor(nphs/2),1)))

%% ------------------------- metric -------------------------
metric = calc_gradEntropy_metric(imgR_uncorrect(:,:,:,:,1), imgR_amoco(:,:,:,:,1), 'imgR uncorrect', 'imgR amoco');

%% ------------------------- phase correction -------------------------
[phscorr_uncorrect,mra_uncorrect] = fun_phaseCorr_2ord_3d(imgR_uncorrect,[0.02,0.08]);
[phscorr_amoco,mra_amoco] = fun_phaseCorr_2ord_3d(imgR_amoco,[0.02,0.08]);

%% ------------------------- save -------------------------
mkdir(fullfile(save_dir,"recon_data"))
mkdir(fullfile(save_dir,"compare_metric"))
mkdir(fullfile(save_dir,"phsCorr_mra"))
% recon data 
save(fullfile(save_dir,"recon_data","uncorrect.mat"),"imgR_uncorrect")
save(fullfile(save_dir,"recon_data","amoco.mat"),"imgR_amoco")
saveas(1,fullfile(save_dir,"recon_data","magnitude_compare.fig"))

% compare metric
save(fullfile(save_dir,"compare_metric","gradient_entropy_metric.mat"),"metric")
saveas(2,fullfile(save_dir,"compare_metric","gradient_entropy_metric.fig"))

% results of phase correction and mra
save(fullfile(save_dir,"phsCorr_mra","phscorr_uncorrect.mat"),"phscorr_uncorrect")
save(fullfile(save_dir,"phsCorr_mra","phscorr_amoco.mat"),"phscorr_amoco")
save(fullfile(save_dir,"phsCorr_mra","mra_uncorrect.mat"),"mra_uncorrect")
save(fullfile(save_dir,"phsCorr_mra","mra_amoco.mat"),"mra_amoco")
saveas(3,fullfile(save_dir,"phsCorr_mra","phscorr_uncorrect_slice.fig"))
saveas(4,fullfile(save_dir,"phsCorr_mra","phscorr_uncorrect_time.fig"))
saveas(5,fullfile(save_dir,"phsCorr_mra","phscorr_uncorrect_mra.fig"))
saveas(6,fullfile(save_dir,"phsCorr_mra","phscorr_amoco_slice.fig"))
saveas(7,fullfile(save_dir,"phsCorr_mra","phscorr_amoco_time.fig"))
saveas(8,fullfile(save_dir,"phsCorr_mra","phscorr_amoco_mra.fig"))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function: splitToNParts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split a sampling pattern into N random subsets.
%
% INPUT:
%   pattern : logical matrix (sampling pattern, non‑zero positions)
%   n       : number of subsets
%
% OUTPUT:
%   patternsCell : 1×n cell array, each element is a logical matrix
%                  containing a subset of the original sampling pattern
%
% NOTE:
%   The union of all subsets equals the original pattern.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function patternsCell = splitToNParts(pattern, n)
    % Check if n is valid
    if n <= 0 || mod(n,1) ~= 0
        error('Number of parts n must be a positive integer');
    end
    % Get linear indices of non-zero entries
    nonZeroIdx = find(pattern);
    % Total number of sampling points
    numOnes = numel(nonZeroIdx);
    % Size of each subset
    partSize = floor(numOnes / n);
    % Random permutation of indices
    randomOrder = randperm(numOnes);
    % Initialize output
    patternsCell = cell(1,n);
    for i = 1:n
        if i < n
            partIdx = randomOrder((i-1)*partSize+1 : i*partSize);
        else
            % Last subset takes remaining points
            partIdx = randomOrder((i-1)*partSize+1 : end);
        end
        tempPattern = false(size(pattern));
        tempPattern(nonZeroIdx(partIdx)) = true;
        patternsCell{i} = tempPattern;
    end
    % Verify that union equals original pattern
    combinedPattern = false(size(pattern));
    for i = 1:n
        combinedPattern = combinedPattern | patternsCell{i};
    end
    assert(isequal(pattern,combinedPattern), ...
        'Output patterns do not match the input pattern.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function: pre_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare reconstruction data and parameters.
%
% INPUT:
%   data : k-space data after preprocessing
%          size = [nx ny nz nt nc resp]
%
%   Info : reconstruction parameter structure
%
% OUTPUT:
%   img0 : initial image estimate
%   mat  : normalized k-space data
%   opt  : reconstruction parameter structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img0,mat,opt] = pre_data(data,Info)
    opt = struct();
    % Initial image estimate (root-sum-of-squares reconstruction)
    img0 = sqrt(sum(abs(ifftc(ifftc(sum(data,6),2),3)).^2,5));
    % Scaling factor
    opt.scale = max(img0(:));
    % Normalize data
    mat = data / opt.scale;
    img0 = img0 / opt.scale;
    % Get dimensions
    [opt.nx,opt.ny,opt.nz,opt.nt,opt.nc,opt.resp] = size(mat);
    opt.size = size(mat);
    %% Reconstruction parameters
    opt.num       = Info.recon_phs_num;   % number of respiratory bins
    opt.itermax   = Info.d4itermax;       % ADMM iterations
    opt.B         = Info.d4B;             % LLR patch size
    opt.lambdaTV  = Info.d4tv;            % TV regularization weight
    opt.lambdaLLR = Info.d4llr;           % LLR regularization weight
    opt.p1 = 1;       % ADMM penalty parameter
    opt.p  = 0.5;     % Schatten-p norm
    opt.gpu = Info.gpu;
    opt.class = 'single';
    % Sampling mask
    opt.Up = abs(mat(floor(opt.nx/2),:,:,:,floor(opt.nc/2),:))>0;

end
