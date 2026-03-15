function [ HFEN_m ] = hfen( U,r )
%HFEN Summary of this function goes here
%This metric has been derived from the implementation in Ravishankar
% et al , IEEE TMI 2011.
% Basically it quantifies the error made in fine features like edges,
% details etc.
% U - reconstructed image time series
% r - ideal image time series


for t = 1:size(r,3)

%% scale factor such that recon and ideal data are at same scale
alpha = sum(dot(U(:,:,t),r(:,:,t)))/(sum(dot(U(:,:,t),U(:,:,t))));
U(:,:,t)=(alpha)*U(:,:,t);



%% the hfen error for one time frame
HFEN_m(t)=norm(imfilter(abs(U(:,:,t)),fspecial('log',15,1.5))-imfilter(abs(r(:,:,t)),fspecial('log',15,1.5)),'fro')./norm(imfilter(abs(r(:,:,t)),fspecial('log',15,1.5)),'fro');

end

%% mean of the hfen error from all time frames
HFEN_m=mean(HFEN_m);