function metric = calc_gradEntropy_metric(img1, img2, name1, name2)

% calc_gradEntropy_metric
% ------------------------------------------------------------
% Compute gradient entropy metric for two 4D images and plot
% the comparison across cardiac phases.
%
% INPUT
%   img1  : 4D image [nx ny nz nphs]
%   img2  : 4D image [nx ny nz nphs]
%   name1 : label for image1 (string)
%   name2 : label for image2 (string)
%
% OUTPUT
%   metric : structure containing mean entropy curves
%
% Example
%   metric = calc_gradEntropy_metric(imgR_uncorrect, imgR_amoco,...
%                                    'uncorrect','amoco');
% ------------------------------------------------------------

[nx,ny,nz,nphs] = size(img1);

grid = [ceil(nx/4), ceil(ny/4), ceil(nz/4), ceil(nphs/4)];

% compute gradient entropy
mmetric1 = gradEntropy(img1,grid);
mmetric2 = gradEntropy(img2,grid);

metric1_mean = zeros(1,nphs);
metric2_mean = zeros(1,nphs);

for i = 1:nphs

    data4metric = mmetric1(:,:,:,i);
    metric1_mean(i) = mean(data4metric(:));

    data4metric = mmetric2(:,:,:,i);
    metric2_mean(i) = mean(data4metric(:));

end

% plot
figure

x = 1:nphs;

hold on

plot(x, metric1_mean,'-s','LineWidth',1.5,'MarkerSize',8)
plot(x, metric2_mean,'-s','LineWidth',1.5,'MarkerSize',8)

xlabel('Cardiac Phase')
ylabel('Gradient Entropy')

title('Gradient Entropy Comparison')

legend(name1,name2)

hold off

% output
metric.metric1_mean = metric1_mean;
metric.metric2_mean = metric2_mean;

end
