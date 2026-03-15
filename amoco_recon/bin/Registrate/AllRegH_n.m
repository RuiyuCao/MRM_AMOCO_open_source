function [output_im] = AllRegH_n(num,im)

output_im = squeeze(im(:,:,:,:,1));

for i = 1:size(im,4)
    for j = 2:num
        regRes = squeeze(im(:,:,:,i,j));
        output_im(:,:,:,i) = output_im(:,:,:,i) + regRes;
    end
end
end

