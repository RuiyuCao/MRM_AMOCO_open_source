function [output_im] = AllRegH(num,im,SRegiResult)

output_im = squeeze(im(:,:,:,:,1));

for i = 1:size(im,4)
    for j = 2:num
        regRes = Registrate(squeeze(im(:,:,:,i,j)),SRegiResult{j-1},0);
        output_im(:,:,:,i) = output_im(:,:,:,i) + regRes;
    end
end
end