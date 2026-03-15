function [output_im] = AllReg(num,im,SRegiResult)
output_im = zeros([size(im),num]);
output_im(:,:,:,:,1) = im;

for i = 1:size(output_im,4)
    for j = 2:num
        output_im(:,:,:,i,j) = Registrate(squeeze(im(:,:,:,i)),SRegiResult{j-1},1);
    end
end

end