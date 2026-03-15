function [output_im] = AllReg_n(num,im)
output_im = zeros([size(im),num]);
for i = 1:num
    output_im(:,:,:,:,i) = im;
end
end