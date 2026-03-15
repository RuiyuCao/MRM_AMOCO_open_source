function output_im = Registrate(im,SRegiResult,ForB)

[nx,ny,nz] = size(im);

output_im = im;

[X,Y,Z] = meshgrid(1:ny,1:nx,1:nz);

if isa(im,'gpuArray')
    X = gpuArray(X);
    Y = gpuArray(Y);
    Z = gpuArray(Z);
end

if ForB == 0

    deform_f(:,:,:,1) = SRegiResult(:,:,:,2,1);
    deform_f(:,:,:,2) = SRegiResult(:,:,:,1,1);
    deform_f(:,:,:,3) = SRegiResult(:,:,:,3,1);

    X_new = X + deform_f(:,:,:,1);
    Y_new = Y + deform_f(:,:,:,2);
    Z_new = Z + deform_f(:,:,:,3);

elseif ForB == 1

    deform_b(:,:,:,1) = SRegiResult(:,:,:,2,2);
    deform_b(:,:,:,2) = SRegiResult(:,:,:,1,2);
    deform_b(:,:,:,3) = SRegiResult(:,:,:,3,2);

    X_new = X + deform_b(:,:,:,1);
    Y_new = Y + deform_b(:,:,:,2);
    Z_new = Z + deform_b(:,:,:,3);

end
edge_cutx = floor(0.04*nx);
edge_cuty = floor(0.04*ny);
edge_cutz = floor(0.04*nz);
xs = edge_cutx:nx-edge_cutx;
ys = edge_cuty:ny-edge_cuty;
zs = edge_cutz:nz-edge_cutz;

output_im_tmp = interpn(Y,X,Z,im,...
                        Y_new(xs,ys,zs),...
                        X_new(xs,ys,zs),...
                        Z_new(xs,ys,zs),...
                        'linear',0);

output_im(xs,ys,zs) = output_im_tmp;

end
