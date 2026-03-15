function ind = patchInd(X,B)

[nx, ny, nz, nt] = size(X);
nx_mod = nx-mod(nx,B);
ny_mod = ny-mod(ny,B);
nz_mod = nz-mod(nz,B);

%Find linear indices within the nx and ny dimension 

%Patch indices in x-dim
%x = repmat([1:B]',[B^2 1]);
x = repmat([1:B]', [ny_mod*nz_mod 1]);
inc = B*ones(numel(x), floor(nx_mod/B)-1);
x = [x inc];
x = cumsum(x,2);
x=x(:);

%Patch indicies in y-dim
y = repmat([1:B]',[1 B])'; y = y(:);
y = repmat(y,[nz_mod,1]);
inc = B*ones(numel(y), floor(ny_mod/B)-1);
y = [y, inc];
y = cumsum(y,2);
y = y(:); y = repmat(y,[nx_mod/B,1]);

%Patch indicies in z-dim
z = ones(B^2, nz_mod);
z = cumsum(z,2);
z = repmat(z(:), [floor(nx_mod/B)*floor(ny_mod/B) 1]);

ind = zeros(nx_mod*ny_mod*nz_mod*nt,1);

% Stack all the x and y indices including shift for
% overlapping patches

	% Convert to linear indicies 
	sub_ind = sub2ind([nx ny nz], x, y, z);
    
    temp_ind = sub_ind+ (0:(nt-1)).*(nx*ny*nz);

	ind(:,1) = temp_ind(:);
