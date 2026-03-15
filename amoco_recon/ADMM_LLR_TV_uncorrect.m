function [imgR,opt] = ADMM_LLR_TV_uncorrect(mat,maps,tranpNtop1,opt,img0,Info)
% recon for uncorrected data
kx = opt.size(1);
ky = opt.size(2);
kz = opt.size(3);
nt = opt.size(4);
nc = opt.size(5);

%% ---------- GPU control ----------
useGPU = opt.gpu;

if useGPU
    mat  = gpuArray(mat);
    maps = gpuArray(maps);
    img0 = gpuArray(img0);
    opt.Up = gpuArray(opt.Up);
end

%% ---------- prepare coil maps ----------
maps = repmat(maps,[1,1,1,1,1,opt.num]);

%% ---------- sampling mask ----------
U = opt.Up;

%% ---------- initial values ----------
[v2,LL] = compTx3d(img0,opt);

v1  = img0;
old = img0;
z_k = img0;

compare = 1000;

if useGPU
    e12 = gpuArray.zeros(LL,1,opt.class);
    e11 = gpuArray.zeros(kx,ky,kz,nt,opt.class);
else
    e12 = zeros(LL,1,opt.class);
    e11 = zeros(kx,ky,kz,nt,opt.class);
end

iter = 1;

%% ---------- ADMM iteration ----------
while iter <= opt.itermax

    tic

    old_compare = compare;

    old_v2 = v2;
    old_v1 = v1;

    %% Z update
    z_k = Z_update_nreg(z_k,U,mat,maps,v1,e11,v2,e12,opt);

    %% convergence
    compare = mat_compare(z_k,old);

    opt.compare = compare;

    if compare > old_compare || compare < 0.001
        break
    end

    old = z_k;

    %% TV update
    v2 = shrink1(compTx3d(z_k,opt) - e12 , opt.lambdaTV/opt.p1);

    %% LLR update
    v1 = LLRp_v3(z_k - e11 , opt.B , opt.lambdaLLR/opt.p1 , opt.p , useGPU);

    %% dual update
    e12 = e12 - compTx3d(z_k,opt) + v2;
    e11 = e11 - z_k + v1;

    %% rho update
    opt = rho_update(v1,old_v1,z_k,v2,old_v2,opt);

    iter = iter + 1;

    elapsedTime = toc;

    disp(['iteration:',num2str(iter-1),...
        '  p:',num2str(gather(opt.p1)),...
        '  Compare:',num2str(gather(compare)),...
        '  elapsedTime:',num2str(elapsedTime),'s'])

end

%% ---------- output ----------
imgR = gather(reshape(z_k,[kx,ky,kz,nt]));

end
