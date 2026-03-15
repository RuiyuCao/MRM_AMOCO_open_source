function At = At_fhp3D(z, sMaps,opt)

At=zeros(opt.size);
At(opt.U)=z;

At = compShx3d(opt.IFT(At),sMaps,opt);

At=At(:);
end

function Shx=compShx3d(Is,sMaps,opt)
% function to calculate Sh*x
% input Is is the multi-coil images
% sMaps is same through nt

% 2013/05/07 Yi Guo

x=conj(repmat(sMaps,[1 1 opt.size(3) 1])).*Is;
Shx=sum(x,4);
end