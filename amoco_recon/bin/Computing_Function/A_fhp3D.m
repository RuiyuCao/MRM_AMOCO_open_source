function A = A_fhp3D(z, sMaps,opt)

z=reshape(z,[opt.size(1),opt.size(2),opt.size(3)]);
A = opt.FT(compSx3d(z,sMaps,opt));

A=A(opt.U);

end

function res=compSx3d(x,sMaps,opt)
res=repmat(sMaps,[1 1 opt.size(3) 1]).*repmat(x,[1 1 1 opt.size(4)]);
end
