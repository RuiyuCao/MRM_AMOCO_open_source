function Z = Z_update_nreg(Z,U,mat,maps,v1,e11,v2,e12,opt)

row = opt.size(1);
col = opt.size(2);
nz  = opt.size(3);
echo = opt.size(4);

right = right1();
grad = left1(Z) - right;

conj1 = -grad;

for i = 1:5

    temp1 = left1(conj1);

    temp = conj(conj1) .* temp1;

    alpha = sum(abs(grad(:)).^2) / abs(sum(temp(:)));

    Z = Z + alpha * conj1;

    grad1 = grad + alpha * temp1;

    beta = sum(abs(grad1(:)).^2) / sum(abs(grad(:)).^2);

    conj1 = -grad1 + beta * conj1;

    grad = grad1;

    if sqrt(sum(abs(grad(:)).^2)) < 1e-9
        break
    end

end


    function l1 = left1(Z)

        tmp = AllReg_n(opt.num,Z);

        tmp = reshape(tmp,[row,col,nz,echo,1,opt.num]);

        tmp = maps .* tmp;

        tmp = fftc(fftc(tmp,3),2);

        tmp = U .* tmp;

        tmp = ifftc(ifftc(tmp,3),2);

        tmp = sum(conj(maps).*tmp,5);

        tmp = AllRegH_n(opt.num,tmp);

        l1 = tmp + opt.p1/2 * Z + opt.p1/2 * compThx3d(compTx3d(Z,opt),opt);

    end


    function r1 = right1()

        tmp = ifftc(ifftc(U.*mat,2),3);

        tmp = sum(conj(maps).*tmp,5);

        tmp = AllRegH_n(opt.num,tmp);

        r1 = tmp + opt.p1/2*(v1+e11) + opt.p1/2*compThx3d(v2+e12,opt);

    end

end
