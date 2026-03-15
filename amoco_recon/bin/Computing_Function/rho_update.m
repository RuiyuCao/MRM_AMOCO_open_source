function opt = rho_update(v1,old_v1,z_k,v2,old_v2,opt)

temp1 = 10;
temp2 = 2;

ans1 = z_k - v1;
ans2 = compTx3d(z_k,opt) - v2;

rZ = sqrt( sum(abs(ans1(:)).^2) + sum(abs(ans2(:)).^2) );

ans1 = old_v1 - v1;
ans2 = compThx3d(v2 - old_v2,opt);

sZ = opt.p1 * sqrt( sum(abs(ans1(:)).^2) + sum(abs(ans2(:)).^2) );

if rZ > temp1*sZ
    opt.p1 = opt.p1 * temp2;

elseif sZ > temp1*rZ
    opt.p1 = opt.p1 / temp2;

else
    opt.p1 = opt.p1;

end

end
