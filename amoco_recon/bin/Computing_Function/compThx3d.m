function Thx=compThx3d(x,opt)
% Calculate 2d spatial transposed Total variation for 3d data set

begin=1;

kx=opt.size(1);
ky=opt.size(2);
kz=opt.size(3);
nt=opt.size(4);
NL=kx*ky*kz*nt;
 
    TV3=reshape(x(begin:begin+NL-1),[kx,ky,kz,nt]);
    TVh3=TV3-circshift(TV3,[0, 0,0, -1]);

    Thx=TVh3;

    

end