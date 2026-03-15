function Tinv=invQT3d_v2(opt,num)
% function to calculate the explicit inverse for TV

p1=opt.p1;

kx=opt.size(1);
ky=opt.size(2);
kz=opt.size(3);
nt=opt.size(4);

   
 T4z=zeros(kx,ky,kz,nt);
 T4z(1,1,1,1:2)=[1;-1];
    
 T4=fftn(T4z);


% make phase shift in Fourier domain to move circushift inside
ph4 = 0:nt-1;

shift4= reshape(exp(1j*2*pi*ph4/nt),[1,1,1,nt]);

Tinv= opt.p1*(-T4.*shift4.*T4)+p1+opt.p3*num;

%Tinv= opt.p1*(-T4.*shift4.*T4)+opt.p3*num;



end
