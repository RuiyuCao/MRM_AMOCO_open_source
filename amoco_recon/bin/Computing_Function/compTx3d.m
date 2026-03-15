function [Tx,LL]=compTx3d(x,opt)
x=reshape(x,[opt.size(1),opt.size(2),opt.size(3),opt.size(4)]);
  
TV3=x-circshift(x,[0, 0,0, 1]);
Tx=[TV3(:)];
LL=numel(x);

end