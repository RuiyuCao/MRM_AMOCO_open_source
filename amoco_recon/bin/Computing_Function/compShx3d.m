function Shx=compShx3d(Is,sMaps)
% function to calculate Sh*x
% input Is is the multi-coil images
% sMaps is same through nt


x=conj(sMaps).*Is;
Shx= squeeze(sum(x,5));
end