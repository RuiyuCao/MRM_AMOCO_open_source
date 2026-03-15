
function Y = svt(X,lambda,p)
[U,S,V] = svd(X,'econ'); s = diag(S);
thres=(lambda).*(s.^(p-1));
s=(s-thres);
s = s.*(s>0);
       
Y = U*diag(s)*V';


end


