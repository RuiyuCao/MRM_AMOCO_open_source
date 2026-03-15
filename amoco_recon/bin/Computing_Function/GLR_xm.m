function [g] = GLR_xm(img, lambda, p)
        
        [n1,n2,n3]=size(img);
        
        [u,sigma,v] = givefastSVD(reshape(img,n1*n2,n3));
        s=diag(sigma);
      
       
        thres=(lambda).*(s.^(p-1));
        
        s=(s-thres);
        s = s.*(s>0);
        g=u*(diag(s))*v';
        g=reshape(g,n1,n2,n3);

end