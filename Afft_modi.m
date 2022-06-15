function [Y]= Afft_modi(X)
global S2 n1 n2 n3 n kk ww gg

    tmp = reshape( fft2( reshape(X, [n1 n2]) ), [n,1]) ;
    Y=zeros(gg(kk),1);
     Y(1:gg(kk))= tmp(S2(1:gg(kk),(ww-1)*n3+kk));
     Y=reshape(Y,[gg(kk),1]);
    
  

