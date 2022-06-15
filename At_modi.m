function x = At_modi(y)
global S2 n1 n2 kk n gg ww n3


w = zeros(n,1);

w(S2(1:gg(kk),(ww-1)*n3+kk))=y(1:gg(kk));


tmp3 = n*reshape( ifft2( reshape(w, [n1 n2]) ), [n,1]) ;
x= tmp3;

