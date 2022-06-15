function x = Att(y)
global S2 n1 n2 n3  m n gg ww 
w = zeros(n,n3);
y_mat=reshape(y, [m,n3]);
for k=1:1:n3
    if ww==1
        w(S2(1:gg(k),((ww-1)*n3)+k),k)=y_mat(1:gg(k),k);
    else
        w(S2(1:gg(k),((ww-1)*n3)+k),k)=y_mat(1:gg(k),k);
    end
    tmp3 = n*reshape( ifft2( reshape(w(:,k), [n1 n2]) ), [n,1]) ;
    x(:,k)= tmp3;
end
%Y1vec=reshape(X1vec,[n*q,1]);