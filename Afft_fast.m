function [tmp]= Afft_fast(X)
global S2 n1 n2 n3 r mk m n ww gg  q kk
for k=1:1:length(X(1,:))
    tmp(:,k) = reshape( fft2( reshape(X(:,k), [n1 n2]) ), [n,1]) ;
end


% Y=zeros(mk(kk),r,1);
% Y(1:mk(kk),:,1)= tmp(S2(1:mk(kk),kk),:);






