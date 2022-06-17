clc;
clear all;close all;
global n1 n2 n m  S2 gg q  jj n3 kk pp ww;
filenames={'lowres_speech.mat'};
ns=[2,4,8,16,32,64];

[fid,msg] = fopen('Subspace_Tracking.txt','wt');
fprintf(fid, '%s & %s    \n','Subspace','AltGDMin MRI');
S = load('lowres_speech.mat');
radial=[16];
X_image=double(cell2mat(struct2cell(S)));
[n1,n2,q]=size(X_image);
n=n1*n2;
X_mat=reshape(X_image,[n,q]);
[mask]=goldencart(n1,n2,q,16);
mask = fftshift(fftshift(mask,1),2);
mask3=reshape(mask,[n1*n2, q]);
mk=[];
for i=1:1:q
    mk(i)=length(find(logical(mask3(:,i))));
    S2(1:mk(i),i)=double(find(logical(mask3(:,i))));
end
m=max(mk);
Y=zeros(m,q);
for k=1:1:q
    ksc = reshape( fft2( reshape(X_mat(:,k), [n1 n2]) ), [n,1]) ;
    Y(1:mk(k),k)=double(ksc(S2(1:mk(k),k)));
end
for ii=1:1:length(ns)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Subspace_Tracking (70,5)%%%%
    n3=q/ns(ii);
    pp=n3;
    tic;
    Zhat=[];
    for ww=1:1:ns(ii)
        gg=mk((ww-1)*n3+1:ww*n3);
        m=max(gg);
        Ys=Y(1:m,(ww-1)*n3+1:ww*n3);
        T=5;
        if ww==1
            T=70;
        end
        [zbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Ys,0,1e-36,10);
        Ytemp=reshape(Afft(zbar_hat),[m,n3]);
        Ybar=Ys-Ytemp;
        if ww==1
            [U0]=initAltGDMin(Ybar);
            Uhat=U0;
        end
        [Uhat, Bhat]=AltGDmin(T,Uhat,Ybar);
        xT=Uhat*Bhat;
        Ymec=Ys-Afft(xT+zbar_hat);
        Ehat=[];
        for kk=1:1:n3
            Ehat(:,kk)=cgls_modi(@Afft_modi,@At_modi, Ymec(1:gg(kk),kk) ,0,1e-36,30);
        end
        Zhat=[Zhat,Ehat+xT+zbar_hat];
    end
    Zhat_Minibatch=reshape(Zhat,[n1, n2,q]);
    Time_Minibatch=  toc;
    Error_Minibatch=RMSE_modi(Zhat_Minibatch,X_image);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subspace_Tracking_Fastest %%%%
    pp=n3;
    Zhat=[];
    tic;
    gg=mk(1:n3);
    m=max(gg);
    Ys=Y(1:m,1:n3);
    ww=1
    T=70;
    [zbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Ys,0,1e-36,10);
    Ytemp=reshape(Afft(zbar_hat),[m,n3]);
    Ybar=Ys-Ytemp;
    [U0]=initAltGDMin(Ybar);
    Uhat=U0;
    [Uhat, Bhat]=AltGDmin(T,Uhat,Ybar);
    xT=Uhat*Bhat;
    Ymec=Ys-Afft(xT+zbar_hat);
    Ehat=[];
    for kk=1:1:n3
        Ehat(:,kk)=cgls_modi(@Afft_modi,@At_modi, Ymec(1:gg(kk),kk) ,0,1e-36,30);
    end
    
    Zhat=Ehat+xT+zbar_hat;
    ww=2;
    gg=mk(n3+1:q);
    tmp=Afft_fast(Uhat);
    tmp1=Afft_fast(zbar_hat);
    for kk=1:1:q-n3
        Yhat=Y(1:mk(n3+kk),n3+kk)-tmp1(S2(1:mk(n3+kk),n3+kk),:);
        AU= tmp(S2(1:mk(n3+kk),n3+kk),:);
        Bhat=AU\Yhat;
        Xhat=Uhat* Bhat;
        Yhathat=Yhat-Afft_modi(Xhat);
        Ehat=cgls_modi(@Afft_modi,@At_modi, Yhathat ,0,1e-6,30);
        Zhat(:,n3+kk)=zbar_hat+Xhat+Ehat;
    end
    Zhat_Online=reshape(Zhat,[n1, n2,q]);
    Time_Online=  toc;
    Error_Online=RMSE_modi(Zhat_Online,X_image);
    
    fprintf(fid, '%d & %8.4f (%5.2f) & %8.4f (%5.2f)\n', ns(ii),Error_Minibatch,Time_Minibatch,Error_Online,Time_Online);
end
fclose(fid);


