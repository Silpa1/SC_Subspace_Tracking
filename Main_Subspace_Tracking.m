clc;
clear all;close all;
global n1 n2 n m  S2 gg q  jj n3 kk pp ww;
filenames={'lowres_speech.mat'};
ns=[4,8,16,32,64];

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
        Xhat=[];
        for ww=1:1:ns(ii)
            gg=mk((ww-1)*n3+1:ww*n3);
            m=max(gg);
            Ys=Y(1:m,(ww-1)*n3+1:ww*n3);
            L=[];
            T=5;
            if ww==1
                T=70;
            end
            [Xbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Ys,0,1e-36,10);
            Ybar_hat=Afft(Xbar_hat);
            Ybar_hat=reshape(Ybar_hat,[m,n3]);
            Yinter=Ys-Ybar_hat;
            if ww==1
                [Uhat]=initAltGDMin(Yinter);
            end
            [Uhat, Bhat]=GDMin_wi(T,Uhat,Yinter);
            xT=Uhat*Bhat;
            L(:,1:n3)=xT+Xbar_hat;
            Ymec=Ys-Afft(L);
            E_mec=[];
            for kk=1:1:n3
                E_mec(:,kk)=cgls_modi(@Afft_modi,@At_modi, Ymec(1:gg(kk),kk) ,0,1e-36,30);
            end
            X_out(:,1:n3)=L+E_mec;
            Xhat=[Xhat,L+E_mec];
        end
        Xhat_mat=reshape(Xhat,[n1, n2,q]);
        Time_ST=  toc;
        Error_ST=RMSE_modi(Xhat_mat,X_image);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subspace_Tracking_Fastest %%%%
        pp=n3;
        Xhat=[];
        tic;
        gg=mk(1:n3);
        m=max(gg);
        Ys=Y(1:m,1:n3);
        ww=1
        L=[];
        T=5;
        if ww==1
            T=70;
        end
        [Xbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Ys,0,1e-36,10);
        Ybar_hat=Afft(Xbar_hat);
        Ybar_hat=reshape(Ybar_hat,[m,n3]);
        Yinter=Ys-Ybar_hat;
        [Uhat]=initAltGDMin(Yinter);
        [Uhat, Bhat]=GDMin_wi(T,Uhat,Yinter);
        xT=Uhat*Bhat;
        L(:,1:n3)=xT+Xbar_hat;
        Ymec=Ys-Afft(L);
        E_mec=[];
        for kk=1:1:n3
            E_mec(:,kk)=cgls_modi(@Afft_modi,@At_modi, Ymec(1:mk(kk),kk) ,0,1e-6,3);
        end
        ww=2;
        Xhat_f=[];
        Xhat_f(:,1:n3)=L+E_mec;
        X_out_mat=reshape(Xhat_f,[n1,n2,n3]);
            RMSE_modi(X_out_mat,X_image(:,:,1:n3));
        gg=mk(n3+1:q);
        tmp=Afft_fast(Uhat);
        tmp1=Afft_fast(Xbar_hat);
        for kk=1:1:q-n3
            Yhat=Y(1:mk(n3+kk),n3+kk)-tmp1(S2(1:mk(n3+kk),n3+kk),:);
            AU= tmp(S2(1:mk(n3+kk),n3+kk),:);
            Bhat=AU\Yhat;
            Xhat=Uhat* Bhat;
            Yhathat=Yhat-Afft_modi(Xhat);
            E_mec=cgls_modi(@Afft_modi,@At_modi, Yhathat ,0,1e-6,30);
            Xhatzst(:,kk)=Xbar_hat+Xhat+E_mec;
        end
        X_hat=[Xhat_f,Xhatzst];
        Xhat_mat=reshape(X_hat,[n1, n2,q]);
        Time_STf=  toc;
        Error_STf=RMSE_modi(Xhat_mat,X_image);
        
        fprintf(fid, '%d & %8.4f (%5.2f) & %8.4f (%5.2f)\n', ns(ii),Error_ST,Time_ST,Error_STf,Time_STf);
    end
    fclose(fid);


