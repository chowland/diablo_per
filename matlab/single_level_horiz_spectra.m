fname = [rundir 'start.h5'];
dname = {'/U1/', '/U2/', '/U3/', '/TH1/'};

J = 726;     % index associated with height in the domain
N = 1024;

for n=1:4
    S1 = h5read(fname,dname{n});
    S1 = reshape(S1(:,J,:),[N,N]);
    
    CS1 = fftn(S1/N^2);
    
    NK = floor(N/3);
    KH = 0:NK;
    
    CS1 = 0.5*CS1.*conj(CS1);
    CS1(NK+2:N-NK,:) = 0;
    CS1(:,NK+2:N-NK) = 0;
    
    K = [0:N/2-1 -N/2:-1];
    K2 = K.^2;
    
    CS2=zeros(NK+1);
    for i=1:N
        for k=1:N
            m=round(sqrt(K2(i)+K2(k)));
            if m<=NK
                CS2(m+1) = CS2(m+1) + CS1(i,k);
            end
        end
    end
    
    Eh{n}=CS2;
    
    disp(['Computed spectra for ' dname{n}]);
end

Eht=Eh{1}+Eh{2}+Eh{3}+Eh{4};
figure;
loglog(KH,Eht,'k'); hold on
loglog(KH,Eh{1},'r--')
loglog(KH,Eh{2},'g--')
loglog(KH,Eh{3},'b--')

loglog(KH,0.1*KH.^(-5/3),'k--')
set(gcf,'Renderer','painters');
xlabel('$\kappa$')
ylabel('$E(\kappa)$')
title(['$z = ' num2str(yvec(J)) ', \ \log_{10}\overline{\epsilon} = ' num2str(log10(epsilon(1,J))) '$ (P)']);