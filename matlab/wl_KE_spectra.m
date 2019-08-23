% Compute kinetic energy wavelet spectrum

grp = {'noF_start','P_force','R_force','V_force'};
Nruns = [2 17 17 15];
dname={'/U1/','/U2/','/U3/','/TH1/'};

N = 1024;   NK=floor(N/3);      K=[0:N/2-1 -N/2:-1];    KY=0:NK;
Ri_t = 1;

for b=2:4
% for b=1:4
%     for a=1:Nruns(b)
    for a=20
        
%         fname = ['/media/cjh225/sdc2/v2_peta4/Feb2019/' ...
%             grp{b} '/t' int2str(a) '/start.h5'];
        fname = ['/local/scratch/public/cjh225/v2_peta4/Feb2019/' ...
            grp{b} '/end.h5'];
        
        WTR=cell(1,4);
%         E_large=cell(1,4);  E_small=cell(1,4);  E=cell(1,4);
%         Fr=cell(1,4);   Fr_shear=cell(1,4);
        
        for n=[2 1 3 4]
            
            % Load variable from h5 file
            S1 = h5read(fname,dname{n});
            disp('Loaded variable')
           
            
            % WAVELET SPECTRUM CALCULATION
            [WAVE,PERIOD,SCALE,COI] = wavelet(S1(1,:,1)/N,1/N);
            WT = zeros(size(WAVE));
            parfor i=1:N
                for k=1:N
                    WAVE = wavelet(S1(i,:,k)/N,1/N);
                    WT = WT + abs(WAVE).^2;
                end
                if mod(i,100)==1; disp(['Completed loop ' int2str(i)]); end
            end
            WTR{n}=WT/N^2;
            Kw = 1./SCALE;
            disp(['Calculated ' dname{n}(2:end-1) ' wavelet spectrum'])
           
            
            % FOURIER SPECTRUM (WITH SCALE SEPARATION)
%             if n==2; Kh2=Ri_t/mean(S1(:).^2); end % Horizontal cutoff for large scale spectrum
%             CS1 = fftn(S1/N^3); clear S1
%             CS1=0.5*CS1.*conj(CS1);     % Calculate 2-sided energy spectrum
%             for i=NK+2:N-NK             % Dealias the high wavenumbers
%                 CS1(i,:,:)=0;
%                 CS1(:,i,:)=0;
%                 CS1(:,:,i)=0;
%             end
%             CS2=zeros(N,NK+1,N);
%             CS2(:,1,:)=CS1(:,1,:);
%             for j=2:NK+1
%                 CS2(:,j,:)=CS1(:,j,:)+CS1(:,N-j+2,:);
%             end
%             if n==4;    CS2=Ri_t*CS2;   end
%             E_large{n}=zeros(1,NK+1); E_small{n}=zeros(1,NK+1);
%             for i=1:N
%                 for k=1:N
%                     if K(i)^2+K(k)^2<Kh2
%                         E_large{n}=E_large{n}+CS2(i,:,k);
%                     else
%                         E_small{n}=E_small{n}+CS2(i,:,k);
%                     end
%                 end
%             end
%             E{n}=sum(sum(CS2,3),1);
%             if mod(n,2)==1
%                 Fr{n}=KY.^2.*E{n}/Ri_t;
%                 Fr_shear{n}=KY.^2.*CS2(1,:,1)/Ri_t;
%             end
%             disp(['Computed spectra for ' dname{n}(2:end-1)])
%             clear CS1 CS2
            
        end
        
%         save(['Spectra/' grp{b}(1) num2str(a,'%02u') '.mat'],'Kw','WTR', ...
%             'PERIOD','SCALE','COI','KY','E_large','E_small','E', ...
%             'Fr','Fr_shear')
        save(['Spectra/no_mean/' grp{b}(1) num2str(a,'%02u') '.mat'], ...
            'Kw','WTR','PERIOD','SCALE','COI')
        disp(['Saved wavelet and Fourier spectra for ' fname(37:end-9)])
        
    end
end