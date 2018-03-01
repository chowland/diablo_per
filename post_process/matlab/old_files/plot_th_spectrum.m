eta=0.111;
Pr=300;
l_b=eta/sqrt(300);
L=0.009;

NBINS=170;
%NBINS=1024;
%NBINS=256;
delta_x=2*pi/512;
clear turbo;
clear kmag;
clear th_spectrum;
turbo=dlmread('th_spectrum.txt');
%turbo=dlmread('energy_spectrum.txt');

nk=floor(length(turbo(:,1))/(NBINS+1))

row=1;

for k=1:nk
  tii_k(k)=turbo(row,1);
  row=row+1;
  for i=1:NBINS
    kmag(i)=turbo(row,1);
    counter(i,k)=turbo(row,2);
    th_spectrum(i,k)=turbo(row,3);
    row=row+1;
  end
end   

%surf(tii_k,kmag,log(th_spectrum));
%h=gca;
%set(h,'Yscale','log');
%figure

for k=1:10:nk
h=loglog(kmag(1:NBINS-1)*eta/sqrt(Pr)/mult,kmag(1:NBINS-1)'.*th_spectrum(1:NBINS-1,k)*mult,'b-');
hold on
color_num=min(0.9,1-k/nk);
set(h,'Color',[1-0.1*color_num color_num color_num ]);
end

%for k=1:nk
%l_test(k)=trapz(kmag,squeeze(th_spectrum(:,k)))/trapz(kmag,kmag.*squeeze(th_spectrum(:,k))');
%l_integral(k)=(pi/2)*trapz(kmag,squeeze(th_spectrum(:,k))')/thvar(k);
%t_integral(k)=l_integral(k)*0.009/V_s*mult^2;
%end





