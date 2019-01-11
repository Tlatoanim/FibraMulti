clear all; close all;  clc;
addpath('../KTP_crystal');
addpath('../KTP_TEMP');
addpath('../20130828_PPKTP');
load('../constants.mat');
load('BS.mat');
%%%%%%%%%%%%%%%%%GRID DE FRECUENCIAS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%PAR�METROS DEL CRISTAL (INICIALES)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lktp = 10.11;    %  periodo
kktp = 2.*pi./Lktp;
L = 10*(1e3);            % long. cristal
 Lp0 = .405;  
FWHM = .001;
c=3.*10.^14;
FWHMw = (2*pi*c/(Lp0^2))*FWHM;
sigma = FWHMw/(sqrt(2*log(2)));
wp0 = 2*pi*c/Lp0;
% sigma=3.*(1e13);
lammin=0.800;
lammax=0.820;
ommin=2.*pi.*c./lammax;
ommax=2.*pi.*c./lammin;
Lp01 = 0.405028;
Lp02 = 0.405291;
Lp03 =0.405544;         %Lp0=.404+0.002*(r-1)/4;%Lp0=0.407; % bombeo
FWHM1 = .000122969;
FWHM2 = .000157202;
FWHM3 = .000198404;     % ancho
c=3.*10.^14;
FWHMw1 = (2.*pi.*c/(Lp01.^2)).*FWHM1;
FWHMw2 = (2.*pi.*c/(Lp02.^2)).*FWHM2;
FWHMw3 = (2.*pi.*c/(Lp03.^2)).*FWHM3;
sigma1 = FWHMw1/(2.*sqrt(2.*log(2)));
sigma2 = FWHMw2/(2.*sqrt(2.*log(2)));
sigma3 = FWHMw3/(2.*sqrt(2.*log(2)));
                        % sigma=2*pi*c./FWHMw;
wp01 = 2.*pi.*c./Lp01;      %wp0=wp0*(1e-15);
wp02 = 2.*pi.*c./Lp02;  
wp03 = 2.*pi.*c./Lp03;

N=600;
omsvec=linspace(ommin,ommax,N);%% genera vector de N entradas igualmente espaciados, incluyendo a ommin y ommax
omivec=linspace(ommin,ommax,N);
[ws,wi]=meshgrid(omsvec,omivec);
ws=ws./1e15;
wi=wi./1e15;
%%%%%%%%%%%%%%%%%%%%C�LCULO DE JSI (Dependencia de la temperatura)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%a=0.6*(1e-6);
alpha=6.7.*(1e-6);
beta=11.*(1e-9);
th=pi/2; ph=0.0;
%T=37;
% T=50; 

M=100;
TA=linspace(23,70,M);
%mov(M) = struct('cdata',[],'colormap',[]);
T=TA
for j=1:M
    T=TA(j);
   j
    %T=38.5;

LkT=Lktp.*(1.0+alpha.*(T-25.0)+beta.*(T-25.0).^2);
kktp=2.*pi./LkT;
LT=L.*(1.0+alpha.*(T-25.0)+beta.*(T-25.0).^2);

%LT=L.*exp(a.*(T));

Sp=[sin(th).*cos(ph); sin(th).*sin(ph);cos(th)];

kp = kwfsT1(1e15.*(ws+wi),Sp,T);     % modo 00
ki = kwfsT1(1e15.*(ws),Sp,T);        % modo 00
ks = kwssT1(1e15.*(wi),Sp,T);        % modo 00

Dk1 = kp - ks - ki - kktp;


PMF1 = (sin(LT.*Dk1./2)./(LT.*Dk1./2)).*exp(-1i.*(LT.*Dk1./2));
%PE = exp(-(((ws + wi) - wp0).^2)./(sigma.^2));
PE1 = exp(-((1e15.*(ws + wi) - wp01).^2)./(2.*sigma1.^2)); %degenerado a 23 grados
PE2=  exp(-((1e15.*(ws + wi) - wp02).^2)./(2.*sigma2.^2)); %degenerado a 37.5
PE3=  exp(-((1e15.*(ws + wi) - wp03).^2)./(2.*sigma3.^2)); %degenerado a 50
%PE2=PE1;
%PE3=PE1;
PE=(0.395916.*PE1+0.975567.*PE2+0.677332.*PE3)./...
    (0.395916+0.975567+0.677332);
%lfi=2.*pi.*c./wp02-0.0001; lff=2.*pi.*c./wp02+0.0001;
%wfi=2.*pi.*c./lff; wff=2.*pi.*c./lfi;
%filtro=(heaviside(1e15.*(ws+wi)-wfi)-heaviside(1e15.*(ws+wi)-wff));
%filtro=1.0
%pcolor(filtro)
%PE=filtro.*PE;

a=0.000005;
lfi=2.*pi.*c./(wp03)-a; lff=2.*pi.*c./(wp03)+a;
wfi=2.*pi.*c./lff; wff=2.*pi.*c./lfi;
wfi=wfi./1e15; wff=wff./1e15;
filtro=(heaviside((ws+wi)-wfi)-heaviside((ws+wi)-wff));
%filtro=1.0;
JSI = PE.*PMF1.*filtro;


str=num2str(T);


% f=sum(0.5*abs(JSI+transpose(JSI)).^2);%%%%Suma incoherente
wii=wi(:,1);

figure(2)
f=sum(JSI);
f=abs(f).^2; %%%%Suma coherente
waterfall((2.*pi.*c./wii).*1e-12,T,f)
hold all;
waterfall((2.*pi.*c./wii).*1e-12,T,abs(sum(transpose(JSI))).^2)
hold all;
xlabel('\lambda (nm)');
ylabel('T(C)');


%f=sum(abs(JSI).^2)
%plot((2.*pi.*c./wi),f);


%domega=wi(2)-(wi(1))
N=300;
NORM=sqrt(sum(sum(abs(JSI).^2)));

G11=JSI.*output11;
G12=JSI.*output12;
%G13=JSI.*output13;
G21=JSI.*output21;
G22=JSI.*output22;

GN=sqrt(sum(sum(abs(G11).^2+abs(G12).^2+abs(G21).^2+abs(G22).^2)));

G11=G11./GN;
G12=G12./GN;
%G13=JSI.*output13;
G21=G21./GN;
G22=G22./GN;

%G11P=(G11+transpose(G11))./sqrt(2);
%%G12P=(G12+transpose(G12))./sqrt(2);
%G21P=(G21+transpose(G21))./sqrt(2);
%G22P=(G22+transpose(G22))./sqrt(2);
figure(1)
% pcolor(wi*(1e3),ws*(1e3), 0.5*abs(JSI+transpose(JSI)).^2);
% pcolor(wi*(1e3),ws*(1e3), abs((JSI)).^2);
pcolor(1e-12.*2.*pi.*c./wi, 1e-12.*2.*pi.*c./ws, abs(G12).^2);
%mov(j)=getframe(gcf);
shading interp
axis square
colorbar
title(['JSI(\lambda_{s},\lambda_{x|x|i}) en ', str,'C']);
xlabel('\lambda_{s} (nm)');
ylabel('\lambda_{i} (nm)');
%mov(j) = getframe;


%G12P=(G12+G21)./(sqrt(2))
%G21P=(G12+G21)./(sqrt(2))
%G21=(G21+G12)./(sqrt(2))
%v(j)=sum(sum(2.*abs(G12).*abs(transpose(G21))))./...
%    (sum(sum(abs(G12).^2))+sum(sum(abs(transpose(G21)).^2)))
%SG=sum(sum(abs(G11+G12+G21+G22).^2));
%SG=sum(sum(abs(G11).^2))+sum(sum(abs(G12).^2))+sum(sum(abs(G21).^2))+...
%    sum(sum(abs(G22).^2))
taumax=5000;
taumin=-taumax;
delta=(taumax-taumin)/N;
for lm=1:N
    %a=2.60e-12;

tau(lm)=taumin+lm*delta;
G=exp(1i.*(ws-wi).*tau(lm));
N11=(abs(G11+G.*transpose(G11)).^2);
N22=(abs(G22+G.*transpose(G22)).^2);
%N33=abs(G33+G.*transpose(G33)).^2;
N12=(abs(G12+G.*transpose(G21)).^2);
%N13=abs(G13+G.*transpose(G31)).^2;
%N23=abs(G23+G.*transpose(G32)).^2;
N21=(abs(G21+G.*transpose(G12)).^2);
%HOM=(1./sqrt(2)).*abs(JSI./NORM+transpose(G.*JSI./NORM)).^2;   
%HOM3(lm)=sum(sum(HOM));

CN11(lm)=sum(sum(N11));
CN22(lm)=sum(sum(N22));
CN12(lm)=sum(sum(N12));
CN21(lm)=sum(sum(N21));
CN(lm)=CN11(lm)+CN12(lm)+CN22(lm)+CN21(lm);output11.^2;

AS(lm)=sum(sum(abs(G12).*abs(transpose(G21)).*cos(tau(lm).*(ws-wi)...
    -(angle(G12)-angle(transpose(G21))))));
HOM(lm)=2.*sum(sum(abs(G12).*abs(transpose(G21)).*cos(tau(lm).*(wi-ws)...
    -((angle(G12)-angle(transpose(G21)))))));
end
AMAX=max(max(AS));
AMIN=min(min(AS));
B=(sum(sum(abs(G12).^2))+sum(sum(abs(transpose(G21)).^2)));
v2(j)=(AMAX-AMIN)./(B+AMAX+AMIN);


%visibilidad
ncmax=max((CN12+CN21)./CN);
ncmin=min((CN12+CN21)./CN);
%ncmax=max(CN12)
%ncmin=min(CN12)
v3(j)=(ncmax-ncmin)./(ncmax+ncmin);


%dw=domega;

figure(4)   
lm=linspace(delta.*1,delta.*N,N);
%plot(1e-15.*c.*tau,(CN12+CN21)./CN)
waterfall(1e-15.*c.*tau,T,(0.5+HOM));
xlabel('\Delta l (um)');
ylabel('T(C)');
zlabel('p_c');
hold all;

figure(5)
waterfall(tau,T,(CN12));
xlabel('\Delta t (fs)');
ylabel('T(C)');
zlabel('p_c');
hold all;

figure(100)
waterfall(1e-15.*c.*tau,T,(CN12));
xlabel('\Delta t (fs)');
ylabel('T(C)');
zlabel('p_c');
hold all;
end

% figure(6)
% plot(TA,v);
% xlabel('T');
% ylabel('visibilidad');
% 
figure(7)
plot(TA,v2);
xlabel('T');
ylabel('visibilidad');

figure(8)
plot(TA,v3);
xlabel('T');
ylabel('visibilidad');
