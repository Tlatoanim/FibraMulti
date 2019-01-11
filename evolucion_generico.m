clear all
close all
clc

rcore=4.5; % core radius
NA=0.14;   % numerical aperture
s=1.0*9;   % separation
lambda=linspace(0.5,1.8,200);  %wavelength

fname=strcat('./coeff9_014','.txt');%%aqu� est�n los coeficientes que se utilizan para evaluar la funci�n CoupCoeff
coeff=load(fname);   %% dispersion coefficients

kappa=CoupCoeff(coeff,rcore,NA,s,lambda); %(1/mm)
beta=polyval(coeff,lambda);%%coeff es el vector de coeficientes de un polinomio de lambda
%%el n�mero de coeficientes determina el grado del polinomio.

% figure
% subplot(2,1,1)
% plot(lambda,kappa)
% xlabel('lambda')
% ylabel('k')
% subplot(2,1,2)
% plot(lambda,beta)
% xlabel('lambda')
% ylabel('beta')

%k = polyval(coeff,lambda);
%figure
%plot(lambda,k)
%xlabel('hola')
c=3.*10.^14; %% velocidad de la luz en Micras
lamc=.810;  %%alguna longitud de onda
omc=2.*pi.*c./lamc; %% frecuencia angular de la luz
lammax=.780+.4; ommin=2.*pi.*c./lammax; %%l�mite m�nimo para la frecuencia
ommax=2.*omc-ommin; %%l�mite m�ximo para la omega
N=500; %% se calculan 200 puntos
omsvec=linspace(ommin,ommax,N);%% genera vector de N entradas igualmente espaciados, incluyendo a ommin y ommax
omivec=linspace(ommin,ommax,N);%% Lo mismo de arriba, no entiendo porque est� declarado de nuevo y se usa nuevo nombre

output11=zeros(N,N); %% Matriz cuadrada de puros ceros A CONTINUACI�N se declaran puras matrices cuadradas de ceros.
output22=output11;output33=output11; %% LOS DOS INDICES corresponden a la correlaci�n entre pares de fotones
output12=output11; output13=output11; output23=output11; %% 3^2=9 formas de tener dos fotones en tres n�cleos (tomando en cuenta que los dos fotones
output21=output11; output31=output11; output32=output11; %% pueden estar en el mismo n�cleo
vec=[0,0,1,0,0,0,0,0,0]; %% vector con 9 entradas

for ll=1:41
  ls=ll-1;  
zv=(ll-1)*250.0

for ii=1:N
    oms=omsvec(ii); lams=2.*pi.*c./oms; %%Define las longitudes de onda en las que se va a evaluar
    ns=n0(lams);
    CCs=ns*CoupCoeff(coeff,rcore,NA,s,lams); %%Calcula los coeficientes de acoplamiento para el signal en la longitud de onda lams
    for jj=1:N
        omi=omivec(jj); lami=2.*pi.*c./omi; %%Define las longitudes de onda en las que se va a evaluar
        ni=n0(lami);
        CCi=ni*CoupCoeff(coeff,rcore,NA,s,lami); %% Calcula los coeficientes de acoplamiento para el idler en la longitud de onda lami
        betap=ns*polyval(coeff,lams)+ni*polyval(coeff,lami); %% Calcula las Betas a partir de los polinomios formados por los coeficientes de dispersi�n en lams y lami
       
       
        A=[-betap,0,0,-CCi,-CCi,0,-CCs,-CCs,0;  %1  a1a1
            0,-betap,0,-CCs,0,-CCi,-CCi,0,-CCs;  %2   a2a2
            0,0,-betap,0,-CCs,-CCs,0,-CCi,-CCi;  %3   a3a3
            -CCi,-CCs,0,-betap,-CCi,0,0,0,-CCs;  %4   a1a2
            -CCi,0,-CCs,-CCi,-betap,-CCs,0,0,0;   %5   a1a3
            0,-CCi,-CCs,0,-CCs,-betap,-CCi,0,0;   %6   a2a3
            -CCs,-CCi,0,0,0,-CCi,-betap,-CCs,0;    %7  a2a1
            -CCs,0,-CCi,0,0,0,-CCs,-betap,-CCi;   %8   a3a1
            0,-CCs,-CCi,-CCs,0,0,0,-CCi,-betap ]; %9   a3a2
        
        B=expm(-i.*A.*zv)*transpose(vec); %Se definen todas las soluciones para el Hamiltoniano de "nearest neighbour"
        output11(ii,jj)=B(1); % Se refiere a dos fotones en el núcleo 1
        output22(ii,jj)=B(2); % Dos fotones en el núcleo 2
        output33(ii,jj)=B(3); % Dos fotones en el núcleo 3
        output12(ii,jj)=B(4); % Un fotón en el núcleo 1 y un fotón en el núcleo 2
        output13(ii,jj)=B(5); % Un fotón en el núcleo 1 y un fotón en el núcleo 3
        output23(ii,jj)=B(6); % Un fotón en el núcleo 2 y otro en el núcleo 3
        output21(ii,jj)=B(7); % Un fotón en el núcleo 2 y otro en el núcleo 1
        output31(ii,jj)=B(8); % Un fotón en el núcleo 3 y otro en el núcleo 1
        output32(ii,jj)=B(9); % Un fotón en el núcleo 3 y otro en el núcleo 2
    end

end

omsvec=omsvec/(1e15)
omivec=omivec/(1e15)
fig=figure(2)
subplot(3,3,1)
pcolor(omsvec,omivec,abs(output11+transpose(output11)).^2)
shading interp
%xlabel('\omega_s')
%ylabel('\omega_i')
axis 'square'
caxis( [0 1] )


subplot(3,3,2)
pcolor(omsvec,omivec,abs(output12+transpose(output21)).^2)
shading interp
%xlabel('\omega_s')
%ylabel('\omega_i')
axis 'square'
caxis( [0 1] )

subplot(3,3,3)
pcolor(omsvec,omivec,abs(output13+transpose(output31)).^2)
shading interp
%xlabel('\omega_s')
%ylabel('\omega_i')
axis 'square'
caxis( [0 1] )

subplot(3,3,4)
pcolor(omsvec,omivec,abs(output21+transpose(output12)).^2)
shading interp
%xlabel('\omega_s')
%ylabel('\omega_i')
axis 'square'
caxis( [0 1] )

subplot(3,3,5)
pcolor(omsvec,omivec,abs(output22+transpose(output22)).^2)
shading interp
%xlabel('\omega_s')
%ylabel('\omega_i')
axis 'square'
caxis( [0 1] )

subplot(3,3,6)
pcolor(omsvec,omivec,abs(output23+transpose(output32)).^2)
shading interp
%xlabel('\omega_s')
%ylabel('\omega_i')
axis 'square'
caxis( [0 1] )

subplot(3,3,7)
pcolor(omsvec,omivec,abs(output31+transpose(output13)).^2)
shading interp
%xlabel('\omega_s')
%ylabel('\omega_i')
axis 'square'
caxis( [0 1] )


subplot(3,3,8)
pcolor(omsvec,omivec,abs(output32+transpose(output23)).^2)
shading interp
%xlabel('\omega_s')
%ylabel('\omega_i')
axis 'square'
caxis( [0 1] )


subplot(3,3,9)
pcolor(omsvec,omivec,abs(output33+transpose(output33)).^2)
shading interp
%xlabel('\omega_s')
%ylabel('\omega_i')
axis 'square'
caxis( [0 1] )
omsvec=omsvec*(1e15)
omivec=omivec*(1e15)
M(ll) = getframe(fig);
fileName= ['evolgen33_' num2str(ls) '.mat']
save (fileName);
end
hf = figure;

movie(hf,M);
mplay(M,1)
movie2avi(M, 'evolgen33.avi', 'compression','None', 'fps',1)

%% COMENTADO POR ALI
%% FIN COMENTADO POR ALI

% 
% figure(3)
% subplot(1,3,1)
% pcolor(
% 
% 
% output11=zeros(1,N);
% output22=output11;output33=output11;
% output12=output11; output13=output11; output23=output11;
% output21=output11; output31=output11; output32=output11;
% 
% N=200;
% 
% c=8.*10.^14;% figure(3)
% subplot(1,3,1)
% pcolor(omsvec,omivec,abs(output12+output21).^2)
% shading interp
% axis 'square'
% caxis( [0 1] )
% 
% 
% subplot(1,3,2)
% pcolor(omsvec,omivec,abs(output13+output31).^2)
% shading interp
% axis 'square'
% caxis( [0 1] )
% 
% subplot(1,3,3)
% pcolor(omsvec,omivec,abs(output23+output32).^2)
% shading interp
% axis 'square'
% caxis( [0 1] )
% 
% figure(4) 
% pcolor(omsvec,omivec,abs(output11).^2+abs(output22).^2+abs(output33).^2+ ...
% abs(output12).^2+abs(output13).^2+abs(output23).^2+ ...
% abs(output21).^2+abs(output31).^2+abs(output32).^2)
% shading interp
% caxis( [0 1] )

% lamc=.810;
% omc=2.*pi.*c./lamc;
% lammax=.810+.5; ommin=2.*pi.*c./lammax;
% ommax=2.*omc-ommin;
% 
% omsvec=linspace(ommin,ommax,N);
% 
% zv=30000;
%         
% for ii=1:N
%     oms=omsvec(ii); lams=2.*pi.*c./oms;
%     CCs=CoupCoeff(coeff,rcore,NA,s,lams);
%     
%         lami=lams;
%         CCi=CoupCoeff(coeff,rcore,NA,s,lami);
%         betap=polyval(coeff,lams)+polyval(coeff,lami);
%         
%         A=[-betap,0,0,-CuntitleCi,-CCi,0,-CCs,-CCs,0;  %1  a1a1
%             0,-betap,0,-CCs,0,-CCi,-CCi,0,-CCs;  %2   a2a2
%             0,0,-betap,0,-CCs,-CCs,0,-CCi,-CCi;  %3   a3a3
%             -CCi,-CCs,0,-betap,-CCi,0,0,0,-CCs;  %4   a1a2
%             -CCi,0,-CCs,-CCi,-betap,-CCs,0,0,0;   %5   a1a3
%             0,-CCi,-CCs,0,-CCs,-betap,-CCi,0,0;   %6   a2a3
%             -CCs,-CCi,0,0,0,-CCi,-betap,-CCs,0;    %7  a2a1
%             -CCs,0,-CCi,0,0,0,-CCs,-betap,-CCi;   %8   a3a1
%             0,-CCs,-CCi,-CCs,0,0,0,-CCi,-betap ]; %9   a3a2
%         
%         B=expm(-i.*A.*zv)*transpose(vec);
%         output11(ii)=B(1);
%         output22(ii)=B(2);
%         output33(ii)=B(3);
%         output12(ii)=B(4);
%         output13(ii)=B(5);
%         output23(ii)=B(6);
%         output21(ii)=B(7);
%         output31(ii)=B(8);
%         output32(ii)=B(9);
%         
%       
%    
% end
% 
% figure(3)
% subplot(3,3,1)
% plot(abs(output11).^2)
% shading interp
% 
% 
% subplot(3,3,2)
% plot(abs(output12).^2)
% shading interp
% 
% subplot(3,3,3)
% plot(abs(output12).^2)
% shading interp
% 
% subplot(3,3,4)
% plot(abs(output12).^2)
% shading interp
% 
% 
% subplot(3,3,5)
% plot(abs(output13).^2)
% shading interp
% 
% 
% subplot(3,3,6)
% plot(abs(output23).^2)
% shading interp
% 
% subplot(3,3,7)
% plot(abs(output21).^2)
% shading interp
% 
% 
% subplot(3,3,8)
% plot(abs(output31).^2)
% shading interp
% 
% 
% subplot(3,3,9)
% plot(abs(output32).^2)
% shading interp



% %%%%%%%%% as a function of z
%         
% 
% N=200;
% 
% 
% zmin=0;
% zmax=10000;
% zvec=linspace(zmin,zmax,N);
% 
% for ii=1:N
%     lams=.810;oms=2.*pi.*c./lams;
%     CCs=CoupCoeff(coeff,rcore,NA,s,lams);
%     
%         lami=lams;
%         CCi=CoupCoeff(coeff,rcore,NA,s,lami);
%         betap=polyval(coeff,lams)+polyval(coeff,lami);
%         zv=zvec(ii);
%         
%         A=[-betap,0,0,-CCi,-CCi,0,-CCs,-CCs,0;  %1  a1a1
%             0,-betap,0,-CCs,0,-CCi,-CCi,0,-CCs;  %2   a2a2
%             0,0,-betap,0,-CCs,-CCs,0,-CCi,-CCi;  %3   a3a3
%             -CCi,-CCs,0,-betap,-CCi,0,0,0,-CCs;  %4   a1a2
%             -CCi,0,-CCs,-CCi,-betap,-CCs,0,0,0;   %5   a1a3
%             0,-CCi,-CCs,0,-CCs,-betap,-CCi,0,0;   %6   a2a3
%             -CCs,-CCi,0,0,0,-CCi,-betap,-CCs,0;    %7  a2a1
%             -CCs,0,-CCi,0,0,0,-CCs,-betap,-CCi;   %8   a3a1
%             0,-CCs,-CCi,-CCs,0,0,0,-CCi,-betap ]; %9   a3a2
%         
%         B=expm(-i.*A.*zv)*transpose(vec);
%         output11(ii)=abs(B(1)).^2;
%         output22(ii)=abs(B(2)).^2;
%         output33(ii)=abs(B(3)).^2;
%         output12(ii)=abs(B(4)).^2;
%         output13(ii)=abs(B(5)).^2;
%         output23(ii)=abs(B(6)).^2;
%         output21(ii)=abs(B(7)).^2;
%         output31(ii)=abs(B(8)).^2;
%         output32(ii)=abs(B(9)).^2;
%    
% end
% 
% figure(4)
% subplot(3,3,1)
% plot(output11)
% shading interp
% 
% 
% subplot(3,3,2)
% plot(output12)
% shading interp
% 
% subplot(3,3,3)
% plot(output12)
% shading interp
% 
% subplot(3,3,4)
% plot(output12)
% shading interp
% 
% 
% subplot(3,3,5)
% plot(output13)
% shading interp
% 
% 
% subplot(3,3,6)
% plot(output23)
% shading interp
% 
% subplot(3,3,7)
% plot(output21)
% shading interp
% 
% 
% subplot(3,3,8)
% plot(output31)
% shading interp
% 
% 
% subplot(3,3,9)
% plot(output32)
% shading interp
%             
%             
% 
