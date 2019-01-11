    clear all 
close all
clc
format long e
% =========================================================================
rcore=4.5; % core radius
NA=0.14;   % numerical aperture
s=1.0*11.0;   % separation
lambda=linspace(0.5,1.8,200);  %wavelength
% =========================================================================
vec=[0,1,0]; %% vector con 3 entradas correspondiente a estados de entrada
lams=0.810;
 fname=strcat('./coeff9_014','.txt');%%aqu� est�n los coeficientes que se utilizan para evaluar la funci�n CoupCoeff
    coeff=load(fname);  
    kappa=CoupCoeff(coeff,rcore,NA,s,lambda); %(1/mm)
    beta=polyval(coeff,lambda);%%coeff es el vector de coeficientes de un polinomio de lambda
    %el n�mero de coeficientes determina el grado del polinomio.   
    CC=CoupCoeff(coeff,rcore,NA,s,lams); %%Calcula los coeficientes de acoplamiento para el signal en la longitud de onda lams
    betap=polyval(coeff,lams);
    ns=n0(lams);  
for lm=1:10001
    ls=lm-1;
    zv=(lm-1)*10^0.0; %(micras) %zv=0;
   
      
        
    A=[-betap,-CC,-CC;
       -CC,-betap,-CC;
       -CC,-CC,-betap];
        
        %% LA B es un vector que representa el estado vec evolucionado a lo largo de z
        %% Las output representan matrices para cada estado de entrada, evolucionado a una distancia zv y evaluado en las frecuencias ii y jj
        B1=expm(-i.*ns.*A.*zv)*transpose(vec); %Se definen todas las soluciones para el Hamiltoniano de "nearest neighbour"
        output(lm,1)=B1(1);
        output(lm,2)=B1(2);
        output(lm,3)=B1(3);
        output(lm,4)=0;
        
end
 x = linspace(1,lm,lm);
 y = linspace(1,4,4);
 [X Y] = meshgrid(x,y);
% surf(X,Y, transpose(abs(output).^2));
pcolor(X,Y, transpose(abs(output).^2));
shading flat
colorbar
%shading interp
set(gca,'YTick',[1:1:3])
set(gca,'FontSize',16, 'FontName', 'Helvetica');
%TeXString = texlabel('lambda12^(3/2)/pi - pi*delta^(2/3)');
xlabel('longitud (\mu m)','FontSize',14,...
       'FontWeight','bold','Color','k')
ylabel('guia','FontSize',14,...
       'FontWeight','bold','Color','k')
%fileName= ['evo_ini' num2str(ls) 'cm_clasp.mat']
print -depsc evo_uno_810.eps
fileName= ['evo_uno_810.mat']
save (fileName);

