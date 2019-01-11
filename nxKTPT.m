function nt = nxKTPT(lw,T)
lw = lw; %ESTO DEBE DE ESTAR EN MICRAS)
%Applied Optics Vol.42 No.33 20/11/2003
%En este programa se considera que x=z (en el art'iculo)
A0=9.9587*(10^-6);
A1=9.9228*(10^-6);
A2=-8.9603*(10^-6);
A3=4.1010*(10^-6);

B0=-1.1882*(10^-6);
B1=10.459*(10^-6);
B2=-9.8136*(10^-6);
B3=3.1481*(10^-6);

TERM11=A1./(lw.^1);
TERM12=A2./(lw.^2);
TERM13=A3./(lw.^3);
n1= A0+TERM11+TERM12+TERM13;
TERM21=B1./(lw.^1);
TERM22=B2./(lw.^2);
TERM23=B3./(lw.^3);
n2=B0+TERM21+TERM22+TERM23;
nt = n1*(T-25)+n2*(T-25).^2;