function nt = nyKTPT(lw,T)
lw = lw; %(DEBE DE ESTAR EN MICRAS)
%Applied Optocs Vol.42 No.33 20/11/2003
%Temperature-dependent dispersion equations
A0=6.2897*(10^-6);
A1=6.3061*(10^-6);
A2=-6.0629*(10^-6);
A3=2.6486*(10^-6);

B0=-0.14445*(10^-8);
B1=2.2244*(10^-8);
B2=-3.5770*(10^-8);
B3=1.3470*(10^-8);
%A0.+
n1=A1./(lw.^1);
n2=n1;

TERM11=A1./(lw.^1);
TERM12=A2./(lw.^2);
TERM13=A3./(lw.^3);
n1= A0+TERM11+TERM12+TERM13;
TERM21=B1./(lw.^1);
TERM22=B2./(lw.^2);
TERM23=B3./(lw.^3);
n2=B0+TERM21+TERM22+TERM23;
nt = n1*(T-25)+n2*(T-25).^2;