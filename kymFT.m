function kw =kymFT(w,m,T)
m=m+1;
load('../constants.mat');
addpath('../KTP_TEMP');
dnF = [0.00238 0.00089 0.00018]';
lw = 2.*pi.*c./w;
kw = (w./c).*(nyKTPT(lw,T)+dnF(m,1));
