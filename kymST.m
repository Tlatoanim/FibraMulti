function kw =kymST(w,m,T)
m=m+1;
load('../constants.mat');
addpath('../KTP_TEMP');
dnS=[0.00515 0.00379 0.00289]';
lw = 2.*pi.*c./w;
kw = (w./c).*(nyKTPT(lw,T)+dnS(m,1));
