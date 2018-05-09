clc
clear
close all
u = 0;
lam = 0;
f=1;
u0=u;
lam0=lam;
ds = 0.1
lam = 0.11;
for i = 1:10
   k = [cos(u), -f; 2*(u-u0),  2*(lam-lam0)] ;
   fg=[lam*f-sin(u) ;...
       (lam-lam0)^2+(u-u0)^2-ds^2 ];
   du = inv(k)*fg;
   u = u + du(1);
   lam = lam+du(2);
end