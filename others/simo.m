clc
clear
close all


L = 1;
EA=0.1; EI=0.1^3;GA=0.1;
nelm = 10;
x = 0:L/nelm:L;
y = x*0;
nnod = length(x); nelm = nnod-1;
for lmn = 1:nelm
   icon(lmn,:) = [lmn lmn+1];
end

u =  zeros(nnod*3,1);
delu =  zeros(nnod*3,1);
niter = 50;
nstep = 1;
for istep = 1:nstep
    istep
for iter = 1:niter
%     iter

[Kg Fg_int Fg]=planar_simo(x, y, icon, u, EA, EI, GA, nnod, nelm);

Kg(1:3, :) = 0;Kg(:,1:3) = 0;
Kg(1,1) = 1; Kg(2,2) = 1; Kg(3,3) = 1; 
           rad = L/pi;   Fg(3*nnod) =EI/rad*2*istep/nstep ;
% Fg(end-1) =EI*istep/nstep *9;
Res = (Fg-Fg_int);
Res(1:3) = 0;



if norm(Res)/norm(Fg) < 1e-6
    [istep iter]
    break
end
u = planar_simo_solve_update(Kg, Res, u);

 plot(x'+u(1:3:end),y'+u(2:3:end),'-o');  axis equal
 
 pause(.2)

 hold on
end

end
