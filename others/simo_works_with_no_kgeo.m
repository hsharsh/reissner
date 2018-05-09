clc
clear
close all
tic
L = 1;
nelm = 10;
nnod = nelm+1;
x = 0:L/nelm:L;
y = x*0;
x2 = x; y2 = y+1; x2(1) = []; y2(1) = [];
x1 = y; y1 = x;
x = [x1 x2]*300; y = [y1 y2]*300;
nnod = length(x); nelm = nnod-1;
x0 = x;
y0 = y;
% EA = .1; EI = .1^3; GA = .1;
EA = 1e6; EI = 2e5; GA = 1e5;
x_gp = [-1 1]/sqrt(3);
w_gp = [1 1];
% xgp = 0;
% w_gp=2;


% x  =  sin(x0/1);
% y  =  1 - cos(x0/1);

u =  zeros(nnod*3,1);
delu =  zeros(nnod*3,1);
theta = u(3:3:end);
niter = 500;
nstep = 10;
for istep = 1:nstep
    istep
for iter = 1:niter
%     iter
Fg_int = u *0; 
Kg = zeros(3*nnod);
Fg = zeros(3*nnod,1);

total_length = 0;
for lmn = 1:nelm
    node1 = lmn; node2 = lmn+1;
    iv = [3*(node1-1)+1 3*(node1-1)+2 3*(node1-1)+3 ...
          3*(node2-1)+1 3*(node2-1)+2 3*(node2-1)+3 ];
    Uloc = u(iv);
    xi = 0;
    [Bmat B0 ENN Rmat Tmat len G_N G_Q Uloc]= Bmats(Uloc, x, y, lmn, theta, xi);
    strain = Rmat * B0 * Uloc + (Rmat -eye(3)) * ENN ;
    Dmat = [EA 0 0;0 GA 0; 0 0 EI];
    stress(lmn).comp = Dmat * strain;
    Fg_int(iv) = Fg_int(iv) +  Tmat'*  Bmat'*stress(lmn).comp * len  ;
    total_length = total_length + len;
end


for lmn = 1:nelm
    node1 = lmn; node2 = lmn+1;
    iv = [3*(node1-1)+1 3*(node1-1)+2 3*(node1-1)+3 ...
          3*(node2-1)+1 3*(node2-1)+2 3*(node2-1)+3 ];
         
    Uloc = u(iv);
    Kl = zeros(6);Kgeo = zeros(6);
    for igp = 1:length(w_gp)
        Dmat = [EA 0 0;0 0 0; 0 0 EI];
        [Bmat B0 ENN Rmat Tmat len G_N G_Q Uloc_junk]= Bmats(Uloc, x, y, lmn, theta, x_gp(igp));
        Kl = Kl + Bmat'*Dmat*Bmat *len/2*w_gp(igp);
        Kgeo = Kgeo + (stress(lmn).comp(1) * G_N) *len/2*w_gp(igp);
    end
    Dmat = [0 0 0;0 GA 0; 0 0 0];
    [Bmat B0 ENN Rmat Tmat len G_N G_Q Uloc_junk] = Bmats(Uloc, x, y, lmn, theta, 0);
    Kl = Kl + Bmat'*Dmat*Bmat *len/2*2;
    Kgeo = Kgeo + (stress(lmn).comp(2) * G_Q) *len/2*2;
                
    Kl = Kl + Kgeo*1;
    Kg(iv,iv) = Kg(iv,iv) +    Tmat'*Kl*Tmat ;
    
end

% Kg(1:3, :) = 0;Kg(:,1:3) = 0;
% Kg(1,1) = 1; Kg(2,2) = 1; Kg(3,3) = 1; 
%            rad = L/pi;   Fg(3*nnod) =EI/rad*2*istep/nstep ;
% % Fg(end-1) =EI*istep/nstep *9;
% Res = (Fg-Fg_int);
% Res(1:3) = 0;

% Kg(1:2, :) = 0;Kg(:,1:2) = 0;
% Kg(1,1) = 1; Kg(2,2) = 1;% Kg(3,3) = 1; 
% Kg(end-2, :) = 0;Kg(:,end-2) = 0;
% Kg(end-2,end-2) = 1; %Kg(3,3) = 1; 
% Fg(end-1) =-EI*9*istep/nstep   ;
% Res = (Fg-Fg_int);
% Res(1:2) = 0;
% Res(end-2) = 0;

Kg(1:2, :) = 0;Kg(:,1:2) = 0; Kg(end-2:end-1, :) = 0;Kg(:,end-2:end-1) = 0; 
Kg(1,1) = 1; Kg(2,2) = 1; Kg(end-2,end-2) = 1;  Kg(end-1,end-1) = 1; 
Fg(3*12-1) = -45*istep/nstep ;
Res = (Fg-Fg_int);
Res([1:2  end-2:end-1]) = 0; 

% disp('normalised residue = ')
% norm(Res)/norm(Fg)
if norm(Res)/norm(Fg) < 1e-6
%     norm(Res)/norm(Fg)
    [istep iter]
    break
end
delu = pinv(Kg)* Res;
% delu'
u = u + delu;
theta = theta + delu(3:3:end);

%  plot(x0'+u(1:3:end),y0'+u(2:3:end),'-o');  axis equal
 
%  pause(.2)

%  hold on
end
plot(-u(3*12-1), -Fg(3*12-1),'o')
% x = x+u(1:3:end)';
% y = y+u(2:3:end)';
% u = u*0;
% u
hold on 

end
toc
wrig_ana= [1.2562041041418723, 6.844106463878319
4.1416809373155665, 13.840304182509485
8.654882496928185, 20.684410646387832
14.787294768718734, 25.703422053231932
24.57918517013516, 29.809885931558924
35.79214195183778, 33.15589353612167
47.409127410288406, 35.893536121672994
59.22735320581663, 38.17490494296578
71.85982836521252, 40.45627376425854
81.84212308555617, 41.977186311787065
93.45291653363522, 43.498098859315576
104.24946061784654, 45.019011406844086
112.60016060526902, 45.93155893536121
123.38122466355131, 44.4106463878327
128.85341382946817, 39.69581749049429
]
plot(wrig_ana(:,1), wrig_ana(:,2))
