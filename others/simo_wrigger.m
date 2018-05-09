clc
clear
close all
tic
u_plot=[];
 f_plot=[];
%  mov = VideoWriter('simo_wrigger.avi'); 
L = 1;
nelm = 10;
nnod = nelm+1;
x = 0:L/nelm:L;
y = x*0;
x2 = x; y2 = y+1; x2(1) = []; y2(1) = [];
x1 = y; y1 = x;
x = [x1 x2]*300; y = [y1 y2]*300;
nnod = length(x); nelm = nnod-1;
for lmn = 1:nelm
   icon(lmn,:) = [lmn lmn+1];
end
x0 = x;
y0 = y;
EA=1e6;EI=2e5;GA=1e5;
u =  zeros(nnod*3,1);
delu =  zeros(nnod*3,1);
theta = u(3:3:end);
niter = 100;
nstep = 10;
for istep = 1:nstep
    istep
for iter = 1:niter
%     iter
[Kg Fg_int Fg]=planar_simo(x, y, icon, u, EA, EI, GA, nnod, nelm);

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

ibc = [1 2 nnod*3-2 nnod*3-1  ]';
bc = [0  0 0 0 ]' ;
% % % ibc = [1 2 3*12-1  nnod*3-2 nnod*3-1  ]';
% % % if iter == 1
% % %     bc = [0 0 -130/nstep 0 0 ]' ;
% % % else
% % %     bc = [0 0 0 0 0 ]' ;
% % % end
% Kg(1:2, :) = 0;Kg(:,1:2) = 0; Kg(end-2:end-1, :) = 0;Kg(:,end-2:end-1) = 0; 
% Kg(1,1) = 1; Kg(2,2) = 1; Kg(end-2,end-2) = 1;  Kg(end-1,end-1) = 1; 
Fg(3*12-1) = -45*istep/nstep ;
Res = (Fg-Fg_int);

Res = Res-Kg(:,ibc)*bc;
Res(ibc) = bc; 
Kg(:,ibc)=0; Kg(ibc,:) = 0; Kg(ibc,ibc) = eye(length(ibc));


% Kg(1:3, :) = 0;Kg(:,1:3) = 0; 
% Kg(1,1) = 1; Kg(2,2) = 1; Kg(3,3) = 1; 
% Fg(end-2) = -pi^2*EI/L^2*istep/nstep/3;
% Res = (Fg-Fg_int);
% Res([1:3  ]) = 0; 
% if min(eig(Kg))<0
%     min(eig(Kg))
%     Fg(end-2)/(pi^2*EI/L^2/4)
%     break
% end


% disp('normalised residue = ')
% norm(Res)/norm(Fg)
if norm(Res)/norm(Fg) < 1e-4
%     norm(Res)/norm(Fg)
    [istep iter]
    break
end
[u ] = planar_simo_solve_update(Kg, Res, u);
 
%  pause(.2)

end
 plot(x0'+u(1:3:end),y0'+u(2:3:end),'k-o');  axis equal
 xlim([-100 300])
%  frame=getframe(gcf);
%  mov = addframe(mov,frame); 
 hold on
u_plot=[u_plot -u(3*12-1)];
f_plot=[f_plot -Fg(3*12-1)];
% plot(-u(3*12-1), -Fg_int(3*12-1),'o')
% hold on
% x = x+u(1:3:end)';
% y = y+u(2:3:end)';
% u = u*0;
% u

end
toc
% mov = close(mov);
%  movie2avi(M,'WaveMovie.avi'); 
% plot(u_plot,f_plot)
% wrig_ana= [0  0  
%     1.2562041041418723, 6.844106463878319
% 4.1416809373155665, 13.840304182509485
% 8.654882496928185, 20.684410646387832
% 14.787294768718734, 25.703422053231932
% 24.57918517013516, 29.809885931558924
% 35.79214195183778, 33.15589353612167
% 47.409127410288406, 35.893536121672994
% 59.22735320581663, 38.17490494296578
% 71.85982836521252, 40.45627376425854
% 81.84212308555617, 41.977186311787065
% 93.45291653363522, 43.498098859315576
% 104.24946061784654, 45.019011406844086
% 112.60016060526902, 45.93155893536121
% 123.38122466355131, 44.4106463878327
% 128.85341382946817, 39.69581749049429
% ]
% plot(wrig_ana(:,1), wrig_ana(:,2))
