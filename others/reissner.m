clc
clear
close all

tic

u_plot=[];
f_plot=[];
%  mov = VideoWriter('simo_wrigger.avi'); 

L = 1;
nelm = 10;

x = 0:L/nelm:L;
y = x*0;
x2 = x; y2 = y+1; x2(1) = []; y2(1) = [];
x1 = y; y1 = x;
x = [x1 x2]*300; y = [y1 y2]*300;

nnod = length(x);
nelm = nnod-1;

icon = zeros(nelm,2);

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
    disp(istep);
    for iter = 1:niter
        [Kg, Fg_int, Fg] = planar_simo(x, y, icon, u, EA, EI, GA, nnod, nelm);

        ibc = [1 2 nnod*3-2 nnod*3-1  ]';
        bc = [0  0 0 0 ]' ;
        Fg(3*12-1) = -45*istep/nstep ;
        Res = (Fg-Fg_int);

        Res = Res-Kg(:,ibc)*bc;
        Res(ibc) = bc; 
        Kg(:,ibc)=0; Kg(ibc,:) = 0; Kg(ibc,ibc) = eye(length(ibc));
    
        if norm(Res)/norm(Fg) < 1e-4
            disp([istep iter]);
            break
        end
        
        [u ] = planar_simo_solve_update(Kg, Res, u);
    end
    plot(x0'+u(1:3:end),y0'+u(2:3:end),'k-o');  axis equal
    xlim([-100 300])
    hold on
    u_plot=[u_plot -u(3*12-1)];
    f_plot=[f_plot -Fg(3*12-1)];
end
toc