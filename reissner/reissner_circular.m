clc
clear
close all

tic

u_plot=[];
f_plot=[];

L = 1;
nelm = 10;

x = (0:L/nelm:L)*300;
y = x*0;

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
niter = 1;
nstep = 100;

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';

for istep = 1:nstep
    istep
    for iter = 1:niter
        [Kg, Fg_int, Fg] = planar_simo(x, y, icon, u, EA, EI, GA, nnod, nelm);

        ibc = [1 2 3]';
        bc = [0  0 0]' ;
        Fg(3*11) = -4000*istep/nstep ;
        Res = (Fg-Fg_int);

        Res = Res-Kg(:,ibc)*bc;
        Res(ibc) = bc; 
        Kg(:,ibc)=0;
        Kg(ibc,:) = 0;
        Kg(ibc,ibc) = eye(length(ibc));
    
        if norm(Res)/norm(Fg) < 1e-4
            disp([istep iter]);
            break
        end
        
        [u ] = planar_simo_solve_update(Kg, Res, u);
    end
    plot(x0'+u(1:3:end),y0'+u(2:3:end),'k-o');  axis equal
    
    
    xlim([-100 300])
    title(['Step ', num2str(istep)]);
    xlabel('x');
    ylabel('y');
    
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if istep == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
        pause
    u_plot=[u_plot -u(3*11-1)];
    f_plot=[f_plot -Fg(3*11-1)];
end
toc