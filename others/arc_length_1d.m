clc
clear
close all
u = [0 0]';
du = [0 0]';
alpha = .1;
for step =1:15
    step
    ds = 0.5;
    du = [0 0]';
    uold = u;
for i = 1:200
%    res = tan(u(1)) - sqrt(2)*sin(u(1))+u(2)
   res = sin(u(1)) - u(2);
   p = (u(1)-uold(1))^2+alpha*(u(2)-uold(2))^2-ds^2;
%   p = u(2)^2 -ds^2;
%    K  = 1+tan(u(1))^2 - sqrt(2)*cos(u(1));
     K = cos(u(1));
%    du = inv(K)*-res;
%    u = u+du;
   alpha = sign(K);
   ktan = [K -1; 2*(u(1)-uold(1)) 2*alpha*ds];
   f = [-res -p]';
   du = inv(ktan)*f;
   u = u+du ;
   
   hold on
end

u-uold;

du;
sqrt((u(1)-uold(1))^2+alpha*(u(2)-uold(2))^2)
   plot(u(1), u(2),'o')
   x = 0:.01:5*pi;
   plot(x,sin(x),'r')
%    pause
    
    
end
