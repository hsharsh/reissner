clc
clear
close all
u = [1e-16 1e-16]';
du = [0 0]';
ds = .1;
for j = 1:5
for ds = .05
for i = 1:10
%    res = tan(u(1)) - sqrt(2)*sin(u(1))+u(2)
   res = sin(u(1)) - u(2)
   p = u(1)^2+u(2)^2-ds^2;
%   p = u(2)^2 -ds^2;
%    K  = 1+tan(u(1))^2 - sqrt(2)*cos(u(1));
     K = cos(u(1));
%    du = inv(K)*-res;
%    u = u+du;
   
   ktan = [K -1; 2*u(1) 2*ds]
   f = [-res -p]';
   du = inv(ktan)*f;
   u = u+du ;
   
   hold on
   end
   plot(u(1), u(2),'o')
   x = 0:.01:pi;
   plot(x,sin(x),'r')
%    pause
    
    
end
end