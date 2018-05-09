clc
clear
close all
u = [0 0]';
du = [0 0]';
alpha = 1;
for step =1:25
    step
    ds = .7;
    du = [0 0]';
    uold = u;
for i = 1:20
   res = sin(u(1)) - u(2);
   p = (u(1)-uold(1))^2+alpha*(u(2)-uold(2))^2-ds^2;
   K = cos(u(1));
   if i ==1 ;  
        du(1) = 1/abs(K);
        du(2) = ds/abs(du(1));
   else
       ktan = [K -1; 2*(u(1)-uold(1)) 2*alpha*(u(2)-uold(2))];
       f = [-res -p]';
       du = inv(ktan)*f;
   end
   u = u+du ;
   
   hold on
end

   plot(u(1), u(2),'o')
   x = 0:.01:5*pi;
   plot(x,sin(x),'r')
%    pause
    
    
end
