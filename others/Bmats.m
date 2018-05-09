function [Bmat B0 ENN Rmat Tmat len G_N G_Q Uloc] = Bmats(Uloc, x, y, lmn, theta,xi)
    len = sqrt((x(lmn+1) - x(lmn))^2+ (y(lmn+1) - y(lmn))^2);
    rotation = (theta(lmn+1) + theta(lmn))/2;
    Rmat = [...
        cos(rotation), sin(rotation), 0; ...
        -sin(rotation), cos(rotation), 0; ...
            0      ,     0 ,    1         ];
    angle = atan2(y(lmn+1) - y(lmn), x(lmn+1) - x(lmn));
    Tmat = [ cos(angle)   sin(angle)     0    0 0      0;...
           -sin(angle)   cos(angle)     0    0 0      0;...
                 0 0 1 0 0 0; ...
              0 0         0       cos(angle)   sin(angle) 0;
              0  0        0      -sin(angle)   cos(angle) 0;
              0 0 0 0 0 1];
    Uloc = Tmat*Uloc  ;
        
    B0 = 1/len*[-1  0  0 1  0  0; ...
                 0 -1  0 0  1  0;...
                 0  0 -1 0  0  1];
    ENN = [1 0 0]';

    
    
    u_e_prime = (Uloc(4) - Uloc(1))/len;
    w_e_prime = (Uloc(5) - Uloc(2))/len;
    DN1 = -1/len; N1 = (1-xi)/2;
    DN2 =  1/len; N2 = (1+xi)/2;
    alpha = -(1+u_e_prime)*sin(rotation) + w_e_prime * cos(rotation);
    beta = -(1+u_e_prime)*cos(rotation) - w_e_prime * sin(rotation); 

      Bmat = [DN1*cos(rotation)  DN1*sin(rotation) alpha*N1   DN2*cos(rotation)  DN2*sin(rotation) alpha*N2; ...
           -DN1*sin(rotation)  DN1*cos(rotation) beta*N1   -DN2*sin(rotation)  DN2*cos(rotation)  beta*N2; ...
                 0            0             DN1            0            0             DN2      ...
           ];
       alpha3 = -(1+u_e_prime)*cos(rotation) - w_e_prime * sin(rotation);
       G_11= GN_IJ(rotation, N1, N1, DN1, DN1, alpha3);
       G_12= GN_IJ(rotation, N1, N2, DN1, DN2, alpha3);
       G_21= GN_IJ(rotation, N2, N1, DN2, DN1, alpha3);
       G_22= GN_IJ(rotation, N2, N2, DN2, DN2, alpha3);
       G_N = [G_11 G_12; G_21 G_22];
       
       alpha4 = (1+u_e_prime)*sin(rotation) - w_e_prime * cos(rotation);
       G_11= GQ_IJ(rotation, N1, N1, DN1, DN1, alpha4);
       G_12= GQ_IJ(rotation, N1, N2, DN1, DN2, alpha4);
       G_21= GQ_IJ(rotation, N2, N1, DN2, DN1, alpha4);
       G_22= GQ_IJ(rotation, N2, N2, DN2, DN2, alpha4);
       G_Q = [G_11 G_12; G_21 G_22];
       
  
end

function G_IJ= GN_IJ(rotation, NI, NJ, DNI, DNJ, alpha3)
   G_IJ = [         0                  0            -DNI*NJ*sin(rotation)  ; ...
                   0                  0              DNI*NJ*cos(rotation)  ; ...
           -NI*DNJ*sin(rotation)   NI*DNJ*cos(rotation)     alpha3*NI*NJ  ;  ...
             ]    ;
end

function G_IJ= GQ_IJ(rotation, NI, NJ, DNI, DNJ, alpha4)
   G_IJ = [         0                  0            -DNI*NJ*cos(rotation)  ; ...
                   0                  0            -DNI*NJ*sin(rotation)  ; ...
           -NI*DNJ*cos(rotation)    -NI*DNJ*sin(rotation)     alpha4*NI*NJ  ;  ...
             ]    ;
end