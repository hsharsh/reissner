function [Kg Fg_int Fg]=planar_simo(x, y, icon, u, EA, EI, GA, nnod, nelm)

x_gp = [-1 1]/sqrt(3);
w_gp = [1 1];
Fg_int = u *0; 
Kg = zeros(3*nnod);
Fg = zeros(3*nnod,1);
theta = u(3:3:end);

total_length = 0;
for lmn = 1:nelm
    node1 = icon(lmn,1); node2 = icon(lmn,2);
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
    node1 = icon(lmn,1); node2 = icon(lmn,2);
    iv = [3*(node1-1)+1 3*(node1-1)+2 3*(node1-1)+3 ...
          3*(node2-1)+1 3*(node2-1)+2 3*(node2-1)+3 ];
         
    Uloc = u(iv);
    Kl = zeros(6);Kgeo = zeros(6);
    for igp = 1:length(w_gp)
        Dmat = [EA 0 0;0 0 0; 0 0 EI];
        [Bmat B0 ENN Rmat Tmat len G_N G_Q Uloc_junk]= Bmats(Uloc, x, y, lmn, theta, x_gp(igp));
        Kl = Kl + Bmat'*Dmat*Bmat *len/2*w_gp(igp);
        Kgeo = Kgeo + (stress(lmn).comp(1) * G_N) *len/2*w_gp(igp);
        Kgeo = Kgeo + (stress(lmn).comp(2) * G_Q) *len/2*w_gp(igp);
    end
    Dmat = [0 0 0;0 GA 0; 0 0 0];
    [Bmat B0 ENN Rmat Tmat len G_N G_Q Uloc_junk] = Bmats(Uloc, x, y, lmn, theta, 0);
    Kl = Kl + Bmat'*Dmat*Bmat *len/2*2;
    Kl = Kl + Kgeo;
    Kg(iv,iv) = Kg(iv,iv) +    Tmat'*Kl*Tmat ;
    
end