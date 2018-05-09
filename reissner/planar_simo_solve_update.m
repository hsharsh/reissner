function u = planar_simo_solve_update(Kg, Res, u)
    delu = (Kg)\ Res;
    u = u + delu;
end