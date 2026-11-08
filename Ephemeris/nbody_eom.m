function dx = nbody_eom(t, x, pars)
%NBODY_EOM  Point-mass N-body equations of motion relative to a central body.
%
%   r'' = -GM_c r/|r|^3 - sum_j GM_j [ (r - r_j)/|r - r_j|^3 + r_j/|r_j|^3 ]
%
% r is the spacecraft position relative to the centre and r_j the position of
% third body j relative to the centre (from SPICE). The second term in the
% bracket is the indirect acceleration of the non-inertial centre.
% Nondimensional units: t in TU since pars.et0, x scaled by pars.sfVec.

    r  = x(1:3);
    et = (pars.et0 + t) * pars.TU;
    a  = zeros(3, 1);
    for k = 1:numel(pars.bodies)
        if strcmp(pars.bodies{k}, pars.center)
            a = a - pars.GM(k) * r / norm(r)^3;
        else
            rj = cspice_spkpos(pars.bodies{k}, et, pars.frame, 'NONE', pars.center) / pars.LU;
            d  = r - rj;
            a  = a - pars.GM(k) * (d / norm(d)^3 + rj / norm(rj)^3);
        end
    end
    dx = [x(4:6); a];
end
