function [value, isterminal, direction] = nbody_impacts(t, x, pars)
%NBODY_IMPACTS  Event: spacecraft reaches the surface of any attracting body.

    et = (pars.et0 + t) * pars.TU;
    n = numel(pars.bodies);
    value = zeros(n, 1);
    for k = 1:n
        if strcmp(pars.bodies{k}, pars.center)
            rk = zeros(3, 1);
        else
            rk = cspice_spkpos(pars.bodies{k}, et, pars.frame, 'NONE', pars.center) / pars.LU;
        end
        value(k) = norm(x(1:3) - rk) - pars.R(k);
    end
    isterminal = ones(n, 1);
    direction  = -ones(n, 1);
end
