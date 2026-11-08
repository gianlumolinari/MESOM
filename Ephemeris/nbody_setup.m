function pars = nbody_setup(et0, center, frame, bodies, LU, TU)
%NBODY_SETUP  Parameters for nbody_eom.
%
%   pars = nbody_setup(et0, center, frame, bodies, LU, TU)
%
%   et0    : reference epoch [s past J2000 TDB]
%   center : SPICE ID of the integration centre, e.g. '399'
%   frame  : inertial frame, e.g. 'ECLIPJ2000'
%   bodies : cell array of SPICE IDs of the attracting bodies (incl. centre)
%   LU, TU : length [km] and time [s] units
%
% Gravitational parameters and mean radii are read from the loaded kernels.

    pars.center = center;
    pars.frame  = frame;
    pars.bodies = bodies;
    pars.LU = LU;
    pars.TU = TU;
    pars.et0 = et0 / TU;                                 % epoch in time units
    pars.sfVec = [LU; LU; LU; LU/TU; LU/TU; LU/TU];      % state scale factors
    n = numel(bodies);
    pars.GM = zeros(1, n);
    pars.R  = zeros(1, n);
    for k = 1:n
        pars.GM(k) = cspice_bodvrd(bodies{k}, 'GM', 1) * TU^2 / LU^3;
        radii = cspice_bodvrd(bodies{k}, 'RADII', 3);
        pars.R(k) = mean(radii) / LU;
    end
end
