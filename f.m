function dXdt = f(X, pars)
% This integrates the EoMs of the CR3BP returning Xdot
    mu = pars.mu;
    x = X(1); y = X(2); z = X(3);
    dx = X(4); dy = X(5); dz = X(6);

    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

    ddx = 2*dy + x - (1 - mu)*(x + mu)/r1^3 - mu*(x - 1 + mu)/r2^3;
    ddy = -2*dx + y - (1 - mu)*y/r1^3 - mu*y/r2^3;
    ddz = -(1 - mu)*z/r1^3 - mu*z/r2^3;

    dXdt = [dx; dy; dz; ddx; ddy; ddz];
end