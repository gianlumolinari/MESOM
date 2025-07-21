function dXdt = cr3bp_eom(~, X, mu)

    x = X(1); y = X(2); z = X(3);
    vx = X(4); vy = X(5); vz = X(6);

    r1 = sqrt((x + mu)^2 + y^2 + z^2);       
    r2 = sqrt((x - (1 - mu))^2 + y^2 + z^2); 

    
    ax = 2*vy + x - (1 - mu)*(x + mu)/r1^3 - mu*(x - (1 - mu))/r2^3;
    ay = -2*vx + y - (1 - mu)*y/r1^3 - mu*y/r2^3;
    az = -(1 - mu)*z/r1^3 - mu*z/r2^3;

    dXdt = [vx; vy; vz; ax; ay; az];
end