function C = Jacobi(X, mu)
    x = X(1); y = X(2); z = X(3);
    vx = X(4); vy = X(5); vz = X(6);

    r1 = sqrt((x + mu)^2 + y^2 + z^2);       
    r2 = sqrt((x - (1 - mu))^2 + y^2 + z^2); 

    C = (x.^2 + y.^2) + 2 * ((1 - mu) ./ r1 + mu ./ r2) - (vx.^2 + vy.^2 + vz.^2);
    
end