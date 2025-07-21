clear all; close all; clc; format longG;
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontsize',14);
tmoon = 0.507987575964444;
%Define Mass ratio parameter for the Sun-Earth system
mu           = 3.0542E-6; 
pars.mu = mu;

% Single-shooting method
Iter         = 100;      % Max. no. of Newton's Method iterations
Tol          = 1e-10;   % Newton's method tolerance


% Initial guess - Build initial guess from 

X0_tilde = [1.00105871003043
-1.11829405105640e-10
-6.86288108338043e-67
5.52930876445842e-09
0.0525941829622717
5.46150883100796e-64];

T_tilde  = 0.126996893991111;



% Predictor / Corrector
Tdesired = 0.126996893991111; 

    % 1. Predictor
    X0 = X0_tilde;
    T  = T_tilde ;


    % 2. Corrector
    for jj = 1:Iter


        Phi0    = eye(6);
        Xaug0   = [X0; reshape(Phi0, 36, 1)];


        opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
        [t, Xaug] = ode113(@(t,Xaug) STM_dot(t, Xaug, pars), [0 T], Xaug0, opts);

        %2a. Collect final point, vector field at final point, STM
        XT           = Xaug(end, 1:6).';
        STM          = reshape(Xaug(end, 7:end), 6, 6);
        f_XT         = f(XT,pars); 


        %2b. Phase Condition
        phs = X0(2);       % y = 0 for Lyapunov Planar Orbits


        %2c. Target period
        prd = T-Tdesired ;


        %2d. Build error vector and error jacobian
        F = [XT - X0; X0(2); prd];

        DF = [STM - Phi0,f_XT;         
              f(X0_tilde,pars)',0;           
              0,0,0,0,0,0,1];  % 8 x 7

        dZ = -DF\F;
        fprintf('|F| = %.6e, |dZ| = %.6e\n',norm(F),norm(dZ));

        if(norm(F) < Tol)
            fprintf('Periodic Orbit has been found!\n\n');

            %2e. Store results
            Xpo = X0;
            Tpo = T;

            %2f. Jacobi integral
            Cpo  = Jacobi(X0, pars.mu);

            % 2i. Break and Repeat!
            break;
        else

            % Newton's update
            X0 = X0 + dZ(1:6);
            T  = T + dZ(end);
        end
    end








