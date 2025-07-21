function dXaug = STM_dot(t, Xaug, pars)
    % Extract state and STM
    X   = Xaug(1:6);                  % State vector
    Phi = reshape(Xaug(7:end), 6, 6); % STM (reshaped from 36x1 to 6x6)

    % Compute state derivative
    dX  = f(X, pars);              % CR3BP dynamics

    % Compute Jacobian at X
    A   = dfdx(X, pars);              % 6x6 Jacobian

    % STM derivative: dPhi/dt = A * Phi
    dPhi = A * Phi;

    % Combine into single column vector
    dXaug = [dX; reshape(dPhi, 36, 1)];
end