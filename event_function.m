% % This event function is implemented in the ode integration of the
% % CR3BP trajectory to determine the times at which the spacecraft is
% % inside the occultation cone caused by the Moon obstructing the Sun.
% % The conditions for the event to trigger are determined geomtrically
% % by computing the angles (beta1 and beta2) between the vector
% % connecting the two vertices at the edge of the occultation zone and
% % the position vector between the S/C and each of the two vertices. For
% % the spacecraft to be inside of the occultation zone, two conditions
% % are required to be respected simultaneously: delta1 = alpha1 - beta1 > 0
% % and delta2 = alpha2 - beta2 > 0.
% 
% % When the spacecraft enters the occultation zone, delta1 and delta2
% % cross the 0 in positive direction, i.e. triggering the event. To
% % ensure that both are positive at the same time, the condition is
% % triggered only when the smallest of the two values is positive, since
% % when that is the case, the other is necessarily going to be positive
% % as well.
% %Similarly, when the spacecraft exits the occultation zone, the values
% %cross the zero in negative direction, i.e. from positive to negative.
% %When the smallest of the two values crosses the zero and becomes
% %negative, the condition for occultation is no longer satisfied,
% %therefore triggering the event and indicating when the spacecraft
% %exits the occultation zone.
% 
% 
function [value, isterminal, direction] = event_function(t, X, mu, aM, aE, Rs, Rs2, Rm, nE, nM, thetaM0)
    % t is already in TU (normalized time units)
    % No conversion needed since all angular velocities are in rad/TU

    pos = X(1:3); % Spacecraft position

    % Compute angles (already in correct units)
    thetaE = nE * t; % Earth angle
    thetaM = nM * t + thetaM0; % Moon angle

    % Direction cosine matrix for Moon position
    DCM = [cos(thetaE - thetaM), sin(thetaE - thetaM), 0;
           -sin(thetaE - thetaM), cos(thetaE - thetaM), 0;
           0, 0, 1];

    % Moon position in Sun-Earth frame
    moon_SE = [1 - mu; 0; 0] + DCM * [aM / aE; 0; 0];
    dSM = norm(moon_SE);
    X_hat = moon_SE / dSM;

    % Occultation geometry
    l1 = (Rm * dSM) / (Rs2 - Rm);
    l2 = (Rm * dSM) / (Rs - Rm);

    alpha1 = atan(Rs / (dSM + l1));
    alpha2 = atan(Rs2 / (dSM + l2));

    % Vertex positions
    V1 = moon_SE + l2 * X_hat - [mu; 0; 0];
    V2 = moon_SE + l1 * X_hat - [mu; 0; 0];

    % Angle calculations
    Vs1 = V1 - pos;
    V21 = V1 - V2;
    Vs2 = V2 - pos;
    V12 = V2 - V1;

    dot1 = max(min(dot(Vs1, V21) / (norm(Vs1) * norm(V21)), 1.0), -1.0);
    dot2 = max(min(dot(Vs2, V12) / (norm(Vs2) * norm(V12)), 1.0), -1.0);

    beta2 = acos(dot1);
    beta1 = acos(dot2);

    delta1 = alpha1 - beta1;
    delta2 = alpha2 - beta2;

    % Combined events
    min_delta = min(delta1, delta2);
    value = [min_delta; min_delta];
    isterminal = [0; 0]; % Don't stop integration
    direction = [1; -1]; % [entry: positive-going, exit: negative-going]
    
end




