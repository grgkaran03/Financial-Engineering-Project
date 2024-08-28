function [Delta, Gamma, Vega, Theta, Rho] = calculate_greeks(U, S, Time, sig, r, k, h)
    [m, n] = size(U);
    
    % Delta
    Delta = zeros(m, n);
    Delta(2:end-1, :) = (U(3:end, :) - U(1:end-2, :)) / (2*h);
    Delta(1, :) = (U(2, :) - U(1, :)) / h;
    Delta(end, :) = (U(end, :) - U(end-1, :)) / h;
    
    % Gamma
    Gamma = zeros(m, n);
    Gamma(2:end-1, :) = (U(3:end, :) - 2*U(2:end-1, :) + U(1:end-2, :)) / (h^2);
    
    % Vega (approximation using finite difference)
    dSig = 0.01;
    U_up = FTCS(Time(end), K, r, sig+dSig, delta, S, Time, h, k, isTerminal, epsilon);
    U_down = FTCS(Time(end), K, r, sig-dSig, delta, S, Time, h, k, isTerminal, epsilon);
    Vega = (U_up - U_down) / (2*dSig);
    
    % Theta
    Theta = zeros(m, n);
    Theta(:, 2:end-1) = -(U(:, 3:end) - U(:, 1:end-2)) / (2*k);
    Theta(:, 1) = -(U(:, 2) - U(:, 1)) / k;
    Theta(:, end) = -(U(:, end) - U(:, end-1)) / k;
    
    % Rho (approximation using finite difference)
    dR = 0.0001;
    U_up = FTCS(Time(end), K, r+dR, sig, delta, S, Time, h, k, isTerminal, epsilon);
    U_down = FTCS(Time(end), K, r-dR, sig, delta, S, Time, h, k, isTerminal, epsilon);
    Rho = (U_up - U_down) / (2*dR);
end