isTerminal = true;

T = 1;
K = 20;
r = 0.07;
sig = 0.2;
delta = 0.01;

% Boundary
S_min = 0;
S_max = 30;

h = 1;
k = 1e-3;
m = (S_max - S_min)/h;
n = ceil(T/k);

S = S_min:h:S_max;
Time = 0:k:T;

% Tau = T - t
% U = FTCS(T, K, r, sig, delta, S, Time, h, k, isTerminal,epsilon);
% U = Crank(T, K, r, sig, delta, S, Time, h, k, isTerminal,epsilon);
% U = BTCS(T, K, r, sig, delta, S, Time, h, k, isTerminal,epsilon);

% figure; plot(S, U(1, :)); hold on; plot(S, U(end, :)); hold off;

epsilons = [1e-2, 1e-4, 1e-6];

for i = 1:length(epsilons)
    figure;
    epsilon = epsilons(i);

    U_FTCS = FTCS(T, K, r, sig, delta, S, Time, h, k, isTerminal, epsilon);
    U_Crank = Crank(T, K, r, sig, delta, S, Time, h, k, isTerminal, epsilon);
    U_BTCS = BTCS(T, K, r, sig, delta, S, Time, h, k, isTerminal, epsilon);
    
    subplot(1, 3, 1);
    surf(Time, S, U_FTCS');
    title(['FTCS, \epsilon = ', num2str(epsilon)]);
    xlabel('Time');
    ylabel('Stock Price');
    zlabel('Option Value (U)');
    shading interp;
    
    % Crank-Nicolson subplot
    subplot(1, 3, 2);
    surf(Time, S, U_Crank');
    title(['Crank-Nicolson, \epsilon = ', num2str(epsilon)]);
    xlabel('Time');
    ylabel('Stock Price');
    zlabel('Option Value (U)');
    shading interp;
    
    % BTCS subplot
    subplot(1, 3, 3);
    surf(Time, S, U_BTCS');
    title(['BTCS, \epsilon = ', num2str(epsilon)]);
    xlabel('Time');
    ylabel('Stock Price');
    zlabel('Option Value (U)');
    shading interp;
end

epsilon = 1e-2;

plot_times = [0, 0.25, 0.5, 0.75, 1];

U_FTCS = FTCS(T, K, r, sig, delta, S, Time, h, k, isTerminal, epsilon);
U_Crank = Crank(T, K, r, sig, delta, S, Time, h, k, isTerminal, epsilon);
U_BTCS = BTCS(T, K, r, sig, delta, S, Time, h, k, isTerminal, epsilon);

figure;
subplot(3, 1, 1);
hold on;
for t = plot_times
    time_index = find(abs(Time - t) < 1e-10, 1);
    plot(S, U_FTCS(time_index, :), 'DisplayName', sprintf('t = %.2f', t));
end
hold off;
title(['FTCS, \epsilon = ', num2str(epsilon)]);
xlabel('Stock Price (S)');
ylabel('Option Value (U)');
legend('Location', 'best');

% Plot for Crank-Nicolson
subplot(3, 1, 2);
hold on;
for t = plot_times
    time_index = find(abs(Time - t) < 1e-10, 1);
    plot(S, U_Crank(time_index, :), 'DisplayName', sprintf('t = %.2f', t));
end
hold off;
title(['Crank-Nicolson, \epsilon = ', num2str(epsilon)]);
xlabel('Stock Price (S)');
ylabel('Option Value (U)');
legend('Location', 'best');

% Plot for BTCS
subplot(3, 1, 3);
hold on;
for t = plot_times
    time_index = find(abs(Time - t) < 1e-10, 1);
    plot(S, U_BTCS(time_index, :), 'DisplayName', sprintf('t = %.2f', t));
end
hold off;
title(['BTCS, \epsilon = ', num2str(epsilon)]);
xlabel('Stock Price (S)');
ylabel('Option Value (U)');
legend('Location', 'best');

sgtitle(['Numerical Solutions at Different Times, \epsilon = ', num2str(epsilon)]);


[Delta, Gamma, Vega, Theta, Rho] = calculate_greeks(U_Crank, S, Time, sig, r, k, h, K, delta, isTerminal, epsilon);
% Plot Delta

figure;
subplot(2, 3, 2);
hold on;
for t = plot_times
    time_index = find(abs(Time - t) < 1e-10, 1);
    plot(S, Delta(time_index, :), 'DisplayName', sprintf('t = %.2f', t));
end
hold off;
title('Delta');
xlabel('Stock Price (S)');
ylabel('Delta');
legend('Location', 'best');

% Plot Gamma
subplot(2, 3, 3);
hold on;
for t = plot_times
    time_index = find(abs(Time - t) < 1e-10, 1);
    plot(S, Gamma(time_index ,:), 'DisplayName', sprintf('t = %.2f', t));
end
hold off;
title('Gamma');
xlabel('Stock Price (S)');
ylabel('Gamma');
legend('Location', 'best');

% Plot Vega
subplot(2, 3, 4);
hold on;
for t = plot_times
    time_index = find(abs(Time - t) < 1e-10, 1);
    plot(S, Vega(time_index, :), 'DisplayName', sprintf('t = %.2f', t));
end
hold off;
title('Vega');
xlabel('Stock Price (S)');
ylabel('Vega');
legend('Location', 'best');

% Plot Theta
subplot(2, 3, 5);
hold on;
for t = plot_times
    time_index = find(abs(Time - t) < 1e-10, 1);
    plot(S, Theta(time_index, :), 'DisplayName', sprintf('t = %.2f', t));
end
hold off;
title('Theta');
xlabel('Stock Price (S)');
ylabel('Theta');
legend('Location', 'best');

% Plot Rho
subplot(2, 3, 6);
hold on;
for t = plot_times
    time_index = find(abs(Time - t) < 1e-10, 1);
    plot(S, Rho(time_index, :), 'DisplayName', sprintf('t = %.2f', t));
end
hold off;
title('Rho');
xlabel('Stock Price (S)');
ylabel('Rho');
legend('Location', 'best');

sgtitle(['Option Value and Greeks, \epsilon = ', num2str(epsilon)]);


function [Delta, Gamma, Vega, Theta, Rho] = calculate_greeks(U, S, Time, sig, r, k, h, K, delta, isTerminal, epsilon)
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