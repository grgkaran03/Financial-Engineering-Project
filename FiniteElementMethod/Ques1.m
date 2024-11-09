function [V, S, t] = Ques1(x_min, x_max, N, M, flag, q_idx)
    quadrature_rules = { @(f, x_min, x_max, h, N, i, j, a, b) trapezoidal(f, x_min, x_max, h, N, i, j, a, b),
                         @(f, x_min, x_max, h, N, i, j, a, b) simpsons(f, x_min, x_max, h, N, i, j, a, b) };
                     
    basis = @(x_min, x_max, h, N, i, x) linear_basis(x_min, x_max, h, N, i, x);
    
    der_b = @(x_min, x_max, h, N, i, x, flag) derivative_linear_basis(x_min, x_max, h, N, i, x, flag);
    
    
    [T, K] = deal(1, 40);
    [r, sigma, delta] = deal(0.05, 0.2, 0);
    
    q = 2 * r / (sigma ^ 2);
    q_delta = 2 * (r - delta) / (sigma ^ 2);

    [x, tau, h, k] = create_grid(x_min, x_max, N, M, T, sigma);
    x_eff = x(2 : length(x)-1); 
    
    quadrature = quadrature_rules{q_idx};
    
    % Construct stiffness matrix A and mass Matrix B
    [A, B] = deal(zeros(length(x_eff)));
    for i = 1 : length(x_eff)
        for j = 1 : length(x_eff)
            if i == j
                for idx = 0 : 1
                    func1 = @(x_min, x_max, h, N, i, j, val, fl) der_b(x_min, x_max, h, N, i, val, fl) * der_b(x_min, x_max, h, N, j, val, fl);
                    A(i, j) = A(i, j) + quadrature(func1, x_min, x_max, h, length(x), i, j, x(i+idx), x(i+idx+1));

                    func3 = @(x_min, x_max, h, N, i, j, val, fl) basis(x_min, x_max, h, N, i, val) * basis(x_min, x_max, h, N, j, val);
                    B(i, j) = B(i, j) + quadrature(func3, x_min, x_max, h, length(x), i, j, x(i+idx), x(i+idx+1));
                end
            elseif abs(i - j) == 1
                idx = min(i, j);
                func1 = @(x_min, x_max, h, N, i, j, val, fl) der_b(x_min, x_max, h, N, i, val, fl) * der_b(x_min, x_max, h, N, j, val, fl);
                A(i, j) = A(i, j) + quadrature(func1, x_min, x_max, h, length(x), i, j, x(idx+1), x(idx+2));

                func3 = @(x_min, x_max, h, N, i, j, val, fl) basis(x_min, x_max, h, N, i, val) * basis(x_min, x_max, h, N, j, val);
                B(i, j) = B(i, j) + quadrature(func3, x_min, x_max, h, length(x), i, j, x(idx+1), x(idx+2));
            end
        end
    end
    
    B_final = B + (0.5 * k) .* A;
    A_final = B - (0.5 * k) .* A;

    [w, b] = deal(zeros(length(x_eff), length(tau)));

    for i = 1 : length(tau)
        b(:, i) = ((x_eff - x_min)/(x_max - x_min)) * h * (0.25 * (q_delta + 1) ^ 2) * beta(x_max, tau(i), q_delta);
    end

    % Construct weights w
    w(: , 1) = IC(x_eff, q_delta) - phi_b(x_eff, 0, q_delta, x_min, x_max);
    for i = 2 : length(tau)
        w(:, i) = B_final \ (A_final * w(:, i - 1) - (k/2) * (b(:, i) + b(:, i-1)));
    end

    % Get option price at nodes
    base = zeros(N + 1, M + 1);
    for i = 1 : length(x)
        for j = 1 : length(tau)
            base(i, j) = phi_b(x(i), tau(j), q_delta, x_min, x_max);
        end
    end

    temp = zeros(1, M + 1);
    nodes = [temp; w; temp];
    nodes = nodes + base;

    V = nodes;
    temp1 = (-0.5) * (q_delta - 1) * x;
    temp2 = -(0.25 * (q_delta - 1)^2 + q) * tau;

    for i = 1 : length(x)
        for j = 1 : length(tau)
            V(i, j) = (K * exp(temp1(i) + temp2(j))) * nodes(i, j);
        end
    end

    S = K * exp(x);
    t = T - tau * 2 / sigma^2;

    % Plot graphs
    if(flag)
        plots(S, t, V, q_idx, "European Call option using");
    end

end


% ***************  HELPER FUNCTIONS  ***************

function y = IC(x, q_delta)
    y = max(0, exp(0.5 * x * (q_delta + 1)) - exp(0.5 * x * (q_delta - 1)));
end

% Boundary condition at x_min
function y = alpha(x_min, tau, q_delta)
    y = zeros(size(tau));
end

% Boundary condition at x_max
function y = beta(x_max, tau, q_delta)
    y = exp(0.5 * (q_delta + 1) * x_max + 0.25 * tau * (q_delta + 1) ^ 2);
end

function val = phi_b(x, tau, q_delta, x_min, x_max)
    val = beta(x_max, tau, q_delta) - alpha(x_min, tau, q_delta);
    val = alpha(x_min, tau, q_delta) + val * (x - x_min) / (x_max - x_min);
end