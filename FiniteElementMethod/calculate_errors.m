function calculate_errors(x_min, x_max, N_array, M_array,q_idx)
    dx_array = (x_max - x_min)./N_array;

    En = [];
    E2n = [];
    for p = 1:length(N_array)
        N1 = N_array(p);
        N2 = 2 * (N1 + 1);
        N3 = 2 * (N2 + 1);
        M = M_array(p);

        [V1, ~, ~] = Ques1(x_min, x_max, N1, M, 0, q_idx);
        [V2, ~, ~] = Ques1(x_min, x_max, N2, M, 0, q_idx);
        [V3, ~, ~] = Ques1(x_min, x_max, N3, M, 0, q_idx);

        Err1 = zeros(size(V1));
        Err2 = zeros(size(V2));
        for i = 1:N1
            for j = 1:M
                Err1(i,j) = abs(V1(i,j) - V2(2*i, j));
            end
        end

        for i = 1:N2
            for j = 1:M
                Err2(i,j) = abs(V2(i,j) - V3(2*i, j));
            end
        end

        En = [En, norm(Err1,2)];
        E2n = [E2n, norm(Err2,2)];
    end

    % Plot the errors in a log-log scale
    figure;
    loglog(dx_array, log2(En./E2n), 'LineWidth', 1);
    xlabel('Δx');
    ylabel('L2-Norm Error');
    if(q_idx == 1)
        title('Convergence Analysis For Trapezoidal quadrature.');
    else
        title('Convergence Analysis For Simpsons quadrature.');
    end
    hold off;
end
