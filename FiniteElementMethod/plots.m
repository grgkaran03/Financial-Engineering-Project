function plots(S, time, V, idx, plot_title)
    if idx == 1
        type = " Trapezoidal Rule";
    else
        type = " Simpson's Rule";
    end
    
    figure;
    surf(S', time', V');
    xlabel("S");
    ylabel("t");
    zlabel("V(S, t)");
    title(plot_title + type);

    figure;
    plot(S', V(:, end), S', V(:, 1), 'LineWidth',1);
    title(strcat(plot_title, type, ' at final time level'));
    legend('t = 0', 't = T');
    xlabel('S');
    ylabel('V(S,0) and V(S,T)');
    hold off;
end

