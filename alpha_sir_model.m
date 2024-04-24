function alpha_sir_model
    % Main function to run the SIR model simulation with different alpha values

    % Parameters
    beta = 1/10;  % Recovery rate
    I0 = 1e-6;  % Initial fraction of infected individuals
    S0 = 1 - I0;  % Initial fraction of susceptible individuals
    t_span = [0 200];  % Time span for simulation

    % Alpha values to plot
    alpha_values = [0.2, 0.25, 0.5, 1.0, 2.0];
    colors = [
        [102, 194, 165] / 255;  % Bluish green
        [252, 141, 98] / 255;   % Vermilion
        [141, 160, 203] / 255;  % Blue Purple
        [231, 138, 195] / 255;  % Pale Red Purple
        [166, 216, 84] / 255    % Yellow Green
    ];
    % Create figure
    figure;
    hold on;
    title('Impact of α on i(t), s(t) in SIR model');
    xlabel('Time (days)');
    ylabel('Proportion of Population');

    lineWidth = 2;

    % Loop over alpha values
    for i = 1:length(alpha_values)
        alpha = alpha_values(i);
        % Solve SIR model ODEs
        [T, Y] = ode45(@(t, y) sir_model(t, y, alpha, beta), t_span, [S0 I0]);
        S = Y(:, 1);
        I = Y(:, 2);
        
        % Plot I(t) with solid line
        plot(T, I, 'Color', colors(i, :), 'LineStyle', '-','LineWidth',lineWidth, 'DisplayName', ['α = ' num2str(alpha)]);
        % Plot S(t) with dashed line
        plot(T, S, 'Color', colors(i, :), 'LineStyle', '--','LineWidth',lineWidth);
    end

    legend('show');
    grid on;
    hold off;

    % Save the figure with high resolution
    saveas(gcf, 'sir_model_impact.png');
end

% Define the SIR model ODEs
function dy = sir_model(t, y, alpha, beta)
    S = y(1);
    I = y(2);
    dS_dt = -alpha * S * I;
    dI_dt = alpha * S * I - beta * I;
    dy = [dS_dt; dI_dt];
end
