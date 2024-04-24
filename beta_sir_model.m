% Parameters
alpha = 0.2;  % Contact rate
I0 = 1e-6;  % Initial fraction of infected individuals
S0 = 1 - I0;  % Initial fraction of susceptible individuals
tspan = [0 200];  % Time span

% Recovery rate values to plot and setup colors
beta_values = [0.4, 0.2, 0.1, 0.05, 0.025];
colors = [
        [102, 194, 165] / 255;  % Bluish green
        [252, 141, 98] / 255;   % Vermilion
        [141, 160, 203] / 255;  % Blue Purple
        [231, 138, 195] / 255;  % Pale Red Purple
        [166, 216, 84] / 255    % Yellow Green
    ];
lineWidth = 2;
figure; % Create new figure
hold on; % Hold on for multiple plots

% Solve SIR model ODEs for different beta values
for i = 1:length(beta_values)
    beta = beta_values(i);
    color = colors(i, :);
    [t, result] = ode45(@(t, y) sir_model(t, y, alpha, beta), tspan, [S0; I0]);
    
    S = result(:, 1);
    I = result(:, 2);
    
    % Plot I(t) with solid line
    plot(t, I, 'Color', color, 'LineStyle', '-', 'LineWidth', lineWidth, 'DisplayName', sprintf('β = %.3f', beta));
    % Plot S(t) with dashed line
    plot(t, S, 'Color', color, 'LineStyle', '--', 'LineWidth', lineWidth);
end

% Adding labels and legend
xlabel('Time (days)');
ylabel('Proportion of Population');
legend show;
grid on;
title('Impact of β on i(t), s(t) in SIR model');

hold off; % Release hold

% Define the ODE system for SIR model at the end of the file
function dy = sir_model(t, y, alpha, beta)
    S = y(1);
    I = y(2);
    dS_dt = -alpha * S * I;
    dI_dt = alpha * S * I - beta * I;
    dy = [dS_dt; dI_dt];
end
