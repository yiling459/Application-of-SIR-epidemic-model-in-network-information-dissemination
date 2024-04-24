% Set model parameters
number = 1e5;
alpha = 0.2;
sigma = 2.5;
beta = alpha / sigma;
tEnd = 200;
t = 0:tEnd;

% Initial infected proportions list
i0_values = [0.0001, 0.001, 0.01, 0.1, 0.25, 0.5];
colors = [
        [102, 194, 165] / 255;  % Bluish green
        [252, 141, 98] / 255;   % Vermilion
        [141, 160, 203] / 255;  % Blue Purple
        [231, 138, 195] / 255;  % Pale Red Purple
        [166, 216, 84] / 255;   % Yellow Green
        [255, 165, 0] / 255     % Orange
    ];
lineWidth = 2;

% Setup plot
figure;
hold on;
title('Impact of i0 on i(t) and s(t) in SIR model');

for i = 1:length(i0_values)
    i0 = i0_values(i);
    s0 = 1 - i0;
    Y0 = [i0; s0];
    [t, ySIR] = ode45(@(t, y) dySIR(t, y, alpha, beta), t, Y0);

    plot(t, ySIR(:, 1), 'Color', colors(i, :), 'LineWidth', lineWidth,'DisplayName', sprintf('i0=%.3f', i0));
    plot(t, ySIR(:, 2), '--', 'Color', colors(i, :), 'LineWidth', lineWidth);
end

xlabel('Time (days)');
ylabel('Proportion of Population');
legend show;
grid on;
hold off;

% Define the ODE system for SIR model at the end of the file
function dy = dySIR(t, y, alpha, beta)
    i = y(1);
    s = y(2);
    di_dt = alpha * s * i - beta * i;
    ds_dt = -alpha * s * i;
    dy = [di_dt; ds_dt];
end
