% Main executable block at the top

% Set model parameters
number = 1e5; % Total population
lamda = 0.2; % Daily contact rate, average number of susceptible individuals effectively contacted per day by an infectious person
sigma = 2.5; % Number of contacts during infectious period
mu = lamda / sigma; % Daily recovery rate, proportion of infected individuals who recover each day
fsig = 1 - 1 / sigma;
fprintf('lamda=%.2f\tmu=%.2f\tsigma=%.2f\t(1-1/sig)=%.2f\n', lamda, mu, sigma, fsig);

% Colors for the plot
colors = [
        [102, 194, 165] / 255;  % Bluish green
        [252, 141, 98] / 255;   % Vermilion
        [141, 160, 203] / 255;  % Blue Purple
        [231, 138, 195] / 255;  % Pale Red Purple
        [166, 216, 84] / 255;   % Yellow Green
        [255, 165, 0] / 255;    % Orange
        [0, 128, 128] / 255;    % Teal
        [255, 192, 203] / 255;  % Pink
        [0, 0, 128] / 255       % Dark Blue
    ];
lineWidth = 2; % Line width for plotting

% Numerical solution using ODE solver
tEnd = 200; % Prediction period length
t = 0:tEnd-1; % (start,stop,step)
s0List = linspace(0.01, 0.9, 9); % Creates a linearly spaced vector with 9 elements
figure;
hold on;
for i = 1:length(s0List) % s0, initial proportion of susceptible individuals
    s0 = s0List(i);
    i0 = 1 - s0; % i0, initial proportion of infected individuals
    Y0 = [i0; s0]; % Initial conditions for the ODE system
    [t, ySIR] = ode45(@(t, y) dySIR(t, y, lamda, mu), t, Y0);
    plot(ySIR(:,2), ySIR(:,1), 'Color', colors(i, :), 'LineWidth', lineWidth);
end

% Plot settings
title('Phase trajectory of SIR models');
axis([0 1 0 1]);
plot([0 1], [1 0], 'b-', 'LineWidth', lineWidth);
plot([1/sigma, 1/sigma], [0, 1-1/sigma], 'b--', 'LineWidth', lineWidth);
xlabel('s(t) - susceptible proportion');
ylabel('i(t) - infected proportion');
text(0.8, 0.9, sprintf('$1/\\sigma$ = %.2f', 1/sigma), 'Color', 'blue');
hold off;

% Function definitions at the end of the file

% Define the SIR model differential equations
function dy = dySIR(t, y, lamda, mu)
    i = y(1);
    s = y(2);
    di_dt = lamda * s * i - mu * i; % di/dt = lamda * s * i - mu * i
    ds_dt = -lamda * s * i; % ds/dt = -lamda * s * i
    dy = [di_dt; ds_dt];
end
