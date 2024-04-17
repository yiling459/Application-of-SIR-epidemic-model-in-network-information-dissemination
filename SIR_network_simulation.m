%% SIR Model
% Initial parameters
S0 = 980; % Initial number of SUSCEPTIBLE
I0 = 20;  % Initial number of INFECTED
R0 = 0;   % Initial number of REMOVED (including both cured and deceased)
alpha = 0.001; % Infection rate coefficient
beta = 0.01;   % Recovery rate coefficient

% Total population
N = S0 + I0 + R0;

% Simulation days
days = 365;

% Initialize the S, I, R arrays
S = zeros(1, days);
I = zeros(1, days);
R = zeros(1, days);

% Set initial conditions
S(1) = S0;
I(1) = I0;
R(1) = R0;

% Simulate the propagation process every day
for t = 2:days
    S(t) = S(t-1) - alpha * S(t-1) * I(t-1);
    I(t) = I(t-1) + alpha * S(t-1) * I(t-1) - beta * I(t-1);
    R(t) = R(t-1) + beta * I(t-1);
    
    % Ensure the values of S, I, R remain non-negative
    S(t) = max(S(t), 0);
    I(t) = max(I(t), 0);
    R(t) = max(R(t), 0);
    
    % Optional: If using a fixed population model (S+I+R=N), re-normalize
    % If not using this approach, comment out the next three lines
    %total = S(t) + I(t) + R(t);
    %S(t) = N * S(t) / total;
    %I(t) = N * I(t) / total;
    %R(t) = N * R(t) / total;
end

% Plot results
plot(1:days, S, 'b', 1:days, I, 'r', 1:days, R, 'g');
xlabel('Days');
ylabel('Population');
legend('Susceptible', 'Infected', 'Removed');
title('SIR Model Simulation');
grid on;

%% Modified SIR Model
% initial parameters
S0 = 980; % initial number of SUSCEPTIBLE
I0 = 20;  % initial number of INFECTED
R0 = 0;   % initial number of REMOVED (including both cured and deceased)
alpha0 = 0.001; % basic infection rate
alpha = 0.00001; % heat propagation rate
beta = 0.01;     % removal rate

% Total population
N = S0 + I0 + R0;

% Simutaltion days
days = 365;

% Initialize the S, I, R array
S = zeros(1, days);
I = zeros(1, days);
R = zeros(1, days);

S(1) = S0;
I(1) = I0;
R(1) = R0;

% Simulate the propagation process every day
for t = 2:days
    % Update the values of S, I, R
    S(t) = S(t-1) - (alpha0 + alpha * I(t-1)) * S(t-1) * I(t-1);
    I(t) = I(t-1) + (alpha0 + alpha * I(t-1)) * S(t-1) * I(t-1) - beta * I(t-1);
    R(t) = R(t-1) + beta * I(t-1);

    % Make sure that S, I, R values do not become negative
    S(t) = max(S(t), 0);
    I(t) = max(I(t), 0);
    R(t) = max(R(t), 0);
end

% Plot Results
plot(1:days, S, 'b', 1:days, I, 'r', 1:days, R, 'g');
xlabel('Days');
ylabel('Population');
legend('Susceptible', 'Infected', 'Removed');
title('Modified SIR Model for Information Spread on a Network');
grid on;
