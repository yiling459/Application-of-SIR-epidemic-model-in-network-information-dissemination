%% SIR Model Initialization
S0 = 5e6;      % Initial number of SUSCEPTIBLE
I0 = 80630;    % Initial number of INFECTED
R0 = 0;        % Initial number of REMOVED
% alpha0 = 1.2e-6; % Base infection rate 
% beta = 1.4e-7;  % Recovery rate 

% alpha0 = 3.2e-8;      
% beta = 0.1228;  
% const = 1.6e-10;

alpha0 = 3.89e-8;
beta = 0.1384;
const = 2.6e-11;
days = 359; % Simulation based on the data length


S = zeros(1, days);
I = zeros(1, days);
R = zeros(1, days);
S(1) = S0;
I(1) = I0;
R(1) = R0;


%% SIR Model Simulation
for t = 2:days
    % Calculate dynamic infection rate alpha(I(t)) based on the given formula
    N = S(t-1) + I(t-1) + R(t-1);
    alphaI_t = (1 - alpha0) * I(t-1) / N * const;
  
    S(t) = S(t-1) - (alpha0 + alphaI_t) * S(t-1) * I(t-1);
    I(t) = I(t-1) + (alpha0 + alphaI_t) * S(t-1) * I(t-1) - beta * I(t-1);
    R(t) = R(t-1) + beta * I(t-1);
    
    % Ensure values don't drop below zero
    S(t) = max(S(t), 0);
    I(t) = max(I(t), 0);
    R(t) = max(R(t), 0);
end

%% Plot Results
figure;
hold on;
plot(1:days, S, 'b', 'LineWidth', 2, 'DisplayName', 'Susceptible');
plot(1:days, I, 'r', 'LineWidth', 2, 'DisplayName', 'Infected');
plot(1:days, R, 'g', 'LineWidth', 2, 'DisplayName', 'Removed');
xlabel('Days');
ylabel('Population');
legend('show');
title('Modified SIR Model with Dynamic Infection Rate');
grid on;

%% Calculate the prediction accuracy
filename = 'Problem_C_Data_Wordle.xlsx';  
dataTable = readtable(filename);  
I_actual = dataTable.Infected;  
I_actual = I_actual(~isnan(I_actual));
days =  359;

% plot the graph
figure;
hold on;
plot(1:days, I, 'r', 'LineWidth', 2, 'DisplayName', 'Predicted Infected');
plot(1:days, I_actual, 'b--', 'LineWidth', 2, 'DisplayName', 'Actual Infected');
xlabel('Days');
ylabel('Number of Infected Individuals');
title('Comparison of Predicted and Actual Infection Data');
legend('show');
grid on;
