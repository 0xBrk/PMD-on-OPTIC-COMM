%% SECTION 1: Setup and Constans Parameters (for 10 Gb/s)
%% clear all;
close all;
clc;

% 10 Gb/s System Parameters (From project_2 document at section 4)
T0_10G = 20;       % Input Pulse Width (ps)
Tb_10G = 100;      % Bit period (ps)
t_pulse_10G = linspace(-Tb_10G/2, Tb_10G/2, 201); 
dt_10G = t_pulse_10G(2) - t_pulse_10G(1); 

% --- Original Project Parameters for PMD and Lenght  ---
DPMD_list = [0.05, 0.1, 0.2]; % ps/sqrt(km)
L_list = [50, 100, 200];     % km

%% SECTION 2: Creating Input Pulse Signal (For Single Pulse Analysis)
% Time Vector (According to project documentation)
t_single = linspace(-100, 100, 401); 
% Gaussian input pulse (from project file section 5)
E0_single = exp(-(t_single/T0_10G).^2); 
I_in_single = abs(E0_single).^2; 

%% SECTION 3: Simulaton Loop (Single Pulse Graphs) 

% To store reslt we created 3x3 blank matris
penalty_results = zeros(length(DPMD_list), length(L_list));

fprintf('----- SECTION 3 & 4: 10G Single Pulse Analysis (With Documentation Variables) -----\n');
fprintf('DPMD (ps/sqrt(km))\t L (km)\t DGD (ps)\t Penalty Ratio (10G)\n');
fprintf('-------------------------------------------------------------\n');

for i = 1:length(DPMD_list)
    for j = 1:length(L_list)
        
        DPMD_exist = DPMD_list(i);
        L_exist = L_list(j);
        
        % --- Step 1 (DGD Calculation) ---
        DGD = DPMD_exist * sqrt(L_exist); 
        
        % --- Step 2 (Broadening Simulation)  ---
        E_fast_single = exp(-((t_single + DGD/2)/T0_10G).^2); 
        E_slow_single = exp(-((t_single - DGD/2)/T0_10G).^2); 
        I_out_single = (abs(E_fast_single).^2 + abs(E_slow_single).^2) / 2; 
        
        % --- Step 3 (Graph Drawing)  ---
        figure; % Pulse graph for each situation (totaly 9 graphs will be creating)
        plot(t_single, I_in_single, 'b-', 'LineWidth', 2, 'DisplayName', 'Input Pulse (I_in)');
        hold on; 
        plot(t_single, I_out_single, 'r--', 'LineWidth', 2, 'DisplayName', 'Output Pulse (I_out)');
        hold off; 
        title(sprintf('10G Single Pulse Broadening : D_PMD = %.2f, L = %d km, DGD = %.3f ps', DPMD_exist, L_exist, DGD));
        xlabel('Time (ps)');
        ylabel('Normalized Intensity');
        legend;
        grid on;
        
        % --- Step 4 (Penalty Ratio Calculation)  ---
        penalty_ratio = (DGD / Tb_10G)^2; 
        fprintf('%.2f\t\t\t %d\t %.3f\t\t %.6f\n', DPMD_exist, L_exist, DGD, penalty_ratio);
        penalty_results(i, j) = penalty_ratio;
        
    end
end
fprintf('-------------------------------------------------------------\n');

%% SECTION 4: Creating  Analysis Results Graphs 

% --- Graph 1: Penalty Ratio vs. Lenght (L) ---
figure; % Figure 10
plot(L_list, penalty_results(1, :), 'b-o', 'LineWidth', 2, 'DisplayName', 'D_{PMD} = 0.05');
hold on;
plot(L_list, penalty_results(2, :), 'g-s', 'LineWidth', 2, 'DisplayName', 'D_{PMD} = 0.1');
plot(L_list, penalty_results(3, :), 'r-^', 'LineWidth', 2, 'DisplayName', 'D_{PMD} = 0.2');
hold off;
title('Analysis 1: Dependence of the Penalty Ratio on Fiber Length (L)');
xlabel('Fiber Length (L) (km)');
ylabel('Penalty Ratio (No Unit)');
legend('show', 'Location', 'northwest');
grid on;

% --- Analysis Graph 2: Penalty vs. PMD Coefficient (DPMD) ---
figure; % Figure 11
plot(DPMD_list, penalty_results(:, 1), 'b-o', 'LineWidth', 2, 'DisplayName', 'L = 50 km');
hold on;
plot(DPMD_list, penalty_results(:, 2), 'g-s', 'LineWidth', 2, 'DisplayName', 'L = 100 km');
plot(DPMD_list, penalty_results(:, 3), 'r-^', 'LineWidth', 2, 'DisplayName', 'L = 200 km');
hold off;

% To fix the MATLAB text error, ‘Interpreter’ was added to ‘none’.
title('Analysis 2: Dependence of the Penalty Rate on the PMD Coefficient (D_{PMD})', 'Interpreter', 'none');
xlabel('PMD Coefficient D_{PMD} (ps/\sqrt{km})', 'Interpreter', 'none');
ylabel('Penalty Rate (Unitless)');
legend('show', 'Location', 'northwest');
grid on;


%% Section 5: EYE DIAGRAM SIMULATION (The Worst Situaiton)

fprintf('----- CHAPTER 5: Eye Diagram Simulation Begins -----\n');

% --- 5.1. Creating Signal Series ---
num_bits = 200; % Simulated bit amount
bits = randi([0 1], 1, num_bits); % Random bit array

samples_per_bit = length(t_pulse_10G);
single_pulse_E = exp(-(t_pulse_10G/T0_10G).^2); % T0=20ps
E_in_train = []; % Empty input signal train
for k = 1:num_bits
    if bits(k) == 1
        E_in_train = [E_in_train, single_pulse_E];
    else
        E_in_train = [E_in_train, zeros(1, samples_per_bit)];
    end
end

% --- 5.2. Applying Worst-Case Scenario Fiber to Signal ---
% Worst-case scenario: DPMD = 0.2, L = 200 (Original parameters)
i = 3; 
j = 3;
DPMD = DPMD_list(i); % DPMD = 0.2
L = L_list(j);       % L = 200
DGD = DPMD * sqrt(L); % DGD = 0.2 * sqrt(200) ≈ 2.828 ps
shift_samples = round(DGD / dt_10G); % Convert the delay into the number of samples

% Apply PMD to the entire signal train
E_slow_train = [zeros(1, shift_samples), E_in_train(1:end-shift_samples)];
E_fast_train = E_in_train; 
I_out_train = (abs(E_fast_train).^2 + abs(E_slow_train).^2) / 2; 

% --- 5.3. Draw the Eye Diagram Manually  ---
figure; % Figure 12 (Worst-Case Scenario Eye Diagram)
hold on;

% Draw an eye with a width of 2 bit periods (2*Tb = 200 ps)
samples_per_eye = 2 * samples_per_bit;
t_eye = (0:samples_per_eye-1) * dt_10G;
num_eyes = floor(length(I_out_train) / samples_per_eye);

for k = 10 : num_eyes-10 % Skip the transitions at the beginning and end
    start_index = (k-1) * samples_per_eye + 1;
    end_index = start_index + samples_per_eye - 1;
    
    % Plot the PMD output signal (distorted)
    plot(t_eye, I_out_train(start_index:end_index), 'r');
    % Also plot the input signal (uncorrupted) as a reference.
    plot(t_eye, abs(E_in_train(start_index:end_index)).^2, 'b', 'LineWidth', 0.1); 
end

title(sprintf('EYE DIAGRAM (Worst Case): D_PMD = %.2f, L = %d km, DGD = %.3f ps', DPMD, L, DGD));
xlabel('Time (2 * T_b = 200 ps)');
ylabel('Normalized Density');
grid on;
hold off;

fprintf('All tasks are done exiting.\n');
