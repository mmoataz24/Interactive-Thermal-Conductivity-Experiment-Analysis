%% Interactive Thermal Conductivity Experiment Analysis
% This program calculates thermal conductivity (k) for multiple materials
% with complete uncertainty propagation and comparative analysis

clear; clc; close all;

%% Display Header
fprintf('=========================================================\n');
fprintf('   THERMAL CONDUCTIVITY EXPERIMENT ANALYZER\n');
fprintf('=========================================================\n\n');

%% Constants
m = 0.3;                    % Mass of water [kg]
cp = 4186;                  % Specific heat of water [J/(kg·°C)]
fprintf('Constants:\n');
fprintf('  Mass of water (m): %.2f kg\n', m);
fprintf('  Specific heat (cp): %.0f J/(kg·°C)\n\n', cp);

%% Get number of materials to test
num_materials = input('How many materials do you want to test? (e.g., 2): ');
fprintf('\n');

% Initialize storage arrays
materials = cell(num_materials, 1);
all_results = cell(num_materials, 1);

%% Loop through each material
for mat = 1:num_materials
    fprintf('=========================================================\n');
    fprintf('MATERIAL #%d\n', mat);
    fprintf('=========================================================\n');
    
    % Get material name
    materials{mat}.name = input(sprintf('Enter material name (e.g., Aluminum, Copper): '), 's');
    fprintf('\n');
    
    %% Get geometry with uncertainties
    fprintf('--- GEOMETRY MEASUREMENTS ---\n');
    
    % Length
    L_cm = input('Enter rod length in cm (e.g., 16): ');
    delta_L_cm = input('Enter uncertainty in length in cm (e.g., 0.5): ');
    L = L_cm / 100;                % Convert to meters
    delta_L = delta_L_cm / 100;    % Convert to meters
    
    % Diameter/Radius
    d_cm = input('Enter rod diameter in cm (e.g., 2): ');
    delta_d_cm = input('Enter uncertainty in diameter in cm (e.g., 0.5): ');
    r_cm = d_cm / 2;
    delta_r_cm = delta_d_cm / 2;
    r = r_cm / 100;                % Convert to meters
    delta_r = delta_r_cm / 100;    % Convert to meters
    
    % Calculate cross-sectional area with uncertainty
    A = pi * r^2;                  % [m²]
    delta_A = 2 * pi * r * delta_r; % Uncertainty propagation
    
    fprintf('\n--- CALCULATED GEOMETRY ---\n');
    fprintf('  Length: L = %.3f ± %.3f m\n', L, delta_L);
    fprintf('  Radius: r = %.4f ± %.4f m\n', r, delta_r);
    fprintf('  Area: A = (%.2e ± %.2e) m²\n\n', A, delta_A);
    
    % Store geometry
    materials{mat}.L = L;
    materials{mat}.delta_L = delta_L;
    materials{mat}.A = A;
    materials{mat}.delta_A = delta_A;
    
    %% Get number of trials
    num_trials = input('How many trials for this material? (e.g., 4): ');
    fprintf('\n');
    
    % Initialize trial arrays
    Q_dot = zeros(num_trials, 1);
    delta_T_rod = zeros(num_trials, 1);
    k_values = zeros(num_trials, 1);
    trial_data = zeros(num_trials, 5);
    
    %% Input trial data
    for trial = 1:num_trials
        fprintf('--- TRIAL %d ---\n', trial);
        
        T_wi = input('  Initial water temp (T_wi) [°C]: ');
        T_wf = input('  Final water temp (T_wf) [°C]: ');
        T_hot = input('  Hot end temp (T_hot) [°C]: ');
        T_cold = input('  Cold end temp (T_cold) [°C]: ');
        time_input = input('  Time [enter in seconds, or "5m" for minutes]: ', 's');
        
        % Parse time input
        if contains(time_input, 'm')
            time_min = str2double(strrep(time_input, 'm', ''));
            t = time_min * 60;
        else
            t = str2double(time_input);
        end
        
        % Store trial data
        trial_data(trial, :) = [T_wi, T_wf, T_hot, T_cold, t];
        
        % Calculate heat transfer rate Q_dot
        delta_T_water = T_wf - T_wi;
        Q_dot(trial) = m * cp * delta_T_water / t;
        
        % Calculate temperature difference along rod
        delta_T_rod(trial) = T_hot - T_cold;
        
        % Calculate thermal conductivity: Q_dot = k * A * delta_T / L
        % Therefore: k = Q_dot * L / (A * delta_T)
        k_values(trial) = Q_dot(trial) * L / (A * delta_T_rod(trial));
        
        fprintf('  → Q_dot = %.2f W\n', Q_dot(trial));
        fprintf('  → ΔT_rod = %.2f °C\n', delta_T_rod(trial));
        fprintf('  → k = %.1f W/(m·K)\n\n', k_values(trial));
    end
    
    %% Calculate averages and statistics
    Q_dot_avg = mean(Q_dot);
    Q_dot_std = std(Q_dot);
    
    delta_T_avg = mean(delta_T_rod);
    delta_T_std = std(delta_T_rod);
    delta_T_uncertainty = 0.71;  % From measurement uncertainty
    delta_T_sem = delta_T_uncertainty / sqrt(num_trials);  % Standard error of mean
    
    k_avg = mean(k_values);
    k_std = std(k_values);
    k_sem = k_std / sqrt(num_trials);  % Standard error of mean
    
    %% Uncertainty propagation for k
    % δk/k = sqrt((δL/L)² + (δA/A)² + (δΔT/ΔT)²)
    rel_uncertainty_L = delta_L / L;
    rel_uncertainty_A = delta_A / A;
    rel_uncertainty_T = delta_T_sem / delta_T_avg;
    
    rel_uncertainty_k = sqrt(rel_uncertainty_L^2 + ...
                             rel_uncertainty_A^2 + ...
                             rel_uncertainty_T^2);
    delta_k = rel_uncertainty_k * k_avg;
    
    %% Display results
    fprintf('=========================================================\n');
    fprintf('RESULTS FOR %s\n', upper(materials{mat}.name));
    fprintf('=========================================================\n');
    fprintf('Average Heat Transfer Rate:\n');
    fprintf('  Q_dot = %.2f ± %.2f W\n\n', Q_dot_avg, Q_dot_std);
    
    fprintf('Average Temperature Gradient:\n');
    fprintf('  ΔT_rod = %.2f ± %.2f °C\n\n', delta_T_avg, delta_T_sem);
    
    fprintf('Thermal Conductivity:\n');
    fprintf('  k = %.1f ± %.1f W/(m·K)  [from statistical variation]\n', k_avg, k_sem);
    fprintf('  k = (%.1f ± %.1f) × 10² W/(m·K)  [with measurement uncertainty]\n\n', ...
            k_avg/100, delta_k/100);
    
    fprintf('Relative uncertainties:\n');
    fprintf('  δL/L = %.3f (%.1f%%)\n', rel_uncertainty_L, rel_uncertainty_L*100);
    fprintf('  δA/A = %.3f (%.1f%%)\n', rel_uncertainty_A, rel_uncertainty_A*100);
    fprintf('  δΔT/ΔT = %.3f (%.1f%%)\n', rel_uncertainty_T, rel_uncertainty_T*100);
    fprintf('  Total δk/k = %.3f (%.1f%%)\n\n', rel_uncertainty_k, rel_uncertainty_k*100);
    
    % Store results
    materials{mat}.num_trials = num_trials;
    materials{mat}.Q_dot = Q_dot;
    materials{mat}.Q_dot_avg = Q_dot_avg;
    materials{mat}.delta_T_rod = delta_T_rod;
    materials{mat}.delta_T_avg = delta_T_avg;
    materials{mat}.k_values = k_values;
    materials{mat}.k_avg = k_avg;
    materials{mat}.k_sem = k_sem;
    materials{mat}.delta_k = delta_k;
    materials{mat}.trial_data = trial_data;
end %loop through material end

%% Generate Comparison Plots
fprintf('=========================================================\n');
fprintf('GENERATING COMPARISON PLOTS...\n');
fprintf('=========================================================\n\n');

% Create figure with multiple subplots
figure('Position', [100, 100, 1200, 900]);

% Colors for different materials
colors = lines(num_materials);

%% Plot 1: Heat Transfer Rate Comparison
subplot(3, 2, 1);
hold on;
for mat = 1:num_materials
    trials = 1:materials{mat}.num_trials;
    plot(trials, materials{mat}.Q_dot, '-o', 'LineWidth', 2, ...
         'MarkerSize', 8, 'Color', colors(mat,:), ...
         'DisplayName', materials{mat}.name);
end
hold off;
xlabel('Trial Number', 'FontSize', 11);
ylabel('Q_{dot} (W)', 'FontSize', 11);
title('Heat Transfer Rate by Trial', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;

%% Plot 2: Average Heat Transfer Rate with Error Bars
subplot(3, 2, 2);
x_pos = 1:num_materials;
Q_avg_vals = zeros(num_materials, 1);
Q_err_vals = zeros(num_materials, 1);
material_names = cell(num_materials, 1);
for mat = 1:num_materials
    Q_avg_vals(mat) = materials{mat}.Q_dot_avg;
    Q_err_vals(mat) = std(materials{mat}.Q_dot);
    material_names{mat} = materials{mat}.name;
end
bar(x_pos, Q_avg_vals, 'FaceColor', 'flat', 'CData', colors);
hold on;
errorbar(x_pos, Q_avg_vals, Q_err_vals, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
hold off;
set(gca, 'XTick', x_pos, 'XTickLabel', material_names);
ylabel('Average Q_{dot} (W)', 'FontSize', 11);
title('Average Heat Transfer Rate', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%% Plot 3: Temperature Gradient by Trial
subplot(3, 2, 3);
hold on;
for mat = 1:num_materials
    trials = 1:materials{mat}.num_trials;
    plot(trials, materials{mat}.delta_T_rod, '-s', 'LineWidth', 2, ...
         'MarkerSize', 8, 'Color', colors(mat,:), ...
         'DisplayName', materials{mat}.name);
end
hold off;
xlabel('Trial Number', 'FontSize', 11);
ylabel('\DeltaT_{rod} (°C)', 'FontSize', 11);
title('Temperature Gradient Along Rod', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;

%% Plot 4: Thermal Conductivity by Trial
subplot(3, 2, 4);
hold on;
for mat = 1:num_materials
    trials = 1:materials{mat}.num_trials;
    plot(trials, materials{mat}.k_values, '-d', 'LineWidth', 2, ...
         'MarkerSize', 8, 'Color', colors(mat,:), ...
         'DisplayName', materials{mat}.name);
end
hold off;
xlabel('Trial Number', 'FontSize', 11);
ylabel('k (W/(m·K))', 'FontSize', 11);
title('Thermal Conductivity by Trial', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;

%% Plot 5: Average Thermal Conductivity Comparison
subplot(3, 2, 5);
k_avg_vals = zeros(num_materials, 1);
k_err_vals = zeros(num_materials, 1);
for mat = 1:num_materials
    k_avg_vals(mat) = materials{mat}.k_avg;
    k_err_vals(mat) = materials{mat}.delta_k;
end
bar(x_pos, k_avg_vals, 'FaceColor', 'flat', 'CData', colors);
hold on;
errorbar(x_pos, k_avg_vals, k_err_vals, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
hold off;
set(gca, 'XTick', x_pos, 'XTickLabel', material_names);
ylabel('k (W/(m·K))', 'FontSize', 11);
title('Average Thermal Conductivity with Uncertainty', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%% Plot 6: Relative Performance (Normalized to first material)
subplot(3, 2, 6);
k_normalized = k_avg_vals / k_avg_vals(1) * 100;
bar(x_pos, k_normalized, 'FaceColor', 'flat', 'CData', colors);
set(gca, 'XTick', x_pos, 'XTickLabel', material_names);
ylabel('Relative Performance (%)', 'FontSize', 11);
title(sprintf('Thermal Conductivity (Relative to %s)', materials{1}.name), ...
      'FontSize', 12, 'FontWeight', 'bold');
grid on;
ylim([0 max(k_normalized)*1.2]);
% Add percentage labels on bars
for i = 1:num_materials
    text(i, k_normalized(i) + 5, sprintf('%.0f%%', k_normalized(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

sgtitle('Thermal Conductivity Experiment: Comparative Analysis', ...
        'FontSize', 14, 'FontWeight', 'bold');

%% Generate Summary Table
fprintf('=========================================================\n');
fprintf('SUMMARY TABLE\n');
fprintf('=========================================================\n');
fprintf('%-15s %12s %15s %15s\n', 'Material', 'Trials', 'k_avg (W/m·K)', 'Uncertainty');
fprintf('---------------------------------------------------------\n');
for mat = 1:num_materials
    fprintf('%-15s %12d %15.1f %15.1f\n', ...
            materials{mat}.name, ...
            materials{mat}.num_trials, ...
            materials{mat}.k_avg, ...
            materials{mat}.delta_k);
end
fprintf('=========================================================\n\n');

%% Additional Statistical Analysis
if num_materials == 2
    fprintf('COMPARISON BETWEEN TWO MATERIALS:\n');
    fprintf('=========================================================\n');
    
    ratio = materials{2}.k_avg / materials{1}.k_avg;
    fprintf('%s conducts heat %.2fx better than %s\n', ...
            materials{2}.name, ratio, materials{1}.name);
    
    diff = abs(materials{2}.k_avg - materials{1}.k_avg);
    fprintf('Absolute difference: %.1f W/(m·K)\n\n', diff);
end

fprintf('Analysis complete! All plots generated.\n');