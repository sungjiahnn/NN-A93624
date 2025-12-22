% SR101 signal correct NADH
% Sung Ji Ahn 2025

clear all
close all

% Load only the 'Mean' values (second column) from CSVs
nadh_data = readmatrix('ValuesN.csv');
sr101_data = readmatrix('ValuesS.csv');

% Extract only second column (mean intensity values)
nadh = nadh_data(:, 2);
sr101 = sr101_data(:, 2);

% Make sure they are column vectors
nadh = nadh(:);
sr101 = sr101(:);

% Parameters
baseline_frames = 11:14;  % Frame indices for baseline
K =1.15;                % Correction factor

% Normalize each to its own baseline (ΔF/F)
nadh_baseline = mean(nadh(baseline_frames));
sr101_baseline = mean(sr101(baseline_frames));

nadh_dff = (nadh - nadh_baseline) / nadh_baseline;
sr101_dff = (sr101 - sr101_baseline) / sr101_baseline;

% Correct NADH signal
corrected_nadh = nadh_dff - K * sr101_dff;

% Time vector (assuming 1.087 sec per frame)
T = length(nadh);
time = (1:T) * 1.087;

% Plot
figure;
plot(time, nadh_dff * 100, 'k-', 'LineWidth', 1.5); hold on;
plot(time, sr101_dff * 100, 'r--', 'LineWidth', 1.5);
plot(time, corrected_nadh * 100, 'b-', 'LineWidth', 2);
xline(15, '--k', 'Stim ON');
xline(45, '--k', 'Stim OFF');
xlabel('Time (s)');
ylabel('ΔF/F (%)');
legend('NADH (raw)', 'SR101', 'Corrected NADH');
title('NADH Correction with SR101');
grid on;

%save('correctedNADH_workspace.mat', 'nadh_dff', 'sr101_dff', 'corrected_nadh', 'time');
%savefig('correctedNADH_plot.fig');

corrected_nadh'