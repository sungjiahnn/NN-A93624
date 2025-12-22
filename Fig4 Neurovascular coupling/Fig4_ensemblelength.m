%% length of ensemble. plot diameter and significant ensemble
% Sung Ji Ahn 2025

clear all
%close all

load('ensemblere.mat','pcntactive', 'shuffmean');
load('VesselA2.mat', 'diameter');
%load('Vesselpostprocess.mat', 'VesselA1');
%diameter=VesselA1';
%T = readtable('vessel.csv');
%diameter = T{:,2};

% Example data
t = 1:length(pcntactive);  % time axis (adjust to your actual time vector)

% resample for 1140 points
%x_old = 1:1140;                 % original time points
%x_new = linspace(1, 1140, 1150); % new resampled grid
%diameter = interp1(x_old, diameter, x_new, 'linear');


Fs = 3.839;               % samples per second
d = designfilt('lowpassiir','FilterOrder',4, ...      % a bit steeper
    'HalfPowerFrequency',0.7,'DesignMethod','butter', ...
    'SampleRate',Fs);
filtdiameter = filtfilt(d, diameter);

% normalize to lower 20%
sortedVals = sort(filtdiameter(:), 'ascend');
nBase = round(0.3 * numel(sortedVals));          % number of points in lower 20%
baseline = median(sortedVals(1:nBase), 'omitnan'); % mean of lowest 20%
norm_diameter = filtdiameter / baseline;   

%filtdiameter = lowpass(diameter,0.7,Fs);
%filtdiameter = sgolayfilt(diameter, 3, 11);

%x_med = smoothdata(diameter,'movmedian',5);       % kill spikes (≈1.3 s window)
%filtdiameter  = lowpass(x_med, 0.7, Fs);   

filtE=sgolayfilt(pcntactive,3,11);
signensemble = filtE > shuffmean;

figure;
yyaxis left
%plot(t, diameter, 'c', 'LineWidth', 0.5);
%ylim([min(diameter) max(diameter)]); 
hold on;
%plot(t, filtdiameter, 'b', 'LineWidth', 1.5);
plot(t, norm_diameter, '-b', 'LineWidth', 1.5);
ylabel('Diameter (\mum)')
xlabel('Time')

% Find start and end indices of consecutive 1s in sighensemble
d = diff([0 signensemble 0]);
start_idx = find(d == 1);
end_idx   = find(d == -1) - 1;

% Shade regions where sighensemble == 1
yl = ylim;
for i = 1:length(start_idx)
    x = [t(start_idx(i)) t(end_idx(i)) t(end_idx(i)) t(start_idx(i))];
    y = [yl(1) yl(1) yl(2) yl(2)];
    fill(x, y, [0.9 0.9 0.9], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
end

yyaxis right
%plot(t, pcntactive, 'r', 'LineWidth', 0.5);
plot(t, filtE, 'r', 'LineWidth', 1);
ylabel('Ensemble activity (0–100)')
ylim([0 100]);

xlabel('Time (s)');
title('Diameter trace with shaded stimulation periods');

% Find transitions and segment lengths
d = diff([0 signensemble 0]);        % pad with zeros to catch edges
start_idx = find(d == 1); % where 0 → 1
end_idx   = find(d == -1) - 1; % where 1 → 0
run_lengths = end_idx - start_idx + 1;

% Optional: show results
%disp('Start indices of 1s:'); disp(start_idx);
%disp('Lengths of 1s:'); disp(run_lengths);

pc  = pcntactive(:); 
thr = shuffmean(:);  

% Sampling info (edit Fs if your sampling rate is known)
dt = 1/Fs;

% Per-bout AUC (area above threshold), plus some handy stats
nBout = numel(start_idx);
auc_above = zeros(nBout,1);
peak_above = zeros(nBout,1);
mean_above = zeros(nBout,1);

for k = 1:nBout
    s = start_idx(k); e = end_idx(k);
    y = pc(s:e) - thr(s:e);    % amount above threshold
    y(y < 0) = 0;              % clip below-threshold to zero
    % AUC in “(pcntactive units) * seconds”
    auc_above(k)  = trapz(y) * dt;
    peak_above(k) = max(y, [], 'omitnan');
    mean_above(k) = mean(y, 'omitnan');
end

start_idx   = start_idx(:);
end_idx     = end_idx(:);
run_lengths = run_lengths(:);
duration_s  = run_lengths * dt;

% summary table
BoutSummary = table( ...
    start_idx, end_idx, run_lengths, duration_s, ...
    auc_above, peak_above, mean_above, ...
    'VariableNames', {'start','stop','length_frames','duration_s','auc_above_thr','peak_above','mean_above'});

disp(BoutSummary);

% save figure and workspace
saveas(gcf, 'couplinglength.fig');  
save('couplinglength.mat')

% Optional: save results
%writetable(BoutSummary, 'ensemble_bout_summary.csv');

% Also report totals if useful
%total_auc = sum(auc_above, 'omitnan');
%fprintf('Total AUC above threshold across all bouts: %.3f (units: pcntactive*sec)\n', total_auc);