% Sung Ji Ahn 2025
% plot intravenous pO2 traces A/V, export first peak and shallow peak values for
% stats

function plot_pO2_traces_WT_vs_PS19()

    % Parameters
    resample_time = 0:0.1:60;  % 10 Hz
    baseline_range = find(resample_time >= 10 & resample_time <= 15);
    vessel_types = {'a', 'v'};
    resample_time = 0:0.1:60;  % 10 Hz smoothing
    %time = [];  % will store time vector from first file

    % Select WT files 
    disp("Select CSV files for WT mice");
    [wt_files, wt_path] = uigetfile('*.csv', 'Select WT CSV files', 'MultiSelect', 'on');
    if isequal(wt_files, 0), return; end
    if ischar(wt_files), wt_files = {wt_files}; end

    % Select PS19 files
    disp("Select CSV files for PS19 mice");
    [ps_files, ps_path] = uigetfile('*.csv', 'Select PS19 CSV files', 'MultiSelect', 'on');
    if isequal(ps_files, 0), return; end
    if ischar(ps_files), ps_files = {ps_files}; end

    % Get time vector from first WT file
    first_file = readcell(fullfile(wt_path, wt_files{1}));
    time = cell2mat(first_file(1,2:end));
    peak1_range = find(time >= 17 & time <= 19);
    peak2_range = find(time >= 35 & time <= 45);

    % Process groups with smoothing 
    data.WT = process_group(wt_files, wt_path, vessel_types, baseline_range, time, resample_time);
    data.PS19 = process_group(ps_files, ps_path, vessel_types, baseline_range, time, resample_time);

    % Plot for each vessel type 
    for v = 1:length(vessel_types)
        vt = vessel_types{v};
        plot_group(data, resample_time, vt, false);
        plot_group(data, resample_time, vt, true);
    end
    
    % Export sO2 and OEF peaks to CSV 
    export_peak_summary(data, time, peak1_range, peak2_range);

    save('pO2_traces_by_mouse.mat', 'data', 'resample_time');
    disp('Saved raw and normalized pO₂ traces to: pO2_traces_by_mouse.mat');
end

function group_data = process_group(file_list, path, vessel_types, baseline_range, time, resample_time)
    group_data.raw = struct();
    group_data.norm = struct();
    group_data.smooth = struct();

    for v = 1:length(vessel_types)
        group_data.raw.(vessel_types{v}) = [];
        group_data.norm.(vessel_types{v}) = [];
        group_data.smooth.(vessel_types{v}) = [];
    end

    for i = 1:length(file_list)
        T = readcell(fullfile(path, file_list{i}));
        labels = string(T(2:end, 1));
        traces = cell2mat(T(2:end, 2:end));

        for v = 1:length(vessel_types)
            vt = vessel_types{v};
            idx = labels == vt;
            if ~any(idx), continue; end

            vessel_traces = traces(idx, :);
            smoothed = nan(size(vessel_traces, 1), length(resample_time));

            for j = 1:size(vessel_traces, 1)
                smoothed(j, :) = interp1(time, vessel_traces(j, :), resample_time, 'pchip');
            end

            mean_raw = mean(vessel_traces, 1, 'omitnan');
            mean_smooth = mean(smoothed, 1, 'omitnan');
            baseline = mean(mean_smooth(baseline_range), 'omitnan');
            mean_norm = mean_smooth / baseline;

            group_data.raw.(vt)(end+1, :) = mean_raw;
            group_data.norm.(vt)(end+1, :) = mean_norm;
            group_data.smooth.(vt)(end+1, :) = mean_smooth;
        end
        %i
    end
end

function plot_group(data, time, vessel_type, normalized)
    % Prepare data 
    groupnames = {'WT', 'PS19'};
    color_map = {[0 0.45 0.74], [0.85 0.33 0.1]};
    fill_map  = {[0.7 0.85 1], [1 0.8 0.7]};
    fig_tag = sprintf('%s - %s', vessel_type, tern(normalized, 'Normalized', 'Raw'));

    figure('Name', fig_tag); hold on;

    for g = 1:2
        name = groupnames{g};
        traces = tern(normalized, data.(name).norm.(vessel_type), data.(name).smooth.(vessel_type));
        if isempty(traces), continue; end

        mean_trace = nanmean(traces, 1);
        sem_trace = nanstd(traces, 0, 1) ./ sqrt(size(traces, 1));

        % Shaded error
        fill([time, fliplr(time)], ...
             [mean_trace + sem_trace, fliplr(mean_trace - sem_trace)], ...
             fill_map{g}, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(time, mean_trace, 'Color', color_map{g}, 'LineWidth', 2, ...
             'DisplayName', name);
    end

    yline(1, '--', 'Color', [0.5 0.5 0.5]);
    xlabel('Time (s)');
    ylabel(tern(normalized, 'Normalized pO₂', 'pO₂ (mmHg)'));
    title(sprintf('%s Traces (%s)', upper(vessel_type), tern(normalized, 'Normalized', 'Raw')));
    legend('Location', 'best');
    grid on; box on; hold off;
end

function export_peak_summary(data, time, peak1_range, peak2_range)
    groups = {'WT', 'PS19'};
    vessels = {'a', 'v'};
    peak_data = {};
    header = {'Group', 'Vessel', 'MouseID', 'pO2_Peak1', 'pO2_Peak2'};

    for g = 1:2
        group = groups{g};
        for v = 1:2
            vessel = vessels{v};
            pO2_mat = data.(group).raw.(vessel);
            for i = 1:size(pO2_mat, 1)
                pO2_1 = mean(pO2_mat(i, peak1_range), 'omitnan');
                pO2_2 = mean(pO2_mat(i, peak2_range), 'omitnan');

                peak_data(end+1, :) = {group, vessel, i, pO2_1, pO2_2};
            end
        end
    end

    % Write to CSV
    cell2csv('pO2_peaks_by_mouse.csv', [header; peak_data]);
    disp('Exported pO2peaks to pO2_peaks_by_mouse.csv');
end
function cell2csv(filename, cellArray)
    fid = fopen(filename, 'w');
    for i = 1:size(cellArray, 1)
        for j = 1:size(cellArray, 2)
            var = cellArray{i, j};
            if isnumeric(var)
                fprintf(fid, '%g', var);
            else
                fprintf(fid, '%s', var);
            end
            if j ~= size(cellArray, 2)
                fprintf(fid, ',');
            end
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end

function out = tern(cond, a, b)
    % Ternary operator
    out = a;
    if ~cond, out = b; end
end
