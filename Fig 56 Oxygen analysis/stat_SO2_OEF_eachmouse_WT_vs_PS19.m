function plot_SO2_traces_WT_vs_PS19()
    % === Parameters ===
    time = [];  % will be set from first file
    vessel_types = {'a', 'v'};

    % === Prompt user to select WT files ===
    disp("Select CSV files for WT mice");
    [wt_files, wt_path] = uigetfile('*.csv', 'Select WT CSV files', 'MultiSelect', 'on');
    if isequal(wt_files, 0), return; end
    if ischar(wt_files), wt_files = {wt_files}; end

    % === Prompt user to select PS19 files ===
    disp("Select CSV files for PS19 mice");
    [ps_files, ps_path] = uigetfile('*.csv', 'Select PS19 CSV files', 'MultiSelect', 'on');
    if isequal(ps_files, 0), return; end
    if ischar(ps_files), ps_files = {ps_files}; end

    % === Get time vector from first WT file ===
    first_file = readcell(fullfile(wt_path, wt_files{1}));
    time = cell2mat(first_file(1,2:end));

    % === Define analysis windows ===
    baseline_range = find(time >= 10 & time <= 15);
    peak1_range = find(time >= 17 & time <= 19);
    peak2_range = find(time >= 40 & time <= 45);
    stim_range = find(time >= 15 & time <= 45);
    recovery_range = find(time >= 45 & time <= 55);

    % === Process groups ===
    data.WT = process_group(wt_files, wt_path, vessel_types,  'WT');
    data.PS19 = process_group(ps_files, ps_path, vessel_types, 'PS19');

    % === Export trace-level sO2 and OEF for nested analysis ===
    export_trace_level_summary(data, baseline_range, peak1_range, peak2_range, stim_range, recovery_range);

    % === Save data ===
    save('SO2_traces_by_mouse.mat', 'data', 'time');
    disp('Saved raw and normalized SO₂ traces to: SO2_traces_by_mouse.mat');
end

function data = process_group(file_list, path, vessel_types, group_label)
    for v = 1:length(vessel_types)
        vt = vessel_types{v};
        raw.(vt) = [];
        norm.(vt) = [];
        sO2.(vt) = [];
    end
    trace_meta = {};

    for i = 1:length(file_list)
        T = readcell(fullfile(path, file_list{i}));
        labels = string(T(2:end, 1));
        traces = cell2mat(T(2:end, 2:end));
        [~, base_name, ~] = fileparts(file_list{i});

        % Collect arterioles and venules in chunks of 7
        for chunk = 1:floor(size(traces,1)/7)
            chunk_idx = (chunk-1)*7 + (1:7);
            chunk_labels = labels(chunk_idx);
            chunk_traces = traces(chunk_idx, :);

            a_idx = find(chunk_labels == "a");
            v_idx = find(chunk_labels == "v");

            if isempty(a_idx) || isempty(v_idx), continue; end

            a_mean = mean(chunk_traces(a_idx, :), 1, 'omitnan');
            v_mean = mean(chunk_traces(v_idx, :), 1, 'omitnan');

            sa = (a_mean.^2.59) ./ (a_mean.^2.59 + 40.2^2.59);
            sv = (v_mean.^2.59) ./ (v_mean.^2.59 + 40.2^2.59);
            oef = (sa - sv) ./ sa;

            trace_meta(end+1, :) = {group_label, base_name, chunk, sa, sv, oef};
        end
    end

    data.trace_meta = trace_meta;
end

function export_trace_level_summary(data, baseline_range, peak1_range, peak2_range, stim_range, recovery_range)
    groups = {'WT', 'PS19'};
    output = {};
    header = {'Group', 'MouseID', 'Chunk','sO2_Peak1', 'sO2_Peak2', 'OEF_baseline', 'OEF_stim', 'OEF_Peak1', 'OEF_Peak2','OEF_recovery'};

    for g = 1:2
        group = groups{g};
        meta = data.(group).trace_meta;

        for i = 1:size(meta, 1)
            group_label = meta{i,1};
            mouse_id = meta{i,2};
            chunk_id = meta{i,3};
            sa = meta{i,4};
            sv = meta{i,5};
            oef = meta{i,6};

            sO2_1 = mean(sa(peak1_range), 'omitnan');
            sO2_2 = mean(sa(peak2_range), 'omitnan');
            oef1 = mean(oef(peak1_range), 'omitnan');
            oef2 = mean(oef(peak2_range), 'omitnan');
            baseline_1=mean(oef(baseline_range), 'omitnan');
            stim_1=mean(oef(stim_range), 'omitnan');
            recovery_1=mean(oef(recovery_range), 'omitnan');

            output(end+1, :) = {group_label, mouse_id, chunk_id,  sO2_1, sO2_2, baseline_1, oef1, oef2, stim_1, recovery_1};
        end
    end

    cell2csv('sO2_OEF_trace_summary.csv', [header; output]);
    disp('Exported per-chunk sO2 and OEF summary to sO2_OEF_trace_summary.csv');
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
