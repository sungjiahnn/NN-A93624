function analyze_pO2_nested_normalized()
    % === Parameters ===
    baseline_range = 11:17;
    dip_range=18:20;
    first_peak_range = 19:22;
    second_peak_range = 41:52;

    % === Select WT files ===
    disp("Select CSV files for WT mice");
    [wt_files, wt_path] = uigetfile('*.csv', 'Select WT CSV files', 'MultiSelect', 'on');
    if isequal(wt_files, 0), return; end
    if ischar(wt_files), wt_files = {wt_files}; end

    % === Select PS19 files ===
    disp("Select CSV files for PS19 mice");
    [ps_files, ps_path] = uigetfile('*.csv', 'Select PS19 CSV files', 'MultiSelect', 'on');
    if isequal(ps_files, 0), return; end
    if ischar(ps_files), ps_files = {ps_files}; end

    % === Process and compile long-format data ===
    all_data = table();
    all_data = process_mouse_group_normalized(wt_files, wt_path, 'WT', all_data, baseline_range, dip_range, first_peak_range, second_peak_range);
    all_data = process_mouse_group_normalized(ps_files, ps_path, 'PS19', all_data, baseline_range,dip_range, first_peak_range, second_peak_range);

    % === Save long-format normalized data ===
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    longname = ['normalized_nested_pO2_peaks_' timestamp '.csv'];
    writetable(all_data, longname);
    disp(['Saved normalized long-format data to: ' longname]);

    % === Run nested t-tests per vessel type ===
    run_nested_tests(all_data);
    % === Export per-mouse summary ===
    export_per_mouse_summary(all_data);
end

function all_data = process_mouse_group_normalized(file_list, path, genotype, all_data, baseline_range, range1, range2, range3)
    for i = 1:length(file_list)
        file = file_list{i};
        mouse_id = erase(file, '.csv');  % use filename as ID

        T = readcell(fullfile(path, file));
        labels = string(T(2:end, 1));           % vessel type
        traces = cell2mat(T(2:end, 2:end));     % trace values

        for r = 1:size(traces, 1)
            vessel = labels(r);
            trace = traces(r, :);

            if length(trace) < max(range2), continue; end

            % === Normalize by baseline ===
            baseline = mean(trace(baseline_range), 'omitnan');
            if isnan(baseline) || baseline == 0, continue; end  % skip bad data
            norm_trace = trace / baseline;

            % === Extract peaks from normalized trace ===
            val1 = min(norm_trace(range1), [], 'omitnan');%mean(norm_trace(range1), 'omitnan'); %%%%%
            val2 = max(norm_trace(range1), [], 'omitnan');
            val3 = mean(norm_trace(range2), 'omitnan');

            temp = table( ...
                repmat(string(mouse_id), 3, 1), ...
                repmat(string(genotype), 3, 1), ...
                repmat(vessel, 3, 1), ...
                ["dip";"first"; "second"], ...
                [val1; val2; val3], ...
                'VariableNames', {'mouse_id', 'genotype', 'vessel_type', 'peak', 'value'} ...
            );

            all_data = [all_data; temp];
        end
    end
end

function run_nested_tests(all_data)
    fprintf('\n=== Running Nested T-Tests per Vessel Type ===\n');

    all_data.mouse_id = categorical(all_data.mouse_id);
    all_data.genotype = categorical(all_data.genotype);
    all_data.vessel_type = categorical(all_data.vessel_type);
    all_data.peak = categorical(all_data.peak);

    vessel_types = categories(all_data.vessel_type);

    for v = 1:length(vessel_types)
        vt = vessel_types{v};
        tbl = all_data(all_data.vessel_type == vt, :);

        fprintf('\n--- Vessel Type: %s ---\n', vt);
        lme = fitlme(tbl, 'value ~ genotype + (1|mouse_id)');
        disp(lme);
    end
end

function export_per_mouse_summary(all_data)
    all_data.mouse_id = string(all_data.mouse_id);
    all_data.genotype = string(all_data.genotype);
    all_data.vessel_type = string(all_data.vessel_type);
    all_data.peak = string(all_data.peak);

    summary = groupsummary(all_data, {'mouse_id', 'genotype', 'vessel_type', 'peak'}, 'mean', 'value');
    wide_tbl = unstack(summary, 'mean_value', 'peak');
    wide_tbl = movevars(wide_tbl, {'mouse_id', 'genotype', 'vessel_type'}, 'Before', 1);

    % Rename
    wide_tbl.Properties.VariableNames{'dip'} = 'dip_min';
    wide_tbl.Properties.VariableNames{'first'} = 'first_peak';
    wide_tbl.Properties.VariableNames{'second'} = 'second_peak';

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = ['normalized_mouse_summary_' timestamp '.csv'];
    writetable(wide_tbl, filename);
    disp(['Per-mouse normalized summary saved to: ' filename]);
end
