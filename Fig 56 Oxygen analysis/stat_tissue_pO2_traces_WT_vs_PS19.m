function analyze_pO2_nested()
    % === Parameters ===
    base_range = 12:17; % 10:15 
    stim_range = 18:52; % 
    first_peak_range = 19:21; % 17:19
    second_peak_range = 41:52; % 35:45

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
    all_data = process_mouse_group(wt_files, wt_path, 'WT', all_data, base_range, stim_range, first_peak_range, second_peak_range);
    all_data = process_mouse_group(ps_files, ps_path, 'PS19', all_data, base_range, stim_range, first_peak_range, second_peak_range);

    % === Save long-format data ===
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    longname = ['nested_pO2_peaks_' timestamp '.csv'];
    writetable(all_data, longname);
    disp(['Saved long-format data to: ' longname]);

    % === Export per-mouse summary ===
    export_per_mouse_summary(all_data);
end

function all_data = process_mouse_group(file_list, path, genotype, all_data, range1, range2, range3, range4)
    dt = 0.8638;

    for i = 1:length(file_list)
        file = file_list{i};
        mouse_id = erase(file, '.csv');  % use filename as mouse_id

        T = readcell(fullfile(path, file));
        labels = string(T(2:end, 1));           % vessel type
        traces = cell2mat(T(2:end, 2:end));     % trace values

        for r = 1:size(traces, 1)
            vessel = labels(r);
            trace = traces(r, :);

            neededLen = max([range1(:); range2(:); range3(:); range4(:)]);
            if length(trace) < neededLen
                continue;
            end
            

            val1 = mean(trace(range1), 'omitnan'); % baseline

            trace_r2= trace(range2);
            %pct_r2 = (trace_r2 - val1) / val1 * 100; % calculate percent
            %pct_r2(pct_r2 < 0) = 0; % above 0
            net_r2 = trace_r2 - val1;
            %net_r2(net_r2 < 0) = 0; % above 0
            val2 = trapz(net_r2) * dt; % stimulation AUC
            
            val3 = max(trace(range3), [], 'omitnan'); % first peak
            val4 = mean(trace(range4), 'omitnan'); % second peak

            peaks  = ["baseline"; "stim"; "first"; "second"];
            values = [val1; val2; val3; val4];


            temp = table( ...
                repmat(string(mouse_id), numel(peaks), 1), ...
                repmat(string(genotype), numel(peaks), 1), ...
                repmat(vessel,          numel(peaks), 1), ...
                peaks, ...
                values, ...
                'VariableNames', {'mouse_id', 'genotype', 'vessel_type', 'peak', 'value'} ...
            );

            all_data = [all_data; temp];
        end
    end
end


function export_per_mouse_summary(all_data)
    % Convert to string for grouping
    all_data.mouse_id = string(all_data.mouse_id);
    all_data.genotype = string(all_data.genotype);
    all_data.vessel_type = string(all_data.vessel_type);
    all_data.peak = string(all_data.peak);

    % Group and average
    summary = groupsummary(all_data, {'mouse_id', 'genotype', 'vessel_type', 'peak'}, 'mean', 'value');

    % Pivot to wide format
    wide_tbl = unstack(summary, 'mean_value', 'peak');

    % Reorder columns
    wide_tbl = movevars(wide_tbl, {'mouse_id', 'genotype', 'vessel_type'}, 'Before', 1);

    % Rename columns
    wide_tbl.Properties.VariableNames{'first'} = 'first_peak';
    wide_tbl.Properties.VariableNames{'second'} = 'second_peak';

    % Save
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = ['mouse_summary_' timestamp '.csv'];
    writetable(wide_tbl, filename);
    disp(['Per-mouse summary saved to: ' filename]);
end
