function plot_sO2_traces_WT_vs_PS19()
    % === Parameters ===
    peak1_range = [];
    peak2_range = [];
    vessel_types = {'a', 'v'};
    resample_time = 0:0.1:60;  % 10 Hz smoothing

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
    peak1_range = find(time >= 17 & time <= 19);
    peak2_range = find(time >= 35 & time <= 45);
    base_range = find(time >= 10 & time <= 15);
    stim_range = find(time >= 15 & time <= 45);
    recovery_range = find(time >= 45 & time <= 55);

    % === Process groups ===
    data.WT = process_group(wt_files, wt_path, vessel_types, resample_time, time);
    data.PS19 = process_group(ps_files, ps_path, vessel_types, resample_time, time);

    % === Plot sO2 and OEF ===
    plot_sO2_and_OEF(data, time);

    % === Export sO2 and OEF peaks to CSV ===
    export_peak_summary(data, time, peak1_range, peak2_range, base_range, stim_range, recovery_range);

    % === Export full OEF traces to CSV ===
    export_oef_traces(data, time);

end

function data = process_group(file_list, path, vessel_types, resample_time, time)
    for v = 1:length(vessel_types)
        vt = vessel_types{v};
        raw.(vt) = [];
        sO2.(vt) = [];
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
            %smoothed = nan(size(vessel_traces, 1), length(resample_time));

            %for j = 1:size(vessel_traces, 1)
            %    smoothed(j, :) = interp1(time, vessel_traces(j, :), resample_time, 'pchip');
            %end

            mean_raw = mean(vessel_traces, 1, 'omitnan');
            sat = (mean_raw.^2.59) ./ (mean_raw.^2.59 + 40.2^2.59);

            raw.(vt)(end+1, :) = mean_raw;
            sO2.(vt)(end+1, :) = sat;
        end
    end

    data.raw = raw;
    data.sO2 = sO2;
end


function plot_sO2_and_OEF(data, time)
    figure('Name', 'sO2 and OEF');
    vessel_types = {'a', 'v'};
    groupnames = {'WT', 'PS19'};
    colors = {[0 0.45 0.74], [0.85 0.33 0.1]};
    
    figure; hold on;
    for g = 1:2
        group = groupnames{g};
        sa = data.(group).sO2.a;
        sv = data.(group).sO2.v;

        m_sa = mean(sa, 1, 'omitnan');
        se_sa = std(sa, 0, 1, 'omitnan') ./ sqrt(size(sa, 1));
        fill([time, fliplr(time)], [m_sa+se_sa, fliplr(m_sa-se_sa)], colors{g}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(time, m_sa, 'Color', colors{g}, 'LineWidth', 2);
        m_sv = mean(sv, 1, 'omitnan');
        se_sv = std(sv, 0, 1, 'omitnan') ./ sqrt(size(sv, 1));
        fill([time, fliplr(time)], [m_sv+se_sv, fliplr(m_sv-se_sv)], colors{g}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(time, m_sv, 'Color', colors{g}, 'LineWidth', 2);
        xlim([5 60]);  
        ylabel('sO₂'); title('sO₂');
    end

    figure; hold on;
    for g = 1:2
        group = groupnames{g};
        sa = data.(group).sO2.a;
        sv = data.(group).sO2.v;
        oef = (sa - sv) ./ sa;
        m_oef = mean(oef, 1, 'omitnan');
        se_oef = std(oef, 0, 1, 'omitnan') ./ sqrt(size(oef, 1));
        fill([time, fliplr(time)], [m_oef+se_oef, fliplr(m_oef-se_oef)], colors{g}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(time, m_oef, 'Color', colors{g}, 'LineWidth', 2);
        xlim([5 60]); 
        ylabel('OEF'); title('OEF');
    end
    xlabel('Time (s)');
end

function export_peak_summary(data, time, peak1_range, peak2_range, base_range, stim_range, recovery_range)
    groups = {'WT', 'PS19'};
    vessels = {'a', 'v'};
    peak_data = {};
    header = {'Group', 'Vessel', 'MouseID', 'sO2_Peak1', 'sO2_Peak2', 'OEF_Peak1', 'OEF_Peak2', 'baseline', 'stim', 'recovery'};

    for g = 1:2
        group = groups{g};
        for v = 1:2
            vessel = vessels{v};
            sO2_mat = data.(group).sO2.(vessel);
            for i = 1:size(sO2_mat, 1)
                sO2_1 = mean(sO2_mat(i, peak1_range), 'omitnan');
                sO2_2 = mean(sO2_mat(i, peak2_range), 'omitnan');
                sO2_3 = mean(sO2_mat(i, base_range), 'omitnan');
                sO2_4 = mean(sO2_mat(i, stim_range), 'omitnan');
                sO2_5 = mean(sO2_mat(i, recovery_range), 'omitnan');

                if strcmp(vessel, 'a') && i <= size(data.(group).sO2.v, 1)
                    sv = data.(group).sO2.v(i, :);
                    sa = sO2_mat(i, :);
                    oef = (sa - sv) ./ sa;
                    oef1 = mean(oef(peak1_range), 'omitnan');
                    oef2 = mean(oef(peak2_range), 'omitnan');
                    oef3 = mean(oef(base_range), 'omitnan');
                    oef4 = mean(oef(stim_range), 'omitnan');
                    oef5 = mean(oef(recovery_range), 'omitnan');
                else
                    oef1 = NaN; oef2 = NaN; oef3 = NaN; oef4 = NaN; oef5 = NaN;
                end

                peak_data(end+1, :) = {group, vessel, i, sO2_1, sO2_2, oef1, oef2, oef3, oef4, oef5};
            end
        end
    end

    % Write to CSV
    cell2csv('sO2_OEF_peaks_by_mouse.csv', [header; peak_data]);
    disp('Exported sO2 and OEF peaks to sO2_OEF_peaks_by_mouse.csv');
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
    out = a;
    if ~cond, out = b; end
end

function export_oef_traces(data, time)
    groups = {'WT', 'PS19'};
    output = {};
    header = [{'Group', 'MouseID'}, arrayfun(@(t) sprintf('t%.1f', t), time, 'UniformOutput', false)];

    for g = 1:2
        group = groups{g};
        sa = data.(group).sO2.a;
        sv = data.(group).sO2.v;
        n = min(size(sa,1), size(sv,1));
        for i = 1:n
            oef = (sa(i,:) - sv(i,:)) ./ sa(i,:);
            output(end+1, :) = [{group, i}, num2cell(oef)];
        end
    end

    cell2csv('OEF_traces_by_mouse.csv', [header; output]);
    disp('Exported full OEF traces to OEF_traces_by_mouse.csv');
end