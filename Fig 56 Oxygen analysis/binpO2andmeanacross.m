function plot_pO2_gradient_WT_vs_PS19_perFileAveraged()
    % -------- Parameters --------
    bin_width = 11;     % µm
    max_dist  = 160;    % µm
    bin_edges   = 0:bin_width:max_dist;
    bin_centers = bin_edges(1:end-1) + bin_width/2;

    % ---- Select WT files ----
    disp("Select WT (wild-type) CSV files");
    [wt_files, wt_path] = uigetfile('*.csv', 'Select WT CSV files', 'MultiSelect', 'on');
    if isequal(wt_files,0), disp('No WT files selected'); return; end
    if ischar(wt_files), wt_files = {wt_files}; end

    % ---- Select PS19 files ----
    disp("Select PS19 CSV files");
    [ps_files, ps_path] = uigetfile('*.csv', 'Select PS19 CSV files', 'MultiSelect', 'on');
    if isequal(ps_files,0), disp('No PS19 files selected'); return; end
    if ischar(ps_files), ps_files = {ps_files}; end

    % ---- Per-file binning (each file -> 1 value per bin) ----
    WT_bins  = per_file_binned_means(wt_files, wt_path, bin_edges);   % nWT x nBins
    PS_bins  = per_file_binned_means(ps_files, ps_path, bin_edges);   % nPS x nBins

    % ---- Average across files (unweighted) ----
    [mean_wt, sem_wt, n_wt] = across_files_mean_sem(WT_bins);
    [mean_ps, sem_ps, n_ps] = across_files_mean_sem(PS_bins);

    % ---- Plot (shaded SEM) ----
    figure; hold on; box on; grid on;
    colWT   = [0 0.45 0.74];
    colPS19 = [0.85 0.33 0.1];

    % WT band + mean
    fill([bin_centers, fliplr(bin_centers)], ...
         [mean_wt+sem_wt, fliplr(mean_wt-sem_wt)], ...
         colWT, 'FaceAlpha',0.25, 'EdgeColor','none');
    plot(bin_centers, mean_wt, '-', 'Color', colWT, 'LineWidth', 2, 'DisplayName', sprintf('WT (n=%d)', n_wt));

    % PS19 band + mean
    fill([bin_centers, fliplr(bin_centers)], ...
         [mean_ps+sem_ps, fliplr(mean_ps-sem_ps)], ...
         colPS19, 'FaceAlpha',0.25, 'EdgeColor','none');
    plot(bin_centers, mean_ps, '-', 'Color', colPS19, 'LineWidth', 2, 'DisplayName', sprintf('PS19 (n=%d)', n_ps));

    xlabel('Distance from arteriole (µm)');
    xticks(0:50:150)
    ylabel('Tissue pO₂ (mmHg)');
    title('Average pO₂ gradient (per-file bin means, mean \pm SEM across files)');
    legend('Location','northeast'); legend boxoff;
    xlim([0 max_dist]);
    % ylim can be set as you wish, or auto:
    % ylim([10 65]);
end

% ===================== Helpers =====================

function M = per_file_binned_means(files, folder, bin_edges)
    % Returns an nFiles x nBins matrix of per-file bin means (NaN if a bin has no points in that file)
    nBins = numel(bin_edges) - 1;
    M = nan(numel(files), nBins);

    for i = 1:numel(files)
        T = readtable(fullfile(folder, files{i}));

        % Center: first row
        x0 = T{1,1}; y0 = T{1,2};
        % Data from row 2..end: columns [x,y,pO2]
        x   = T{2:end,1};
        y   = T{2:end,2};
        pO2 = T{2:end,3};

        valid = ~isnan(x) & ~isnan(y) & ~isnan(pO2);
        x = x(valid); y = y(valid); pO2 = pO2(valid);

        r = sqrt((x - x0).^2 + (y - y0).^2);

        % Assign distances to bins
        binID = discretize(r, bin_edges, 'IncludedEdge','left'); % [edge(i), edge(i+1))
        % Compute mean per bin for THIS file
        % Use accumarray; default fill NaN if no values
        for b = 1:nBins
            vals = pO2(binID == b);
            if ~isempty(vals)
                M(i,b) = mean(vals, 'omitnan');
            else
                M(i,b) = NaN;
            end
        end
    end
end

function [mu, se, nFiles] = across_files_mean_sem(M)
    % M: nFiles x nBins, each row = one file’s bin means
    nFiles = size(M,1);
    mu = mean(M, 1, 'omitnan');
    % SEM across files (each file contributes equally)
    n  = sum(~isnan(M), 1);
    se = std(M, 0, 1, 'omitnan') ./ max(sqrt(n),1);
end
