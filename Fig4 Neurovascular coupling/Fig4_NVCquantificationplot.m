%% plot pairing eventsFig4
% 2025 Sung Ji Ahn

clearvars
%close all

filename = 'couplingquantification.xlsx';   % change to your Excel filename
T = readtable(filename);

% Ensure all key columns exist
requiredCols = {'genotype','filenum','pair','puff','auc'};
if ~all(ismember(requiredCols, T.Properties.VariableNames))
    error('Missing one or more required columns: genotype, filenum, pair, puff, auc');
end

% Pairing ratio per file (and per puff condition)
% Paired count per (genotype, filenum, puff)
Gsum = varfun(@nansum, T, ...
    'InputVariables', 'pair', ...
    'GroupingVariables', {'genotype','filenum','puff'});

% Total events per (genotype, filenum, puff)
Gcnt = varfun(@numel, T, ...
    'InputVariables', 'pair', ...
    'GroupingVariables', {'genotype','filenum','puff'});

% Align tables (they are in the same group order)
paired_count = Gsum.nansum_pair;
total_count  = Gcnt.numel_pair;

pairing_ratio = paired_count ./ total_count;

PairTable = table( ...
    Gsum.genotype, Gsum.filenum, Gsum.puff, pairing_ratio, ...
    'VariableNames', {'genotype','filenum','puff','pairing_ratio'});

% Separate by condition
WT_air  = PairTable.pairing_ratio(strcmpi(PairTable.genotype,'WT')   & PairTable.puff==1);
PS_air  = PairTable.pairing_ratio(strcmpi(PairTable.genotype,'PS19') & PairTable.puff==1);
WT_spont  = PairTable.pairing_ratio(strcmpi(PairTable.genotype,'WT')   & PairTable.puff==0);
PS_spont  = PairTable.pairing_ratio(strcmpi(PairTable.genotype,'PS19') & PairTable.puff==0);

% Boxplots: pairing ratio 
figure('Position',[100 100 1050 480]);

subplot(1,2,1); hold on; box on; grid on;
x_air = categorical([repmat("WT",numel(WT_air),1); repmat("PS19",numel(PS_air),1)]);
y_air = [WT_air; PS_air];
boxchart(x_air, y_air);

subplot(1,2,2); hold on; box on; grid on;
x_sp = categorical([repmat("WT",numel(WT_spont),1); repmat("PS19",numel(PS_spont),1)]);
y_sp = [WT_spont; PS_spont];
boxchart(x_sp, y_sp);

sgtitle('Pairing ratio by genotype and condition');

% AUC for spontaneous paired events only (WT vs PS19) 
spont_paired = T(T.pair==1, :);   % rows that are paired & spontaneous
WT_auc = spont_paired.duration_s(strcmpi(spont_paired.genotype,'WT'));
PS_auc = spont_paired.duration_s(strcmpi(spont_paired.genotype,'PS19'));
%
figure('Position',[200 200 460 480]); hold on; box on; grid on;
x_auc = categorical([repmat("WT",numel(WT_auc),1); repmat("PS19",numel(PS_auc),1)]);
y_auc = [WT_auc; PS_auc];
boxchart(x_auc, y_auc);
ylabel('AUC (spontaneous paired ensembles)');
title('Spontaneous paired duration%: WT vs PS19');
%{
figure('Position',[200 200 460 480]);
hold on; box on; grid on;

% choose bin edges so both groups use the same bins
edges = linspace(min([WT_auc; PS_auc],[],'omitnan'), ...
                 max([WT_auc; PS_auc],[],'omitnan'), 20);

colWT   = [0.20 0.44 0.84];   % blue
colPS19 = [0.80 0.20 0.20];   % red

% WT histogram
h1 = histogram(WT_auc, edges, 'Normalization', 'probability', ...
               'FaceColor', colWT, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

% PS19 histogram
h2 = histogram(PS_auc, edges, 'Normalization', 'probability', ...
               'FaceColor', colPS19, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

xlabel('AUC (spontaneous paired ensembles)');
ylabel('Fraction of events');
title('Spontaneous paired AUC distribution: WT vs PS19');
legend({'WT','PS19'}, 'Location','best'); legend boxoff;

%
%% Spontaneous events: compare duration_s for paired vs unpaired within each genotype
spont = T(T.puff==0, :);  % spontaneous only

% Helper to extract durations by pairing within a genotype
get_durs = @(tbl, g) deal( ...
    tbl.duration_s(strcmpi(tbl.genotype,g) & tbl.pair==1), ... % paired
    tbl.duration_s(strcmpi(tbl.genotype,g) & tbl.pair==0));    % unpaired

[ps19_paired, ps19_unpaired] = get_durs(spont, 'PS19');
[wt_paired,   wt_unpaired]   = get_durs(spont, 'WT');

% ---- Boxplots with stats (rank-sum) ----
figure('Position',[100 100 1080 460])

% PS19
subplot(1,2,1); hold on; box on; grid on
x_ps = categorical([repmat("Paired",numel(ps19_paired),1); repmat("Unpaired",numel(ps19_unpaired),1)]);
y_ps = [ps19_paired; ps19_unpaired];
boxchart(x_ps, y_ps);
ylabel('Duration (s)');
title('PS19: spontaneous events');

% stats
if ~isempty(ps19_paired) && ~isempty(ps19_unpaired)
    p_ps19 = ranksum(ps19_paired, ps19_unpaired);  % Wilcoxon rank-sum
    text(0.5, max(y_ps)*1.02, sprintf('ranksum p=%.3g', p_ps19), 'Units','data','HorizontalAlignment','center')
else
    text(0.5, 0.5, 'Insufficient data', 'Units','normalized','HorizontalAlignment','center')
end

% WT
subplot(1,2,2); hold on; box on; grid on
x_wt = categorical([repmat("Paired",numel(wt_paired),1); repmat("Unpaired",numel(wt_unpaired),1)]);
y_wt = [wt_paired; wt_unpaired];
boxchart(x_wt, y_wt);
ylabel('Duration (s)');
title('WT: spontaneous events');

% stats
if ~isempty(wt_paired) && ~isempty(wt_unpaired)
    p_wt = ranksum(wt_paired, wt_unpaired);  % Wilcoxon rank-sum
    text(0.5, max(y_wt)*1.02, sprintf('ranksum p=%.3g', p_wt), 'Units','data','HorizontalAlignment','center')
else
    text(0.5, 0.5, 'Insufficient data', 'Units','normalized','HorizontalAlignment','center')
end

sgtitle('Spontaneous ensemble duration: Paired vs Unpaired (by genotype)')
%}