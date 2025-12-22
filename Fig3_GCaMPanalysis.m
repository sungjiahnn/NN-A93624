% Sung Ji Ahn 2023

clear all 
% input calcium signals per manually drawn roi on CalciumDX
load('matlab.mat','tr','ncells','len','tmpres')

% tr      : ncells x len matrix of raw fluorescence time series
% ncells  : number of cells
% len     : number of time points
% tmpres  : temporal resolution (seconds per frame)
%close all

%% 1) Compute dF/F using a sliding baseline (f0 = mean of lower intensities in window)

% Length of sliding window in frames (here corresponding to 15 seconds)
% comment notes: 15 for 5-sec data; 50 for 30-sec data, depending on tmpres

windowL=round(15/tmpres); 
deltaF=NaN(ncells,len); % dF/F for each cell x time
f0matrix=NaN(ncells,len); % baseline f0 used at each time
dettr=NaN(1,len); % detrended trace (single neuron at a time)

activity=zeros(ncells,len); % binary matrix: 1 = active, 0 = not active
% Low-pass filter design for smoothing dF/F traces (2nd-order Butterworth)
d1 = designfilt('lowpassiir','FilterOrder',2,'HalfPowerFrequency',0.2,'DesignMethod','butter');%zero-phase filter
activityF=zeros(ncells,len); % dF/F masked by activity (0 when inactive)
astd=zeros(1,ncells); % (later used as global STD, but preallocated as vector)

for i=1:ncells
    % Detrend the trace over time (remove slow drift) and anchor to first value
    dettr(1,:)=detrend(tr(i,:))+tr(i,1); 

    % Compute sliding baseline f0 and dF/F
    for j=1:len-1
        % Build the window of past data points up to current time j
        if j < (windowL +1)   % smaller window for early time points
            window=NaN(1,j);
            for k=1:j
                window(1,k)=dettr(1,k);
            end
        else    % full running window of length windowL
            window=NaN(1,windowL);
            for k=1:windowL
                window(1,k)=dettr(1,j-k);
            end
        end

        % Use only the lower-intensity values in the window to define baseline
        % 1) Find 25th percentile of values in the window
        w25=prctile(window,25);     % 25 for 5sec 10 for 30sec
        % 2) Keep only values <= 25th percentile, set the rest to NaN
        window(window>w25)=NaN;
        % 3) Baseline f0 = mean of these low values
        f0=nanmean(window);
        
        % dF/F = (signal - f0) / f0
        deltaF(i,j)=(dettr(1,j)-f0)/f0;
        f0matrix(i,j)=f0;
    end

    % For the last time point, copy the previous dF/F value
    deltaF(i,len)=deltaF(i,len-1);
end

% Global standard deviation of dF/F across all cells and time points
astd=std2(deltaF);% std of whole matrix

%% 2) Detect activity events per cell using thresholds on filtered dF/F
for i=1:ncells
    % Smooth dF/F trace with zero-phase low-pass filter
    y = filtfilt(d1,deltaF(i,:)); %zero phase filter
    flag  = 0;              % activity flag: 0 = inactive, 1 = currently in an event
    ydiff = diff(y,1);      % first derivative of filtered trace
    A     = std(ydiff);     % STD of derivative (used to detect onset)
    
    % B marks large positive changes in derivative (onsets)
    B = zeros(1, length(ydiff));
    B(ydiff > A/2) = 1;
    
    % Baseline level for this cell: 15th percentile of filtered trace
    w25 = prctile(y,15);    % effectively a low baseline
     for j = 2:len-1
        % ---- Event onset detection ----
        % If not currently in an event (flag == 0) and signal rises above
        % (baseline + 2 * global STD), start an event.
        if (flag == 0) && (y(1,j) > astd*2 + w25)   % 2 SD above baseline
            flag        = 1;
            activity(i,j) = 1;   % mark this time point as active
            idx         = j;
            
            % Walk backward to include earlier samples leading up to onset
            % while no new large jump in derivative is detected.
            while B(1,idx) - B(1,idx-1) < 1
                idx = idx - 1;
                activity(i,idx) = 1;
                if idx == 1
                    break
                end
            end
        
        % ---- Event continuation ----
        % If already in an event (flag == 1) and signal stays above
        % (baseline + 0.5 * global STD), continue marking as active.
        elseif (flag == 1) && (y(1,j) > astd*0.5 + w25)  % 0.5 SD above
            activity(i,j) = 1;
        
        % ---- Event end ----
        % If signal drops below (baseline + 0.5 * global STD), end event.
        elseif (y(1,j) < astd*0.5 + w25) && (flag == 1)
            flag = 0;
        end
    end
end

% dF/F masked by activity: non-active times set to 0
activityF = deltaF .* activity;

%% 3) Remove cells that are never active

for i = 1:ncells
    if sum(activityF(i,:)) == 0
        % Mark cell as invalid by inserting NaN (later we remove rows with NaN)
        activityF(i,1) = NaN;
        activity(i,1)  = NaN;
    end
end

% Remove all rows that contain any NaN (cells with no activity)
activityF = activityF(all(~isnan(activityF),2),:);
activity  = activity(all(~isnan(activity),2),:);

% Recompute ncells and len after removing inactive cells
[ncells,len] = size(activityF);

% Percent of cells active at each time point
pcntactive(1,:) = sum(activity) / ncells * 100;   % 0–100%

%% 4) Compute pairwise Spearman's correlation across cells

% Sort cells by their activity at time index 161 (arbitrary alignment point)
sorinf = sortrows(activityF, 161);

% Spearman correlation (rho) between cells (across time)
[rho,pval] = corr(sorinf','type','Spearman');

figure(5)
imagesc(rho)
title('correlation matrix (Spearmans''s rho)', 'FontSize', 10);
colormap('hot'); 
caxis([-0.2 0.8])
colorbar 

% savefig(figure(5),'correlationre.fig')

%% 5) Generate null distribution for ensemble activity by time-shuffling

shuffpcntactive = NaN(1000,len);  % percent-active traces from shuffled data

for i = 1:1000
    % Random circular time offset for each cell
    offset   = round(rand(ncells,1) * len);
    shuffleM = NaN(ncells,len);   % shuffled activity matrix
    
    for ii = 1:ncells
        for j = 1:len
            k = offset(ii,1);
            % Circularly shift activity(ii,:) by k
            if j + k > len
                shuffleM(ii,j + k - len) = activity(ii,j);
            else
                shuffleM(ii,j + k) = activity(ii,j);
            end
        end
    end
    
    % Percent of cells active at each time point in shuffled data
    shuffpcntactive(i,:) = sum(shuffleM) / ncells * 100;
end

% Time vector in seconds
matrixx = tmpres:tmpres:tmpres*len;

% Mean and variability of shuffled percent-active traces
shuffmean = mean(shuffpcntactive);                     % across 1000 shuffles
shuffSEM  = std(shuffpcntactive) / sqrt(ncells);       % SEM (uses ncells as N)

% 99.9% confidence intervals from t-distribution
CI99 = tinv([0.001 0.999], ncells-1);                  
yCI99 = bsxfun(@times, shuffSEM, CI99(:));             % upper/lower envelopes

% Ensemble threshold: baseline (shuffled mean) + lower CI limit
% (Conceptually: anything above this is unlikely under time-shuffled null)
thresh = shuffmean + CI99(1,1);  

%% 6) Define significant ensembles based on threshold

signensemble       = zeros(1,len);      % binary: 1 = ensemble at that time
cellproporensemble = pcntactive;        % percent of cells in ensemble at each time

% Mark time points where observed percent-active exceeds threshold
signensemble(pcntactive > thresh)  = 1;
cellproporensemble(pcntactive < thresh) = NaN;   % keep only ensemble frames

% Ensemble frequency over the entire recording (in % of time points)
ensemblefreq = sum(signensemble) / len * 100;

% Mean size of ensembles (average % of cells active during ensemble frames)
proportionensemble = nanmean(cellproporensemble);

% Vector of ensemble sizes for all ensemble frames
ensembledata = cellproporensemble(~isnan(cellproporensemble))';

%% 7) Plot raster, activity statistics, and ensembles

figure(3)
tiledlayout(11,1)

% (1) dF/F during active periods
nexttile([3 1])
imagesc(activityF)
ylabel('cells')
colormap('jet'); 
colorbar
caxis([0 1.5])

% (2) Binary activity raster
nexttile([3 1])
imagesc(activity)
ylabel('cells (binary)')

% (3) Percent of cells active vs shuffle-based null
nexttile([3 1])
hold on
plot(matrixx, pcntactive)             % real data
plot(matrixx, shuffmean)              % shuffle mean
plot(matrixx, yCI99 + shuffmean)      % 99.9% CI around shuffle mean
ylabel('% cells active')
xticks([])
hold off
grid

% (4) Ensemble time points (1 = ensemble, 0 = none)
nexttile
imagesc(signensemble)
ylabel('ensembles (99.9% CI)')
xticks([])
yticks([])

% whisker puff trace (air-puff timing) 
%{
nexttile
imagesc(whiskerpuff)
ylabel('air')
yticks([])
%}

% Save figure and variables
savefig(figure(3),'ensemblere.fig')
save('ensemblere')