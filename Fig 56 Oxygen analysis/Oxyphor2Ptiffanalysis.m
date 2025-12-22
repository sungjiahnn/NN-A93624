% Sung Ji Ahn 2024

close all
clear all

% Ask user to select a single .tif movie (stack)
% channel 1 (Oxyphor 2P channel) tif, separated using imageJ macro
disp('Select a single .tif stack from the current directory'); disp(' ')
filename = uigetfile('*.tif', 'Multiselect', 'off');

% Get info about the TIF (number of frames, sizes, etc.)
movieInfo = imfinfo(filename);
disp(['Loading: ' filename '...']); disp(' ')

% Read first frame (or page) to get dimensions
[X,a] = imread(filename);

% X is a 2D image; R = number of rows (e.g. positions along vessel),
% T = number of columns (e.g. time * repeats concatenated)
[R,T] = size(X);

% series = number of frames (time points) in the tif stack
[series,a] = size(movieInfo);

% Each repeat is 150 columns wide; infer number of repeats in the image
repeat = T/150;

% Preallocate: averagetrace has size [time(150) x rows x frames]
% For each tif frame (series), for each row, we compute the mean over repeats
% and subtract the overall median of that row's matrix.
averagetrace = NaN(150,R,series);    

for l = 0:series-1
    % Read frame (slice) l+1 of the tif stack
    [X,map] = imread(filename, l+1);
    
    for i = 1:R
        % For each row i, reshape columns into [150 x repeat] matrix
        % matrix(k, j) = intensity at time k, repeat j
        matrix = NaN(150, repeat);
        
        for j = 0:repeat-1
            for k = 1:150
                matrix(k, j+1) = X(i, j*150 + k);
            end
        end
        
        % Average across repeats at each time point, then subtract median
        % (baseline correction)
        averagetrace(:, i, l+1) = mean(matrix,2) - median(matrix,'all');
    end
end

%{
% extracting a subset of time points (8:147) to analysistrace
analysistrace = NaN(R,140,series);
x = [0:2:278];

for k = 1:series
    for j = 1:R
        for i = 8:147
            analysistrace(j, i-7, k) = averagetrace(i, j, k);
        end
    end
end
%}

%% Time-averaging across series (frames) and extraction of time window

time = 2;               % Number of frames to average together (temporal bin size)
newT = ceil(series/time);  % New "time length" after averaging

% x corresponds to time points used for fitting (140 points, step 2 ms?)
x = [0:2:278];

temp          = NaN(150, time);         % temporary buffer for averaging
analysistrace = NaN(R,140,newT);        % final data to fit: [rows x 140 time pts x newT]

for i = 1:R
    for j = 0:newT-1
        % Average over "time" consecutive frames (if they exist)
        for k = 1:time
            if (j*time + k) <= series
                temp(:,k) = averagetrace(:, i, j*time + k);
            end
        end
        
        % Mean across the 'time' frames at each of 150 points
        temp1 = nanmean(temp, 2);
        
        % Keep only time points 8 to 147 (140 points total) for fitting
        for k = 8:147
            analysistrace(i, k-7, j+1) = temp1(k,1);
        end
    end
end

%% Single exponential curve fit to each row & time segment

% fitABC(i, :, j) stores parameters [a, b, c] of fit a*exp(b*x) + c
fitABC = NaN(R, 3, newT);

% gofcoefficients(i, :, j) = [SSE, adjusted R^2, RMSE]
gofcoefficients = NaN(R, 3, newT);

for i = 1:R
    for j = 1:newT
        
        % y is the decay curve for one row (i) and time bin (j)
        y = analysistrace(i,:,j);
        
        % Prepare xData and yData for the Curve Fitting Toolbox
        [xData, yData] = prepareCurveData( x, y );

        % Define exponential fit model: y = a * exp(b*x) + c
        ft = fittype( 'a*exp(b*x)+c' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';

        % Initial guess for parameters [a, b, c]
        % Here: a=100, b=0, c=-5 (e.g., example tuned for a particular vessel)
        opts.StartPoint = [100 0 -5];

        % Fit model to data
        [fitresult, gof] = fit( xData, yData, ft, opts );

        % If fit is poor (RMSE > 2), randomly adjust starting point
        % and refit until RMSE < 2
        while gof.rmse > 2
            adj = 100 .* rand(1,1);
            opts.StartPoint = [adj+20 0 -5 ];
            [fitresult, gof] = fit( xData, yData, ft, opts );
        end

        % Store fit parameters and goodness-of-fit metrics
        fitABC(i,:,j)          = coeffvalues(fitresult);
        gofcoefficients(i,1,j) = gof.sse;
        gofcoefficients(i,2,j) = gof.adjrsquare;
        gofcoefficients(i,3,j) = gof.rmse;
        
        % Simple progress printouts
        i
        j
    end
end

%% Stern–Volmer calibration: convert lifetime -> pO2

% pO2results(i, j) = oxygen value for row i and time bin j
pO2results = NaN(R, newT);

for i = 1:R
    for j = 1:newT
        % For an exponential y = a*exp(b*x)+c, b is related to decay constant
        % Here b ~ -1/tau, so lifetime (tau) ~ -1/b (scale unit conversion)
        inverselambda = 1/fitABC(i,2,j);
        lifetime      = -1 * inverselambda * 0.000001;  % convert to seconds
        
        % Empirical/experimental calibration curve: pO2 = f(tau)
        pO2results(i,j) = ...
            560.89597  * exp(-1 * lifetime / 0.00000439093) + ...
            158.72502  * exp(-1 * lifetime / 0.0000181168)  - ...
            18.83227;
    end
end

%% Visualization: pO2 map and time courses

figure(77)
% Alternative visualization (commented out) using a grid of x,y coordinates:
% scatter(grid.xcoo(:),grid.ycoo(:),50,averagetrace20(5,:),'filled')
% scatter(grid.xcoo(:),grid.ycoo(:),100,pO2results(:),'s','filled')
% set(gca,'Ydir', 'reverse')

% Show pO2results as an image (rows x time bins)
imagesc(pO2results)
colorbar; 
% caxis([15 200]);    % optional: fix color scale
% meanpO2 = mean(pO2results);
% minpO2  = min(pO2results);
% maxpO2  = max(pO2results);

savefig(figure(77), 'experiment_2.fig')

% Another figure: pO2 traces over time for each row
figure(88)
plot(pO2results','DisplayName','pO2results','LineWidth',2)
savefig(figure(88), 'experiment_2line.fig')

% Save all variables to experiment_2.mat
save('experiment_2')

% clean up outliers
%}
