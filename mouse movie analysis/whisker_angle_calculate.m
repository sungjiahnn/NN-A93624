%% Whisker tracking (Radon-based) + event detection
% Sung Ji Ahn 2019
%  select mp4
%  draw ROI box around whisker
%  Radon transform over theta = -40:90 and variance-based angle estimate
%  low-pass filter + downsample
%  acceleration + thresholded whisking events
%  save figures + .mat results

clearvars; close all; clc;

%% BLOCK 0: Select video
disp('Select mouse .mp4 stack from the current directory'); disp(' ');
[mousemovie, filepath] = uigetfile('*.mp4', 'Select .mp4 file', 'MultiSelect', 'off');
if isequal(mousemovie, 0)
    disp('No file selected. Exiting.'); 
    return;
end

videoPath = fullfile(filepath, mousemovie);
v = VideoReader(videoPath);

% Safer frame count (VideoReader.NumFrames can be unreliable for some codecs)
try
    nFrames = v.NumFrames;
catch
    nFrames = floor(v.Duration * v.FrameRate);
end

%% BLOCK 1: ROI selection 
firstFrameRGB = readFrame(v);
firstFrameGray = sum(firstFrameRGB, 3);

figure(9); clf;
imagesc(firstFrameGray); axis image; colormap gray;
title('Draw a box around whisker, double-click to finalize');
thebox = imrect(gca);
wait(thebox); % blocks until user finalizes ROI

roiPos = thebox.getPosition(); % [x y w h]
roiPos = round(roiPos);        % integer pixel coords

% Clamp ROI to frame bounds
frameH = size(firstFrameGray, 1);
frameW = size(firstFrameGray, 2);
roiPos(1) = max(1, roiPos(1));
roiPos(2) = max(1, roiPos(2));
roiPos(3) = max(1, min(roiPos(3), frameW - roiPos(1)));
roiPos(4) = max(1, min(roiPos(4), frameH - roiPos(2)));

theta = -40:90; % preserve key element
whiskerAngle = NaN(1, nFrames);

%% BLOCK 1b: Radon tracking loop 
disp('Tracking whisker angle (Radon)...');
radonStart = tic;

v.CurrentTime = 0; % rewind to start
for n = 1:nFrames
    % Read next frame
    if hasFrame(v)
        frameRGB = readFrame(v);
    else
        break;
    end

    frameGray = sum(frameRGB, 3);

    % Crop ROI
    x1 = roiPos(1); y1 = roiPos(2); w = roiPos(3); h = roiPos(4);
    crop_slice = frameGray(y1:(y1+h), x1:(x1+w));

    % Edge emphasis (preserve logic)
    J = imgaussfilt(crop_slice, 0.3);
    edgeImg = imabsdiff(crop_slice, J);
    bcrop = imbinarize(edgeImg, 0.2);

    % Radon transform + variance selection (preserve logic)
    [R, ~] = radon(bcrop, theta);
    colVar = var(R, 0, 1);          % variance over rows for each theta
    ordVar = sort(colVar);
    threshIdx = round(numel(ordVar) * 0.9); % top 10% by variance (preserve)
    sieve = colVar > ordVar(threshIdx);

    angles = nonzeros(theta(:) .* sieve(:)); % corresponding theta values
    whiskerAngle(n) = mean(angles);

    % Carry forward previous value if NaN (preserve behavior)
    if isnan(whiskerAngle(n))
        if n > 1
            whiskerAngle(n) = whiskerAngle(n-1);
        else
            whiskerAngle(n) = 0;
        end
    end
end

radonTime = toc(radonStart);
disp(['Whisker tracking time was ' num2str(radonTime) ' seconds.']); disp(' ');

% Remove any trailing NaNs (e.g., if frame count estimate was high)
whiskerAngle = whiskerAngle(~isnan(whiskerAngle));

%% BLOCK 2: Filter + event detection
samplingRate    = 120; % Hz
downSampleRate  = 30;  % Hz
filterThreshold = 20;  % Hz
filterOrder     = 2;

[z, p, k] = butter(filterOrder, filterThreshold / (samplingRate / 2), 'low');
[sos, g]  = zp2sos(z, p, k);

% Filtered angle (demeaned) + downsample (preserve)
filteredWhiskerAngle = filtfilt(sos, g, whiskerAngle - mean(whiskerAngle));
resampledWhiskers    = resample(filteredWhiskerAngle, downSampleRate, samplingRate);

% Acceleration + threshold metric (preserve)
whiskerAcceleration = diff(whiskerAngle, 2);
whiskingThreshold   = abs(diff(whiskerAngle, 2)) * downSampleRate^2;

threshoffset = 20000; % change if needed (preserve)
threshLine   = ones(1, length(whiskingThreshold)) * threshoffset;
whiskerbinary = whiskingThreshold > threshoffset;

%% BLOCK 3: Plotting 
figure(77); clf;

ax1 = subplot(4,1,1);
plot((1:length(whiskerAngle)) / samplingRate, whiskerAngle, 'k');
title('Whisker Angle');
ylabel('Angle (degrees)'); xlabel('Time (seconds)');

ax2 = subplot(4,1,2);
plot((1:length(whiskerAcceleration)) / samplingRate, whiskerAcceleration, 'k');
title('Whisker Acceleration');
ylabel('Acceleration (degrees/sec^2)'); xlabel('Time (seconds)');

ax3 = subplot(4,1,3);
plot((1:length(whiskingThreshold)) / samplingRate, whiskingThreshold, 'k'); hold on;
hThr = plot((1:length(whiskingThreshold)) / samplingRate, threshLine, 'r', 'LineWidth', 2);
title('Setting Whisking Event Threshold');
ylabel('a.u.'); xlabel('Time (seconds)');
legend(hThr, {'Whisking Threshold'});

ax4 = subplot(4,1,4);
plot((1:length(whiskerbinary)) / samplingRate, whiskerbinary, 'k');
title('Whisking Event');
ylabel('off             on'); xlabel('Time (seconds)');

linkaxes([ax1 ax2 ax3 ax4], 'x');

%% BLOCK 4: Save outputs
% Remove ".mp4" safely
[~, baseFileName, ~] = fileparts(mousemovie);

saveas(figure(77), fullfile(filepath, [baseFileName, '_whisker.fig']));
saveas(figure(9),  fullfile(filepath, [baseFileName, '_whiskerbox.fig']));

% Save results (preserve intent, but make explicit what's saved)
save(fullfile(filepath, [baseFileName, '_Wresults.mat']), ...
    'mousemovie', 'videoPath', 'roiPos', 'theta', ...
    'whiskerAngle', 'filteredWhiskerAngle', 'resampledWhiskers', ...
    'whiskerAcceleration', 'whiskingThreshold', 'threshoffset', 'whiskerbinary', ...
    'samplingRate', 'downSampleRate', 'filterThreshold', 'filterOrder', 'radonTime');
