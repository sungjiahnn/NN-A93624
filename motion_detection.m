%% Leg motion tracking using ROI-based frame differencing
% Sung Ji Ahn 2019
%   mp4 selection
%   Interactive ROI selection
%   frame-to-frame displacement (abs difference)
%   low-pass filtering + downsampling
%   acceleration-based motion detection
%   optional robust thresholding and debouncing
%   figure + .mat saving

clearvars; close all; clc;

%% BLOCK 0: Select video 
disp('Select mouse .mp4 stack from the current directory'); disp(' ');
[mousemovie, filepath] = uigetfile('*.mp4', 'Select .mp4 file', 'MultiSelect', 'off');
if isequal(mousemovie,0)
    disp('No file selected. Exiting.');
    return;
end

videoPath = fullfile(filepath, mousemovie);
v = VideoReader(videoPath);

%% BLOCK 1: ROI selection 
v.CurrentTime = 0;
firstFrameRGB = readFrame(v);
firstFrameGray = sum(firstFrameRGB, 3);

figure(9); clf;
imagesc(firstFrameGray); axis image; colormap gray;
title('Draw a box around leg, double-click to finalize');

thebox = imrect(gca);
wait(thebox);
roiPos = round(thebox.getPosition());   % [x y w h]

% Clamp ROI to frame bounds
frameH = size(firstFrameGray,1);
frameW = size(firstFrameGray,2);
roiPos(1) = max(1, roiPos(1));
roiPos(2) = max(1, roiPos(2));
roiPos(3) = max(1, min(roiPos(3), frameW - roiPos(1)));
roiPos(4) = max(1, min(roiPos(4), frameH - roiPos(2)));

%% BLOCK 2: Motion tracking loop 
disp('Tracking leg motion (frame-to-frame displacement)...');

v.CurrentTime = 0;

% Preallocate conservatively
nFramesEst = floor(v.Duration * v.FrameRate);
displacement = NaN(1, nFramesEst);

% Read first frame as reference
refRGB  = readFrame(v);
refGray = sum(refRGB, 3);
x1 = roiPos(1); y1 = roiPos(2); w = roiPos(3); h = roiPos(4);
refCrop = refGray(y1:(y1+h), x1:(x1+w));

radonStart = tic;
idx = 1;

while hasFrame(v)
    idx = idx + 1;

    frameRGB  = readFrame(v);
    frameGray = sum(frameRGB, 3);
    curCrop   = frameGray(y1:(y1+h), x1:(x1+w));

    % Mean absolute frame-to-frame difference
    diffImg = imabsdiff(refCrop, curCrop);
    displacement(idx) = mean(diffImg(:));

    % Update reference
    refCrop = curCrop;
end

radonTime = toc(radonStart);
disp(['Movement tracking time was ' num2str(radonTime) ' seconds.']); disp(' ');

% Trim unused preallocation
displacement = displacement(~isnan(displacement));

%% BLOCK 3: Filtering + motion detection
samplingRate    = 120;   % Hz
downSampleRate  = 30;    % Hz
filterThreshold = 20;    % Hz
filterOrder     = 2;

[z, p, k] = butter(filterOrder, filterThreshold / (samplingRate / 2), 'low');
[sos, g]  = zp2sos(z, p, k);

filteredDisplacement = filtfilt(sos, g, displacement - mean(displacement));
resampledDisplacement = resample(filteredDisplacement, downSampleRate, samplingRate);

% Acceleration-based metric
displacementAcceleration = diff(displacement, 2);
displacementThreshold = abs(diff(displacement, 2)) * downSampleRate^2;

%% BLOCK 4: Thresholding 
% Fixed threshold (kept for consistency with older analysis)
fixedThresh = 20000;
threshLine  = ones(1, length(displacementThreshold)) * fixedThresh;
displacementBinary = displacementThreshold > fixedThresh;

% Robust threshold (recommended)
medDT = median(displacementThreshold, 'omitnan');
madDT = mad(displacementThreshold, 1);
kRobust = 7;  % tweak if needed
robustThresh = medDT + kRobust * madDT;

motion_bin = displacementThreshold > robustThresh;

% Debounce / clean motion events
Fs = samplingRate;
min_bout_s = 0.10;   % 100 ms
min_gap_s  = 1.0;    % 1 s

min_bout_n = max(1, round(min_bout_s * Fs));
min_gap_n  = max(1, round(min_gap_s  * Fs));

motion_bin = imclose(motion_bin, true(1, min_gap_n));
motion_bin = bwareaopen(motion_bin, min_bout_n);

percent_moved = 100 * sum(motion_bin) / numel(motion_bin);
fprintf('Percent time moved: %.2f%%\n', percent_moved);

%% BLOCK 5: Figure generation
figure(77); clf;

ax1 = subplot(4,1,1);
plot((1:length(displacement)) / samplingRate, displacement, 'k');
title('Frame-to-frame displacement');
ylabel('a.u.'); xlabel('Time (seconds)');

ax2 = subplot(4,1,2);
plot((1:length(displacementAcceleration)) / samplingRate, displacementAcceleration, 'k');
title('Displacement acceleration');
ylabel('a.u.'); xlabel('Time (seconds)');

ax3 = subplot(4,1,3);
plot((1:length(displacementThreshold)) / samplingRate, displacementThreshold, 'k'); hold on;
plot((1:length(displacementThreshold)) / samplingRate, threshLine, 'r', 'LineWidth', 2);
title('Motion threshold');
ylabel('a.u.'); xlabel('Time (seconds)');
legend({'Signal','Fixed threshold'});

ax4 = subplot(4,1,4);
plot((1:length(displacementBinary)) / samplingRate, displacementBinary, 'k');
title('Motion event');
ylabel('off             on'); xlabel('Time (seconds)');

linkaxes([ax1 ax2 ax3 ax4], 'x');

%% BLOCK 6: Save outputs 
[~, baseFileName, ~] = fileparts(mousemovie);

saveas(figure(77), fullfile(filepath, [baseFileName, '_motion.fig']));
saveas(figure(9),  fullfile(filepath, [baseFileName, '_motionbox.fig']));

save(fullfile(filepath, [baseFileName, '_Mresults.mat']), ...
    'videoPath', 'roiPos', 'displacement', 'filteredDisplacement', ...
    'resampledDisplacement', 'displacementAcceleration', ...
    'displacementThreshold', 'displacementBinary', 'motion_bin', ...
    'percent_moved', 'samplingRate', 'downSampleRate', 'radonTime');
