%% Surface vessel diameter extraction from 2P TIFF projection movie
% Sung Ji Ahn adopted and customized. 
%
% Original code credits:
%   Kevin L. Turner (Penn State) + Patrick J. Drew (DrewLab)
%   https://github.com/KL-Turner   https://github.com/DrewLab

clearvars; close all; clc;

%% SETTINGS (edit here)
MScanData = struct();
MScanData.notes.frameRate  = 1.081;        % Hz
micronsPerPixel = 0.9944;             % um/pixel calibration

doLowpassPlot = true;                  % show filtered trace subplot
lowpassHz     = 1;                     % lowpass cutoff for plot (Hz)
saveName      = 'Vessel1.mat';         % output mat file

%% Select TIFF stack 
disp('Select a single .tif stack from the current directory'); disp(' ');
[tifStack, tifPath] = uigetfile('*.tif', 'Select TIFF stack', 'MultiSelect', 'off');
if isequal(tifStack, 0)
    disp('No file selected. Exiting.');
    return;
end
tifFile = fullfile(tifPath, tifStack);

movieInfo = imfinfo(tifFile);
nFrames   = numel(movieInfo);
disp(['Loading: ' tifStack '...']); disp(' ');

%% Basic metadata
MScanData.notes.header.fileName       = movieInfo(1).Filename;
MScanData.notes.header.frameWidth     = movieInfo(1).Width;
MScanData.notes.header.frameHeight    = movieInfo(1).Height;
MScanData.notes.header.numberOfFrames = nFrames;

MScanData.notes.xSize = movieInfo(1).Width;
MScanData.notes.ySize = movieInfo(1).Height;

MScanData.notes.startframe = 1;
MScanData.notes.endframe   = nFrames;

MScanData.notes.micronsPerPixel = micronsPerPixel;
MScanData.notes.xFactor         = micronsPerPixel;
MScanData.notes.header.timePerLine = 1/(MScanData.notes.frameRate * MScanData.notes.ySize);

%% Draw vessel ROI + axis line 
firstFrame = imread(tifFile, 'TIFF', 'Index', 1);

figure('Name','Select vessel ROI'); clf;
imagesc(double(firstFrame)); axis image; colormap(gray);
xlabel('pixels'); ylabel('pixels');
title('Draw ROI polygon around vessel (double-click to finish)');

roiPoly = impoly(gca, [1 1; 1 20; 20 20; 20 1]); 
wait(roiPoly);  % user double-clicks to finish
roiAPI = iptgetapi(roiPoly);
MScanData.notes.vesselROI.boxPosition.xy = roiAPI.getPosition();

title('Draw a line along the vessel center-axis (double-click to finish)');
axisLine = imline(gca, round(MScanData.notes.xSize*[.25 .75]), round(MScanData.notes.ySize*[.25 .75])); %#ok<IMLINE>
wait(axisLine); % user double-clicks to finish
lineAPI = iptgetapi(axisLine);
MScanData.notes.vesselROI.vesselLine.position.xy = lineAPI.getPosition();

%% Radon projections (with registration) 
disp('Analyzing vessel projections from ROI...'); disp(' ');
MScanData = GetDiameterFromMovie(MScanData, tifFile);

%% Diameter from FWHM 
try
    MScanData = FWHM_MovieProjection(MScanData, [MScanData.notes.startframe MScanData.notes.endframe]);
catch ME
    warning('FWHM calculation failed: %s', ME.message);
end

%% Remove motion artifacts
rawThresh = 0.5;

MScanData.data.vesselDiameter = RemoveMotion( ...
    MScanData.data.tempVesselDiameter, ...
    MScanData.notes.vesselROI.modalFixedDiameter, ...
    2, rawThresh);

%% Plot
Fs = MScanData.notes.frameRate;           % sampling rate (Hz)
nyquistHz = Fs / 2;

t = (1:MScanData.notes.endframe) / Fs;

figure('Name','Vessel diameter'); clf;

if doLowpassPlot
    % Ensure cutoff is valid for Butterworth filter
    cutoffHz = min(lowpassHz, 0.95 * nyquistHz);   % stay below Nyquist
    Wn = cutoffHz / nyquistHz;

    [B, A] = butter(4, Wn, 'low');
    filtVesselDiam = filtfilt(B, A, MScanData.data.vesselDiameter);

    subplot(2,1,1);
    plot(t, MScanData.data.vesselDiameter, 'k');
    title('Vessel diameter over time');
    xlabel('Time (sec)');
    ylabel('Diameter (\mum)');
    axis tight;

    subplot(2,1,2);
    plot(t, filtVesselDiam, 'k');
    title(sprintf('Lowpass filtered (%.2f Hz cutoff)', cutoffHz));
    xlabel('Time (sec)');
    ylabel('Diameter (\mum)');
    axis tight;
else
    plot(t, MScanData.data.vesselDiameter, 'k');
    title('Vessel diameter over time');
    xlabel('Time (sec)');
    ylabel('Diameter (\mum)');
    axis tight;
end

%% Save 
save(saveName, 'MScanData');
disp('Diameter calculation complete.');
disp(['Saved: ' saveName]); disp(' ');


%% Local functions
function MScanData = GetDiameterFromMovie(MScanData, fileID)
refFrame = imread(fileID, 'TIFF', 'Index', 1);
fftRef   = fft2(double(refFrame));

X = repmat(1:MScanData.notes.xSize, MScanData.notes.ySize, 1);
Y = repmat((1:MScanData.notes.ySize)', 1, MScanData.notes.xSize);

xy = MScanData.notes.vesselROI.vesselLine.position.xy;
MScanData.notes.vesselROI.projectionAngle = atand( diff(xy(:,1)) / diff(xy(:,2)) );

for f = MScanData.notes.startframe:MScanData.notes.endframe
    rawFrame = imread(fileID, 'TIFF', 'Index', f);
    fftFrame = fft2(double(rawFrame));
    [MScanData.notes.pixelShift(:, f), ~] = DftRegistration(fftRef, fftFrame, 1);

    mask = inpolygon( ...
        X + MScanData.notes.pixelShift(3, f), ...
        Y + MScanData.notes.pixelShift(4, f), ...
        MScanData.notes.vesselROI.boxPosition.xy(:,1), ...
        MScanData.notes.vesselROI.boxPosition.xy(:,2));

    boundedFrame = rawFrame .* uint16(mask);
    MScanData.notes.vesselROI.projection(f, :) = radon(boundedFrame, MScanData.notes.vesselROI.projectionAngle); %#ok<AGROW>
end
end

function MScanData = FWHM_MovieProjection(MScanData, theFrames)
for f = min(theFrames):max(theFrames)
    proj = MScanData.notes.vesselROI.projection(f, :);
    proj = medfilt1(proj, 5);
    MScanData.data.rawVesselDiameter(f) = CalcFWHM(proj);
end

MScanData.data.tempVesselDiameter = MScanData.data.rawVesselDiameter * MScanData.notes.xFactor;

[counts, bins] = hist(MScanData.data.tempVesselDiameter, 0:0.25:100);
[~, idx] = max(counts);
MScanData.notes.vesselROI.modalFixedDiameter = bins(idx);
end

function width = CalcFWHM(data, smoothing, threshold)
data = double(data(:));

if nargin < 2, smoothing = 1; end
if smoothing > 1
    data = conv2(data, rectwin(smoothing)./smoothing, 'valid');
end

if nargin < 3
    offset    = min(data);
    threshold = max(data - offset)/2 + offset;
end

aboveI = find(data > threshold);
if isempty(aboveI)
    width = 0; return;
end

firstI = aboveI(1);
lastI  = aboveI(end);

if (firstI-1 < 1) || (lastI+1 > numel(data))
    width = 0; return;
end

point1offset = (threshold - data(firstI-1)) / (data(firstI) - data(firstI-1));
point2offset = (threshold - data(lastI))   / (data(lastI+1) - data(lastI));

point1 = (firstI-1) + point1offset;
point2 = lastI + point2offset;

width = point2 - point1;
end

function newVesselDiameter = RemoveMotion(vesselDiameter, baseline, diffThresh, rawThresh)
jumpIdx = find(diff(vesselDiameter) > diffThresh);
devIdx  = find(abs((vesselDiameter - baseline)/baseline) > rawThresh);
badIdx  = union(jumpIdx + 1, devIdx);

goodIdx = 1:numel(vesselDiameter);
goodIdx(badIdx) = [];

if isempty(badIdx)
    newVesselDiameter = vesselDiameter;
    return;
end

if goodIdx(1) ~= 1
    goodIdx = [1:goodIdx(1)-1, goodIdx];
end

newVesselDiameter = zeros(size(vesselDiameter));
count = 1;

for a = 1:numel(goodIdx)-1
    step = goodIdx(a+1) - goodIdx(a);

    if step == 1
        newVesselDiameter(count) = vesselDiameter(count);
    else
        newVesselDiameter(count) = vesselDiameter(count);
        newVesselDiameter(count+1:count+step-1) = mean([vesselDiameter(goodIdx(a)), vesselDiameter(goodIdx(a+1))]);
    end
    count = count + step;
end

newVesselDiameter(count) = vesselDiameter(goodIdx(end));
if badIdx(end) == numel(vesselDiameter)
    newVesselDiameter(goodIdx(end)+1:end) = vesselDiameter(goodIdx(end));
end
end

function [output, Greg] = DftRegistration(buf1ft, buf2ft, usfac)
% Efficient subpixel image registration by cross-correlation (Guizar-Sicairos et al., 2008)

if exist('usfac','var')~=1, usfac=1; end 

if usfac == 0
    CCmax = sum(sum(buf1ft.*conj(buf2ft)));
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2);
    err = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero);
    err = sqrt(abs(err));
    diffphase = atan2(imag(CCmax),real(CCmax));
    output = [err,diffphase];

elseif usfac == 1
    [m,n]=size(buf1ft);
    CC = ifft2(buf1ft.*conj(buf2ft));
    [max1,loc1] = max(CC);
    [~,loc2] = max(max1);
    rloc=loc1(loc2); cloc=loc2;
    CCmax=CC(rloc,cloc);
    rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n);
    err = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1));
    err = sqrt(abs(err));
    diffphase=atan2(imag(CCmax),real(CCmax));

    md2 = fix(m/2);
    nd2 = fix(n/2);

    if rloc > md2, row_shift = rloc - m - 1; else, row_shift = rloc - 1; end
    if cloc > nd2, col_shift = cloc - n - 1; else, col_shift = cloc - 1; end

    output=[err,diffphase,row_shift,col_shift];

else
    [m,n]=size(buf1ft);
    mlarge=m*2; nlarge=n*2;
    CC=zeros(mlarge,nlarge);
    CC(m+1-fix(m/2):m+1+fix((m-1)/2), n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
        fftshift(buf1ft).*conj(fftshift(buf2ft));

    CC = ifft2(ifftshift(CC));
    [max1,loc1] = max(CC);
    [~,loc2] = max(max1);
    rloc=loc1(loc2); cloc=loc2;
    CCmax=CC(rloc,cloc);

    [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
    if rloc > md2, row_shift = rloc - m - 1; else, row_shift = rloc - 1; end
    if cloc > nd2, col_shift = cloc - n - 1; else, col_shift = cloc - 1; end
    row_shift=row_shift/2; col_shift=col_shift/2;

    if usfac > 2
        row_shift = round(row_shift*usfac)/usfac;
        col_shift = round(col_shift*usfac)/usfac;
        dftshift = fix(ceil(usfac*1.5)/2);

        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);

        [max1,loc1] = max(CC);
        [~,loc2] = max(max1);
        rloc = loc1(loc2); cloc = loc2;
        CCmax = CC(rloc,cloc);

        rg00 = dftups(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
        rf00 = dftups(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);

        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;
    else
        rg00 = sum(sum(buf1ft.*conj(buf1ft)))/m/n;
        rf00 = sum(sum(buf2ft.*conj(buf2ft)))/m/n;
    end

    err = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
    err = sqrt(abs(err));
    diffphase=atan2(imag(CCmax),real(CCmax));

    if md2 == 1, row_shift = 0; end
    if nd2 == 1, col_shift = 0; end

    output=[err,diffphase,row_shift,col_shift];
end

if (nargout > 1)&&(usfac > 0)
    [nr,nc]=size(buf2ft);
    Nr = ifftshift((-fix(nr/2):ceil(nr/2)-1));
    Nc = ifftshift((-fix(nc/2):ceil(nc/2)-1));
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(1i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(1i*diffphase);
end

    function out=dftups(in,nor,noc,usfac,roff,coff)
        [nr,nc]=size(in);
        if exist('roff','var')~=1, roff=0; end
        if exist('coff','var')~=1, coff=0; end
        if exist('usfac','var')~=1, usfac=1; end
        if exist('noc','var')~=1, noc=nc; end
        if exist('nor','var')~=1, nor=nr; end

        kernc = exp((-1i*2*pi/(nc*usfac))*(ifftshift((0:nc-1)).' - floor(nc/2))*((0:noc-1) - coff));
        kernr = exp((-1i*2*pi/(nr*usfac))*(((0:nor-1).' - roff))*(ifftshift((0:nr-1)) - floor(nr/2)));
        out = kernr*in*kernc;
    end
end
