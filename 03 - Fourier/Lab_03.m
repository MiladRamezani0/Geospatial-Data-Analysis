clc
clear
close all

addpath('SAR_Data');
addpath('Matlab_Functions');

if ispc
    par = '\';
elseif ismac
    par = '/';
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% |                  FOURIER ANALYSIS FOR 1D TIME-SERIES                  |
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


%% STEP 1: Data Import
% This section has the purpose of importing the given environmetal and
% displacement data.

% Displacement data -------------------------------------------------------
% The variable t0_SAR is initialized to the day  of the first observation, 
% while the time is the one of the given environmental data.
filename_SAR = 'SAR_data/S1_Nocelle_4y_ASC.xlsx';

[data_SAR, time_SAR, coord_SAR] = SAR_DataImport(filename_SAR);
t0_SAR = datetime('07-Jul-2019', 'Format', 'dd-MMM-uuuu');
t_fin = 153;
time_SAR = time_SAR(1:t_fin);
data_SAR = data_SAR(:, 1:t_fin);

% PSs definition ----------------------------------------------------------
PS_shift = 21;   % starting PS id
PS_id = 3;       % chosen PS

% Coordinates' conversion -------------------------------------------------
% Convesion of PS coordinates from WGS84 to UTM 33N.
[x_SAR, y_SAR, ~] = deg2utm(coord_SAR(:,2), coord_SAR(:,1));

% The software geoSplinter requires the data to be given ordered for 
% increasing x values.
[~, sort_idx] = sort(x_SAR);
x_SAR1 = x_SAR(sort_idx);
y_SAR1 = y_SAR(sort_idx);

% Import a GeoTIFF for the dam area ---------------------------------------
dam_tif = 'DAM_Nocelle.tif';
[area_dam, R] = readgeoraster(dam_tif);

% Create a grid with the coordinates of the dam area
[rows, cols] = meshgrid(1:size(area_dam,2), 1:size(area_dam,1));
[x_dam, y_dam] = R.intrinsicToWorld(rows,cols);



%% STEP 1: Matlab Functions for Fourier Analysis/Synthesis
% This step is meant to understand how the Matlab built-in functions for
% the Fourier analysis and synthesis work.

% Create a simple vector with 5 values
x = ([1 2 1 3 1.5]');

% Compute the DFT (Discrete Fourier Transform)
% The fft function estimates the Fourier coefficients with the constant
% term as first value of the estimation vector (X1). So, the first half of
% the spectrum will be made by the first (n-1)/2 elements.
X0 = fft(x);
disp(X0);

% For plotting and interpretation reasons it is easier and better to put
% the constant term (corresponding to zero frequency) to the middle of the
% spectrum so that it is symmetric. To reach this goal, the last (n-1)/2 
% elements are mirrored to becomes the first (n-1)/2 elements.
% Example.
% X1(1) -> X((n-1)/2)     = X(3)
% X1(2) -> X((n-1)/2 + 1) = X(4)
% X1(3) -> X((n-1)/2 + 2) = X(5)
% X1(4) -> X((n-1)/2 - 2) = X(1)
% X1(5) -> X((n-1)/2 - 1) = X(2)
% Therefore, we also use the fftshift function as follows.
X = fftshift(fft(x));
disp(X);

% The matlab function ifft can be used to perform the Fourier synthesis,
% which means retriveing the signal from the values of the coefficients of
% the harmonics. Of course, if the order of the parameters has been altered
% through the fftshift function, than the same re-organization of the
% elements inside the vector has to be made.
x0 = ifft(X0);
disp(x0);

x1 = ifft(ifftshift(X));
disp(x1);



%% STEP 2: White Noise
% This section is meant to understand what is the meaning of the white
% noise and how its Fourier spectrum looks like.

% Creation of a random vector with num elements
num = 1001;
wn = randn(1, num);

% Compute the frequency of each element
% The frequencies are symmetric with respect to the 0 (origin of the axis)
% and both positive and negative side will be half of the total number of
% data (this is valid only if the data is odd).
freq_wn = 1/(1*num) * ((-num+1)/2 : 1 : (num-1)/2)';

% Fourier analysis: compute the Fourier coefficients of the white noise
% vector
cF_wn = fftshift(fft(wn));

% Plot the Fourier spectrum of the white noise
figure
plot(freq_wn, abs(cF_wn), 'LineWidth', 2);
xlabel('Frequency')
ylabel('Amplitude')
title('Fourier spectrum of the white noise', 'FontSize', 20)
set(gca, 'FontSize', 15)



%% STEP 3: Data Manipulation
% This step is meant to prepare the SAR displacement data for the Fourier
% analysis. The fundamental requirement is that the chosen dataset has a
% constant sampling rate. Our SAR displacements do not have a constant
% sampling rate over the entire time span. In fact, up to the end of 2021
% there are observations each 6 days (almost), while after that they become
% each 12 days (due to the loss of one of the two satellites).
% Therefore we decide to interpolate each SAR time series with linear
% splines each 6 days. This means that the displacement estimation will be
% exact when we have actual observations, while will be interpolated in the
% holes.

% Creation of a regularized time vector (constant sampling rate)
time_reg = 0:6:time_SAR(end);
nObs = length(time_reg);

% Interpolation of the SAR time series
data_SAR_spl = zeros(length(time_reg),size(data_SAR,1));
for i = 1:size(data_SAR,1)
    data_SAR_spl(:,i) = interp1(time_SAR, data_SAR(i,:), time_reg, 'linear')';
end



%% STEP 4: Computation of the Fourier Spectrum
% This step is meant to compute the Fourier spectrum of the SAR
% displacement time series, in order to estimate the dominating period.
% Indeed, f = 1/T.

% Compute the frequencies for the spectrum
% As the number of observations is the same for all the PSs, it remains
% always the same.
freq_PS = 1/(6*nObs) * ((-nObs+1)/2 : 1 : (nObs-1)/2)';

% PS 24 -------------------------------------------------------------------
PS_24 = 3;
cF_1 = fftshift(fft(data_SAR_spl(:,PS_24)));

% Plot
figure
plot(freq_PS, abs(cF_1), 'LineWidth', 1.5);
xlabel('Frequency')
ylabel('Amplitude')
title('Fourier spectrum of SAR displacement', 'FontSize', 20)
set(gca, 'FontSize', 15)

% PS 39 -------------------------------------------------------------------
PS_39 = 18;
cF_2 = fftshift(fft(data_SAR_spl(:,PS_39)));

% Plot
hold on
plot(freq_PS, abs(cF_2), 'LineWidth', 1.5);

% PS 33 -------------------------------------------------------------------
PS_33 = 12;
cF_3 = fftshift(fft(data_SAR_spl(:,PS_33)));

% Plot
hold on
plot(freq_PS, abs(cF_3), 'LineWidth', 1.5);
legend('PS 24', 'PS 39', 'PS 33');



%% STEP 5: Evaluation of the Peak(s)
% This step is meant to find the frequency corresponding to the peak of the
% Fourier spectrum. This corresponds to finding the dominant (or
% fundamental) period of oscillation for the dam.
% We decide to extract the the peak from the PS where the oscillatory
% behavior can be seen, and then maintain the same period for all the PSs,
% to limit the possible influence of outlier PSs. Of course, the Fourier
% analysis will estimate the coefficents appropriately for each PS.
% Therefore, we choose PS 24 (which peak frequency corresponds also to PS
% 39).

% Find the index of the peak frequency (= fundamental frequency)
[f_max, idx_max] = max(abs(cF_1));
idx_cc = (size(cF_1,1)-1)/2+1;   % index of the constant term (central frequency (=0))

% Creation of an empty vector (all frequencies amplitudes are zero) apart
% from the peak(s). Remember that the spectrum is symmetric.
cF_peak = zeros(size(cF_1));
cF_peak(96) = cF_1(96);           % negative peak
cF_peak(102) = cF_1(102);         % positive peak
cF_peak(idx_cc) = cF_1(idx_cc);   % constant

% Plot the "reduced" spectrum
figure
plot(freq_PS, abs(cF_peak), 'LineWidth', 2)
xlabel('Frequency')
ylabel('Amplitude')
title('Reduced Fourier spectrum of SAR displacement', 'FontSize', 20)
set(gca, 'FontSize', 15)



%% STEP 6: Fourier Syntesis
% This step is meant to perform the Fourier synthesis on the reduced
% spectrum. This means that a single sinusoidal harmonic will decribe the
% seasonality of our signal.

% Perform the synthesis (IFFT)
data_SAR_PS_harm = ifft(ifftshift(cF_peak));

% Plot the synthesis and the raw time series
figure
plot(time_reg, data_SAR_spl(:,PS_24), 'LineWidth', 2)
hold on
plot(time_reg, data_SAR_PS_harm, 'LineWidth', 2)
xlabel('Time [days]')
ylabel('Displacement [mm]')
title('SAR displacement', 'FontSize', 20)
set(gca, 'FontSize', 15)
legend('Raw time series', 'Fourier syntheis', 'FontSize', 15)



%% STEP 7: Fourier Spectrum for all PSs
% This step is meant to replicate the Fourier spectrum for all the PSs. As
% stated before, we will keep the same fundamental frequencies for all the
% PSs, and we will do the synthesis from the corresponding amplitude for
% each time series.

% Creation of the reduced spectrum
cF_PS_peak = zeros(size(data_SAR,1), size(cF_1,1));

for i = 1:size(data_SAR,1)
    cF_i = fftshift(fft(data_SAR_spl(:,i)));
    cF_PS_peak(i,96) = cF_i(96);
    cF_PS_peak(i,102) = cF_i(102);
    cF_PS_peak(i,idx_cc) = cF_i(idx_cc);
end

% Extract the real and imaginery part
rF_PS_np = real(cF_PS_peak(:,96));       % negative peak, real part = cos
iF_PS_np = imag(cF_PS_peak(:,96));       % negative peak, imaginary part = sin

rF_PS_pp = real(cF_PS_peak(:,102));      % positive peak, real part = cos
iF_PS_pp = imag(cF_PS_peak(:,102));      % positive peak, imaginary part = sin

rF_PS_cp = real(cF_PS_peak(:,idx_cc));   % consant peak, only real part

% Compute the residuals of each PS time series
common_epochs = ismember(time_reg, time_SAR);
v_SAR = zeros(size(data_SAR));
data_SAR_Fourier = zeros(size(data_SAR));
for i = 1:size(data_SAR,1)
    data_SAR_PS_harm_tmp = ifft(ifftshift(cF_PS_peak(i,:)));
    data_SAR_Fourier(i,:) = data_SAR_PS_harm_tmp(common_epochs);
    v_SAR(i,:) = data_SAR(i,:) - data_SAR_PS_harm_tmp(common_epochs);
end


%% STEP 8: Spatial Interpolation of the Fourier Coefficients
% This step is meant to perform the spatial interpolation of the Fourier
% coefficients of each PS. We will do this with bilinear splines.
% As a result we will obtain a surface for all the parameters. The result
% is then resampled on a grid of 1 x 1 m over the dam.

% Creation of the necessary folders
if ~exist('./data_input', 'dir')    % data_input folder
    mkdir('./data_input');
end
if ~exist('./data_output', 'dir')   % data_output folder
    mkdir('./data_output');
end
if ~exist('./job', 'dir')           % job folder
    mkdir('./job');
end

% Splines' parameters
n_row = 22;      % number of rows
n_col = 13;      % number of columns
lambda = 0.05;   % regularization parameter

% Negative peak, real part (cosine) ---------------------------------------
filename = 'cosNP';
[raster_rn, parm_rn, sigmaParm_rn, x_nodes, y_nodes] = bil_interp_gS ...
    (x_SAR1, y_SAR1, rF_PS_np, sort_idx, filename, n_row, n_col, lambda);

% Grid resampling
[x_dam1, y_dam1, bil_rn_grid] = gridResampling ...
    (x_nodes, y_nodes, raster_rn, x_dam, y_dam, area_dam);

% Negative peak, imaginary part (sine) ------------------------------------
filename = 'sinNP';
[raster_in, parm_in, sigmaParm_in, ~, ~] = bil_interp_gS ...
    (x_SAR1, y_SAR1, iF_PS_np, sort_idx, filename, n_row, n_col, lambda);

% Grid resampling
[~, ~, bil_in_grid] = gridResampling ...
    (x_nodes, y_nodes, raster_in, x_dam, y_dam, area_dam);

% Positive peak, real part (cosine) ---------------------------------------
filename = 'cosPP';
[raster_rp, parm_rp, sigmaParm_rp, ~, ~] = bil_interp_gS ...
    (x_SAR1, y_SAR1, rF_PS_pp, sort_idx, filename, n_row, n_col, lambda);

% Grid resampling
[~, ~, bil_rp_grid] = gridResampling ...
    (x_nodes, y_nodes, raster_rp, x_dam, y_dam, area_dam);

% Positive peak, imaginary part (sine) ------------------------------------
filename = 'sinPP';
[raster_ip, parm_ip, sigmaParm_ip, ~, ~] = bil_interp_gS ...
    (x_SAR1, y_SAR1, iF_PS_pp, sort_idx, filename, n_row, n_col, lambda);

% Grid resampling
[~, ~, bil_ip_grid] = gridResampling ...
    (x_nodes, y_nodes, raster_ip, x_dam, y_dam, area_dam);

% Constant peak, real part ------------------------------------------------
filename = 'cosPP';
[raster_rc, parm_rc, sigmaParm_rc, ~, ~] = bil_interp_gS ...
    (x_SAR1, y_SAR1, rF_PS_cp, sort_idx, filename, n_row, n_col, lambda);

% Grid resampling
[~, ~, bil_rc_grid] = gridResampling ...
    (x_nodes, y_nodes, raster_rc, x_dam, y_dam, area_dam);



%% STEP 9: Displacement Time Series Estimation for all Coordintes
% This step is meant to estimate the harmonic displacement time series for
% all the coordinates of the dam by performing a Fourier synthesis on the
% spatially interpolated coefficients.

% Re-estimation of the complex numbers (i.e., Fourier coefficients) for
% each coordinate over the dam
bic_n_grid = complex(bil_rn_grid, bil_in_grid);   % negative peak
bic_p_grid = complex(bil_rp_grid, bil_ip_grid);   % positive peak
bic_c_grid = complex(bil_rc_grid);                % constant peak

% Fourier synthesis for each coordinate
dam_SAR_harm = zeros(size(area_dam,1), size(area_dam,2), length(time_reg));

for i = 1:size(bic_n_grid,1)
    for j = 1:size(bic_n_grid,2)

        cF_ij = zeros(size(cF_peak,1),1);
        cF_ij(96) = bic_n_grid(i,j);
        cF_ij(102) = bic_p_grid(i,j);
        cF_ij(idx_cc) = bic_c_grid(i,j);
        if ~all(cF_ij == 0)
            dam_SAR_harm(i,j,:) = ifft(ifftshift(cF_ij));
        end

    end
end

dam_SAR_harm_tobs = dam_SAR_harm(:,:,common_epochs);
    
% Plot --------------------------------------------------------------------
plot_gif = 0;
if plot_gif == 1

    im = cell(length(time_reg), 1);
    
    for i = 1:length(time_reg)
    
        fig = figure('units','normalized','outerposition',[0 0 1 1], 'Color', 'w');
        surf(x_dam1, y_dam1, dam_SAR_harm(:,:,i));
        shading interp;
        title(sprintf('LOS deformation model at epoch %d - Fourier', i), 'FontSize', 35);
        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Displacement [mm]');
        set(gca, 'FontSize', 15);
        colorbar;
        clim([-5 5]);
        zlim([-5 5]);
        % view(2);
        pause(1)
        drawnow
        frame = getframe(fig);
        im{i} = frame2im(frame);
        close;
    
    end
    
    % Save it as a GIF
    filename = 'Dam_LOS_displ_ASC_Fourier_3D.gif';
    for idx = 1:length(time_reg)
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',1,'DelayTime',0.3);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
        end
    end

end



%% STEP 10: Residuals' Computation
% This step is meant to compute the residuals for each SAR displacement
% time series (so, for each PS).

% Extract the modelled harmonic time series for each PS from the 4D
% function
data_SAR_harm = zeros(size(data_SAR));
data_SAR_harm_tmp = zeros(size(data_SAR_spl'));

for i = 1:size(x_SAR,1)

    % Calculate the Euclidean distance
    distances = sqrt((x_dam - x_SAR(i)).^2 + (y_dam - y_SAR(i)).^2);
    
    % Find the row and column indices of the point with the minimum distance
    [minDistance, linearIndex] = min(distances(:));
    [rowIndex, colIndex] = ind2sub(size(distances), linearIndex);

    % Store the result in a matrix
    data_SAR_harm_tmp(i,:) = reshape(dam_SAR_harm(rowIndex, colIndex,:), 1, size(dam_SAR_harm,3));

    % Take only the epochs in which we had the observations
    data_SAR_harm(i,:) = data_SAR_harm_tmp(common_epochs);

end

% Computation of residuals as difference between the observed displacement
% and the estimated harmonic
v_SAR_interp = data_SAR - data_SAR_harm;



%% STEP 12: Temporal and Spatial Interpolation of Residuals
% This step is meant to perform the temporal and spatial interpolation of
% the residuals with bilinear splines.
% [SAME CODE AS IN LAB 02]

[res_dam] = splineInterp4D ...
    (v_SAR_interp, time_SAR, PS_shift, sort_idx, x_SAR1, y_SAR1, dam_tif);

% Plot
filename = 'Dam_LOS_res_ASC_Fspl_3D.gif';
plotGIF(x_dam1, y_dam1, res_dam, filename, plot_gif);



%% STEP 13: Variables for Collocation
% This step is meant to save the required variables for the collocation
% analysis.
if ~exist('../Lab 04 - Collocation/', 'dir')
    mkdir('../Lab 04 - Collocation/');
end

save('../Lab 04 - Collocation/Fourier_out.mat', 'v_SAR', 'data_SAR_Fourier', 'freq_PS', 'x_dam1', 'y_dam1', 'dam_SAR_harm_tobs');


fprintf('\nScript execution completed!\n');

