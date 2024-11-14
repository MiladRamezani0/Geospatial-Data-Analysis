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
% |                SPLINE INTERPOLATION FOR 1D AND 2D DATA                |
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
epoch = 10;       % chosen epoch

% Coordinates' conversion -------------------------------------------------
% Convesion of PS coordinates from WGS84 to UTM 33N
[x_SAR, y_SAR, ~] = deg2utm(coord_SAR(:,2), coord_SAR(:,1));



%% STEP 2: Data Preparation
% This step is meant to create the dataset .txt files to be used during the
% spline analysis.

% Directories -------------------------------------------------------------
% data_input
if ~exist('data_input', 'dir')
    mkdir('data_input')
end

% data_output
if ~exist('data_output', 'dir')
    mkdir('data_output')
end

% 1D dataset --------------------------------------------------------------
% Extraction of the displacement time-series for one PS from the file. It
% is then created a .txt file with time of observation as first column and
% corresponding displacement as the second.
if ~exist(strcat('./data_input/S1_PS', num2str(PS_id+PS_shift), '_ASC.txt'), 'file')
    writematrix([time_SAR', data_SAR(PS_id,:)'], ...
        strcat('./data_input/S1_PS', num2str(PS_id+PS_shift), '_ASC.txt'), 'Delimiter', 'space')
end

% 1D dataset with hole ----------------------------------------------------
% The previously defined time series is modified so that it now containes a
% hole (= gap in the data).
h_in = 41;
h_fin = 60;
data_SAR_hole = data_SAR(PS_id,:);
time_SAR_hole = time_SAR;
time_SAR_hole(h_in:h_fin) = [];
data_SAR_hole(h_in:h_fin) = [];

if ~exist(strcat('./data_input/S1_PS', num2str(PS_id+PS_shift), '_ASC_hole.txt'), 'file')
    writematrix([time_SAR', data_SAR(PS_id,:)'], ...
        strcat('./data_input/S1_PS', num2str(PS_id+PS_shift), '_ASC_hole.txt'), 'Delimiter', 'space')
end

% 2D dataset --------------------------------------------------------------
% Extraction of the displacement time-series for one PS from the file. It
% is then created a .txt file with time of observation as first column and
% corresponding displacement as the second.
% The software geoSplinter requires the data to be given ordered for 
% increasing x values.
[~, sort_idx] = sort(x_SAR);
x_SAR1 = x_SAR(sort_idx);
y_SAR1 = y_SAR(sort_idx);
data_SAR_epoch = data_SAR(:,epoch);
data_SAR_epoch = data_SAR_epoch(sort_idx);

if ~exist(strcat('./data_input/S1_t', num2str(epoch), '_ASC.txt'), 'file')
    writematrix([x_SAR1, y_SAR1, data_SAR_epoch], ...
        strcat('./data_input/S1_t', num2str(epoch), '_ASC.txt'), 'Delimiter', 'space')
end



%% STEP 3: 1D Linear Spline (16), without Regularization λ = 0)
% - Job file: S1_PS24_ASC_lin0.job
% Parameters:
% Dimension of dataset (1D/2D):         1
% Type of splines (linear/cubic):       1
% Delta discretization (1delta/2delta): 1
% Input file:                           ./data_input/S1_PS24_ASC.txt
% Output file:                          ./data_output/S1_PS24_ASC_lin0
% Number of observations:               153
% Number of nodes (= splines):          16
% First abscissa (= time):              0
% Last abscissa (= time):               1176
% Regularization parameter (λ):         0
% Number of significant digits:         8

suffix_0 = 'lin0';

[data_lin0, raster_lin0, par_lin0, sigmaPar_lin0] = ...
  geoSplinter (strcat('.', par, 'data_output', par, 'S1_PS', num2str(PS_id+PS_shift), '_ASC_', suffix_0), 'lin');



%% STEP 4: 1D Linear Spline (31), without Regularization λ = 0)
% - Job file: S1_PS24_ASC_lin1.job
% Parameters:
% Dimension of dataset (1D/2D):         1
% Type of splines (linear/cubic):       1
% Delta discretization (1delta/2delta): 1
% Input file:                           ./data_input/S1_PS24_ASC.txt
% Output file:                          ./data_output/S1_PS24_ASC_lin1
% Number of observations:               153
% Number of nodes (= splines):          31
% First abscissa (= time):              0
% Last abscissa (= time):               1176
% Regularization parameter (λ):         0
% Number of significant digits:         8

suffix_1 = 'lin1';

[data_lin1, raster_lin1, par_lin1, sigmaPar_lin1] = ...
    geoSplinter (strcat('.', par, 'data_output', par, 'S1_PS', num2str(PS_id+PS_shift), '_ASC_', suffix_1), 'lin');



%% STEP 5: 1D Linear Spline (31), with Regularization λ = 0.1)
% - Job file: S1_PS24_ASC_lin2.job
% Parameters:
% Dimension of dataset (1D/2D):         1
% Type of splines (linear/cubic):       1
% Delta discretization (1delta/2delta): 1
% Input file:                           ./data_input/S1_PS24_ASC.txt
% Output file:                          ./data_output/S1_PS24_ASC_lin2
% Number of observations:               153
% Number of nodes (= splines):          31
% First abscissa (= time):              0
% Last abscissa (= time):               1176
% Regularization parameter (λ):         0.1
% Number of significant digits:         8

suffix_2 = 'lin2';

[data_lin2, raster_lin2, par_lin2, sigmaPar_lin2] = ...
    geoSplinter (strcat('.', par, 'data_output', par, 'S1_PS', num2str(PS_id+PS_shift), '_ASC_', suffix_2), 'lin');



%% STEP 6: 1D Linear Spline (31), with Regularization λ = 10)
% - Job file: S1_PS24_ASC_lin3.job
% Parameters:
% Dimension of dataset (1D/2D):         1
% Type of splines (linear/cubic):       1
% Delta discretization (1delta/2delta): 1
% Input file:                           ./data_input/S1_PS24_ASC.txt
% Output file:                          ./data_output/S1_PS24_ASC_lin3
% Number of observations:               153
% Number of nodes (= splines):          31
% First abscissa (= time):              0
% Last abscissa (= time):               1176
% Regularization parameter (λ):         10
% Number of significant digits:         8

suffix_3 = 'lin3';

[data_lin3, raster_lin3, par_lin3, sigmaPar_lin3] = ...
    geoSplinter (strcat('.', par, 'data_output', par, 'S1_PS', num2str(PS_id+PS_shift), '_ASC_', suffix_3), 'lin');



%% STEP 7: 1D Linear Spline (31), with Regularization λ = 1000)
% - Job file: S1_PS24_ASC_lin4.job
% Parameters:
% Dimension of dataset (1D/2D):         1
% Type of splines (linear/cubic):       1
% Delta discretization (1delta/2delta): 1
% Input file:                           ./data_input/S1_PS24_ASC.txt
% Output file:                          ./data_output/S1_PS24_ASC_lin4
% Number of observations:               153
% Number of nodes (= splines):          31
% First abscissa (= time):              0
% Last abscissa (= time):               1176
% Regularization parameter (λ):         1000
% Number of significant digits:         8

suffix_4 = 'lin4';

[data_lin4, raster_lin4, par_lin4, sigmaPar_lin4] = ...
    geoSplinter (strcat('.', par, 'data_output', par, 'S1_PS', num2str(PS_id+PS_shift), '_ASC_', suffix_4), 'lin');



%% STEP 8: 1D Linear Spline (31), with Regularization λ = 10000)
% - Job file: S1_PS24_ASC_lin5.job
% Parameters:
% Dimension of dataset (1D/2D):         1
% Type of splines (linear/cubic):       1
% Delta discretization (1delta/2delta): 1
% Input file:                           ./data_input/S1_PS24_ASC.txt
% Output file:                          ./data_output/S1_PS24_ASC_lin5
% Number of observations:               153
% Number of nodes (= splines):          31
% First abscissa (= time):              0
% Last abscissa (= time):               1176
% Regularization parameter (λ):         10000
% Number of significant digits:         8

suffix_5 = 'lin5';

[data_lin5, raster_lin5, par_lin5, sigmaPar_lin5] = ...
    geoSplinter (strcat('.', par, 'data_output', par, 'S1_PS', num2str(PS_id+PS_shift), '_ASC_', suffix_5), 'lin');



%% STEP 9: 1D Linear Spline (31), with Regularization λ = 5, Hole)
% - Job file: S1_PS24_ASC_lin_hole0.job
% Parameters:
% Dimension of dataset (1D/2D):         1
% Type of splines (linear/cubic):       1
% Delta discretization (1delta/2delta): 1
% Input file:                           ./data_input/S1_PS24_ASC_hole.txt
% Output file:                          ./data_output/S1_PS24_ASC_lin_hole0
% Number of observations:               153
% Number of nodes (= splines):          31
% First abscissa (= time):              0
% Last abscissa (= time):               1176
% Regularization parameter (λ):         5
% Number of significant digits:         8

suffix_6 = 'lin_hole0';

[data_lin_hole0, raster_lin_hole0, par_lin_hole0, sigmaPar_lin_hole0] = ...
    geoSplinter (strcat('.', par, 'data_output', par, 'S1_PS', num2str(PS_id+PS_shift), '_ASC_', suffix_6), 'lin');



%% STEP 10: 1D Cubic Spline (26), without Regularization λ = 0)
% - Job file: S1_PS24_ASC_cub0.job
% Parameters:
% Dimension of dataset (1D/2D):         1
% Type of splines (linear/cubic):       2
% Input file:                           ./data_input/S1_PS24_ASC.txt
% Output file:                          ./data_output/S1_PS24_ASC_cub0
% Number of observations:               153
% Number of nodes (= splines):          26
% First abscissa (= time):              0
% Last abscissa (= time):               1176
% Regularization parameter (λ):         0
% Number of significant digits:         8

suffix_7 = 'cub0';

[data_cub0, raster_cub0, par_cub0, sigmaPar_cub0] = ...
    geoSplinter (strcat('.', par, 'data_output', par, 'S1_PS', num2str(PS_id+PS_shift), '_ASC_', suffix_7), 'cub');



%% STEP 11: 1D Cubic Spline (26), with Regularization λ = 5, Hole)
% - Job file: S1_PS24_ASC_cub_hole0.job
% Parameters:
% Dimension of dataset (1D/2D):         1
% Type of splines (linear/cubic):       2
% Input file:                           ./data_input/S1_PS24_ASC_hole.txt
% Output file:                          ./data_output/S1_PS24_ASC_cub_hole0
% Number of observations:               153
% Number of nodes (= splines):          26
% First abscissa (= time):              0
% Last abscissa (= time):               1176
% Regularization parameter (λ):         5
% Number of significant digits:         8

suffix_8 = 'cub_hole0';

[data_cub_hole0, raster_cub_hole0, par_cub_hole0, sigmaPar_cub_hole0] = ...
    geoSplinter (strcat('.', par, 'data_output', par, 'S1_PS', num2str(PS_id+PS_shift), '_ASC_', suffix_8), 'cub');



%% STEP 11: 1D Cubic Spline (26), with Regularization λ = 1, Hole)
% - Job file: S1_PS24_ASC_cub_hole1.job
% Parameters:
% Dimension of dataset (1D/2D):         1
% Type of splines (linear/cubic):       2
% Input file:                           ./data_input/S1_PS24_ASC_hole.txt
% Output file:                          ./data_output/S1_PS24_ASC_cub_hole1
% Number of observations:               153
% Number of nodes (= splines):          26
% First abscissa (= time):              0
% Last abscissa (= time):               1176
% Regularization parameter (λ):         1
% Number of significant digits:         8

suffix_9 = 'cub_hole1';

[data_cub_hole1, raster_cub_hole1, par_cub_hole1, sigmaPar_cub_hole1] = ...
    geoSplinter (strcat('.', par, 'data_output', par, 'S1_PS', num2str(PS_id+PS_shift), '_ASC_', suffix_9), 'cub');



%% STEP 11: 2D Bilinear Spline (15x13), with Regularization λ = 1, Holes)
% - Job file: S1_t10_ASC_bil0.job
% Parameters:
% Dimension of dataset (1D/2D):         2
% Type of splines (linear/cubic):       1
% Delta discretization (1delta/2delta): 1
% Input file:                           ./data_input/S1_t10_ASC.txt
% Output file:                          ./data_output/S1_t10_ASC_bil0
% Number of observations (=PS):         29
% Number of rows:                       15
% Number of columns:                    13
% Minimum X:                            633346.153342255
% Maximum X:                            633446.788004437
% Minimum Y:                            4345054.03968278
% Maximum Y:                            4345220.06722414
% Regularization parameter (λ):         5
% Number of significant digits:         8

suffix_10 = 'bil0';

[data_bil0, raster_bil0, par_bil0, sigmaPar_bil0] = ...
    geoSplinter (strcat('.', par, 'data_output', par, 'S1_t', num2str(epoch), '_ASC_', suffix_10), 'bil');




%% STEP 12: 1D and 2D Combination
% APPLICATION OF THESE CONCEPTS TO THE CREATION OF A CONTINOUS SPATIAL-
% TEMPORAL DEFORMATION MODEL FOR THE DAM.

% 1D linear splines  ------------------------------------------------------
% For the modelling of all PSs deformation time-series.

% Creation of the input .txt files
for i = 1:size(data_SAR,1)
    writematrix([time_SAR', data_SAR(i,:)'], ...
        strcat('./data_input/S1_PS', num2str(i+PS_shift), '_ASC_ts.txt'), 'Delimiter', 'space')
end

% Creation of the .job file (S1_PSi_ASC_ts_lin.job)
% Parameters:
% Dimension of dataset (1D/2D):         1
% Type of splines (linear/cubic):       1
% Delta discretization (1delta/2delta): 1
% Input file:                           ./data_input/S1_PSi_ASC_ts.txt
% Output file:                          ./data_output/S1_PSi_ASC_ts_lin
% Number of observations:               153
% Number of nodes (= splines):          31
% First abscissa (= time):              0
% Last abscissa (= time):               1176
% Regularization parameter (λ):         5
% Number of significant digits:         8

n_spl = 31;   % number of splines
par_ts = zeros(n_spl, size(data_SAR,1));
sigmaPar_ts = zeros(n_spl, size(data_SAR,1));

for i = 1:size(data_SAR,1)

    job_file = strcat('.', par, 'job', par, 'S1_PS', num2str(i+PS_shift), '_ASC_ts_lin.job');
    fid = fopen(job_file, 'w');
    
    fprintf(fid, '\n%s\n', '1');
    fprintf(fid, '%s\n', '1');
    fprintf(fid, '%s\n', '1');
    fprintf(fid, '%s\n', strcat('.', par, 'data_input', par, 'S1_PS', num2str(i+PS_shift), '_ASC_ts.txt'));
    fprintf(fid, '%s\n', strcat('.', par, 'data_output', par, 'S1_PS', num2str(i+PS_shift), '_ASC_ts_lin'));
    fprintf(fid, '%s\n', '153');
    fprintf(fid, '%s\n', num2str(n_spl));
    fprintf(fid, '%s\n', '0');
    fprintf(fid, '%s\n', '1176');
    fprintf(fid, '%s\n', '5');
    fprintf(fid, '%s\n', '8');
    fprintf(fid, '%s', '.');
    
    fclose(fid);
    
    if ismac
        job_execution = strcat('./geoSplinter_analysis_macOS < ./job/S1_PS', num2str(i+PS_shift), '_ASC_ts_lin.job');
    elseif ispc
        job_execution = strcat('.\geoSplinter_analysis_win < .\job\S1_PS', num2str(i+PS_shift), '_ASC_ts_lin.job');
    end
    
    system(job_execution)
    system('exit')

    [data_lin_tmp, raster_lin_tmp, par_lin_tmp, sigmaPar_lin_tmp, time_1Dlin] = ...
        geoSplinter_noFig (strcat('./data_output/S1_PS', num2str(i+PS_shift), '_ASC_ts_lin'), 'lin');

    par_ts(:,i) = par_lin_tmp;
    sigmaPar_ts(:,i) = sigmaPar_lin_tmp;

end


% -------------------------------------------------------------------------


% 2D blinear splines  -----------------------------------------------------
% To spatially interpolate all the splines' parameters.

for i = 1:size(par_ts,1)

    par_ts_i = par_ts(i,:);
    par_ts_i = par_ts_i(sort_idx);
    
    writematrix([x_SAR1, y_SAR1, par_ts_i'], ...
        strcat('.', par, 'data_input', par, 'S1_ASC_parm', num2str(i), '.txt'), 'Delimiter', 'space')


end

% Creation of the .job file (S1_ASC_parmi_bil.job)
% Parameters:
% Dimension of dataset (1D/2D):         2
% Type of splines (linear/cubic):       1
% Delta discretization (1delta/2delta): 1
% Input file:                           ./data_input/S1_ASC_parmi.txt
% Output file:                          ./data_output/S1_ASC_parmi_bil
% Number of observations (=PS):         29
% Number of rows:                       15
% Number of columns:                    13
% Minimum X:                            633346.153342255
% Maximum X:                            633446.788004437
% Minimum Y:                            4345054.03968278
% Maximum Y:                            4345220.06722414
% Regularization parameter (λ):         0.05
% Number of significant digits:         8

n_row = 15;   % number of rows
n_col = 13;   % number of columns
par_ST = zeros(n_row, n_col, size(par_ts,1));
sigmaPar_ST = zeros(n_row, n_col, size(par_ts,1));
x_ST = zeros(n_row, n_col, size(par_ts,1));
y_ST = zeros(n_row, n_col, size(par_ts,1));

for i = 1:size(par_ts,1)

    job_file = strcat('.', par, 'job', par, 'S1_ASC_parm', num2str(i), '_bil.job');
    fid = fopen(job_file, 'w');
    
    fprintf(fid, '\n%s\n', '2');
    fprintf(fid, '%s\n', '1');
    fprintf(fid, '%s\n', '1');
    fprintf(fid, '%s\n', strcat('.', par, 'data_input', par, 'S1_ASC_parm', num2str(i), '.txt'));
    fprintf(fid, '%s\n', strcat('.', par, 'data_output', par, 'S1_ASC_parm', num2str(i), '_bil'));
    fprintf(fid, '%s\n', '29');
    fprintf(fid, '%s\n', num2str(n_row));
    fprintf(fid, '%s\n', num2str(n_col));
    fprintf(fid, '%s\n', num2str(min(x_SAR1)));
    fprintf(fid, '%s\n', num2str(max(x_SAR1)));
    fprintf(fid, '%s\n', num2str(min(y_SAR1)));
    fprintf(fid, '%s\n', num2str(max(y_SAR1)));
    fprintf(fid, '%s\n', '0.05');
    fprintf(fid, '%s\n', '8');
    fprintf(fid, '%s', '.');
    
    fclose(fid);
    
    if ismac
        job_execution = strcat('./geoSplinter_analysis_macOS < ./job/S1_ASC_parm', num2str(i), '_bil.job');
    elseif ispc
        job_execution = strcat('.\geoSplinter_analysis_win < .\job\S1_ASC_parm', num2str(i), '_bil.job');
    end
    
    system(job_execution)
    system('exit')

    [data_ST_tmp, raster_ST_tmp, par_ST_tmp, sigmaPar_ST_tmp, nodes_coord_ST_tmp] = ...
        geoSplinter_noFig (strcat('.', par, 'data_output', par, 'S1_ASC_parm', num2str(i), '_bil'), 'bil');

    par_ST(:,:,i) = par_ST_tmp;
    sigmaPar_ST(:,:,i) = sigmaPar_ST_tmp;
    x_ST(:,:,i) = nodes_coord_ST_tmp(:,:,1);   % X coordinates of the bilinear splines' nodes
    y_ST(:,:,i) = nodes_coord_ST_tmp(:,:,2);   % X coordinates of the bilinear splines' nodes

end



%% STEP 13: Deformation Time-Series for Every Point on the Dam
% The estimated parameters of the bilinear splines interpolating all the
% coefficients are now used to estimate the displacement time-series all
% over the dam.
plot_gif = 1;

% Import a GeoTIFF for the dam area ---------------------------------------
dam_tif = 'DAM_Nocelle.tif';
[area_dam, R] = readgeoraster(dam_tif);

% Create a grid with the coordinates of the dam area
[rows, cols] = meshgrid(1:size(area_dam,2), 1:size(area_dam,1));
[x_dam, y_dam] = R.intrinsicToWorld(rows,cols);

% Compute the displacement time-series for each point ---------------------
% The variable displ_dam_S will contain the displacement over the dam for
% each epoch where there is a spline (31) and, resampled over a spatial
% grid with the same resolution as for Laboratory 01.
displ_dam_S = zeros(size(x_dam,1), size(x_dam,2), size(par_ts,1));
for i = 1:size(par_ts,1)

    par_ST_tmp = par_ST(:,:,i);   % one column below each other
    x_ST_tmp = x_ST(:,:,i);
    y_ST_tmp = y_ST(:,:,i);
    F = scatteredInterpolant(x_ST_tmp(:), y_ST_tmp(:), par_ST_tmp(:), 'linear');
    displ_dam_tmp = F(x_dam(:), y_dam(:));

    displ_dam_S(:,:,i) = reshape(displ_dam_tmp, size(x_dam,1), size(x_dam,2));

    clear displ_dam_tmp

    displ_dam_S(:,:,i) = (displ_dam_S(:,:,i)) .* area_dam;
    displ_dam_S(area_dam == 0) = nan;
    x_dam1 = x_dam .* area_dam;
    x_dam1(area_dam == 0) = nan;
    y_dam1 = y_dam .* area_dam;
    y_dam1(area_dam == 0) = nan;

end

% Compute the displacement time-series for each point ---------------------
% Linear interpolation for the times where there are no splines. We put 31
% linear splines so we have 31 displacements, while for the other times we
% need to obtain the displacement by linear interpolation.
displ_dam = nan(size(x_dam,1), size(x_dam,2), length(time_SAR));
for i = 1:size(x_dam,1)
    for j = 1:size(x_dam,2)

        displ_dam(i,j,:) = interp1(time_1Dlin', reshape(displ_dam_S(i,j,:), size(displ_dam_S,3), 1, 1), time_SAR, 'linear');
    
    end
end


% Plot --------------------------------------------------------------------
if plot_gif == 1

    im = cell(length(time_SAR), 1);
    
    for i = 1:length(time_SAR)
    
        fig = figure('units','normalized','outerposition',[0 0 1 1], 'Color', 'w');
        surf(x_dam1, y_dam1, displ_dam(:,:,i));
        shading interp;
        title(sprintf('LOS deformation model at epoch %d - Splines', i), 'FontSize', 35);
        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Displacement [mm]');
        set(gca, 'FontSize', 15);
        colorbar;
        clim([-5 5]);
        zlim([-10 5]);
        %view(2);
        pause(1)
        drawnow
        frame = getframe(fig);
        im{i} = frame2im(frame);
        close;
    
    end
    
    % Save it as a GIF
    filename = 'Dam_LOS_displ_ASC_spline_3D.gif';
    for idx = 1:length(time_SAR)
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',1,'DelayTime',0.3);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
        end
    end

end

fprintf('\nScript execution completed!\n');

