function [res_dam] = splineInterp4D(v_SAR, time_SAR, PS_shift, sort_idx, x_SAR1, y_SAR1, dam_tif)

% INPUT VARIABLES:
% - residuals time series;
% - time of observations;
% - PS_shift;
% - index to sort the coordinates with increasing X;
% - sorted X coordinates of PSs;
% - sorted Y coordinates of PSs;
% - tif file of the dam area;



if ispc
    par = '\';
elseif ismac
    par = '/';
end

% 1D linear splines  ------------------------------------------------------
% For the modelling of all PSs residuals time-series.

% Creation of the input .txt files
for i = 1:size(v_SAR,1)
    writematrix([time_SAR', v_SAR(i,:)'], ...
        strcat('./data_input/S1_PS', num2str(i+PS_shift), '_ASC_v.txt'), 'Delimiter', 'space')
end

% Creation of the .job file (S1_PSi_ASC_v_lin.job)
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
% Regularization parameter (λ):         1
% Number of significant digits:         8

n_spl = 31;   % number of splines
par_v = zeros(n_spl, size(v_SAR,1));
sigmaPar_v = zeros(n_spl, size(v_SAR,1));

for i = 1:size(v_SAR,1)

    job_file = strcat('.', par, 'job', par, 'S1_PS', num2str(i+PS_shift), '_ASC_v_lin.job');
    fid = fopen(job_file, 'w');
    
    fprintf(fid, '\n%s\n', '1');
    fprintf(fid, '%s\n', '1');
    fprintf(fid, '%s\n', '1');
    fprintf(fid, '%s\n', strcat('.', par, 'data_input', par, 'S1_PS', num2str(i+PS_shift), '_ASC_v.txt'));
    fprintf(fid, '%s\n', strcat('.', par, 'data_output', par, 'S1_PS', num2str(i+PS_shift), '_ASC_v_lin'));
    fprintf(fid, '%s\n', '153');
    fprintf(fid, '%s\n', num2str(n_spl));
    fprintf(fid, '%s\n', '0');
    fprintf(fid, '%s\n', '1176');
    fprintf(fid, '%s\n', '1');
    fprintf(fid, '%s\n', '8');
    fprintf(fid, '%s', '.');
    
    fclose(fid);
    
    if ismac
        job_execution = strcat('./geoSplinter_analysis_macOS < ./job/S1_PS', num2str(i+PS_shift), '_ASC_v_lin.job');
    elseif ispc
        job_execution = strcat('.\geoSplinter_analysis_win < .\job\S1_PS', num2str(i+PS_shift), '_ASC_v_lin.job');
    end
    
    system(job_execution)
    system('exit')

    [data_lin_tmp, raster_lin_tmp, par_lin_tmp, sigmaPar_lin_tmp, time_1Dlin] = ...
        geoSplinter_noFig (strcat('./data_output/S1_PS', num2str(i+PS_shift), '_ASC_v_lin'), 'lin');

    par_v(:,i) = par_lin_tmp;
    sigmaPar_v(:,i) = sigmaPar_lin_tmp;

end


% -------------------------------------------------------------------------


% 2D blinear splines  -----------------------------------------------------
% To spatially interpolate all the splines' parameters.

for i = 1:size(par_v,1)

    par_v_i = par_v(i,:);
    par_v_i = par_v_i(sort_idx);
    
    writematrix([x_SAR1, y_SAR1, par_v_i'], ...
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
par_ST = zeros(n_row, n_col, size(par_v,1));
sigmaPar_ST = zeros(n_row, n_col, size(par_v,1));
x_ST = zeros(n_row, n_col, size(par_v,1));
y_ST = zeros(n_row, n_col, size(par_v,1));

for i = 1:size(par_v,1)

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

% Compute the displacement time-series for each point ---------------------
% The variable displ_dam_S will contain the displacement over the dam for
% each epoch where there is a spline (31) and, resampled over a spatial
% grid with the same resolution as for Laboratory 01.

[area_dam, R] = readgeoraster(dam_tif);

% Create a grid with the coordinates of the dam area
[rows, cols] = meshgrid(1:size(area_dam,2), 1:size(area_dam,1));
[x_dam, y_dam] = R.intrinsicToWorld(rows,cols);

res_dam_S = zeros(size(x_dam,1), size(x_dam,2), size(par_v,1));
for i = 1:size(par_v,1)

    par_ST_tmp = par_ST(:,:,i);   % one column below each other
    x_ST_tmp = x_ST(:,:,i);
    y_ST_tmp = y_ST(:,:,i);
    F = scatteredInterpolant(x_ST_tmp(:), y_ST_tmp(:), par_ST_tmp(:), 'linear');
    displ_dam_tmp = F(x_dam(:), y_dam(:));

    res_dam_S(:,:,i) = reshape(displ_dam_tmp, size(x_dam,1), size(x_dam,2));

    clear displ_dam_tmp

    res_dam_S(:,:,i) = (res_dam_S(:,:,i)) .* area_dam;
    res_dam_S(area_dam == 0) = nan;

end

% Compute the displacement time-series for each point ---------------------
% Linear interpolation for the times where there are no splines. We put 31
% linear splines so we have 31 displacements, while for the other times we
% need to obtain the displacement by linear interpolation.
res_dam = nan(size(x_dam,1), size(x_dam,2), length(time_SAR));
for i = 1:size(x_dam,1)
    for j = 1:size(x_dam,2)

        res_dam(i,j,:) = interp1(time_1Dlin', reshape(res_dam_S(i,j,:), size(res_dam_S,3), 1, 1), time_SAR, 'linear');
    
    end
end