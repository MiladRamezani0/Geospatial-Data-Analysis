function [raster, parm, sigmaParm, x_nodes, y_nodes] = bic_interp_gS ...
    (x_coord, y_coord, z_value, sort_idx, filename, n_row, n_col, lambda)


% Definition the parenthesis
if ispc
    par = '\';
elseif ismac
    par = '/';
end

% Creation of the input file
z_value = z_value(sort_idx);
writematrix([x_coord, y_coord, z_value], ...
    strcat('.', par, 'data_input', par, filename, '.txt'), 'Delimiter', 'space')


% Creation of the .job file
% Parameters:
% Dimension of dataset (1D/2D):         2
% Type of splines (linear/cubic):       2
% Input file:                           ./data_input/filenamei.txt
% Output file:                          ./data_output/filenamei_bic
% Number of observations (=PS):         29
% Number of rows:                       15
% Number of columns:                    13
% Minimum X:                            633346.153342255
% Maximum X:                            633446.788004437
% Minimum Y:                            4345054.03968278
% Maximum Y:                            4345220.06722414
% Regularization parameter (λ):         0.05
% Number of significant digits:         8

job_file = strcat('.', par, 'job', par, filename, '_bic.job');
fid = fopen(job_file, 'w');

fprintf(fid, '\n%s\n', '2');
fprintf(fid, '%s\n', '2');
fprintf(fid, '%s\n', strcat('.', par, 'data_input', par, filename, '.txt'));
fprintf(fid, '%s\n', strcat('.', par, 'data_output', par, filename, '_bic'));
fprintf(fid, '%s\n', num2str(size(z_value,1)));
fprintf(fid, '%s\n', num2str(n_row));
fprintf(fid, '%s\n', num2str(n_col));
fprintf(fid, '%s\n', num2str(min(x_coord)));
fprintf(fid, '%s\n', num2str(max(x_coord)));
fprintf(fid, '%s\n', num2str(min(y_coord)));
fprintf(fid, '%s\n', num2str(max(y_coord)));
fprintf(fid, '%s\n', num2str(lambda));
fprintf(fid, '%s\n', '8');
fprintf(fid, '%s', '.');

fclose(fid);

if ismac
    job_execution = strcat('./geoSplinter_analysis_macOS < ./job/', filename, '_bic.job');
elseif ispc
    job_execution = strcat('.\geoSplinter_analysis_win < .\job\', filename, '_bic.job');
end

system(job_execution)
system('exit')

[~, raster_tmp, parm_tmp, sigmaParm_tmp, nodes_coord_tmp] = ...
    geoSplinter_noFig (strcat('.', par, 'data_output', par, filename, '_bic'), 'bic');

raster = raster_tmp;
parm = parm_tmp;
sigmaParm = sigmaParm_tmp;
x_nodes = nodes_coord_tmp(:,:,1);   % X coordinates of the bilinear splines' nodes
y_nodes = nodes_coord_tmp(:,:,2);   % X coordinates of the bilinear splines' nodes
