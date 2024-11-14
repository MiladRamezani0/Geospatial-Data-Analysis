function [SAR_ts, SAR_time, SAR_coord] = SAR_DataImport(SAR_filename)

% This function imports the SAR data from an Excel file obtained as output
% from the PSI processing.
%
% INPUT VARIABLE:
% - SAR_filename: path to the .xlsx file to be imported
%
% OUTPUT VARIABLES:
% - SAR_ts: SAR time series for all the Persistent Scatterers in the area
% - SAR_time: dates of acquisition of the SAR images
% - SAR_coord: coordinates of the SAR PSs
%
%
% (c) Roberto Monti, version 1.0


% Check if the file exists
if ~isfile(SAR_filename)
    error('File does not exist: %s', SAR_filename);
end

% Check if the input file is in the correct format
[SAR_file, ~, ~] = xlsread(SAR_filename);
if size(SAR_file, 1) < 3 || size(SAR_file, 2) < 5
    error('Input file is not in the correct format.');
end

% Prepare the data to be exported
SAR_file = readmatrix(SAR_filename);
SAR_ts = SAR_file(2:end, 5:end);
SAR_time = SAR_file(1, 5:end);
SAR_coord = SAR_file(2:end, 2:3);

disp('SAR data import successfully completed!')
