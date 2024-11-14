clc
clear
close all

addpath('SAR_Data');
addpath('Matlab_Functions');
addpath('Environmental_Data');


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% | REGRESSION ANALYSIS OF ENVIRONMENTAL AND DISPLACEMENT DATA FOR A DAM  |
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


%% STEP 1: Data Import
% This section has the purpose of importing the given environmetal and
% displacement data.

% Displacement data -------------------------------------------------------
% The variable t0_SAR is initialized to the day  of the first observation, 
% while the time is the one of the given environmental data.
filename_SAR = 'S1_Nocelle_4y_ASC.xlsx';

[data_raw_SAR, time_SAR, coord_SAR] = SAR_DataImport(filename_SAR);
t0_SAR = datetime('07-Jul-2019 08:00', 'Format', 'dd-MMM-uuuu HH:mm');
PS_shift = 21;   % starting PS id

% Environmental data ------------------------------------------------------
filename_Twater = 'TEMPERARVO_ACQUA_RAW.DAT';
filename_Hwater = 'IDRO_RAW.DAT';

table_Twater = readtable(filename_Twater);
time_Twater = days(table2array(table_Twater(:,1)) - t0_SAR)';
data_raw_Twater = table2array(table_Twater(:,2));

table_Hwater = readtable(filename_Hwater);
time_Hwater = days(table2array(table_Hwater(:,1)) - t0_SAR)';
data_raw_Hwater = table2array(table_Hwater(:,2));

% Estimation of the minimum final day between the three time series, to be
% set as the last considered day for the analysis.
tfin = min([time_Hwater(end-1), time_Twater(end-1), time_SAR(end)]);

% Flags for processing choices --------------------------------------------
plot_gifs = 0;         % to plot (1) or not (0) the gifs images
data_SAR_choice = 1;   % to use raw (0) or smoothed (1) SAR time-series



%% STEP 2: Data Manipulation
% This step is meant to smooth the SAR displacement data with a moving
% average to partially remove the noise in the observations.

data_raw_smt_SAR = zeros(size(data_raw_SAR,1), size(data_raw_SAR,2));
window_size = 11;   % number of epochs to be used to calculate the moving mean

for i = 1:size(data_raw_SAR,1)
    data_raw_smt_SAR(i,:) = movmean(data_raw_SAR(i,:), window_size, 'Endpoints', 'shrink');
end



%% STEP 3: Data Synchronization
% This step is meant to extract the environmental data at the same time
% stamps as the displacement data.

time_sync = time_SAR(time_SAR < tfin);

% Water level -------------------------------------------------------------
data_Hwater = zeros(1, length(time_sync));
for j = 1:length(time_sync)
    match = find(time_Hwater == time_sync(j));
    if ~isempty(match)
        data_Hwater(j) = data_raw_Hwater(match);
    else
        % if there's no correspondence, it takes the closest value
        [~, closest_idx] = min(abs(time_Hwater - time_sync(j)));
        data_Hwater(j) = data_raw_Hwater(closest_idx);
    end
end
clear match closest_idx

% Surface water temperature -----------------------------------------------
data_Twater = zeros(1, length(time_sync));
for j = 1:length(time_sync)
    match = find(time_Twater == time_sync(j));
    if ~isempty(match)
        data_Twater(j) = data_raw_Twater(match);
    else
        % if there's no correspondence, it takes the closest value
        [~, closest_idx] = min(abs(time_Twater - time_sync(j)));
        data_Twater(j) = data_raw_Twater(closest_idx);
    end
end
clear match closest_idx

% SAR displacement --------------------------------------------------------
data_orig_SAR = data_raw_SAR(:, 1:length(time_sync));
data_smt_SAR = data_raw_smt_SAR(:, 1:length(time_sync));

% Selection of the time-series to be used for the processing
if data_SAR_choice == 0
    data_SAR = data_orig_SAR;   % use the raw time-series
elseif data_SAR_choice == 1
    data_SAR = data_smt_SAR;    % use the smoothed time-series
end

% Plot --------------------------------------------------------------------
% Environmental data
figure
subplot(2,1,1)
plot(time_Twater + t0_SAR, data_raw_Twater, 'k')
hold on
plot(time_sync + t0_SAR, data_Twater, 'o', 'Color', [1 0 0]);
set(gca, 'fontsize', 12);
xlabel('Time [date]')
ylabel('Temperature [°C]')
title('Surface water temperature', 'FontSize', 18)
subplot(2,1,2)
plot(time_Hwater + t0_SAR, data_raw_Hwater, 'k')
hold on
plot(time_sync + t0_SAR, data_Hwater, 'o', 'Color', [0 0 1]);
set(gca, 'fontsize', 12);
xlabel('Time [date]')
ylabel('Height [m a.s.l.]')
title('Water level', 'FontSize', 18)

% SAR data
for i = 1:round(size(data_orig_SAR, 1)/10)
    figure
    k = 0;
    for j = ((i-1)*round(size(data_orig_SAR, 1)/round(size(data_orig_SAR, 1)/10))+1):min(i*round(size(data_orig_SAR, 1)/round(size(data_orig_SAR, 1)/10)),size(data_orig_SAR, 1))
        k = k + 1;
        subplot(5,2,k)
        plot(time_sync + t0_SAR, data_orig_SAR(j, :), 'k', 'LineWidth', 1.5);
        hold on
        plot(time_sync + t0_SAR, data_orig_SAR(j,:), 'o', 'Color', [205 41 144]/255, 'MarkerSize', 3.5, 'LineWidth', 1.3)
        plot(time_sync + t0_SAR, data_smt_SAR(j,:), '-', 'Color', [173 216 230]/255, 'MarkerSize', 3.5, 'LineWidth', 1)
        xlabel('Time [days]')
        ylabel('Displacement [mm]')
        ylim([-10 10])
        title(sprintf('Ascending LOS displacement PS %i', j+PS_shift), 'FontSize', 18, 'FontWeight', 'bold')
        set(gca, 'FontSize', 12);
    end
end



%% STEP 4: TIN
% This step is meant to interpolate with a TIN mesh the SAR displacement
% raw data and to plot the chosen one in the same plot as the environmental
% data.

% Convesion of PS coordinates from WGS84 to UTM 33N
[x_SAR, y_SAR, ~] = deg2utm(coord_SAR(:,2), coord_SAR(:,1));

% Creation of a matrix with the number of the closest PS vertices
T = delaunay(x_SAR, y_SAR);

% Creation of the TIN for all epochs
im = cell(length(time_sync), 1);

if plot_gifs == 1

    for i = 1:length(time_sync)
        
        fig = figure('units','normalized','outerposition',[0 0 1 1], 'Color', 'w');
        TO = triangulation(T, x_SAR, y_SAR, data_SAR(:,i));
        trisurf(TO)
        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Displacement [mm]');
        title(sprintf('TIN at epoch %d - Raw data', i), 'FontSize', 35);
        set(gca, 'FontSize', 15);
        colorbar;
        clim([-10 10]);
        zlim([-10 10]);
        pause(1)
        drawnow
        frame = getframe(fig);
        im{i} = frame2im(frame);
        close;
    
    end
    
    % Save it as a GIF
    filename = 'TIN_rawData_ASC.gif';
    for idx = 1:length(time_sync)
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',1,'DelayTime',0.3);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
        end
    end

end



%% STEP 5: Linear Correlation Coefficients
% This step is meant to explore and visualize the data we are working with,
% in order to assess the possible correlation between the two datasets:
% - environmental data = predictors.
% - SAR displacement data = response.

% Linear correlation coefficient ----------------------------------------------
% Water level - SAR displacement
corr_Hwater_SAR(1,:) = (1:size(data_SAR,1)) + PS_shift;
corr_Hwater_SAR(2,:) = corr(data_Hwater', data_SAR');
fprintf('\nLinear Correlation Coefficient between water level and SAR displacement: \n')
fprintf('PS %i, ρ: %.4f\n', corr_Hwater_SAR);

% Water surface temperature - SAR displacement
corr_Twater_SAR(1,:) = (1:size(data_SAR,1)) + PS_shift;
corr_Twater_SAR(2,:) = corr(data_Twater', data_SAR');
fprintf('\n\nLinear Correlation Coefficient between water surface temperature and SAR displacement: \n')
fprintf('PS %i, ρ: %.4f\n', corr_Twater_SAR);

% Water surface temperature - Water level
corr_Hwater_Twater = corr(data_Twater', data_Hwater');
fprintf('\n\nLinear Correlation Coefficient between water surface temperature and water level: \n')
fprintf('ρ: %.4f\n', corr_Hwater_Twater);



%% STEP 6: Choice of PS
% This step is meant to choose the PS for the next steps. For the selected 
% PS, the relation with the environmental data is assessed with 
% scatterplots. The procedure is eventually repeated to all the other PSs.

% Choice of the PS for the next steps -------------------------------------
% It is selected the point (PS) with the highest correlation coefficient.
PS_id = 1;   % corresponding to PS 24

% Plot --------------------------------------------------------------------
figure
subplot(3,1,1)
plot(time_sync + t0_SAR, data_orig_SAR(PS_id, :), '-ok', 'LineWidth', 1.0)
hold on
plot(time_sync + t0_SAR, data_smt_SAR(PS_id, :), '-m', 'LineWidth', 1.0)
xlabel('Time [days]')
ylabel('Displacement [mm]')
title(sprintf('Ascending LOS displacement PS %i', PS_id + PS_shift), 'FontSize', 20, 'FontWeight', 'bold')
set(gca, 'FontSize', 15);
legend('Raw data', 'Smooth data', 'FontSize', 15)
subplot(3,1,2)
plot(time_Twater + t0_SAR, data_raw_Twater, 'k')
hold on
plot(time_sync + t0_SAR, data_Twater, 'o', 'Color', [1 0 0]);
set(gca, 'fontsize', 12);
xlabel('Time [date]')
ylabel('Temperature [°C]')
title('Surface water temperature', 'FontSize', 18)
subplot(3,1,3)
plot(time_Hwater + t0_SAR, data_raw_Hwater, 'k')
hold on
plot(time_sync + t0_SAR, data_Hwater, 'o', 'Color', [0 0 1]);
set(gca, 'fontsize', 12);
xlabel('Time [date]')
ylabel('Height [m a.s.l.]')
title('Water level', 'FontSize', 18)

% Scatter plot
figure
[~, ~, BigAx, ~, ~] = plotmatrix([data_Twater', data_Hwater', data_SAR(PS_id,:)']);
xlabel(BigAx, 'Surface water temperature [°C]                                                               Water level [m]                                                                    SAR displacement [mm]', 'FontSize', 12);
ylabel(BigAx, 'SAR displacement [mm]                                   Water level [m]                              Surface water temperature [°C]', 'FontSize', 12);
title(BigAx,'Scatter plots between the variables', 'FontSize', 18)



%% STEP 7: Multiple Linear Regression 
% This step is meant to perform the multiple linear regression. We compare
% the results obtained between the manual implementation of 
% the problem and the matlab function 'regress'.

% Definition of linear regression model -----------------------------------
% We subtract the mean to every predictor to avoid the numerical
% instability and ease the invesrion of the design matrix A.
% By subtracing the mean of the predictor columns, we make them 
% stochastically independent from the consant term.

% Design matrix
A = [ones(length(time_sync),1), data_Hwater', data_Twater'];
A_norm = A - [0, mean(A(:,2)), mean(A(:,3))];

% Response vector
y_norm = (data_SAR(PS_id,:) - mean(data_SAR(PS_id,:)))';

% Normal matrix
N = A_norm' * A_norm;

% Useful matrices
N_inv = inv(N);
T_n = A_norm' * y_norm;

% Solution of the problem
x_est = N \ T_n;   % = N_inv * T_n
y_est = A_norm * x_est + mean(data_SAR(PS_id,:))';
x_est_ct = mean(data_SAR(PS_id,:)) - x_est(2)*mean(A(:,2)) - x_est(3)*mean(A(:,3));

% Covariance matrices
Cxx_est = N_inv;
Cyy_est = A_norm / N * A_norm';

% T-student test on the estimated parameters with confidence level alpha
alpha = 0.05;

t_x_est1 = x_est(1) / sqrt(Cxx_est(1,1));   % constant term
t_x_est2 = x_est(2) / sqrt(Cxx_est(2,2));   % Hwater term
t_x_est3 = x_est(3) / sqrt(Cxx_est(3,3));   % Twater term

t_lim = tinv(1 - alpha/2, size(A_norm,1) - size(A_norm,2));
fprintf('\n\nt-student test on the significance of the estimated parameters: \n')
fprintf('Constant term, t_obs: %.4f (t_lim: %.4f)\n', t_x_est1, t_lim);
fprintf('Hwater term, t_obs: %.4f (t_lim: %.4f)\n', t_x_est2, t_lim);
fprintf('Twater term, t_obs: %.4f (t_lim: %.4f)\n', t_x_est3, t_lim);

% Re-compensation of the least squares adjustment
A_norm_1 = A_norm(:, 2:3);
N_1 = A_norm_1' * A_norm_1;
x_est_1 = N_1 \ (A_norm_1' * y_norm);
y_est_1 = A_norm_1 * x_est_1 + mean(data_SAR(PS_id,:))';
Cxx_est_1 = inv(N_1);
Cyy_est_1 = A_norm_1 / N_1 * A_norm_1';

t_x_est1 = x_est_1(1) / sqrt(Cxx_est_1(1,1));   % Hwater term
t_x_est2 = x_est_1(2) / sqrt(Cxx_est_1(2,2));   % Twater term

t_lim = tinv(1 - alpha/2, size(A_norm_1,1) - size(A_norm_1,2));
fprintf('\n\nt-student test on the significance of the estimated parameters: \n')
fprintf('Hwater term, t_obs: %.4f (t_lim: %.4f)\n', t_x_est1, t_lim);
fprintf('Twater term, t_obs: %.4f (t_lim: %.4f)\n', t_x_est2, t_lim);

% Matlab function 'regress' -----------------------------------------------
[coeffs, coeffs_int, res, res_int, stats] = regress(data_SAR(PS_id,:)', [ones(length(data_Hwater),1), data_Hwater', data_Twater']);
data_pred_SAR = coeffs(1) + coeffs(2) * data_Hwater + 18.4*coeffs(3) * data_Twater;%18.4

% Plot --------------------------------------------------------------------
% The time series obtained with the two approaches are the same even if the
% estimated parameters are not. In particular, the constant term c is null
% in the manual implementation while equal to -608.37 from matlab regress
% function. This difference is due to the fact that we have defined the
% matrix A and vector y subtracting the mean to the predictors and to the 
% response; therefore, the vector of estimation y_est is obtained using 
% that design matrix (see the provided pdf for further explainations).
figure
plot(time_sync + t0_SAR, data_SAR(PS_id, :), 'k', 'LineWidth', 1.5);
hold on
plot(time_sync + t0_SAR, data_pred_SAR, 'b', 'LineWidth', 2);
plot(time_sync + t0_SAR, y_est_1, 'om', 'LineWidth', 1.5);
xlabel('Time [days]')
ylabel('Displacement [mm]')
title(sprintf('Ascending LOS displacement PS %i', PS_id + PS_shift), 'FontSize', 20, 'FontWeight', 'bold')
set(gca, 'FontSize', 15);
legend('SAR smoothed data', 'Matlab function', 'Manual implementation', 'FontSize', 20);

% Comparison between constant term from the two approaches
fprintf('\nTrue constant term from manual implementation: %.2f', x_est_ct);
fprintf('\nConstant term from Matlab function: %.2f', coeffs(1));



%% STEP 8: Auto Regression Preparation
% This step is meant to prepare the data to perform a multiple 
% linear regression with auto-regression terms. This means evaluating 
% the correlation length between predictors and response (in other words, 
% the time length over which there is no more correlation between 
% predictors and response).

% Correlation index as function of time lag -------------------------------
% It is plotted the correlation index as function of the time lag to decide
% how many previous instants have to be included in the auto regression 
% analysis.

% Vectors/matrices initialization
data_shift_Hwater = zeros(length(time_sync), size(data_Hwater,2));
data_shift_Twater = zeros(length(time_sync), size(data_Twater,2));
corr_Hwater_SAR_ts = zeros(1, length(time_sync));
corr_Twater_SAR_ts = zeros(1, length(time_sync));

% Correlation indices estimation
for j = 1:size(data_SAR,1)
    for i = 1:length(time_sync)
        data_shift_Hwater(i, i:end) = data_Hwater(1, 1:end-(i-1));
        data_shift_Twater(i, i:end) = data_Twater(1, 1:end-(i-1));
        corr_Hwater_SAR_ts(j,i) = corr(data_shift_Hwater(i,:)', data_SAR(j,:)');
        corr_Twater_SAR_ts(j,i) = corr(data_shift_Twater(i,:)', data_SAR(j,:)');
    end
end

% Average the correlation indices for all PSs
corr_Hwater_SAR_ts = mean(corr_Hwater_SAR_ts, 1);
corr_Twater_SAR_ts = mean(corr_Twater_SAR_ts, 1);

% First change of sign in the correlation index
time_lag = 1:length(time_sync);
tolerance = 1e-3;
if corr_Hwater_SAR_ts(1) > 0   % case first term positive
    first_corr_HS = find(corr_Hwater_SAR_ts < -tolerance , 1, 'first');
else                           % case first term negative
    first_corr_HS = find(corr_Hwater_SAR_ts > -tolerance , 1, 'first');
end
if corr_Twater_SAR_ts(1) > 0   % case first term positive
    first_corr_TS = find(corr_Twater_SAR_ts < -tolerance , 1, 'first');
else                           % case first term negative
    first_corr_TS = find(corr_Twater_SAR_ts > -tolerance , 1, 'first');
end

% Plot
figure
subplot(2,1,1)
plot(time_lag, corr_Hwater_SAR_ts, 'b', 'LineWidth', 1.5);
xline(time_lag(first_corr_HS), 'r', 'LineWidth', 1.5)
text(first_corr_HS, 0.21, sprintf('Time lag: %.0f', first_corr_HS-1), 'VerticalAlignment', 'top', 'Color', 'r', 'FontSize', 15);
yline(0, 'k', 'LineWidth', 1.8);
title('Correlation index evolution - Hwater vs SAR', 'FontSize', 20)
xlabel('Time lag [n]', 'FontSize', 15);
ylabel('Correlation index \rho [-]', 'FontSize', 15);
set(gca, 'FontSize', 15);
subplot(2,1,2)
plot(time_lag, corr_Twater_SAR_ts, 'b', 'LineWidth', 1.5);
xline(time_lag(first_corr_TS), 'r', 'LineWidth', 1.5);
text(first_corr_TS, 0.21, sprintf('Time lag: %.0f', first_corr_TS-1), 'VerticalAlignment', 'top', 'Color', 'r', 'FontSize', 15);
yline(0, 'k', 'LineWidth', 1.8);
title('Correlation index evolution - Twater vs SAR', 'FontSize', 20)
xlabel('Time lag [n]', 'FontSize', 15);
ylabel('Correlation index \rho [-]', 'FontSize', 15);
set(gca, 'FontSize', 15);

% Predictors' matrices for the auto regression ----------------------------
% We select a maximum correlation length of 2 epochs due to expected
% physics of the analyzed phenomenon.
% Extract the shifted time series up to the correlation length.
data_ar_Hwater = data_shift_Hwater(1:min([first_corr_HS, first_corr_TS, 3]), :);
data_ar_Twater = data_shift_Twater(1:min([first_corr_HS, first_corr_TS, 3]), :);

% Substitute zeros with NaNs
data_ar_Hwater(data_ar_Hwater == 0) = nan;
data_ar_Twater(data_ar_Twater == 0) = nan;

% Plot
figure
subplot(2,1,1)
plot(time_sync + t0_SAR, data_SAR(PS_id, :)-mean(data_SAR(PS_id, :)), 'k', 'LineWidth', 1.5);
hold on
plot(time_sync + t0_SAR, data_Hwater-mean(data_Hwater), 'b', 'LineWidth', 1.5);
plot(time_sync + t0_SAR, data_ar_Hwater-mean(data_Hwater), 'b');
title('SAR displacement - Hwater shifted', 'FontSize', 20)
xlabel('Time [date]', 'FontSize', 15);
ylabel('Zero mean displacement [mm]', 'FontSize', 15);
set(gca, 'FontSize', 15);
legend('SAR', 'Watel level')
subplot(2,1,2)
plot(time_sync + t0_SAR, data_SAR(PS_id, :)-mean(data_SAR(PS_id, :)), 'k', 'LineWidth', 1.5);
hold on
plot(time_sync + t0_SAR, data_Twater-mean(data_Twater), 'r', 'LineWidth', 1.5);
plot(time_sync + t0_SAR, data_ar_Twater-mean(data_Twater), 'r');
title('SAR displacement - Twater shifted', 'FontSize', 20)
xlabel('Time [date]', 'FontSize', 15);
ylabel('Zero mean variation [mm | °C]', 'FontSize', 15);
set(gca, 'FontSize', 15);
legend('SAR', 'Watel surface temperature')



%% STEP 9: Auto Regression
% This step is meant to perform a multiple linear regression, adding as 
% further predictors the values of water level and water surface
% temperature at the previous time stamps.

% Estimation of the parameters --------------------------------------------
% First, store the predictors and the response observations in a table
tbl_ar = table(data_ar_Hwater(1,:)', data_ar_Twater(1,:)', ...   % lag = 0 (original data)
       data_ar_Hwater(2,:)', data_ar_Twater(2,:)', ...           % lag = 1
       data_ar_Hwater(3,:)', data_ar_Twater(3,:)', ...           % lag = 2
       data_SAR(PS_id, :)', ...                                  % response
       'VariableNames',{'Hwater', 'Twater', ...                  % lag = 0 (original data)
       'ar_Hwater_1', 'ar_Twater_1', ...                         % lag = 1
       'ar_Hwater_2', 'ar_Twater_2', ...                         % lag = 2
       'SAR_displacement'});                                     % response

% Fit the linear regression model with all predictors
mdl_ar = fitlm(tbl_ar,['SAR_displacement ~ Hwater + Twater + ' ...
       'ar_Hwater_1 + ar_Twater_1 + ar_Hwater_2 + ar_Twater_2']);

% Display the results
% 'mdl' will provide you with the results of the:
% - standard deviation of the parameters' estimation (SE).
% - t-tests on the parameters (significance) of the auto regression (as for
%   the simple linear regression (tStat).
% - p-value for the t-statistic of the two-sided hypothesis test. The 
%   p-value is the likelihood of observing the current results or more 
%   extreme ones if the null hypothesis (the coefficient is equal to zero) 
%   is correct (pValue). If it is larger than 0.05 it therefore means that 
%   parameter is not significant and can be removed (i.e., there is a high
%   probability of observing other values larger than the threshold).
% - global F-test on the overall model computed comparing the chosen model
%   and a constant model (i.e., the mean of observations) (F-statistic).
disp(mdl_ar)

% t-student test on the parameters' significance --------------------------
% The t-value is given as output from the fitlm function.
t_obs_ar = mdl_ar.Coefficients.tStat;
t_lim_ar = tinv(1 - alpha/2, mdl_ar.DFE);
i = 1;

while any(abs(t_obs_ar(2:end)) < t_lim_ar)

    % Find the index of the predictor with the smallest t-statistic
    [~, min_idx] = min(abs(mdl_ar.Coefficients.tStat(2:end)));
    var_names_ar = mdl_ar.Formula.TermNames(2:end);
    var_names_ar(min_idx) = [];

    % Remove the predictor with the smallest t-statistic
    mdl_reduced = fitlm(tbl_ar, ['SAR_displacement ~ ', strjoin(var_names_ar, ' + ')]);
    t_obs_ar = mdl_reduced.Coefficients.tStat;
    mdl_ar = mdl_reduced;
    fprintf('\nIteration %.0f -----------------------------------------\n', i)
    disp(mdl_ar)
    i = i + 1;

end

% Re-computation of the prediction time series ----------------------------
data_pred_ar_SAR_PS = mdl_ar.Fitted';

% Plot --------------------------------------------------------------------
% Time-series comparison
figure
subplot(2,1,1)
plot(time_sync + t0_SAR, data_SAR(PS_id, :), 'k', 'LineWidth', 1.5);
hold on
plot(time_sync + t0_SAR, data_pred_SAR, 'g', 'LineWidth', 1.5);
plot(time_sync + t0_SAR, data_pred_ar_SAR_PS, 'm', 'LineWidth', 1.5);
title('LOS displacement time series comparison', 'FontSize', 20)
xlabel('Time [date]', 'FontSize', 15);
ylabel('Displacement [mm]', 'FontSize', 15);
set(gca, 'FontSize', 15);
legend('SAR raw', 'Linear regression', 'Auto regression')
subplot(2,1,2)
plot(time_sync + t0_SAR, data_pred_SAR - data_SAR(PS_id, :), 'g', 'LineWidth', 1.5);
hold on
plot(time_sync + t0_SAR, data_pred_ar_SAR_PS - data_SAR(PS_id, :), 'm', 'LineWidth', 1.5);
title('Residuals comparison', 'FontSize', 20)
xlabel('Time [date]', 'FontSize', 15);
ylabel('Residuals [mm]', 'FontSize', 15);
set(gca, 'FontSize', 15);
legend('Linear regression', 'Auto regression')

% Scatter plot
figure
scatter(data_SAR(PS_id, :), data_pred_ar_SAR_PS, 30, time_sync, 'filled')
title('Raw vs Predicted time series displacement', 'FontSize', 25)
xlabel('Raw displacement [mm]', 'FontSize', 15);
ylabel('Predicted displacement [mm]', 'FontSize', 15);
set(gca, 'FontSize', 15);
colorbar
title(colorbar, 'Time [days]');



%% STEP 10: Application of the linear regression model to all the other PSs
% This step is meant to replicate the procedure exploited for the single PS
% to all the others over the dam.

data_pred_ar_SAR = zeros(size(data_SAR,1), size(data_SAR,2));
coeffs_ar = zeros(size(tbl_ar,2), size(data_SAR,1));
R2_ar = zeros(2, size(data_SAR,1));

for k = 1:size(data_SAR,1)

    % Estimation of the parameters --------------------------------------------
    % First, store the predictors and the response observations in a table
    tbl_ar = table(data_ar_Hwater(1,:)', data_ar_Twater(1,:)', ...   % lag = 0 (original data)
           data_ar_Hwater(2,:)', data_ar_Twater(2,:)', ...           % lag = 1
           data_ar_Hwater(3,:)', data_ar_Twater(3,:)', ...           % lag = 2
           data_SAR(k, :)', ...                                  % response
           'VariableNames',{'Hwater', 'Twater', ...                  % lag = 0 (original data)
           'ar_Hwater_1', 'ar_Twater_1', ...                         % lag = 1
           'ar_Hwater_2', 'ar_Twater_2', ...                         % lag = 2
           'SAR_displacement'});                                     % response
    
    % Fit the linear regression model with all predictors
    mdl_ar = fitlm(tbl_ar,['SAR_displacement ~ Hwater + Twater + ' ...
           'ar_Hwater_1 + ar_Twater_1 + ar_Hwater_2 + ar_Twater_2']);
    
    % t-student test on the parameters' significance --------------------------
    % The t-value is given as output from the fitlm function.
    t_obs_ar = mdl_ar.Coefficients.tStat;
    t_lim_ar = tinv(1 - alpha/2, mdl_ar.DFE);
    
    if all(abs(t_obs_ar) > t_lim_ar)

        data_pred_ar_SAR(k,:) = zeros(1, length(time_sync));
        fprintf('\nAll the auto regression parameters for PS %i are not significant -> the parameter is removed!\n', k + PS_shift);

    else
        
        status = 0;

        % Perform the test: if the observed value is larger than the threshold one
        % the parameters is kept, otherwise its is put equal to zero
        while any(abs(t_obs_ar(2:end)) < t_lim_ar)
    
            % Find the index of the predictor with the smallest t-statistic
            [~, min_idx] = min(abs(mdl_ar.Coefficients.tStat(2:end)));
            var_names_ar = mdl_ar.Formula.TermNames(2:end);
            var_names_ar(min_idx) = [];

            if isempty(var_names_ar)

                data_pred_ar_SAR(k,:) = zeros(1, length(time_sync));
                fprintf('\nAll the auto regression parameters for PS %i are not significant -> the parameter is removed!\n', k + PS_shift);
                status = 1;
                break
                
            end
        
            % Remove the predictor with the smallest t-statistic
            mdl_reduced = fitlm(tbl_ar, ['SAR_displacement ~ ', strjoin(var_names_ar, ' + ')]);
            t_obs_ar = mdl_reduced.Coefficients.tStat;
            mdl_ar = mdl_reduced;
        
        end
        
        if status ~= 1

            data_pred_ar_SAR(k,:) = mdl_ar.Fitted';
            fieldname = strcat('PS', num2str(k + PS_shift));
            coeffs_nv_ar.(fieldname){:,1} = var_names_ar;
            coeffs_nv_ar.(fieldname){:,2} = mdl_ar.Coefficients.Estimate(2:end);
            coeffs_ar(1,k) = mdl_ar.Coefficients.Estimate(1);
            R2_ar(:,k) = [k + PS_shift; mdl_ar.Rsquared.Adjusted];

        end

    end

end

% Plot --------------------------------------------------------------------
for i = 1:round(size(data_pred_ar_SAR, 1)/10)
    figure
    k = 0;
    for j = ((i-1)*round(size(data_pred_ar_SAR, 1)/round(size(data_pred_ar_SAR, 1)/10))+1):min(i*round(size(data_SAR, 1)/round(size(data_SAR, 1)/10)),size(data_SAR, 1))
        k = k + 1;
        subplot(5,2,k)
        plot(time_sync + t0_SAR, data_SAR(j, :), 'k', 'LineWidth', 1.5);
        hold on
        plot(time_sync + t0_SAR, data_pred_ar_SAR(j,:), '-o', 'Color', [205 41 144]/255, 'MarkerSize', 2.5, 'LineWidth', 1.3)
        xlabel('Time [days]')
        ylabel('Displacement [mm]')
        ylim([-5 5])
        title(sprintf('Ascending LOS displacement PS %i', j+PS_shift), 'FontSize', 18, 'FontWeight', 'bold')
        set(gca, 'FontSize', 12);
    end
end

% Indices of the correlated PSs -------------------------------------------
cor_idx_SAR = find(all(data_pred_ar_SAR ~= 0, 2));

% The rows/colums of the PSs for which no significant parameters are found
% are removed from the corresponding matrices/vectors.
x_cor_SAR = x_SAR(cor_idx_SAR);                        % coordinate X
y_cor_SAR = y_SAR(cor_idx_SAR);                        % coordinate Y
data_pred_ar_SAR = data_pred_ar_SAR(cor_idx_SAR, :);   % auto regression time series
coeffs_ar = coeffs_ar(:, cor_idx_SAR);                 % auto regression coefficients
R2_ar = R2_ar(:, cor_idx_SAR);                         % adjusted R-squared

% Display the adjusted R-squared value
fprintf('\n\nAdjusted R-squared value after auto regression: \n')
fprintf('PS %i, R-squared: %.4f\n', R2_ar);


%% STEP 11: TIN from auto regression results
% This step is meant to interpolate with a TIN mesh the SAR displacement
% raw data and to plot the chosen one in the same plot as the environmental
% data.

% Creation of the TIN for all epochs
T = delaunay(x_cor_SAR, y_cor_SAR);
im = cell(length(time_sync)-(size(data_ar_Hwater,1)-1), 1);

if plot_gifs == 1

    for i = size(data_ar_Hwater,1):length(time_sync)
        
        fig = figure('units','normalized','outerposition',[0 0 1 1], 'Color', 'w');
        TO = triangulation(T, x_cor_SAR, y_cor_SAR, data_pred_ar_SAR(:,i));
        trisurf(TO)
        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Displacement [mm]');
        title(sprintf('TIN at epoch %d - AR prediction', i), 'FontSize', 35);
        set(gca, 'FontSize', 15);
        colorbar;
        clim([-5 5]);
        zlim([-5 5]);
        pause(1)
        drawnow
        frame = getframe(fig);
        im{i-(size(data_ar_Hwater,1)-1)} = frame2im(frame);
        close;
    
    end
    
    % Save it as a GIF
    filename = 'TIN_ar_ASC.gif';
    for idx = 1:length(time_sync)-(size(data_ar_Hwater,1)-1)
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',1,'DelayTime',0.3);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
        end
    end

end



%% STEP 12: Spatial interpolation of displacements for all PSs at each epoch
% This step is meant to spatially interpolate the estimated parameters at 
% each epoch with a Lagrange function.

% Data preparation --------------------------------------------------------
ref_coeffs_names_ar = mdl_ar.VariableNames(1:end-1);
coeffs_nv_ar_fn = fieldnames(coeffs_nv_ar);

% Extrapolation of the parameters of all the PS time series from auto
% regression analysis. They are store in the variable coeffs_ar.
for i = 1:size(coeffs_ar,2)

    fieldname = coeffs_nv_ar_fn{i};
    
    for j = 1:size(coeffs_nv_ar.(fieldname){1,1},1)
    
        idx = find(strcmp(coeffs_nv_ar.(fieldname){1,1}{j,1}, ref_coeffs_names_ar));
        idx = idx + 1;
        coeffs_ar(idx,i) = coeffs_nv_ar.(fieldname){1,2}(j);

    end

end

% Import a GeoTIFF for the dam area ---------------------------------------
dam_tif = 'DAM_Nocelle.tif';
[area_dam, R] = readgeoraster(dam_tif);

% Create a grid with the coordinates of the dam area
[rows, cols] = meshgrid(1:size(area_dam,2), 1:size(area_dam,1));
[x_dam, y_dam] = R.intrinsicToWorld(rows,cols);

% Interpolation of the parameters over the dam area -----------------------
% It is used a TIN mesh. In the variable coeffs_ar there are stored the
% interpolated values of a parameter over the area of the dam. Each 'page'
% of the variable correspond to a different parameter, in the same order as
% defined in the auto regression model: Intercept, Hwater, Twater,
% Hwater_1, Twater_1, ... and so on. When a parameter is not significant,
% it si set equal to zero.
coeffs_ar_TIN = zeros(size(x_dam,1), size(x_dam,2), size(coeffs_ar,1));
ref_coeffs_names_ar_full = cell(size(ref_coeffs_names_ar,1)+1,1);
ref_coeffs_names_ar_full(1,1) = cellstr('Intercept');
ref_coeffs_names_ar_full(2:end,1) = ref_coeffs_names_ar;

for i = 1:size(coeffs_ar,1)

    F = scatteredInterpolant(x_cor_SAR, y_cor_SAR, coeffs_ar(i,:)', 'linear');
    coeffs_ar_tmp = F(x_dam(:), y_dam(:));
    coeffs_ar_TIN(:,:,i) = reshape(coeffs_ar_tmp, size(x_dam,1), size(x_dam,2));
    clear coeffs_ar_tmp

    coeffs_ar_TIN(:,:,i) = (coeffs_ar_TIN(:,:,i)) .* area_dam;
    coeffs_ar_TIN(area_dam == 0) = nan;
    x_dam1 = x_dam .* area_dam;
    x_dam1(area_dam == 0) = nan;
    y_dam1 = y_dam .* area_dam;
    y_dam1(area_dam == 0) = nan;

    parameter_name = ref_coeffs_names_ar_full{i,1};
    parameter_name = regexprep(parameter_name, '_', ' ');

    figure
    plot3(x_cor_SAR, y_cor_SAR, coeffs_ar(i,:)','mo')
    hold on
    mesh(x_dam1, y_dam1, coeffs_ar_TIN(:,:,i), 'FaceColor', 'interp', 'FaceAlpha', 0.8)
    xlabel('X [m]', 'FontSize', 20)
    ylabel('Y [m]', 'FontSize', 20)
    zlabel('Value [-]', 'FontSize', 20)
    title(sprintf('Interpolation of parameter "%s" over the dam', parameter_name), 'FontSize', 25)
    legend('Sample Points', 'Interpolated Surface', 'FontSize', 15, 'Location', 'NorthEast')
    colorbar
    set(gca, 'FontSize', 15)

end

% LOS displacement at each coordinate of the grid -------------------------
% The displacement is obtained multiplying each coefficient for the
% corresponding prediction time series (i.e., at the different time lags).
displ_dam = zeros(size(x_dam,1), size(x_dam,2), length(time_sync));
for i = 1:length(time_sync)
    displ_dam(:,:,i) = coeffs_ar_TIN(:,:,1) +  ...
        coeffs_ar_TIN(:,:,2) .* data_ar_Hwater(1,i) + coeffs_ar_TIN(:,:,3) .* data_ar_Twater(1,i) + ...
        coeffs_ar_TIN(:,:,4) .* data_ar_Hwater(2,i) + coeffs_ar_TIN(:,:,5) .* data_ar_Twater(2,i) + ...
        coeffs_ar_TIN(:,:,6) .* data_ar_Hwater(3,i) + coeffs_ar_TIN(:,:,7) .* data_ar_Twater(3,i);
end

% Plot --------------------------------------------------------------------
if plot_gifs == 1

    im = cell(length(time_sync)-(size(data_ar_Hwater,1)-1), 1);
    
    for i = size(data_ar_Hwater,1):length(time_sync)
    
        fig = figure('units','normalized','outerposition',[0 0 1 1], 'Color', 'w');
        surf(x_dam, y_dam, reshape(displ_dam(:,:,i), size(displ_dam,1), size(displ_dam,2), 1));
        shading interp;
        title(sprintf('LOS deformation model at epoch %d', i), 'FontSize', 35);
        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Displacement [mm]');
        set(gca, 'FontSize', 15);
        colorbar;
        clim([-5 5]);
        zlim([-15 10]);
        view(2);
        pause(1)
        drawnow
        frame = getframe(fig);
        im{i-(size(data_ar_Hwater,1)-1)} = frame2im(frame);
        close;
    
    end
    
    % Save it as a GIF
    filename = 'Dam_LOS_displ_ASC.gif';
    for idx = 1:(size(data_ar_Hwater,1)-1)
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',1,'DelayTime',0.3);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
        end
    end

end



%% STEP 13: Final comparison on the chosen PS
% This step is meant to perform the comparison of the displacement 
% time-series obtained through at the different steps:
% - raw data
% - smoothed data
% - multivariate linear regression
% - multivariate auto regression
% - spatial-temporal multivariate auto regression

% Extract PS_id coordinates -----------------------------------------------
xy_PS_id = [x_SAR(PS_id), y_SAR(PS_id)];

% Calculate the Euclidean distance
distances = sqrt((x_dam - xy_PS_id(1)).^2 + (y_dam - xy_PS_id(2)).^2);

% Find the row and column indices of the point with the minimum distance
[minDistance, linearIndex] = min(distances(:));
[rowIndex, colIndex] = ind2sub(size(distances), linearIndex);

% Plot --------------------------------------------------------------------
figure
plot(time_sync + t0_SAR, data_orig_SAR(PS_id, :), 'k', 'LineWidth', 1);
hold on
plot(time_sync + t0_SAR, data_smt_SAR(PS_id, :), 'k', 'LineWidth', 2);
plot(time_sync + t0_SAR, data_pred_SAR, '--g', 'LineWidth', 1.5);
plot(time_sync + t0_SAR, data_pred_ar_SAR(PS_id, :), '--m', 'LineWidth', 1.5);
plot(time_sync + t0_SAR, reshape(displ_dam(rowIndex, colIndex,:), 1, length(time_sync)), '--b', 'LineWidth', 1.5);
title(sprintf('LOS Displacement Comparison for PS %i', PS_id + PS_shift), 'FontSize', 20);
ylabel('Displacement [mm]', 'FontSize', 15);
xlabel('Time [date]', 'FontSize', 15);
set(gca, 'FontSize', 15)
legend('Raw data', 'Smoothed data', 'MV linear regression', 'MV auto regression', 'Space-Time interp.')

% The variable displ_dam contains the LOS displacement values for all the
% points (coordinates) of the dam. So, you can simply evaluate this matrix
% in any point to see its corresponding displacement time-series.

fprintf('\nScript execution completed!\n');












close all


%%
% % Assuming you have radar data, manual regression (y_est_1), and automatic regression (data_pred_SAR)
% radar_data = data_SAR(PS_id, :);
% manual_regression = y_est_1;
% auto_regression = data_pred_SAR;
% 
% % Calculate residuals
% residuals_manual = radar_data - manual_regression;
% residuals_auto = radar_data - auto_regression;
% 
% % Calculate Bhattacharyya distance
% bc_distance_manual = bhattacharyya_distance(residuals_manual);
% bc_distance_auto = bhattacharyya_distance(residuals_auto);
% 
% fprintf('Bhattacharyya Distance between Radar Data and Manual Regression: %.4f\n', bc_distance_manual);
% fprintf('Bhattacharyya Distance between Radar Data and Auto Regression: %.4f\n', bc_distance_auto);
% 
% % Function to calculate Bhattacharyya Distance
% function bc_distance = bhattacharyya_distance(residuals)
%     % Assuming residuals follow a Gaussian distribution
%     mean_residuals = mean(residuals);
%     cov_residuals = cov(residuals);
% 
%     % Bhattacharyya Coefficient
%     BC = 1/8 * mean_residuals * inv(cov_residuals) * mean_residuals' + 0.5 * log(det(cov_residuals) / sqrt(det(cov_residuals)));
% 
%     % Bhattacharyya Distance
%     bc_distance = -log(BC);
% end


%%
% Assume you have radar data, Twater and Hwater
radar_data = data_SAR(PS_id, :);
manual_regression = y_est_1;
automatic_regression = data_pred_SAR;
% manual_regression = data_Twater;
% automatic_regression = data_Hwater;

% Calculate Bhattacharyya distance
bc_distance_radar_data_Twater = abs(bhattacharyya_distance(radar_data(:), data_Twater(:)));
bc_distance_radar_data_Hwater = abs(bhattacharyya_distance(radar_data(:), data_Hwater(:)));
bc_distance_radar_manual = abs(bhattacharyya_distance(radar_data(:), manual_regression(:)));
bc_distance_radar_auto = abs(bhattacharyya_distance(radar_data(:), automatic_regression(:)));

fprintf('Bhattacharyya Distance between Radar Data and Tempreture water: %.4f\n', bc_distance_radar_data_Twater);
fprintf('Bhattacharyya Distance between Radar Data and Hight of the water: %.4f\n', bc_distance_radar_data_Hwater);
fprintf('Bhattacharyya Distance between Radar Data and Manual Regression: %.4f\n', bc_distance_radar_manual);
fprintf('Bhattacharyya Distance between Radar Data and Auto Regression: %.4f\n', bc_distance_radar_auto);

% Function to calculate Bhattacharyya Distance
function bc_distance = bhattacharyya_distance(P, Q)
    % P and Q are vectors of residuals or probability distributions
    % Calculate mean and covariance matrix
    mean_P = mean(P);
    mean_Q = mean(Q);
    cov_P = cov(P);
    cov_Q = cov(Q);
    
    % Calculate Bhattacharyya Coefficient
    BC = 1/8 * (mean_Q - mean_P) * inv(0.5 * (cov_P + cov_Q)) * (mean_Q - mean_P)' + 0.5 * log(det(0.5 * (cov_P + cov_Q)) / sqrt(det(cov_P * cov_Q)));
    
    % Bhattacharyya Distance
    bc_distance = -log(BC);
end




