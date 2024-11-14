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
% |                              COLLOCATION                              |
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

% Residuals ---------------------------------------------------------------
% These are the residuals computed at the end of Lab_03.m, obtained as
% difference between the observations and the estimated harmonic
% (considering just the fundamental frequency of vibration).
load('Fourier_out.mat');

% ! The residuals have a mean that is different from zero, because the
% ! Fourier analysis has been performed on the regularized time series 
% ! (which means constant time step of 6 days) while for collocation there
% ! are considered only the residuals from the observed dates.

% PSs definition ----------------------------------------------------------
PS_shift = 21;   % starting PS id
PS_id = 3;       % chosen PS

% Coordinates' conversion -------------------------------------------------
% Convesion of PS coordinates from WGS84 to UTM 33N.
[x_SAR, y_SAR, ~] = deg2utm(coord_SAR(:,2), coord_SAR(:,1));

% Import a GeoTIFF for the dam area ---------------------------------------
dam_tif = 'DAM_Nocelle.tif';
[area_dam, R] = readgeoraster(dam_tif);

% Create a grid with the coordinates of the dam area
[rows, cols] = meshgrid(1:size(area_dam,2), 1:size(area_dam,1));
[x_dam, y_dam] = R.intrinsicToWorld(rows,cols);



%% STEP 2: Empirical Covariance
% This step is meant to compute the empirical covariance point cloud and
% average for one displacement time series. This means computing the
% empirical covariance for all the couples of distances (in time). We can
% do this using the f1DEmpCovEst.m matlab function.
[tGrid, eCovF, Cecf, h] = f1DEmpCovEst(v_SAR(PS_id,:)', time_SAR', 6, 2);

% Plot
figure
plot(tGrid, eCovF, '.-', 'LineWidth', 1.3, 'MarkerSize', 12);
xlabel('Time lag [days]', 'FontSize', 15); ylabel('Covariance [mm^2]', 'FontSize', 15); 
title(sprintf('Covariance modelling for PS %i', PS_id+PS_shift), 'FontSize', 20);
xlim([0 tGrid(end)/2]);
set(gca, 'FontSize', 15);



%% STEP 3: Covariance Modelling
% This step is meant to interpolate the empirical covariance with a
% mathematical function.
% Two models are proposed:
%
%              C(tau) = B * exp(-b * tau)  % exponential model
%                          or
%              C(tau) = B * exp(-b * tau^2)  % gaussian model
%
% The parameters B and b are estimated by applying linear least squares
% inverting the following relationship:
%
%              log(C) = log(B) - b * tau  % exponential model
%                           or
%              log(C) = log(B) - b * tau^2  % gaussian model
%

% Exponential model -------------------------------------------------------
idx = find(eCovF < 0, 1, 'first') - 1;              % find last non-negative non element in the empirical covariance
A  = [ones(size(tGrid(2:idx))) -tGrid(2:idx)];      % design matrix
yo = log(eCovF(2:idx));                             % observed covariance (log of the empirical one)
pe = inv(A'*A) * A'*yo;                             % parameters estimation
pe(1) = exp(pe(1));                                 % amplitude (B)
pe(2) = pe(2);                                      % exponent scale (b)
fcove = @(tau) pe(1) * exp(-pe(2) * tau);           % declare a function to compute the covariance model given the time lag

% Add the model to the empirical covariance plot
hold on;
plot(tGrid, fcove(tGrid), '-', 'LineWidth', 2);

% Gaussian model ----------------------------------------------------------
idx = find(eCovF < 0, 1, 'first') - 1;              % find last non-negative non element in the empirical covariance
A  = [ones(size(tGrid(2:idx))) -tGrid(2:idx).^2];   % design matrix
yo = log(eCovF(2:idx));                             % observed covariance (log of the empirical one)
pg = inv(A'*A) * A'*yo;                             % parameters estimation
pg(1) = exp(pg(1));                                 % amplitude (B)
pg(2) = pg(2);                                      % exponent scale (b)
fcovg = @(tau) pg(1) * exp(-pg(2) * tau.^2);        % declare a function to compute the covariance model given the time lag

% Add the model to the empirical covariance plot
hold on;
plot(tGrid, fcovg(tGrid), '-', 'LineWidth', 2);

% Check the "nugget" amplitude --------------------------------------------
% To check if it has the same order of magnitude of the observation error.
% Exponential
sv_exp = sqrt(eCovF(1) - fcove(0));
fprintf('\nNugget from Exponential modelling %.4f\n', sv_exp);

% Gaussian
sv_gau = sqrt(eCovF(1) - fcovg(0));
fprintf('Nugget from Gaussian modelling %.4f\n', sv_gau);

% A-priori (2 mm)
s2v_ap = 2;
fprintf('A-priori noise amplitude (std) %.4f\n', sqrt(s2v_ap));

% Choice of the covariance model
fcov = fcovg;
s2v = sv_exp^2;


%% STEP 4: Understanding the Model's Parameters
% This step is meant to understand which is the physical meaning of the
% estimated parameters, considering a certain mathematical modelling of the
% empirical covariance function. We consider an exponential model.

% Exponential model - Amplitude x 10
pe_1(1) = pe(1) * 10;                                   % amplitude (B)
pe_1(2) = pe(2);                                        % exponent scale (b)
fcove_1 = @(tau) pe_1(1) * exp(-pe_1(2) * tau);         % declare a function to compute the covariance model given the time lag

% Plot
figure
plot(tGrid, eCovF, '.-', 'LineWidth', 1.3, 'MarkerSize', 12);
hold on; plot(tGrid, fcove(tGrid), '-', 'LineWidth', 2);
xlabel('Time lag [days]', 'FontSize', 15); ylabel('Covariance [mm^2]', 'FontSize', 15); 
title(sprintf('Different modelling for PS %i', PS_id+PS_shift), 'FontSize', 20);
xlim([0 tGrid(end)/2]);
set(gca, 'FontSize', 15);

hold on;
plot(tGrid, fcove_1(tGrid), '-', 'LineWidth', 2);


% Exponential model - Amplitude / 10
pe_2(1) = pe(1) / 10;                     
pe_2(2) = pe(2);                                       
fcove_2 = @(tau) pe_2(1) * exp(-pe_2(2) * tau);

hold on;
plot(tGrid, fcove_2(tGrid), '-', 'LineWidth', 2);


% Exponential model - Exponential scale * 5
pe_3(1) = pe(1);                     
pe_3(2) = pe(2) * 5;                                       
fcove_3 = @(tau) pe_3(1) * exp(-pe_3(2) * tau);

hold on;
plot(tGrid, fcove_3(tGrid), '-', 'LineWidth', 2);


% Exponential model - Exponential scale / 5
pe_4(1) = pe(1);                     
pe_4(2) = pe(2) / 5;                                       
fcove_4 = @(tau) pe_4(1) * exp(-pe_4(2) * tau);

hold on;
plot(tGrid, fcove_4(tGrid), '-', 'LineWidth', 2);
legend('Emp. cov.', 'LS sol.', 'A x 10', 'A / 10', 'a x 5', 'a / 5')



%% STEP 5: Semi-Variogram
% This step is meant to compute the semi-variogram. It is obtained as:
% γ(τ) = c(0) - c(τ) + nugget
gamma = (fcov(0) - fcov(tGrid)) + sv_exp^2;

% Plot
figure
plot(tGrid, gamma, '.-', 'LineWidth', 1.3, 'MarkerSize', 12);
xlabel('Time lag [days]', 'FontSize', 15); ylabel('Semivariance [mm^2]', 'FontSize', 15); 
title(sprintf('Semivariogram for PS %i', PS_id+PS_shift), 'FontSize', 20);
xlim([0 tGrid(end)/2]);
set(gca, 'FontSize', 15);



%% STEP 6: Collocation
% This step is meant to compute the covariance matrices needed to aplly the
% collocation formula to estimate the residuals.

% Time of estimation
% We want to estimate the residuals in correspondence of the observations,
% so t_est is equal to time_SAR.
t_est = time_SAR';

% Estimation-observation covariance matrix
t1 = repmat(t_est, 1, length(time_SAR));
t2 = repmat(time_SAR, length(t_est), 1);
Cv_vo = fcov(abs(t1 - t2));

% Observation-observation covariance matrix
t3 = repmat(time_SAR', 1, length(time_SAR));
t4 = repmat(time_SAR, length(time_SAR), 1);
Cvo_vo = fcov(abs(t3 - t4));

% Collocation estimate
v_est =  Cv_vo * inv(Cvo_vo + s2v .* eye(length(time_SAR))) * v_SAR(PS_id,:)';
data_SAR_est  = data_SAR_Fourier(PS_id,:) + v_est'; 

% Plot
figure; plot(time_SAR, v_SAR(PS_id,:), '.-');
hold on; plot(t_est, v_est, 'LineWidth', 3);
title('Estimated residuals [mm]');



%% STEP 7: Error of the Prediction
% This step is meant to estimate the prediction/estimation error of the
% collocation approach.
v_var = fcov(0) - diag(Cv_vo * (Cvo_vo + s2v * eye(length(time_SAR)))^(-1) * Cv_vo');
v_std = sqrt(v_var);


%% STEP 8: Final Signal
% This step is meant to obtain the final estimated signal, by summing up
% the contribution due to the deterministic trend (Fourier harmonic) and
% the stochastic one, from collocation.

figure
plot(time_SAR, data_SAR(PS_id,:), 'k', 'LineWidth', 1)
hold on
errorbar(time_SAR, data_SAR_est, v_std, 'b', 'LineWidth', 0.7)
plot(time_SAR, data_SAR_est, 'b', 'LineWidth', 2)
xlabel('Time [days]')
ylabel('Displacement [mm]')
title('SAR displacement', 'FontSize', 20)
set(gca, 'FontSize', 15)
legend('Observations', 'Estimations', 'FontSize', 15)



%% STEP 9: Prediction - Hole in the Data
% This step is meant to show how the collocation works when there is a hole
% in the data. If the hole is shorter than the correlation length, the
% signal will be predicted for the entire hole. If, instead, the hole is
% beigger than the correlation length, the signal will be predicted only up
% to the correlation length.
t_obs = time_SAR';
t_est  = time_SAR';
n_est = length(t_est);
y_obs = data_SAR(PS_id,:)';
n_obs = length(t_obs);
trend = data_SAR_Fourier(PS_id,:)';

% Introduce a hole in the data
idx2 = [1:50 100:n_obs]';
% idx2 = [1:50 65 80:n_obs]';
y_obs2 = y_obs(idx2);
t_obs2 = t_obs(idx2);
n_obs2 = length(t_obs2);
trend2 = data_SAR_Fourier(PS_id,idx2)';
yres_obs2 = y_obs2 - trend2;

T_est  = repmat(t_est, 1, n_obs2);
T_obs = repmat(t_obs2', n_est, 1); 
C_est_obs = fcov(abs(T_est - T_obs));

T_obs = repmat(t_obs2', n_obs2, 1); 
C_obs_obs = fcov(abs(T_obs' - T_obs));

v_est2 = C_est_obs * (C_obs_obs + s2v * eye(n_obs2))^(-1) * yres_obs2;  % filtered signal
v_var2 = fcov(0) - diag(C_est_obs * (C_obs_obs + s2v * eye(n_obs2))^(-1) * C_est_obs');

% Plot
figure; hold on;
plot(t_obs, y_obs, '--.k', 'LineWidth', 1.5);
plot(t_obs2, y_obs2, '.', 'MarkerSize', 15);
plot(t_est, trend, '-', 'LineWidth', 1.5);
plot(t_est, v_est2 + trend, '-', 'LineWidth', 1.5)
scatter(t_est, v_est2 + trend, 40, sqrt(v_var2), 'filled');
legend('All Observations', 'Used Observations', 'Trend', 'Estimation');
title('Prediction');
hc = colorbar;
title(hc, '[mm]');
xlabel('Time [days]');
ylabel('Displacement [mm]');
set(gca, 'FontSize', 20);



%% STEP 10: Extrapolation - After the end of observations
% This step is meant to use collocation to predict the signal after the
% last observation. The prediction will last only for the correlation
% length.
t_est  = time_SAR';
n_est = length(t_est);

% Remove the last n observations
idx3 = (1:n_obs-20)';
y_obs3 = y_obs(idx3);
t_obs3 = t_obs(idx3);
n_obs3 = length(t_obs3);
trend3 = data_SAR_Fourier(PS_id,idx3)';
yres_obs3 = y_obs3 - trend3;

T_est  = repmat(t_est, 1, n_obs3);
T_obs = repmat(t_obs3', n_est, 1); 
C_est_obs = fcov(abs(T_est - T_obs));

T_obs = repmat(t_obs3', n_obs3, 1); 
C_obs_obs = fcov(abs(T_obs' - T_obs));

v_est3 = C_est_obs * (C_obs_obs + s2v * eye(n_obs3))^(-1) * yres_obs3;  % filtered signal
v_var3 = fcov(0) - diag(C_est_obs * (C_obs_obs + s2v * eye(n_obs3))^(-1) * C_est_obs');

% Plot
figure; hold on;
plot(t_obs, y_obs, '--.k', 'LineWidth', 1.5);
plot(t_obs3, y_obs3, '.', 'MarkerSize', 15);
plot(t_est, trend, '-', 'LineWidth', 1.5);
plot(t_est, v_est3 + trend, '-', 'LineWidth', 1.5)
scatter(t_est, v_est3 + trend, 40, sqrt(v_var3), 'filled');
legend('All Observations', 'Used Observations', 'Trend', 'Estimation');
title('Prediction');
hc = colorbar;
title(hc, '[mm]');
xlabel('Time [days]');
ylabel('Displacement [mm]');
set(gca, 'FontSize', 20);



%% STEP 11: Least Squares Collocation
% This step is meant to perform the least squares collocation, which means
% estimating the parameters of the deterministic trend (our harmonic)
% considering as covariance matrix of the observations Q + Cvo_vo instead
% of just Q. This means introducing a relathionship between the
% observations (they are no longer independent - Q is not diagonal).

% Least squares for harmonic coefficients estimation ----------------------
A = [ones(size(time_SAR'))  2*cos(2*pi*freq_PS(102)*time_SAR')  2*sin(2*pi*freq_PS(102)*time_SAR')];

x_est = (A' * inv(Cvo_vo+s2v*eye(length(time_SAR))) * A)^(-1) * A' * inv(Cvo_vo+s2v*eye(length(time_SAR))) * data_SAR(PS_id,:)';
y_est = A * x_est;
v_SAR_lsc = data_SAR(PS_id,:) - y_est';

figure
plot(time_SAR, data_SAR(PS_id,:))
hold on; plot(time_SAR, y_est)
title('Predicted signal')

figure
plot(v_SAR_lsc)
hold on; plot(v_SAR(PS_id,:))
title('Residuals comparison')
legend('LSC', 'Fourier')


% Covariance --------------------------------------------------------------
[tGrid, eCovF, Cecf, h] = f1DEmpCovEst(v_SAR(PS_id,:)', time_SAR', 6, 2);

% Exponential model
idx = find(eCovF < 0, 1, 'first') - 1;              % find last non-negative non element in the empirical covariance
A  = [ones(size(tGrid(2:idx))) -tGrid(2:idx)];      % design matrix
yo = log(eCovF(2:idx));                             % observed covariance (log of the empirical one)
pe = inv(A'*A) * A'*yo;                             % parameters estimation
pe(1) = exp(pe(1));                                 % amplitude (B)
pe(2) = pe(2);                                      % exponent scale (b)
fcov1 = @(tau) pe(1) * exp(-pe(2) * tau);           % declare a function to compute the covariance model given the time lag
s2v1 = eCovF(1) - fcov1(0);                         % nugget amplitude

% Plot
figure
plot(tGrid, eCovF, '.-', 'LineWidth', 1.3, 'MarkerSize', 12);
hold on; plot(tGrid, fcov1(tGrid), '-', 'LineWidth', 2);
xlabel('Time lag [days]', 'FontSize', 15); ylabel('Covariance [mm^2]', 'FontSize', 15); 
title(sprintf('Covariance modelling for PS %i', PS_id+PS_shift), 'FontSize', 20);
xlim([0 tGrid(end)/2]);
set(gca, 'FontSize', 15);


% Collocation -------------------------------------------------------------
% Time of estimation
t_est = time_SAR';

% Estimation-observation covariance matrix
t1 = repmat(t_est, 1, length(time_SAR));
t2 = repmat(time_SAR, length(t_est), 1);
Cv_vo1 = fcov1(abs(t1 - t2));

% Observation-observation covariance matrix
t3 = repmat(time_SAR', 1, length(time_SAR));
t4 = repmat(time_SAR, length(time_SAR), 1);
Cvo_vo1 = fcov1(abs(t3 - t4));

% Collocation estimate
v_est1 =  Cv_vo1 * inv(Cvo_vo1 + s2v1 .* eye(length(time_SAR))) * v_SAR(PS_id,:)';
data_SAR_est1  = data_SAR_Fourier(PS_id,:) + v_est1'; 

% Plot
figure; plot(time_SAR, v_SAR(PS_id,:), '.-');
hold on; plot(t_est, v_est1, 'LineWidth', 3);
title('Estimated residuals [mm]');

% Prediction error
v_var1 = fcov1(0) - diag(Cv_vo1 * (Cvo_vo1 + s2v1 * eye(length(time_SAR)))^(-1) * Cv_vo1');
v_std1 = sqrt(v_var1);

% Plot
figure
plot(time_SAR, data_SAR(PS_id,:), 'k', 'LineWidth', 1)
hold on
errorbar(time_SAR, data_SAR_est1, v_std1, 'b', 'LineWidth', 0.7)
plot(time_SAR, data_SAR_est1, 'b', 'LineWidth', 2)
plot(time_SAR, data_SAR_est, '--r', 'LineWidth', 1)
xlabel('Time [days]')
ylabel('Displacement [mm]')
title('SAR displacement', 'FontSize', 20)
set(gca, 'FontSize', 15)
legend('Observations', 'Estimations', 'FontSize', 15)



%% STEP 12: Geospatial Analysis - 2D Covariance Function
% This step is meant to compute the 2D spatial-temporal covariance
% function, as the product between the temporal and spatial one. All the PS
% are now considered.
%
% C(τ,d) = C₁(τ) x C₂(d)

% Temporal covariance -----------------------------------------------------
% This means creating an empirical covariance function considering all the
% possible couples of PS displacements, neglecting from which PS it comes
% from, just caring about the temporal distance (lag).
for i = 1:size(v_SAR,1)

    [tGrid, eCovF_i, ~, ~] = f1DEmpCovEst(v_SAR(i,:)', time_SAR', 6, 0);
    eCovF_t(:,i) = eCovF_i;

end
eCovF_t_avg = mean(eCovF_t,2);

% Exponential model
idx = find(eCovF_t_avg < 0, 1, 'first') - 1;       
A  = [ones(size(tGrid(2:idx))) -tGrid(2:idx)];     
yo = log(eCovF_t_avg(2:idx));                      
peT = inv(A'*A) * A'*yo;                            
peT(1) = exp(peT(1));                                
peT(2) = peT(2);                                     
fcovT = @(tau) peT(1) * exp(-peT(2) * tau);
s2vT = eCovF_t_avg(1) - fcovT(0);

% Plot
figure
plot(tGrid, eCovF_t_avg, '.-', 'LineWidth', 1.3, 'MarkerSize', 12);
hold on; plot(tGrid, fcovT(tGrid), '-', 'LineWidth', 2);
xlabel('Time lag [days]', 'FontSize', 15); ylabel('Covariance [mm^2]', 'FontSize', 15); 
title('Temporal covariance modelling', 'FontSize', 20);
xlim([0 tGrid(end)/2]);
set(gca, 'FontSize', 15);


% Spatial covariance ------------------------------------------------------
% This means creating an empirical covariance function considering all the
% possible spatial distances, considering a fixed epoch in time.

[~, ~, sigmaGrid, eCovF_i, ~, ~] = f2DCovEmpEst(v_SAR(:), repmat(x_SAR,size(v_SAR,2),1), repmat(y_SAR,size(v_SAR,2),1), 20, 'cartesian', 2);
eCovF_s = eCovF_i;

% Exponential model
idx = find(eCovF_s < 0, 1, 'first') - 1;       
A  = [ones(size(sigmaGrid(2:idx))) -sigmaGrid(2:idx)];     
yo = log(eCovF_s(2:idx));                      
peS_tmp = inv(A'*A) * A'*yo;                            
peS(1) = exp(peS_tmp(1));                                
peS(2) = peS_tmp(2);                                     
fcovS = @(tau) peS(1) * exp(-peS(2) * tau);
s2vS = eCovF_s(1) - fcovS(0);

% Plot
figure
plot(sigmaGrid, eCovF_s, '.-', 'LineWidth', 1.3, 'MarkerSize', 12);
hold on; plot(sigmaGrid, fcovS(sigmaGrid), '-', 'LineWidth', 2);
xlabel('Time lag [days]', 'FontSize', 15); ylabel('Covariance [mm^2]', 'FontSize', 15); 
title('Temporal covariance modelling', 'FontSize', 20);
xlim([0 sigmaGrid(end)/2]);
set(gca, 'FontSize', 15);


% Collocation -------------------------------------------------------------
% Time of estimation
t_est = time_SAR';

% Observations - Observations covariance matrix
Cvo_voST = zeros(size(v_SAR,1)*size(v_SAR,2), size(v_SAR,1)*size(v_SAR,2));
for i = 1:size(v_SAR,2)
    for m = 1:size(v_SAR,2)
        for j = 1:size(v_SAR,1)
            for k = 1:size(v_SAR,1)

                % Calculate temporal and spatial distances
                dt = abs(time_SAR(i)-time_SAR(m));
                ds = sqrt((x_SAR(j)-x_SAR(k))^2 + (y_SAR(j)-y_SAR(k))^2);

                % Calculate indices for matrix assignment
                row_index = (i - 1) * size(v_SAR, 1) + j;
                col_index = (m - 1) * size(v_SAR, 1) + k;

                Cvo_voST(row_index,col_index) = fcovT(dt) * fcovS(ds);


            end
        end
    end
end

% Estimations - Observations covariance matrix
x_lowres = downsample(x_dam(1,:), 10);
y_lowres = downsample(y_dam(:,1), 10);
[x_dam_lowres, y_dam_lowres] = meshgrid(x_lowres, y_lowres);
x_dam_rs = x_dam_lowres(:);
y_dam_rs = y_dam_lowres(:);

Cv_voST = zeros(size(x_dam_rs,1)*size(v_SAR,2), size(v_SAR,1)*size(v_SAR,2));
for i = 1:size(v_SAR,2)
    for m = 1:size(v_SAR,2)
        for j = 1:size(x_dam_rs,1)
            for k = 1:size(v_SAR,1)

                % Calculate temporal and spatial distances
                dt = abs(time_SAR(i)-time_SAR(m));
                ds = sqrt((x_dam_rs(j)-x_SAR(k))^2 + (y_dam_rs(j)-y_SAR(k))^2);

                % Calculate indices for matrix assignment
                row_index = (i - 1) * size(x_dam_rs, 1) + j;
                col_index = (m - 1) * size(v_SAR, 1) + k;

                Cv_voST(row_index,col_index) = fcovT(dt) * fcovS(ds);


            end
        end
    end
end

% Collocation formula
s2vST = (s2vS + s2vT)/2;
v_ST_tmp =  Cv_voST * inv(Cvo_voST + s2vST * eye(size(Cvo_voST))) * v_SAR(:);
v_ST = reshape(v_ST_tmp, size(v_SAR,2), size(x_dam_rs,1))';


% Grid resampling
for i = 1:size(v_ST,2)

    [x_dam1, y_dam1, v_ST_grid_i] = gridResampling ...
        (x_dam_lowres, y_dam_lowres, v_ST(:,i), x_dam, y_dam, area_dam);
    v_ST_grid(:,:,i) = v_ST_grid_i;

end


% Final model -------------------------------------------------------------
% Result from Fourier analysis/synthesis and collocation, spatially
% resampled over the dam area.
d4D_dam = dam_SAR_harm_tobs + v_ST_grid;


% Plot
plot_gif = 1;
filename = 'Dam_LOS_ASC_fin_3D.gif';
plotGIF(x_dam1, y_dam1, d4D_dam, filename, plot_gif);


fprintf('\nScript execution completed!\n');

