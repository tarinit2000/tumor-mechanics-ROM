%% Optimization + Tradeoff Analysis: Full vs Subsampled vs ROM
% Self-contained script that:
%  - loads example data
%  - builds mechanics matrices and ROM operators
%  - computes von Mises for FOM and subsampled runs
%  - computes ROM reduced solve and ROM error
%  - produces per-step von Mises error plot and a tradeoff plot (avg von Mises error vs speedup)
%
% Requirements: helper functions must be on the MATLAB path:
% mech_matrix_build_2D, grad_matrix, buildBoundaries_2D, get_damper, get_damper_reduced,
% augmentCellMaps_2D, getMechanicsMaps_2D, getProjectionMatrix, getDisplacementProjection_2D,
% getStrainProjection_2D, getStressProjection_2D, buildStrainMat

clear all; clc; close all;
disp('--- Optimization + Tradeoff Analysis (Full script) ---');

%% -------------------------
%% 1) Load data and setup
%% -------------------------
addpath(genpath(pwd));

location = fullfile(pwd,'data','Ex5_patient.mat');
if ~isfile(location)
    error('Data file not found: %s', location);
end
tumor = load(location);

% Assemble cell maps
tumor.N(:,:,1) = tumor.image_data.NTC1;
tumor.N(:,:,2) = tumor.image_data.NTC2;
tumor.N(:,:,3) = tumor.image_data.NTC3;

% Inputs
N0      = tumor.N(:,:,1);                % initial cell map
N_true  = tumor.N(:,:,2:end);            % time series for ROM/training
t       = tumor.schedule_info.times;     % imaging times
h       = tumor.schedule_info.imagedims(1); % grid spacing (assumed square)
bcs     = buildBoundaries_2D(tumor.image_data.BreastMask);
Tissues = tumor.image_data.Tissues;

%% -------------------------
%% 2) Mechanics matrices
%% -------------------------
[M,E,nu] = mech_matrix_build_2D(h, Tissues, bcs);
[d_dX, d_dY] = grad_matrix(h, bcs);

%% -------------------------
%% 3) Full-order reference von Mises (single call)
%% -------------------------
% We get sig_vm here just for checking ROM single-point error later
[~, sig_vm] = get_damper(d_dX, d_dY, N0, M, E, nu);

%% -------------------------
%% 4) Build augmented snapshots and ROM bases
%% -------------------------
N_augmented = augmentCellMaps_2D(cat(3,N0,N_true), t(2:end), 4); % augmentation depth = 4

% Collect mechanics outputs across augmented time series (full solves) to build basis
% We call this once just to get the data for POD
[~, ~, Ux_aug, Uy_aug, Exx_aug, Eyy_aug, Exy_aug, Sxx_aug, Syy_aug, Sxy_aug] = ...
    getMechanicsMaps_2D(N_augmented, M, E, nu, d_dX, d_dY, 'full', 1);

% Build projection matrices / bases
[V, k] = getProjectionMatrix(N_augmented, 0);
[V_u, V_Ux, V_Uy] = getDisplacementProjection_2D(Ux_aug, Uy_aug, k);
V_e = getStrainProjection_2D(Exx_aug, Eyy_aug, Exy_aug, k);
V_s = getStressProjection_2D(Sxx_aug, Syy_aug, Sxy_aug, k);

% Stress operator and reduced operators
S_mat = buildStrainMat(N0, E, nu);

M_r = V_u' * M * V_u;

% Compose gradient/projection operators
Vu_gradXY_V = V_u' * [d_dX, zeros(size(d_dX)); zeros(size(d_dX)), d_dY] * ...
              [V, zeros(size(V)); zeros(size(V)), V];

Ve_gradXYY_Vu = V_e' * [d_dX, zeros(size(d_dX)), zeros(size(d_dX)); ...
                        zeros(size(d_dX)), d_dY, zeros(size(d_dX)); ...
                        zeros(size(d_dX)), zeros(size(d_dX)), d_dY] * ...
                        [V_Ux, zeros(size(V_Uy)), zeros(size(V_Ux)); ...
                         zeros(size(V_Ux)), V_Uy, zeros(size(V_Ux)); ...
                         zeros(size(V_Ux)), zeros(size(V_Uy)), V_Ux];

Vs_SMAT_Ve = V_s' * S_mat * V_e;

% Reduced initial state
N0_r = V' * N0(:);

%% -------------------------
%% 5) Baseline ROM run (single reduced solve)
%% -------------------------
% Use timeit for stable measurement
rom_func = @() get_damper_reduced(Vu_gradXY_V, Ve_gradXYY_Vu, Vs_SMAT_Ve, N0_r, M_r, V_s);
t_ROM_mech = timeit(rom_func);

% Actual call to get data
[~, sig_vm_ROM] = rom_func();

fprintf('(1) Comparing baseline ROM, full FOM, and subsampled FOM:\n')
fprintf('- ROM solve time: %.6f s\n', t_ROM_mech);

%% -------------------------
%% 6) Baseline full and subsampled mechanics (default stride)
%% -------------------------
stride_default = 4;

% --- Full mechanics ---
full_mech_func = @() getMechanicsMaps_2D(N_augmented, M, E, nu, d_dX, d_dY, 'full', 1);

% Measure time (averaged)
t_full_mech = timeit(full_mech_func); 

% Get data (Capture VM directly as 2nd output!)
[~, VM_full, Ux_full, Uy_full, Exx_full, Eyy_full, Exy_full, ...
    Sxx_full, Syy_full, Sxy_full] = full_mech_func();

fprintf('- Full mechanics run time: %.6f s\n', t_full_mech);

% --- Subsampled mechanics (default stride) ---
sub_mech_func = @() getMechanicsMaps_2D(N_augmented, M, E, nu, d_dX, d_dY, 'subsample', stride_default);

% Measure time (averaged)
t_sub_mech = timeit(sub_mech_func);

% Get data (Capture VM directly)
[~, VM_sub, Ux_sub, Uy_sub, Exx_sub, Eyy_sub, Exy_sub, ...
    Sxx_sub, Syy_sub, Sxy_sub] = sub_mech_func();

fprintf('- Subsampled (stride=%d) run time: %.6f s\n', stride_default, t_sub_mech);

%% -------------------------
%% 7) Compute von Mises per time step (full and subsampled)
%% -------------------------
% Reshape pre-calculated VM outputs to (npts x nt)
sig_vm_full_t = reshape(VM_full, [], size(VM_full,3));
sig_vm_sub_t  = reshape(VM_sub,  [], size(VM_sub,3));

% (No manual calculation loop needed anymore!)
nt = size(VM_full,3);

%% -------------------------
%% 8) Per-time-step relative errors (von Mises + optional mechanics)
%% -------------------------
err_vm_t  = zeros(nt,1);
err_Ux_t  = zeros(nt,1);
err_Uy_t  = zeros(nt,1);
err_Sxx_t = zeros(nt,1);
err_Syy_t = zeros(nt,1);
err_Sxy_t = zeros(nt,1);

for i = 1:nt
    % Displacement/stress errors
    err_Ux_t(i)  = norm(Ux_full(:,:,i) - Ux_sub(:,:,i))  / max(norm(Ux_full(:,:,i)), 1e-12);
    err_Uy_t(i)  = norm(Uy_full(:,:,i) - Uy_sub(:,:,i))  / max(norm(Uy_full(:,:,i)), 1e-12);
    err_Sxx_t(i) = norm(Sxx_full(:,:,i) - Sxx_sub(:,:,i)) / max(norm(Sxx_full(:,:,i)), 1e-12);
    err_Syy_t(i) = norm(Syy_full(:,:,i) - Syy_sub(:,:,i)) / max(norm(Syy_full(:,:,i)), 1e-12);
    err_Sxy_t(i) = norm(Sxy_full(:,:,i) - Sxy_sub(:,:,i)) / max(norm(Sxy_full(:,:,i)), 1e-12);

    % von Mises error
    den_vm = max(norm(sig_vm_full_t(:,i)), 1e-12);
    err_vm_t(i) = norm(sig_vm_full_t(:,i) - sig_vm_sub_t(:,i)) / den_vm;
end

%% Per-step plot: von Mises (primary) and optional other mechanics curves
figure('Color','w','Position',[200 200 900 420]);

plot(1:nt, 100*err_vm_t, '-ko', 'LineWidth',1.6, 'MarkerSize',6, 'DisplayName','von Mises (subsample vs FOM)');
hold on;
plot(1:nt, 100*err_Ux_t, '--b', 'LineWidth',1.2, 'DisplayName','Ux (subsample vs FOM)');
plot(1:nt, 100*err_Uy_t, '--c', 'LineWidth',1.2, 'DisplayName','Uy (subsample vs FOM)');
plot(1:nt, 100*err_Sxx_t, ':r', 'LineWidth',1.2, 'DisplayName','Sxx (subsample vs FOM)');
plot(1:nt, 100*err_Syy_t, ':m', 'LineWidth',1.2, 'DisplayName','Syy (subsample vs FOM)');
plot(1:nt, 100*err_Sxy_t, ':g', 'LineWidth',1.2, 'DisplayName','Sxy (subsample vs FOM)');

xlabel('Time step','FontWeight','bold');
ylabel('Relative percent error vs FOM (%)','FontWeight','bold');
title(['Per-step relative error (stride=', num2str(stride_default), ')'],'FontWeight','bold');
legend('Location','northeastoutside');
grid on;
set(gca,'FontSize',11);

%% Summary metrics (von Mises)
fprintf('(2) Summary metrics:\n')

% Calculate averages (convert to percent)
avg_err_Ux_pct  = 100 * mean(err_Ux_t);
avg_err_Uy_pct  = 100 * mean(err_Uy_t);
avg_err_Sxx_pct = 100 * mean(err_Sxx_t);
avg_err_Syy_pct = 100 * mean(err_Syy_t);
avg_err_Sxy_pct = 100 * mean(err_Sxy_t);
avg_err_vm_pct  = 100 * mean(err_vm_t);

% Print results
fprintf('- Avg Ux  error (stride=%d): %.4f %%\n', stride_default, avg_err_Ux_pct);
fprintf('- Avg Uy  error (stride=%d): %.4f %%\n', stride_default, avg_err_Uy_pct);
fprintf('- Avg Sxx error (stride=%d): %.4f %%\n', stride_default, avg_err_Sxx_pct);
fprintf('- Avg Syy error (stride=%d): %.4f %%\n', stride_default, avg_err_Syy_pct);
fprintf('- Avg Sxy error (stride=%d): %.4f %%\n', stride_default, avg_err_Sxy_pct);
fprintf('- Avg VM  error (stride=%d): %.4f %%\n', stride_default, avg_err_vm_pct);
%% Speedup factors (default stride + ROM)
speedup_sub = t_full_mech / t_sub_mech;
speedup_rom = t_full_mech / t_ROM_mech;

fprintf('(3) Speedup factors:\n')
fprintf('Subsampled speedup factor (stride=%d): %.2f\n', stride_default, speedup_sub);
fprintf('ROM speedup factor: %.2f\n', speedup_rom);

%% -------------------------
%% 9) Tradeoff sweep across strides (Fixed Duration)
%% -------------------------
stride_values = [1, 2, 4, 8, 16];
runtime_sub   = zeros(length(stride_values),1);
avg_vm_error_sub = zeros(length(stride_values),1);

fprintf('(4) Tradeoff btwn full & subsampled FOM: sweeping stride values\n');

% Set how many times to repeat for averaging. 
% 5 is enough to smooth OS noise without taking forever.
n_repeats = 5; 

for k = 1:length(stride_values)
    stride_k = stride_values(k);

    % 1. Measure Runtime & Capture Data Simultaneously
    % We loop manually to average time, but keep the data from the last run
    tic;
    for r = 1:n_repeats
        [~, VM_sub_k] = getMechanicsMaps_2D(N_augmented, M, E, nu, d_dX, d_dY, 'subsample', stride_k);
    end
    t_total = toc;
    
    runtime_sub(k) = t_total / n_repeats;

    % 2. Compute von Mises error (using data from the loop)
    nt_k = size(VM_sub_k,3);
    err_vm_k = nan(nt_k,1);
    
    for i = 1:nt_k
        vm_f_vec = sig_vm_full_t(:,i);
        vm_s_vec = reshape(VM_sub_k(:,:,i), [], 1);

        den_vm = max(norm(vm_f_vec), 1e-12);
        err_vm_k(i) = norm(vm_f_vec - vm_s_vec) / den_vm;
    end

    avg_vm_error_sub(k) = mean(err_vm_k);
    fprintf('stride=%d: runtime=%.6fs, avg von Mises err=%.4f (%.3f%%)\n', ...
        stride_k, runtime_sub(k), avg_vm_error_sub(k), 100*avg_vm_error_sub(k));
end
%% ROM single-point error (ensure alignment)
rom_error = norm(sig_vm_ROM(:) - sig_vm(:)) / max(norm(sig_vm(:)), 1e-12);
rom_speedup = t_full_mech / t_ROM_mech;
fprintf('ROM error (von Mises) = %.4f (%.3f%%)\n', rom_error, 100*rom_error);

%% -------------------------
%% 10) Tradeoff plot (categorical) using von Mises errors
%% -------------------------
err_pct_sub = 100 * avg_vm_error_sub(:);          
speedup_sub = t_full_mech ./ runtime_sub(:);      
rom_err_pct  = 100 * rom_error;
rom_spd      = rom_speedup;

cats = [arrayfun(@(s) sprintf('s=%d',s), stride_values, 'UniformOutput', false), {'ROM'}];
xcat = categorical(cats);
xcat = reordercats(xcat, cats);

figure('Color','w','Position',[200 200 900 420]);

% Left axis: von Mises error (%)
yyaxis left
h_err = plot(xcat(1:length(stride_values)), err_pct_sub, '-o', ...
    'LineWidth',1.6, 'MarkerSize',8, 'Color',[0 0.4470 0.7410]);
hold on
h_rom_err_main = plot(xcat(end), rom_err_pct, 'kx', 'MarkerSize',10, 'LineWidth',1.8);
ylabel('Average von Mises error vs FOM (%)','FontWeight','bold');
ylim([0, max([err_pct_sub; rom_err_pct])*1.25]);

% Right axis: speedup
yyaxis right
h_spd = plot(xcat(1:length(stride_values)), speedup_sub, '-s', ...
    'LineWidth',1.6, 'MarkerSize',8, 'Color',[0.85 0.33 0.10]);
hold on
h_rom_spd_main = plot(xcat(end), rom_spd, 'k*', 'MarkerSize',10, 'LineWidth',1.8);
ylabel('Speedup factor vs FOM','FontWeight','bold');
ylim([0, max([speedup_sub; rom_spd])*1.25]);

xlabel('Method','FontWeight','bold');
title('von Mises Error vs Speedup','FontWeight','bold');

handles = [h_err, h_spd, h_rom_err_main, h_rom_spd_main];
valid = arrayfun(@(h) isgraphics(h), handles);
if all(valid)
    legendLabels  = {'Subsample avg von Mises error (%)', 'Subsample speedup', 'ROM error (%)', 'ROM speedup'};
    lg = legend(handles, legendLabels, 'Location', 'northwest');
    lg.Box = 'off';
    set(lg, 'FontSize', 10);
else
    legend([h_err, h_spd], {'Subsample avg von Mises error (%)', 'Subsample speedup'}, 'Location','northwest');
end

grid on;
set(gca,'FontSize',11);

disp('--- Script complete ---');