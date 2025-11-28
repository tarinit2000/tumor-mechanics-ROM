function [Damper, VM, Ux, Uy, Exx, Eyy, Exy, Sxx, Syy, Sxy, timings] = ...
    getMechanicsMaps_2D_LUonce(N, M, E, nu, matX, matY, mode, stride)

    % goal: loop through time steps and call mechanics solver (get_damper) 
    % to compute displ, strain, and stresses for each snapshot of N

%     if nargin < 8
%         mode = 'full'; % default
%         % full: compute mechanics every step (baseline)
%         % subsample: compute mechanics only at stride steps, reuse last values otherwise
%     end
%     if nargin < 9
%         stride = 1; % default
%         % stride: how often to recompute mechanics in subsample mode
%     end
%     
    [sy,sx,nt] = size(N);

    % preallocate
    Damper = zeros(sy,sx,nt);
    VM     = zeros(sy,sx,nt);
    Ux     = zeros(sy,sx,nt);
    Uy     = zeros(sy,sx,nt);
    Exx    = zeros(sy,sx,nt);
    Eyy    = zeros(sy,sx,nt);
    Exy    = zeros(sy,sx,nt);
    Sxx    = zeros(sy,sx,nt);
    Syy    = zeros(sy,sx,nt);
    Sxy    = zeros(sy,sx,nt);

    timings = zeros(nt,1); % store per-step runtime

    for i = 1:nt
        tic; % start timing

        if strcmp(mode,'full') || mod(i,stride) == 0 || i == 1
            % compute mechanics fresh
            [temp_damper, temp_vm, temp_x, temp_y, ...
             temp_Exx, temp_Eyy, temp_Exy, ...
             temp_Sxx, temp_Syy, temp_Sxy] = ...
                get_damper_LUonce(matX, matY, N(:,:,i), M, E, nu);

            % reshape results back into 2D grids
            Damper(:,:,i) = reshape(temp_damper, sy,sx);
            VM(:,:,i)     = reshape(temp_vm, sy,sx);
            Ux(:,:,i)     = reshape(temp_x, sy,sx);
            Uy(:,:,i)     = reshape(temp_y, sy,sx);
            Exx(:,:,i)    = reshape(temp_Exx, sy,sx);
            Eyy(:,:,i)    = reshape(temp_Eyy, sy,sx);
            Exy(:,:,i)    = reshape(temp_Exy, sy,sx);
            Sxx(:,:,i)    = reshape(temp_Sxx, sy,sx);
            Syy(:,:,i)    = reshape(temp_Syy, sy,sx);
            Sxy(:,:,i)    = reshape(temp_Sxy, sy,sx);
        else
            % reuse last computed values
            Damper(:,:,i) = Damper(:,:,i-1);
            VM(:,:,i)     = VM(:,:,i-1);
            Ux(:,:,i)     = Ux(:,:,i-1);
            Uy(:,:,i)     = Uy(:,:,i-1);
            Exx(:,:,i)    = Exx(:,:,i-1);
            Eyy(:,:,i)    = Eyy(:,:,i-1);
            Exy(:,:,i)    = Exy(:,:,i-1);
            Sxx(:,:,i)    = Sxx(:,:,i-1);
            Syy(:,:,i)    = Syy(:,:,i-1);
            Sxy(:,:,i)    = Sxy(:,:,i-1);
        end

        timings(i) = toc; % record runtime
    end
end
