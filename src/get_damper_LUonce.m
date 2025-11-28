function [damper, s_vm, Ux, Uy, e_xx, e_yy, e_xy, s_xx, s_yy, s_xy] = get_damper_LUonce(matX, matY, N, M, E, nu)
    [sy,sx] = size(N);

    lambda1 = 2.5e-3;
    lambda2 = 2.5e-3;
%     lambda2 =  1000;
    
    n_vect = reshape(N,[],1);
    % calculate the gradient using matX and matY
    gradX = matX * n_vect;
    gradY = matY * n_vect;
    grad_N = [gradX; gradY];

    % Calculate displacement
    %Uv = M\(lambda1*grad_N);

    % Wrap top-level iteration in try/catch and save a minimal crash dump 
    % (inputs + stack) for postmortem.
    try
        persistent L U P Q
        if isempty(L)
            [L,U,P,Q] = lu(M);   % original 4-output LU
        end

        rhs = lambda1 * [matX * N(:); matY * N(:)];
        y   = U \ (L \ (P * rhs));
        Uv  = Q * y;
    catch ME
        timestamp = datestr(now,'yyyymmdd_HHMMSS');
        save(fullfile('crash_dumps', ['crash_' timestamp '.mat']), 'ME', 'M', 'rhs', 'N', 'matX', 'matY', '-v7.3');
        rethrow(ME);
    end



    % Reshape U into matrix
%     Um = zeros(sy,sx,2);
%     Um(:,:,1) = reshape(Uv(1:sx*sy), sy,sx);
%     Um(:,:,2) = reshape(Uv(sx*sy+1:end),sy,sx);

    Um = zeros(sy*sx,2);
    Um(:,1) = Uv(1:sx*sy); Ux = Um(:,1);
    Um(:,2) = Uv(sx*sy+1:end); Uy = Um(:,2);

%     [e_xx,e_yy,e_xy] = strains(Um, h, h);
    
    e_xx = matX*Um(:,1);
    e_yy = matY*Um(:,2);
    e_xy = matY*Um(:,1);
    
    
    [s_xx, s_yy, s_xy, s_vm] = stresses(e_xx,e_yy, e_xy, E(:), nu);
    
    damper = exp(-1*s_vm * lambda2);
end