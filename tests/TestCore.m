classdef TestCore < matlab.unittest.TestCase
    properties
        sy = 8;
        sx = 8;
        n
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.n = testCase.sy * testCase.sx;
            % deterministic RNG for reproducibility
            rng(0, 'twister');
            % dependency check
            required = {'mech_matrix_build_2D','grad_matrix','buildBoundaries_2D','get_damper'};
            for k = 1:numel(required)
                if exist(required{k}, 'file') ~= 2
                    testCase.assumeFail(sprintf('Missing required helper: %s', required{k}));
                end
            end
        end
    end

    methods(Test)
        function smokeTest(testCase)
            N_small = zeros(testCase.sy, testCase.sx, 2);
            M_small = speye(2*testCase.n);
            E_small = ones(testCase.sy, testCase.sx);
            nu_small = 0.3;
            [d_dX, d_dY] = grad_matrix(1, buildBoundaries_2D(ones(testCase.sy,testCase.sx)));
            [~, vm] = get_damper(d_dX, d_dY, N_small(:,:,1), M_small, E_small, nu_small);
            testCase.verifyTrue(isnumeric(vm) && all(isfinite(vm(:))));
            testCase.verifyGreaterThanOrEqual(min(vm(:)), 0); % von Mises nonnegative
        end

        function gradMatrixShape(testCase)
            [d_dX, d_dY] = grad_matrix(1, buildBoundaries_2D(ones(4,4)));
            testCase.verifySize(d_dX, [16,16]);
            testCase.verifySize(d_dY, [16,16]);
            % check sparsity pattern: expect mostly zeros
            testCase.verifyLessThan(nnz(d_dX)/numel(d_dX), 0.2);
        end
    end
end
