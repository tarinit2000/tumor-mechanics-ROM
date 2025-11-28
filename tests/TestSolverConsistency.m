classdef TestSolverConsistency < matlab.unittest.TestCase
    methods(TestMethodSetup)
        function setup(testCase)
            rng(0,'twister');
        end
    end

    methods(Test)
        function luEqualsBackslash(testCase)
            n = 20;
            A = randn(n); A = A'*A + eye(n); % SPD
            b = randn(n,1);
            x1 = A \ b;
            [L,U,P] = lu(A);
            x2 = U \ (L \ (P*b));
            relerr = norm(x1-x2)/max(norm(x1),1e-16);
            testCase.verifyLessThan(relerr, 1e-12);
        end

        function solveSparseSPD(testCase)
            % test sparse SPD solve path used in mechanics
            n = 50;
            R = sprandsym(n, 0.05, 1e-1, 1); % sparse SPD-ish
            A = R + speye(n);
            b = randn(n,1);
            x1 = A \ b;
            % check residual
            testCase.verifyLessThan(norm(A*x1 - b)/norm(b), 1e-10);
        end
    end
end
