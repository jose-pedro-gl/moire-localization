classdef fin_diff_matrix_sparse_2_tests < matlab.unittest.TestCase
    properties (TestParameter)
        A = {
            eye(2)
            pi * eye(2)
            with_root(@ () comp2mat(exp(j * pi / 3)))
            with_root(@ () comp2vec([1, exp(j * pi / 4)]))
            with_root(@ () comp2mat(exp(j * pi / 3)) ...
                * comp2vec([1, exp(j * pi / 4)]))
        };
        bc = {
            struct("type", "bloch", "k", [0, 0])
            struct("type", "bloch", "k", [0, 0.2])
            struct("type", "bloch", "k", [0.3, 0.4])
            struct("type", "bloch", "k", [0.1, 0])
            struct("type", "zero")
        };
        V_impl = {
            @(UU, VV) 0 * UU
            @(UU, VV) V_2016(UU, VV, UU, VV, 1, 0)
            @(UU, VV) V_2016(UU, VV, UU, VV, 10, 0)
            @(UU, VV) V_2016(UU, VV, UU, VV, 100, 0)
            @(UU, VV) V_2016(UU, VV, UU, VV, 1000, 0)
        };
    end

    methods (TestClassSetup)
        function add_paths(testCase)
            addpath("..");
        end
    end

    methods (TestClassTeardown)
        function rm_paths(testCase)
            rmpath("..");
        end
    end

    methods (Test)
        function fin_diff_matrix_ref_test(testCase, A, bc, V_impl)
            TU = linspace_ex(0, 1, 21);
            TV = linspace_ex(0, 1, 21);

            UV = space_grid_2(TU, TV);
            UU = UV(:, :, 1);
            VV = UV(:, :, 2);

            V = V_impl(UU, VV);

            actual = fin_diff_matrix_sparse_2(A, V, bc);
            expected = fin_diff_matrix_2(A, V, bc);

            expected(abs(expected) < 1e-12) = 0;
            expected = sparse(expected);

            testCase.verifyEqual(actual, expected, "RelTol", 1e-9);
            testCase.verifyTrue(ishermitian(actual));
        end
    end
end
