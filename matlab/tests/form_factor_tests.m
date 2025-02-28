classdef form_factor_tests < matlab.unittest.TestCase
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
        function form_factor_standard_basis_test(testCase)
            q = [1, 2; 3, 4];
            A = eye(2);

            actual = form_factor(q, A, abs(det(A)));
            expected = 2.5087;

            testCase.verifyEqual(actual, expected, "RelTol", 1e-3);
        end

        function form_factor_rotated_basis_test(testCase)
            e = exp(1);

            q = [1, 2; 3, 4];
            A = comp2mat(e^(j * pi / 3));

            actual = form_factor(q, A, abs(det(A)));
            expected = 2.5087;

            testCase.verifyEqual(actual, expected, "RelTol", 1e-3);
        end

        function form_factor_scaled_basis_test(testCase)
            q = [1, 2; 3, 4];
            A = [2, 0; 0, 3];

            actual = form_factor(q, A, abs(det(A)));
            expected = 2.5087 / sqrt(2 * 3);

            testCase.verifyEqual(actual, expected, "RelTol", 1e-3);
        end

        function form_factor_non_orthogonal_basis_test(testCase)
            q = [1, 2; 3, 4];

            a1 = [1; 0];
            a2 = [sqrt(2) / 2; sqrt(2) / 2];

            A = [a1, a2];

            actual = form_factor(q, A, abs(det(A)));
            expected = 2.5087 / sqrt(sqrt(2) / 2);

            testCase.verifyEqual(actual, expected, "RelTol", 1e-3);
        end

        function form_factor_non_orthogonal_rotated_basis_test(testCase)
            e = exp(1);

            q = [1, 2; 3, 4];

            a1 = [1; 0];
            a2 = [sqrt(2) / 2; sqrt(2) / 2];

            A = comp2mat(e^(j * pi / 3)) * [a1, a2];

            actual = form_factor(q, A, abs(det(A)));
            expected = 2.5087 / sqrt(sqrt(2) / 2);

            testCase.verifyEqual(actual, expected, "RelTol", 1e-3);
        end

        function form_factor_non_orthogonal_scaled_basis_test(testCase)
            q = [1, 2; 3, 4];

            a1 = 2 * [1; 0];
            a2 = 3 * [sqrt(2) / 2; sqrt(2) / 2];

            A = [a1, a2];

            actual = form_factor(q, A, abs(det(A)));
            expected = 2.5087 / sqrt(2 * 3 * sqrt(2) / 2);

            testCase.verifyEqual(actual, expected, "RelTol", 1e-3);
        end

        function form_factor_dims_1_test(testCase)
            q = [1, 2; 3, 4];
            q = repmat(q, [1, 1, 2]);

            a1 = 2 * [1; 0];
            a2 = 3 * [sqrt(2) / 2; sqrt(2) / 2];

            A = [a1, a2];

            actual = form_factor(q, A, abs(det(A)));

            expected = 2.5087 / sqrt(2 * 3 * sqrt(2) / 2);
            expected = repmat(expected, [2, 1]);

            testCase.verifyEqual(actual, expected, "RelTol", 1e-3);
        end

        function form_factor_dims_2_test(testCase)
            q = [1, 2; 3, 4];
            q = repmat(q, [1, 1, 2, 2]);

            a1 = 2 * [1; 0];
            a2 = 3 * [sqrt(2) / 2; sqrt(2) / 2];

            A = [a1, a2];

            actual = form_factor(q, A, abs(det(A)));

            expected = 2.5087 / sqrt(2 * 3 * sqrt(2) / 2);
            expected = repmat(expected, [2, 2]);

            testCase.verifyEqual(actual, expected, "RelTol", 1e-3);
        end
    end
end
