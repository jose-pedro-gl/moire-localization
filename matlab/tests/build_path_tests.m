classdef build_path_tests < matlab.unittest.TestCase
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
        function build_path_test(testCase)
            actual = build_path(1 / 2 * [0, 0; 1, 0; 1, 1; 0, 0], 2);
            
            expected = [
                0, 0;
                0.25, 0;
                0.5, 0;
                0.5, 0.25;
                0.5, 0.5;
                0.25, 0.25;
                0, 0
            ];

            testCase.verifyEqual(actual, expected, "RelTol", 1e-3);
        end
    end
end
