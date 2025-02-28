classdef pde_fn_2_tests < matlab.unittest.TestCase
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
            struct("type", "bloch", "k", [0.3, 0.4])
        };
        Nb = {
            30
        };
        V_impl = {
            @(UU, VV) 0 * UU
            @(UU, VV) V_2016(UU, VV, UU, VV, 1, 0)
            @(UU, VV) V_2016(UU, VV, UU, VV, 10, 0)
            @(UU, VV) V_2016(UU, VV, UU, VV, 100, 0)
            @(UU, VV) V_2016(UU, VV, UU, VV, 1000, 0)
        };
        method = {
            struct("Nu", 79, "Nv", 77, "pde_fn_2", @(A, V, bc, Nb, eigs_impl) ...
                fin_diff_fn_2(A, V, bc, Nb, eigs_impl))
            struct("Nu", 49, "Nv", 47, "pde_fn_2", @plane_wave_fn_2)
        };
        eigs_impl = {
            @(M, N) eigs_matlab(M, N, "largestreal", 1e-6)
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
        function pde_test(testCase, A, bc, Nb, V_impl, method, eigs_impl)
            k = bc.k;

            Nu = method.Nu;
            Nv = method.Nv;
            pde_fn_2 = method.pde_fn_2;

            TU = linspace_ex(0, 1, Nu);
            TV = linspace_ex(0, 1, Nv);
            [UU, VV] = meshgrid(TU, TV);

            V = V_impl(UU, VV);

            [bands, q_largest_b] = pde_fn_2(A, V, bc, Nb, eigs_impl);

            testCase.verifyEqual(size(bands), [Nb, 1]);
            testCase.verifyEqual(size(q_largest_b), [Nv, Nu, Nb]);

            b = bands(1, 1);
            q = q_largest_b(:, :, 1);

            testCase.verifyEqual(q(:)' * q(:), 1, "RelTol", 1e-6);

            L1 = b * q;
            L2 = 1 / 2 * laplacian(A, k, q, 1) + V .* q;

            L1(abs(L1) < 1e-9) = 0;
            L2(abs(L2) < 1e-9) = 0;

            L1_max = max(abs(L1), [], "all");
            L2_max = max(abs(L2), [], "all");

            testCase.verifyEqual(L1_max, L2_max, "RelTol", 1e-1);
            testCase.verifyEqual(L1, L2, "AbsTol", abs(L1_max + L2_max) / 2);
        end

        function eigenvalue_stability_test(testCase, A, bc, Nb, V_impl, method, eigs_impl)
            Nu = method.Nu;
            Nv = method.Nv;
            pde_fn_2 = method.pde_fn_2;

            TU = linspace_ex(-0.5, 0.5, Nu);
            TV = linspace_ex(-0.5, 0.5, Nv);
            [UU, VV] = meshgrid(TU, TV);

            V = V_impl(UU, VV);

            [bands, q_largest_b] = pde_fn_2(A, V, bc, Nb, eigs_impl);

            testCase.verifyEqual(size(bands), [Nb, 1]);
            testCase.verifyEqual(size(q_largest_b), [Nv, Nu, Nb]);

            Nu2 = Nu + 4;
            Nv2 = Nv + 4;

            TU2 = linspace_ex(0, 1, Nu2);
            TV2 = linspace_ex(0, 1, Nv2);
            [UU2, VV2] = meshgrid(TU2, TV2);

            V2 = V_impl(UU2, VV2);

            [bands2, q_largest_b2] = pde_fn_2(A, V2, bc, Nb, eigs_impl);

            testCase.verifyEqual(size(bands2), [Nb, 1]);
            testCase.verifyEqual(size(q_largest_b2), [Nv2, Nu2, Nb]);

            bands(abs(bands) < 1e-9) = 0;
            bands2(abs(bands2) < 1e-9) = 0;

            testCase.verifyEqual(bands, bands2, "RelTol", 5e-2);
        end

        function translation_invariance_test(testCase, A, bc, Nb, V_impl, method, eigs_impl)
            Nu = method.Nu;
            Nv = method.Nv;
            pde_fn_2 = method.pde_fn_2;

            TU = linspace_ex(0, 1, Nu);
            TV = linspace_ex(0, 1, Nv);
            [UU, VV] = meshgrid(TU, TV);

            V = V_impl(UU, VV);

            [bands, q_largest_b] = pde_fn_2(A, V, bc, Nb, eigs_impl);

            testCase.verifyEqual(size(bands), [Nb, 1]);
            testCase.verifyEqual(size(q_largest_b), [Nv, Nu, Nb]);

            TU2 = linspace_ex(-0.25, 0.75, Nu);
            TV2 = linspace_ex(-0.5, 0.5, Nv);
            [UU2, VV2] = meshgrid(TU2, TV2);

            V2 = V_impl(UU2, VV2);

            [bands2, q_largest_b2] = pde_fn_2(A, V2, bc, Nb, eigs_impl);

            testCase.verifyEqual(size(bands2), [Nb, 1]);
            testCase.verifyEqual(size(q_largest_b2), [Nv, Nu, Nb]);

            bands(abs(bands) < 1e-9) = 0;
            bands2(abs(bands2) < 1e-9) = 0;

            testCase.verifyEqual(bands, bands2, "RelTol", 1e-6);
        end

        function pde_fft_test(testCase, A, bc, Nb, V_impl, method, eigs_impl)
            e = exp(1);

            k = bc.k;

            Nu = method.Nu;
            Nv = method.Nv;
            pde_fn_2 = method.pde_fn_2;

            N = Nv * Nu;

            TU = linspace_ex(0, 1, Nu);
            TV = linspace_ex(0, 1, Nv);
            [UU, VV] = meshgrid(TU, TV);

            V = V_impl(UU, VV);

            [bands, q_largest_b] = pde_fn_2(A, V, bc, Nb, eigs_impl);

            testCase.verifyEqual(size(bands), [Nb, 1]);
            testCase.verifyEqual(size(q_largest_b), [Nv, Nu, Nb]);

            UV = cat(3, UU, VV);
            KK = tensorprod_shim(UV, 2 * pi * k, 3, 2);

            b = bands(1, 1);
            q = q_largest_b(:, :, 1) .* e.^(-j * KK);

            VF = 1 / N * fftshift(fft2(V));
            Q = 1 / N * fftshift(fft2(q));

            KU = k(1) + ((-(Nu - 1) / 2):((Nu - 1) / 2));
            KV = k(2) + ((-(Nv - 1) / 2):((Nv - 1) / 2));
            [KUU, KVV] = meshgrid(KU, KV);
            KUV = cat(3, KUU, KVV);
            BT = 2 * pi * A^(-1);
            KXY = change_basis(KUV, BT.');
            kxy_sq = sum(KXY .* KXY, 3);

            L1 = b * Q;
            L2 = -1 / 2 * kxy_sq .* Q + conv2(VF, Q, "same");

            L1(abs(L1) < 1e-12) = 0;
            L2(abs(L2) < 1e-12) = 0;

            L1_max = max(abs(L1), [], "all");
            L2_max = max(abs(L2), [], "all");

            testCase.verifyEqual(L1_max, L2_max, "RelTol", 1e-1);
            testCase.verifyEqual(L1, L2, "AbsTol", abs(L1_max + L2_max) / 2);
        end
    end
end
