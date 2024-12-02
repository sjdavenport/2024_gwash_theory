nvec = 200:200:1000;
mvec = 2*nvec;
nref = nvec;

rho = [0.2, 0.9];
nsims = 10000;

for h2 = [0.2,0.5]
    for do_standardize = [0,1]

        saveloc = '/Users/sdavenport/Documents/Code/MATLAB/Papers/2024_gwash_theory/MixtureResults/';

        distbn = 'norm';
        store = cell(1, length(nvec));
        for J = 1:length(nvec)
            J
            [~,~, OG_h2_sum, ref_h2_sum] = ...
                h2sim2( nvec(J), mvec(J), nref(J), h2, rho, 'ar1mix', nsims, do_standardize, distbn );
            store{J}.OG_h2_sum = OG_h2_sum;
            store{J}.ref_h2_sum = ref_h2_sum;
        end

        save(['./store_armix_gauss_rho29_h2_', num2str(10*h2), '_dostd_', num2str(do_standardize), '_nsims_', num2str(nsims)], 'store', 'nvec', 'rho', 'h2', 'nsims', 'do_standardize', 'distbn')
    end
end

%% Plot the results
for h2 = [0.2, 0.5]
    for do_std = [0,1]

        % if do_std == 1
        %     load(['/Users/sdavenport/Documents/Code/MATLAB/Papers/2024_gwash_theory/MixtureResults/store_armix_gauss_rho29_h2_', num2str(10*h2)]);
        % else
        %     load(['/Users/sdavenport/Documents/Code/MATLAB/Papers/2024_gwash_theory/MixtureResults/store_armix_gauss_rho29_h2_', num2str(10*h2), '_dostd_0']);
        % end
        load(['/Users/sdavenport/Documents/Code/MATLAB/Papers/2024_gwash_theory/MixtureResults/store_armix_gauss_rho29_h2_', num2str(10*h2), '_dostd_', num2str(do_std), '_nsims_10000']);

        ldsc_fixed = zeros(1,length(store));
        ldsc_fixedW = zeros(1,length(store));
        ldsc_fixedW1 = zeros(1,length(store));
        ldsc_free = zeros(1,length(store));
        gwash = zeros(1,length(store));
        gwashW = zeros(1,length(store));
        gwashW1 = zeros(1,length(store));


        ldsc_fixed_std = zeros(1,length(store));
        ldsc_fixedW_std = zeros(1,length(store));
        ldsc_fixedW1_std = zeros(1,length(store));
        ldsc_free_std = zeros(1,length(store));
        gwash_std = zeros(1,length(store));
        gwashW_std = zeros(1,length(store));
        gwashW1_std = zeros(1,length(store));

        ldsc_fixed_mse = zeros(1,length(store));
        ldsc_fixedW_mse = zeros(1,length(store));
        ldsc_fixedW1_mse = zeros(1,length(store));
        ldsc_free_mse = zeros(1,length(store));
        gwash_mse = zeros(1,length(store));
        gwashW_mse = zeros(1,length(store));
        gwashW1_mse = zeros(1,length(store));

        do_ref = 1;

        options = {'ldsc_free', 'ldsc_fixed_intercept', 'ld_ratio'};
        if do_ref == 1
            storetype = 'ref_h2_sum';
        else
            storetype = 'OG_h2_sum';
        end

        for J = 1:length(store)
            ldsc_fixed(J) = store{J}.(storetype).ldsc_fixed_intercept.mean;
            ldsc_fixedW(J) = store{J}.(storetype).ldsc_fixed_interceptW.mean;
            ldsc_fixedW1(J) = store{J}.(storetype).ldsc_fixed_interceptW1.mean;
            ldsc_free(J) = store{J}.(storetype).ldsc_free.mean;
            gwash(J) = store{J}.(storetype).gwash.mean;
            gwashW(J) = store{J}.(storetype).gwashW.mean;
            gwashW1(J) = store{J}.(storetype).gwashW1.mean;

            ldsc_fixed_std(J) = store{J}.(storetype).ldsc_fixed_intercept.std;
            ldsc_fixedW_std(J) = store{J}.(storetype).ldsc_fixed_interceptW.std;
            ldsc_fixedW1_std(J) = store{J}.(storetype).ldsc_fixed_interceptW1.std;
            ldsc_free_std(J) = store{J}.(storetype).ldsc_free.std;
            gwash_std(J) = store{J}.(storetype).gwash.std;
            gwashW_std(J) = store{J}.(storetype).gwashW.std;
            gwashW1_std(J) = store{J}.(storetype).gwashW1.std;

            ldsc_fixed_mse(J) = store{J}.(storetype).ldsc_fixed_intercept.std.^2 + (ldsc_fixed(J) - h2).^2;
            ldsc_fixedW_mse(J) = store{J}.(storetype).ldsc_fixed_interceptW.std.^2 + (ldsc_fixedW(J) - h2).^2;
            ldsc_fixedW1_mse(J) = store{J}.(storetype).ldsc_fixed_interceptW1.std.^2 + (ldsc_fixedW1(J) - h2).^2;
            ldsc_free_mse(J) = store{J}.(storetype).ldsc_free.std.^2 + (ldsc_free(J) - h2).^2;
            gwash_mse(J) = store{J}.(storetype).gwash.std.^2 + (gwash(J) - h2).^2;
            gwashW_mse(J) = store{J}.(storetype).gwashW.std.^2 + (gwashW(J) - h2).^2;
            gwashW1_mse(J) = store{J}.(storetype).gwashW1.std.^2 + (gwashW1(J) - h2).^2;
        end

        % STD plot
        clf
        f = plot(nvec, ldsc_fixedW1_std, 'LineWidth', 4, 'Color', [1,0.65,0.65])
        hold on
        plot(nvec, ldsc_fixed_std, 'LineWidth', 4, 'Color', [1,0.4,0.4])
        % plot(nvec, ldsc_fixedW_std, 'LineWidth', 4, 'Color', [1,0,0])

        % plot(nvec, gwashW1_std, 'LineWidth', 4, 'Color', [0.65,0.65,1])
        plot(nvec, gwash_std, 'LineWidth', 4, 'Color', [0.4,0.4,1])
        plot(nvec, gwashW_std, '--', 'LineWidth', 4, 'Color', [0,0,1])

        matniceplot
        % fullscreen
        BigFont(22)
        if h2 == 0.2
            legend('LDSC', 'LDSC-W', 'GWASH', 'GWASH-W', 'Location', 'NorthEast')
        end
        ylim([0,0.35])
        setylim
        xlabel('n')
        ylabel('Standard Deviation')
        title(['h^2 = ', num2str(h2)])
        if do_std == 1
            saveim(['std_', num2str(10*h2), '.pdf'])
        else
            saveim(['std_', num2str(10*h2), '_std_0.pdf'])
        end

        % Bias plot
        scalefactor = 1;
        plot(nvec, (ldsc_fixedW1 - h2)*scalefactor, 'LineWidth', 4, 'Color', [1,0.65,0.65])
        hold on
        plot(nvec, (ldsc_fixed - h2)*scalefactor, 'LineWidth', 4, 'Color', [1,0.4,0.4])
        % plot(nvec, (ldsc_fixedW - h2)*scalefactor, 'LineWidth', 4, 'Color', [1,0,0])

        % plot(nvec, (gwashW1 - h2)*scalefactor, 'LineWidth', 4, 'Color', [0.65,0.65,1])
        plot(nvec, (gwash - h2)*scalefactor,  'LineWidth', 4, 'Color', [0.4,0.4,1])
        plot(nvec, (gwashW - h2)*scalefactor, '--', 'LineWidth', 4, 'Color', [0,0,1])

        plot(nvec, 0*nvec, '--', 'LineWidth', 2, 'Color', ones(1,3)*0.4)

        ylim([-0.05, 0.04])
        matniceplot
        % fullscreen
        BigFont(22)
        if h2 == 0.2
            % legend('LDSC W1', 'LDSC', 'LDSC W', 'GWASH W1', 'GWASH', 'GWASH W', 'Location', 'SouthEast')
            legend('LDSC', 'LDSC-W', 'GWASH', 'GWASH-W', 'Location', 'SouthEast')
        end
        xlabel('n')
        ylabel('Bias')
        title(['h^2 = ', num2str(h2)])
        set(gca, 'position', [0.1728    0.1385    0.7322    0.7830])
        % ytickformat('%.0e')
        if do_std == 1
            saveim(['bias_', num2str(10*h2), '.pdf'])
        else
            saveim(['bias_', num2str(10*h2), '_std_0.pdf'])
        end

    end
end
% ylim([0,0.2])
%% MSE plot
plot(nvec, ldsc_fixedW1_mse, 'LineWidth', 4, 'Color', [1,0.65,0.65])
hold on
plot(nvec, ldsc_fixed_mse, 'LineWidth', 4, 'Color', [1,0.4,0.4])
plot(nvec, ldsc_fixedW_mse, 'LineWidth', 4, 'Color', [1,0,0])

plot(nvec, gwashW1_mse, 'LineWidth', 4, 'Color', [0.65,0.65,1])
plot(nvec, gwash_mse, 'LineWidth', 4, 'Color', [0.4,0.4,1])
plot(nvec, gwashW_mse, '--', 'LineWidth', 4, 'Color', [0,0,1])

matniceplot
fullscreen
BigFont(25)
legend('LDSC W1', 'LDSC', 'LDSC W', 'GWASH W1', 'GWASH', 'GWASH W')
xlabel('n')
ylabel('MSE')
title(['rho = [0.2, 0.9], h2 = ', num2str(h2)])
