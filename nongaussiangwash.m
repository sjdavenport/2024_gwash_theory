nvec = 200:200:1000;
mvec = 2*nvec;
nref = nvec;

for rho = [0.2, 0.9]
    % dfvec = {1,2,3,5, 'norm'};
    dfvec = {2,2.3, 2.5, 2.8, 3,5, 'norm'};
    h2 = 0.2;
    nsims = 10000;
    for do_standardize = [0,1]
        % saveloc = '/Users/sdavenport/Documents/Code/Servers/jobs/Gen2';

        for I = 1:length(dfvec)
            I
            if ~strcmp(dfvec{I}, 'norm')
                distbn = ['t', num2str(dfvec{I})];
            else
                distbn = 'norm';
            end
            store = cell(1, length(nvec));
            for J = 1:length(nvec)
                J
                [~,~, OG_h2_sum, ref_h2_sum] = ...
                    h2sim2( nvec(J), mvec(J), nref(J), h2, rho, 'ar1', nsims, do_standardize, distbn );
                store{J}.OG_h2_sum = OG_h2_sum;
                store{J}.ref_h2_sum = ref_h2_sum;
            end
            save(['./store_ar_distbn_', distbn, '_rho', num2str(10*rho), '_h2_', num2str(10*h2), '_dostd_', num2str(do_standardize),'.mat'], 'store', 'nvec', 'rho', 'h2', 'nsims', 'do_standardize', 'distbn')
        end
    end
end

%%
save('./store_nongauss', 'store', 'nvec', 'rho', 'h2', 'dfvec', 'nsims')

%% Normal ones
norm_store = cell(1,length(nvec));
for J = 1:length(nvec)
    J
    [~,~, OG_h2_sum, ref_h2_sum] = ...
        h2sim( nvec(J), mvec(J), nref(J), h2, rho, 'ar1', nsims, do_standardize, 'norm' );
    norm_store{1,J}.OG_h2_sum = OG_h2_sum;
    norm_store{1,J}.ref_h2_sum = ref_h2_sum;
end
%%
save('./store_gauss', 'norm_store', 'nvec', 'rho', 'h2', 'dfvec', 'nsims')

%%
load('./store_nongauss')
load('./store_gauss')
for I = 1:length(norm_store)
    store{5,I} = norm_store{1,I};
end

%%
ldsc_fixed = zeros(size(store,1), size(store,2));
ldsc_free = zeros(size(store,1), size(store,2));
ld_ratio = zeros(size(store,1), size(store,2));

ldsc_fixed_ss = zeros(size(store,1), size(store,2));
ldsc_free_ss = zeros(size(store,1), size(store,2));
ld_ratio_ss = zeros(size(store,1), size(store,2));

do_ref = 1;

options = {'ldsc_free', 'ldsc_fixed_intercept', 'ld_ratio'};

for I = 1:size(store,1)
    for J = 1:size(store,2)
        ldsc_fixed(I,J) = store{I,J}.ref_h2_sum.ldsc_fixed_intercept.mean;
        ldsc_free(I,J) = store{I,J}.ref_h2_sum.ldsc_free.mean;
        ld_ratio(I,J) = store{I,J}.ref_h2_sum.ld_ratio.mean;

        ldsc_fixed_ss(I,J) = store{I,J}.ref_h2_sum.ldsc_fixed_intercept.ss;
        ldsc_free_ss(I,J) = store{I,J}.ref_h2_sum.ldsc_free.ss;
        ld_ratio_ss(I,J) = store{I,J}.ref_h2_sum.ld_ratio.ss;
    end
end

ldsc_fixed_std = sqrt(ldsc_fixed_ss/nsims - ldsc_fixed.^2);
ldsc_free_std = sqrt(ldsc_free_ss/nsims - ldsc_free.^2);
ld_ratio_std = sqrt(ld_ratio_ss/nsims - ld_ratio.^2);

%%
for I = 1:5
    plot(nvec, gwash_std(I,:))
    hold on
end
legend('3', '5', '10', '20', 'norm')
matniceplot
ylim([0,0.12])
BigFont(25)
fullscreen

%%
do_ref = 1;

options = {'ldsc_free', 'ldsc_fixed_intercept', 'ld_ratio'};

dfvec = {2,2.3, 2.5, 2.8, 3,5, 'norm'};
nvec = 200:200:1000;

ldsc_fixed = zeros(length(dfvec), length(nvec));
ldsc_free = zeros(length(dfvec), length(nvec));
gwash = zeros(length(dfvec), length(nvec));

ldsc_fixed_std = zeros(length(dfvec), length(nvec));
ldsc_free_std = zeros(length(dfvec), length(nvec));
gwash_std = zeros(length(dfvec), length(nvec));

h2 = 0.2;

for rho = [0.2, 0.9]
    for do_standardize = [0,1]

        for I = 1:length(dfvec)
            if isnumeric(dfvec{I})
                distbn = ['t', num2str(dfvec{I})];
            else
                distbn = num2str(dfvec{I});
            end
            if do_standardize == 1
                a = load(['store_ar_distbn_', distbn, '_rho', num2str(10*rho), '_h2_2.mat']);
            else
                a = load(['store_ar_distbn_', distbn, '_rho', num2str(10*rho), '_h2_2_dostd_0.mat']);
            end
            for J = 1:length(nvec)
                store = a.store{J};
                ldsc_fixed(I,J) = store.ref_h2_sum.ldsc_fixed_intercept.mean;
                ldsc_free(I,J) = store.ref_h2_sum.ldsc_free.mean;
                gwash(I,J) = store.ref_h2_sum.gwash.mean;

                ldsc_fixed_std(I,J) = store.ref_h2_sum.ldsc_fixed_intercept.std;
                ldsc_free_std(I,J) = store.ref_h2_sum.ldsc_free.std;
                gwash_std(I,J) = store.ref_h2_sum.gwash.std;
            end
        end

        % gwash_mse = gwash_std.^2 + (gwash-0.2).^2;
        for do_gwash = [0,1]
            if do_gwash == 1
                std_mate = gwash_std;
                values_mate = gwash;
                start_col = [0,0,1];
                inc_col = [1,1,0];
            else
                std_mate = ldsc_fixed_std;
                values_mate = ldsc_fixed;
                start_col = [1, 0, 0];
                inc_col = [0,1,1];
            end

            % gradvalues = fliplr([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]);
            gradvalues = fliplr([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]);


            for I = 1:length(dfvec)
                % plot(nvec, gwash_std(I,:), 'Color', gradvalues(I)*[1,1,0] + [0,0,1], 'LineWidth', 4)
                if I == 7
                    plot(nvec, std_mate(I,:), '--', 'Color', [0,0,0], 'LineWidth', 4)
                else
                    plot(nvec, std_mate(I,:), 'Color', gradvalues(I+1)*inc_col + start_col, 'LineWidth', 4)
                end
                hold on
                BigFont(22)
                xlabel('n')
                ylabel('Standard Deviation')
            end
            if do_standardize == 1
                ylim([0,0.18])
            else
                ylim([0,1])
            end

            title(['h^2 = 0.2, \rho = ', num2str(rho)])
            % hLegend = legend('df = 2','df = 2.3', 'df = 2.5', 'df = 2.8', 'df = 3', 'df = 5', 'Normal', 'Location', 'southoutside', 'Orientation', 'Horizontal');
            if do_standardize == 0 && rho == 0.9
                hLegend = legend('df = 2','df = 2.3', 'df = 2.5', 'df = 2.8', 'df = 3', 'df = 5', 'Normal', 'Location', 'NorthWest');
            end
            % legendPos = get(hLegend, 'Position');
            % % Calculate new position for the legend to be underneath the plot
            % newLegendPos = legendPos;
            % newLegendPos(2) = -1;
            % set(gcf, 'Color', 'white');

            matniceplot
            if do_gwash
                saveim(['varydf_std_rho', num2str(10*rho), '_dostd_', num2str(do_standardize), '.pdf'])
            else
                saveim(['varydf_std_rho', num2str(10*rho), '_dostd_', num2str(do_standardize), '_ldsc.pdf'])
            end

            clf
            % gradvalues = fliplr([0.4, 0.45,0.5,0.55, 0.6, 0.65, 0.9]);
            for I = 1:length(dfvec)
                if I == 7
                    plot(nvec, values_mate(I,:) - h2, '--', 'Color', [0,0,0], 'LineWidth', 4)
                else
                    plot(nvec, values_mate(I,:) - h2, 'Color', gradvalues(I+1)*inc_col + start_col, 'LineWidth', 4)
                end
                % hold on
                % plot(nvec, 0*nvec, '--', 'LineWidth', 2, 'Color', ones(1,3)*0.4)
                hold on
                BigFont(22)
                xlabel('n')
                % ylabel('Expectation')
                ylabel('Bias')
            end
            % if rho == 0.9
            %     legend('df = 2','df = 2.3', 'df = 2.5', 'df = 2.8', 'df = 3', 'df = 5', 'Normal', 'Location', 'South')
            % end
            title(['h^2 = 0.2, \rho = ', num2str(rho)])
            matniceplot
            if do_standardize == 1
                ylim([-0.25, 0.1])
            else
                ylim([-0.6, 0.1])
            end
            if do_gwash
                saveim(['varydf_bias_rho', num2str(10*rho), '_dostd_', num2str(do_standardize), '.pdf'])
            else
                saveim(['varydf_bias_rho', num2str(10*rho), '_dostd_', num2str(do_standardize), '_ldsc.pdf'])
            end
            clf
        end
    end
end