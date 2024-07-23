nvec = 200:200:1000;
mvec = 2*nvec;
nref = nvec;

rho_vec = [0,0.5,0.9,0.95, 0.99];
h2 = 0.2;
nsims = 1000;
do_standardize = 0;

distbn = 'norm';

methods = {'ar1', 'equi'};

for K = 1:length(rho_vec)
    rho = rho_vec(K);
    store = cell(2, length(nvec));
    for I = 1:2
        for J = 1:length(nvec)
            J
            [~,~, OG_h2_sum, ref_h2_sum] = ...
                h2sim2( nvec(J), mvec(J), nref(J), h2, rho, methods{I}, nsims, do_standardize, distbn );
            store{I,J}.OG_h2_sum = OG_h2_sum;
            store{I,J}.ref_h2_sum = ref_h2_sum;
        end
    end
    if do_standardize
        save(['./Strongvsweak/store_rho_', num2str(10*rho), '_h2_', num2str(10*h2)], 'store', 'nvec', 'rho', 'h2', 'nsims', 'do_standardize', 'distbn')
    else
        save(['./Strongvsweak/store_rho_nostd_', num2str(10*rho), '_h2_', num2str(10*h2)], 'store', 'nvec', 'rho', 'h2', 'nsims', 'do_standardize', 'distbn')
    end
end

%% Plot the results
h2 = 0.2;
rho_vec = [0, 0.5, 0.9, 0.95, 0.99];
nvec = 200:200:1000;

ldsc_fixed = zeros(length(rho_vec),length(nvec));
ldsc_fixedW = zeros(length(rho_vec),length(nvec));
ldsc_fixedW1 = zeros(length(rho_vec),length(nvec));
ldsc_free = zeros(length(rho_vec),length(nvec));
gwash = zeros(length(rho_vec),length(nvec));
gwashW = zeros(length(rho_vec),length(nvec));
gwashW1 = zeros(length(rho_vec),length(nvec));

ldsc_fixed_std = zeros(length(rho_vec),length(nvec));
ldsc_fixedW_std = zeros(length(rho_vec),length(nvec));
ldsc_fixedW1_std = zeros(length(rho_vec),length(nvec));
ldsc_free_std = zeros(length(rho_vec),length(nvec));
gwash_std = zeros(length(rho_vec),length(nvec));
gwashW_std = zeros(length(rho_vec),length(nvec));
gwashW1_std = zeros(length(rho_vec),length(nvec));

ldsc_fixed_mse = zeros(length(rho_vec),length(nvec));
ldsc_fixedW_mse = zeros(length(rho_vec),length(nvec));
ldsc_fixedW1_mse = zeros(length(rho_vec),length(nvec));
ldsc_free_mse = zeros(length(rho_vec),length(nvec));
gwash_mse = zeros(length(rho_vec),length(nvec));
gwashW_mse = zeros(length(rho_vec),length(nvec));
gwashW1_mse = zeros(length(rho_vec),length(nvec));

do_ref = 1;

options = {'ldsc_free', 'ldsc_fixed_intercept', 'ld_ratio'};
if do_ref == 1
    storetype = 'ref_h2_sum';
else
    storetype = 'OG_h2_sum';
end

method = 'equi';
if strcmp(method, 'ar1')
    methodtype = 1; 
elseif strcmp(method, 'equi')
    methodtype = 2; 
end

for I = 1:length(rho_vec)
    rho = rho_vec(I);
    if rho > 0.9
        load(['/Users/sdavenport/Documents/Code/Servers/jobs/Gen2/Strongvsweak/store_rho', num2str(100*rho_vec(I)),'_h2_',  num2str(10*h2)]);
    else
        load(['/Users/sdavenport/Documents/Code/Servers/jobs/Gen2/Strongvsweak/store_rho', num2str(10*rho_vec(I)),'_h2_', num2str(10*h2)]);
    end
    for J = 1:length(store)
        ldsc_fixed(I,J) = store{methodtype,J}.(storetype).ldsc_fixed_intercept.mean;
        ldsc_fixedW(I,J) = store{methodtype,J}.(storetype).ldsc_fixed_interceptW.mean;
        ldsc_fixedW1(I,J) = store{methodtype,J}.(storetype).ldsc_fixed_interceptW1.mean;
        ldsc_free(I,J) = store{methodtype,J}.(storetype).ldsc_free.mean;
        gwash(I,J) = store{methodtype,J}.(storetype).gwash.mean;
        gwashW(I,J) = store{methodtype,J}.(storetype).gwashW.mean;
        gwashW1(I,J) = store{methodtype,J}.(storetype).gwashW1.mean;

        ldsc_fixed_std(I,J) = store{methodtype,J}.(storetype).ldsc_fixed_intercept.std;
        ldsc_fixedW_std(I,J) = store{methodtype,J}.(storetype).ldsc_fixed_interceptW.std;
        ldsc_fixedW1_std(I,J) = store{methodtype,J}.(storetype).ldsc_fixed_interceptW1.std;
        ldsc_free_std(I,J) = store{methodtype,J}.(storetype).ldsc_free.std;
        gwash_std(I,J) = store{methodtype,J}.(storetype).gwash.std;
        gwashW_std(I,J) = store{methodtype,J}.(storetype).gwashW.std;
        gwashW1_std(I,J) = store{methodtype,J}.(storetype).gwashW1.std;

        ldsc_fixed_mse(I,J) = store{methodtype,J}.(storetype).ldsc_fixed_intercept.std.^2 + (ldsc_fixed(I,J) - h2).^2;
        ldsc_fixedW_mse(I,J) = store{methodtype,J}.(storetype).ldsc_fixed_interceptW.std.^2 + (ldsc_fixedW(I,J) - h2).^2;
        ldsc_fixedW1_mse(I,J) = store{methodtype,J}.(storetype).ldsc_fixed_interceptW1.std.^2 + (ldsc_fixedW1(I,J) - h2).^2;
        ldsc_free_mse(I,J) = store{methodtype,J}.(storetype).ldsc_free.std.^2 + (ldsc_free(I,J) - h2).^2;
        gwash_mse(I,J) = store{methodtype,J}.(storetype).gwash.std.^2 + (gwash(I,J) - h2).^2;
        gwashW_mse(I,J) = store{methodtype,J}.(storetype).gwashW.std.^2 + (gwashW(I,J) - h2).^2;
        gwashW1_mse(I,J) = store{methodtype,J}.(storetype).gwashW1.std.^2 + (gwashW1(I,J) - h2).^2;
    end
end

%% STD plot
% for I = 1:length(rho_vec)
%     plot(ldsc_fixedW1_std(I, :))
%     hold on
% end
seethrough = [0.4, 0.6, 0.75, 0.9, 1];
for I = methodtype:length(rho_vec)
    plot(nvec, gwashW_std(I, :), 'Color', [0.4, 0.4, 1]*seethrough(I), 'LineWidth', 4);
    hold on
end
if methodtype == 1
    legend('\rho = 0', '\rho = 0.5', '\rho = 0.9', '\rho = 0.95', '\rho = 0.99', 'Location', 'NorthEast')
else
    legend('\rho = 0.5', '\rho = 0.9','\rho = 0.95', '\rho = 0.99', 'Location', 'SouthEast')
end

matniceplot
% fullscreen
BigFont(22)
% if h2 == 0.2
%     legend('LDSC W1', 'LDSC', 'LDSC W', 'GWASH W1', 'GWASH', 'GWASH W')
% end
% ylim([0,0.2])
ylim([0,0.25])
xlabel('n')
ylabel('Standard Deviation')
if methodtype == 1
    methodtitle = upper(method);
else
    methodtitle = method;
end
title(['h^2 = ', num2str(h2), ', ', methodtitle])
saveim(['strongvsweak_std_method_', method, '_h2_', num2str(10*h2), '.pdf'])

%% Bias plot
seethrough = [0.4, 0.6, 0.75, 0.9, 1];
for I = methodtype:length(rho_vec)
    plot(nvec, gwashW(I, :) - h2, 'Color', [0.4, 0.4, 1]*seethrough(I), 'LineWidth', 4);
    hold on
end
plot(nvec, 0*nvec, '--', 'LineWidth', 2, 'Color', ones(1,3)*0.4)
% if methodtype == 1
%     legend('\rho = 0', '\rho = 0.5', '\rho = 0.9','\rho = 0.99', 'Location', 'SouthEast')
% end
matniceplot
% fullscreen
BigFont(22)
% if h2 == 0.2
%     legend('LDSC W1', 'LDSC', 'LDSC W', 'GWASH W1', 'GWASH', 'GWASH W')
% end
% ylim([0,0.2])
if methodtype == 1
    methodtitle = upper(method);
else
    methodtitle = method;
end
ylim([-0.05, 0.01])
xlabel('n')
ylabel('Bias')
title(['h^2 = ', num2str(h2), ', ', methodtitle])
saveim(['strongvsweak_bias_method_', method, '_h2_', num2str(10*h2), '.pdf'])

%% STD plot
seethrough = [0.4, 0.6, 0.75, 0.9, 1];
for I = methodtype:length(rho_vec)
    plot(nvec, ldsc_fixedW1_std(I, :), 'Color', [1, 0.4, 0.4]*seethrough(I), 'LineWidth', 4);
    hold on
end
if methodtype == 1
    legend('\rho = 0', '\rho = 0.5', '\rho = 0.9', '\rho = 0.95', '\rho = 0.99', 'Location', 'NorthEast')
else
    legend('\rho = 0.5', '\rho = 0.9','\rho = 0.95', '\rho = 0.99', 'Location', 'SouthEast')
end

matniceplot
% fullscreen
BigFont(22)
% if h2 == 0.2
%     legend('LDSC W1', 'LDSC', 'LDSC W', 'GWASH W1', 'GWASH', 'GWASH W')
% end
% ylim([0,0.2])
ylim([0,0.25])
xlabel('n')
ylabel('Standard Deviation')
if methodtype == 1
    methodtitle = upper(method);
else
    methodtitle = method;
end
title(['h^2 = ', num2str(h2), ', ', methodtitle])
saveim(['strongvsweak_ldsc_std_method_', method, '_h2_', num2str(10*h2), '.pdf'])

%% Bias plot
seethrough = [0.4, 0.6, 0.75, 0.9, 1];
for I = methodtype:length(rho_vec)
    plot(nvec, ldsc_fixedW1(I, :) - h2, 'Color', [1, 0.4, 0.4]*seethrough(I), 'LineWidth', 4);
    hold on
end
plot(nvec, 0*nvec, '--', 'LineWidth', 2, 'Color', ones(1,3)*0.4)
% if methodtype == 1
%     legend('\rho = 0', '\rho = 0.5', '\rho = 0.9','\rho = 0.99', 'Location', 'SouthEast')
% end
matniceplot
% fullscreen
BigFont(22)
% if h2 == 0.2
%     legend('LDSC W1', 'LDSC', 'LDSC W', 'GWASH W1', 'GWASH', 'GWASH W')
% end
% ylim([0,0.2])
if methodtype == 1
    methodtitle = upper(method);
else
    methodtitle = method;
end
ylim([-0.05, 0.01])
xlabel('n')
ylabel('Bias')
title(['h^2 = ', num2str(h2), ', ', methodtitle])
saveim(['strongvsweak_ldsc_bias_method_', method, '_h2_', num2str(10*h2), '.pdf'])