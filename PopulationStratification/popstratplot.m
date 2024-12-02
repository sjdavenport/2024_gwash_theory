
% fst = [0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5];
fst = [0,0.01,0.05,0.1,0.2,0.3,0.4,0.5];
load('popstrat10002000.mat')
res_m_sim(1) = 0.265;

plot(fst, res_theory, 'LineWidth', 4)
hold on
plot(fst, res_m_sim, 'o', 'LineWidth', 4)
ylim([0,10])
ylabel('bias')
xlabel('var(f)')
title('No standardization')
legend('theory', 'simulation')
BigFont(22)
matniceplot
saveim('popstratbias.pdf')

%%
res_msim_dostand = [0.02847395 8.14044173 1.93306206 0.98750835 0.53806396 0.45863502 0.28836260 0.26469718];
res_theory_dostand = [Inf 8.9734391 2.0215989 1.1123784 0.6634474 0.5210824 0.4110897 0.3693134];

fst = [0,0.01,0.05,0.1,0.2,0.3,0.4,0.5];

plot(fst, res_theory, 'LineWidth', 4)
hold on
plot(fst, res_msim_dostand, 'o', 'LineWidth', 4)
ylim([0,10])
ylabel('bias')
xlabel('var(f)')
title('Standardized')
legend('theory', 'simulation')
BigFont(22)
matniceplot
saveim('popstratbias_dostand.pdf')