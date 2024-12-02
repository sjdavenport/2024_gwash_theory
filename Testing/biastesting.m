[OG_h2, ref_h2, OG_h2_sum, ref_h2_sum] = h2sim2( 100, 200, 100, 0.2, 0, 'ar1', 1000, 0, 'norm' );
          
%%
histogram(ref_h2.ldsc_fixed_intercept)

%%
n = 100; m = 200; h2 = 0.2; rho = 0; method = 'ar1'; do_standardize = 0; distbn = 'norm';
[ chi2, X ] = gengenmodel( n, m, h2, rho, method, do_standardize, distbn);