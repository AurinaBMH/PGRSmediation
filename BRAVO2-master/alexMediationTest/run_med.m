clear all

% The BRAVO toolbox offers two methods for running inference on a mediation
% model (1) permutation and (2) bootstrapping. 
%
% The permutation model randomly and independently shuffles, without 
% replacement, the X, Y and M vectors. It is good for testing the 
% statistical significance of the effect; i.e., is mediation significant 
% relative to the empirical null? Because the shuffling is done
% independently, the null should be centred on zero.
%
% The bootstrap applies the same reshuffling to X, Y and M. The reshffule
% is done with replacement. That is, we generate a random vector of length
% N, with values ranging between 1 and N, such that any single number can
% appear more than once. We extract those values from the original vectors
% and run the model (i.e., we can pick the same value more than once to
% generate our new vector). This shuffling preserves the correlation
% structure between X, Y and M. As a result, the bootstrapped distribution
% will be centred on the observed point estimate of the effect being tested
% (e.g., on the actual mediation coefficient of the ab path) rather than on
% zero. This means that you cant really use this distribution to test the
% significance of an effect. It is however, useful to construct CIs around
% the observed point estimate.
%
% So it seems like the following is sensible:
%
% 1 - use the output of permutation_mediation.m to compute a pval for a
% given effect (e.g., is mediation significant or not?)
% 
% 2 - use the output of bootstrap_mediation.m to construct CIs around the
% point estimate of the effect.
%
% Note that pvals and CIs can be estimated using either ci_percentile.m or 
% ci_bca.m. 
%
% The function ci_percentile.m computes the pval using a normal 
% distribution with a mean and variance take from the empirical null. The
% advantage of using the theoretical (normal) distribution is that you can
% derive more accurate pvals in the tails of the distribution. However, it
% will be inaccurate if the null is not normal. NB: the original code had
% an error such that the pval was computed in the wrong direction. a '1-'
% was added to correct this, so use the function ci_permutation_BF.m
%
% The function ci_bca.m computes the pval and CIs using the empirical nulls
% and a method that corrects for biases in the distribution. This is
% probably a better one to use.
%
% Note that in output:
% a = Relationship between X and each mediator variable, M
%
% b = Relationship between mediator variable and Y. Bravo outputs B', which
% is the path controlling for c.
%
% c = Relationship between X and Y (the direct pathway). 
%
% c' = the c path controlling for b.
%
% ab = Relationship between X and Y via each mediator M (ie.,the indirect
% path).


%--------------------------------------------------------------------------
% Prep data
%--------------------------------------------------------------------------

% load data
load tarrant

% combine into single matrix
data = [icv caars snp age sex];

% get row indices of missing data
inds = [];
for i = 1:size(data,2)
    inds = [inds; find(data(:,i)==-999)];
end
inds = unique(inds);

% remove all subjects with missing data
data_clean = data;
data_clean(inds,:) = [];

% column names of data
data_names = {'icv', 'caars', 'snp', 'age', 'sex'};

% delete orig variables
clear icv caars snp age sex

%--------------------------------------------------------------------------
% Define variables
%--------------------------------------------------------------------------

% independent variable (snps)
X = data_clean(:,3);
% dependent variable (ADHD Sx measured with CAARS)
Y = data_clean(:,2);
% mediator variable (RT variability, or ICV)
M = data_clean(:,1);
% covariates (age and sex)
C = data_clean(:,4:5);
% no moderators
W = [];

%--------------------------------------------------------------------------
% Run mediation using permutations
%--------------------------------------------------------------------------

% run mediation model
[coeffs_perm, sims_perm] = permutation_mediation(X,Y,M,W,C);

% Get CIs and pvals using percentile - assumes normal distribution?? 
[CI_perc_perm p_perc_perm] = ci_percentile_BF(coeffs_perm.ab, sims_perm.ab, .05);

% Get CIs and pvals using bias-corrected non-parametric approach
[CI_bca_perm p_bca_perm] = ci_bca(coeffs_perm.ab, sims_perm.ab, .05);

%--------------------------------------------------------------------------
% Run mediation using bootstrapping
%--------------------------------------------------------------------------

% run mediation model
[coeffs_boot, sims_boot] = bootstrap_mediation(X,Y,M,W,C);

% Get CIs and pvals using percentile - assumes normal distribution?? 
[CI_perc_boot p_perc_boot] = ci_percentile_BF(coeffs_boot.ab, sims_boot.ab, .05);

% Get CIs and pvals using bias-corrected non-parametric approach
[CI_bca_boot p_bca_boot] = ci_bca(coeffs_boot.ab, sims_boot.ab, .05);





