function    stats = connectomeRegression(conn, design, contrast, sides)

% stats = connectomeRegression(conn, design, contrast)
%
% This function will fit a GLM to each edge of a connectivity matrix. The
% model is formulated as used in the NBS. P-vals are one-tailed.
%
% -------
% INPUTS:
% -------
% conn  - N*N*M matrix where N is number of nodes and M is number of
%       subjects. Note that NaN values will be set to zero. This may or may
%       not be desirable.
%
% design - M*K matrix, where K is number of regressors. Should always
%        include a constant term
%
% contrast - 1*K vector of contrast weights.
%
% sides - 1 = one-tailed test; 2 = two-tailed test. 1-tailed test is
%         performed using a t-stat. 2-tailed is performed using an F-stat,
%         computed as t^2, distributed with K-1 and M-K df. This will work
%         fine for simple contrasts but may not be appropriate for
%         interactions. Default = 2;
%
% -------
% OUTPUTS:
% -------
% stats.t - N*N matrix of t-statistics for the contrast of interest
%
% stats.r2 - N*N matrix of r^2 value for contrast. convert from t-stat as
%           per https://afni.nimh.nih.gov/sscc/gangc/tr.html
%
% stats.p - N*N matrix of p-values.
% 
%==========================================================================

% set default to two-tailed.
if nargin<4
    sides = 2;
end

% number of nodes
nNodes = size(conn,1);

% number of subjects
nSubs = size(conn,3);

% upper triangle indices
upperInds = find(triu(ones(nNodes),1));

% total number of covariates
nReg = size(design,2);

% vectorize each subject's matrix
for sub = 1:nSubs
    a = conn(:,:,sub);
    y(sub,:) = a(upperInds);
end

% get NaN stats and set to zero
nans = isnan(y);
nans = sum(nans);
nansMat = zeros(nNodes);
nansMat(upperInds) = nans;
nansMat = nansMat + nansMat';
y(isnan(y)) = 0;

% number of edges
nEdges = size(y,2);

% degrees of freedom
df = nSubs - nReg;

% estimate beta weights    
beta = design\y;

% get residuals
resid = zeros(nSubs, nEdges);
resid = y - design*beta;

% get mean-squared error
mse = zeros(nSubs, nEdges);
mse = sum(resid.^2)/df;

% get standard error
se = sqrt(mse*(contrast*inv(design'*design)*contrast'));

% get t-statistic
t = (contrast*beta)./se;

% compute r^2 statistic
r2 = (t.^2)./((t.^2) + df);                      % conversion based on https://afni.nimh.nih.gov/sscc/gangc/tr.html 

if sides == 1

    % get pval for each t-stat
    df = repmat(df, size(t,1), size(t,2));     % generate repeated df vector
    p = tcdf(-abs(t),df);
    stat = t;
    TorF = 'Tstat';
    
elseif sides == 2
    
    stat = abs(t).^2;
    df1 = repmat(nReg-1, size(stat,1), size(stat,2));
    df2 = repmat(nSubs-nReg, size(stat,1), size(stat,2));
    p = 1- fcdf(abs(stat), df1, df2);
    TorF = 'Fstat';
    
end

% store stats and pvals in matrix form.
% initialize
stats.stat = zeros(nNodes);
stats.r2 = zeros(nNodes);
stats.p = zeros(nNodes);

% store
stats.stat(upperInds) = t;
stats.r2(upperInds) = r2;
stats.p(upperInds) = p;

% symmetrize
stats.stat = stats.stat+stats.stat';
stats.r2 = stats.r2+stats.r2';

stats.p = stats.p+stats.p';
stats.p(1:nNodes+1:end) = 1;            % set diagonal to 1 for pval matrix

stats.df  = df(1);
stats.nans = nansMat;
stats.testStat = TorF;
stats.sides = sides;
end