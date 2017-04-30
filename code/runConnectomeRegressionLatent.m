clear all
close all

%==========================================================================
% INPUTS
%==========================================================================

% base directory
baseDir = '/Users/Aurina/GoogleDrive/Genetics_connectome/PGRSmediation/';
inDir = [baseDir,'results/matchedData/']; 

% Input mat files. These will be the output of matchData.m.
% Since only latent variable and connectivity data are analysed here, the
% disorder of type of PGRS processing does not matter. The key thing will
% be to filter based on desired parcellation, tractography algorithm,
% sifting and connection weight.
inFiles = dir([inDir,'*SCZ*FACT*aparcANDfirst*SIFT2*standard.mat']);

% names of latent variables to use in output
latentNames = {'pFact', 'imp', 'sczPos','sczNeg','aff'};

% contrast weight - do you want to look at positive or negative
% associations?
conWeight = -1 ;

% sides = 1 for 1-tailed (t-test); sides = 2 for 2-tailed (F-test) test
sides = 1;

% directory for output
outDir = [baseDir,'results/regression/psych-conn/'];

% filename for output
outName = 'regResults_connVlatent_1sided_';

%==========================================================================
% Analysis
%==========================================================================

for f = 1:length(inFiles)

    
    % Preliminaries
    %----------------------------------------------------------------------
    
    % load data
    load([inDir,inFiles(f).name]);
    
    % prep connectivity data
    conn = matchedData.mri.conn;
    conn(isnan(conn)) = 0;                      % set NaNs to zero

    % prep design matrix
    cov = matchedData.mri.cov(:,2:end);
    A = cell2mat(cov(:,2));
    A = str2num(A);
    S = cell2mat(cov(:,1)); 
    cov = [S,A]; 
    
    latent = matchedData.mri.latent(:,2:end);
    design = [latent cov];

    % get numbers and indices
    nSubs = size(conn,3);                       % number of subjects
    nNodes = size(conn,1);                      % number of nodes
    nRegs = size(design,2);                     % number of predictors
    upperInds = find(triu(ones(nNodes),1));     % indices of upper triangle


    for lv = 1:length(latentNames)

        tic;
        % define contrast vector
        contrast = zeros(1, nRegs+1);                  %add one for constant
        contrast(lv) = conWeight;
        stats = connectomeRegression(conn, [design ones(nSubs,1)], contrast, sides);

        % force NaNs to 0 or 1. Arise if connectivity vector contains all
        % zeros.
        stats.stat(isnan(stats.stat)) = 0;
        stats.r2(isnan(stats.r2)) = 0;
        stats.p(isnan(stats.p)) = 1;                        % pval to max (i.e., p=1) for NaN

        % store results
        corrConnLatent(lv).latentNames = latentNames{lv};
        corrConnLatent(lv).r2 = stats.r2;
        corrConnLatent(lv).stat = stats.stat;
        corrConnLatent(lv).p = stats.p;
        corrConnLatent(lv).sides = sides;

        fprintf('Variable %d of %d done \n',lv, length(latentNames)); toc;

    end

    save([outDir,outName,inFiles(1).name], 'corrConnLatent', 'matchedData'); 

end

%========================================================================== 
% Check results with regstats
%==========================================================================

% tic;
% for e = 1:length(upperInds)
%     
%         [u v] = ind2sub([nNodes nNodes],upperInds(e));
%         y = squeeze(conn(u,v,:));
%         stats2 = regstats(y,design);
%         t(e,:) = stats2.tstat.t(2:end);
%         df(e) = stats2.tstat.dfe;
%         p(e,:) = stats2.tstat.pval(2:end);
%         
% end
% 
% % force NaNs to 0 or 1
% t(isnan(t)) = 0;
% p(isnan(p)) = 1;
% 
% fprintf('regstats done \n'); toc;
% 
% for lv = 1:length(latentNames)
%     
%     tMat = zeros(size(conn,1));
%     pMat = zeros(size(conn,1));
%     
%     tMat(upperInds) = t(:,lv);
%     tMat = tMat+tMat';
%     
%     pMat(upperInds) = p(:,lv);
%     pMat = pMat+pMat';
%     
%     eval(['regMatlab.',latentNames{lv},'.t = tMat;']);
%     eval(['regMatlab.',latentNames{lv},'.p = pMat;']);        % pval is 2-sided & is double the p above.
%     eval(['regMatlab.',latentNames{lv},'.df = df;']);
%     
% end
% 
% 
