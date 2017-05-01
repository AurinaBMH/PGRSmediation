
% this script will run mediation analysis using the output of
% runRegStats.m

%==========================================================================
% Inputs
%==========================================================================

pgrsRThr = '03';
disorder = 'SCZ'; 
tract = 'FACT';
parcel = 'aparcANDfirst'; 
sift = 'SIFT2'; 
weights = 'standard';
tTest = '1sided'; 

% significance threshold for mediation model
medThr = .05;
% pval threshold - only edges with p<pThr in original regression will be 
% analysed with mediation
pThr = .05;
% run non-parametric BCA pval estimation? 0 for no, 1 for yes, 2 for both
% (parameteric and non-parametric)
runBca = 0;
% number of permutations for mediation analysis
nPerms = 1000;
%==========================================================================

% path to bravo toolbox
addpath(genpath('/Users/Aurina/GoogleDrive/Genetics_connectome/PGRSmediation/BRAVO2-master/'));

% path to Alex's helper functions
addpath(genpath('/Users/alexfornito/Alex_Docs/Matlab_functions/preproc_and_basic/'));

% base working directory
baseDir = '/Users/Aurina/GoogleDrive/Genetics_connectome/PGRSmediation/';

% where the input files are
matchDir = [baseDir,'results/matchedData/'];

% get list of input .mat files
fileName = sprintf('*%s*%s*%s*%s*%s*', disorder, tract, parcel, sift, weights); 
matchFiles = dir([matchDir,fileName]);

% directory where conn v latent regression results are for filtering edges
regDir = [baseDir,'results/regression/psych-conn/'];
regFile = dir([regDir,sprintf('*%s*', tTest)]);

% Output directory
outDir = [baseDir,sprintf('results/mediation/PGRSoptimal/Prune%s/',pgrsRThr)];

% prefix for output name
outPref = sprintf('FiltP05_P%sABOnly_%s_perm%d_bca0_', parcel, tTest, nPerms);

% names of latent variables to use in output
latentNames = {'pFact', 'imp', 'sczPos','sczNeg','aff'};

% names of mediation model parameters to consider
paramNames = {'a','b','c','c_prime','ab'};



%==========================================================================
% run analysis
%==========================================================================

% load regression results
load([regDir,regFile(1).name],'corrConnLatent');

for f = 1:length(matchFiles)
    
    load([matchDir,matchFiles(f).name]);
    conn = matchedData.mri.conn;
    conn(isnan(conn)) = 0;                              % set NaNs to zero
    nNodes = size(conn, 1);
    C = matchedData.mri.cov(:,2:end) ;  % covariates
   % make C into double - was an error concatonating files
    C1 = C(:,1); C1 = cell2mat(C1);
    C2 = C(:,2); C2 = cell2mat(C2); C2 = str2num(C2);
    C3(:,1)=C1; C3(:,2)=C2;
    C=C3;
    
    for lv = 1:length(corrConnLatent)
      
        tic;
        % get indicies of edges to analyse
        pVals = corrConnLatent(lv).p;
        %pVals(linspace(1,numel(pVals),length(pVals))) = 1; 
        pVals = zerodiag(pVals,1);                       % set diagonal to 1
       
        % set all pvals in lower triangle to 1, so that only upper triangle
        % is indexed
        lowInds = find(tril(ones(size(pVals))));         % lower triangle inds
        pVals(lowInds) = 1;                              % force lower triangle to 1 
        inds = find(pVals<pThr);                         % find pvals in upper triangle that pass threshold

        % extract PGRS, latent and covariate data for mediation model
        X = matchedData.mri.pgrs{2, lv}(:,2);                 % IV
        Y = matchedData.mri.latent(:,lv+1);                   % DV
        W = [];                                               % moderators
        
        % store variable names
        med(lv).latentName = latentNames{lv};
        
        % initalize output fields
        for n = 1:length(paramNames) 
            eval(['med(lv).p.',paramNames{n},' = nan(nNodes);']);
            eval(['med(lv).coeffs.',paramNames{n},' = nan(nNodes);']);
            eval(['med(lv).sims.',paramNames{n},'.rowInd = nan(length(inds),1);']);
            eval(['med(lv).sims.',paramNames{n},'.colInd = nan(length(inds),1);']);
            eval(['med(lv).sims.',paramNames{n},'.data = cell(length(inds),1);']);
        end

        cnt=1;
        for idx = 1:length(inds)

            % convert linear index of edge to subscript
            [indRow indCol] = ind2sub(size(pVals),inds(idx));
            
            % extract connectivity weights for edge as mediator
            % variable
            M = [];
            M = squeeze(conn(indRow,indCol,:));
            
 
           
            % run mediation model
            [permCoeffs, permSims] = permutation_mediation(X,Y,M,W,C,'n_iter',nPerms);
            
            for n = 1:length(paramNames)
                
                coeffs = []; sims = [];
                eval(['coeffs = permCoeffs.',paramNames{n},';']);
                eval(['sims = permSims.',paramNames{n},';']);
                
                if runBca == 0 
                    [ci p] = ci_percentile_BF(coeffs, sims, medThr);
                elseif runBca == 1
                    [ci p] = ci_bca(coeffs, sims, medThr);
                end

                eval(['med(lv).p.',paramNames{n},'(indRow,indCol) = p;']);
                eval(['med(lv).coeffs.',paramNames{n},'(indRow,indCol) = coeffs;']);
                eval(['med(lv).sims.',paramNames{n},'.rowInd(cnt,1) = indRow;']);
                eval(['med(lv).sims.',paramNames{n},'.colInd(cnt,1) = indCol;']);
                eval(['med(lv).sims.',paramNames{n},'.data{cnt} = sims;']);
                
            end          
            cnt = cnt+1;
        end
        
        fprintf('Variable %d of %d for Dataset %d of %d done \n', lv, length(corrConnLatent), f, length(matchFiles)); toc;
        dataFiles = {matchFiles(f).name; regFile(1).name};
        outName = [outPref,matchFiles(f).name(1:end-4)];
        save([outDir,outName], 'med', 'dataFiles', 'runBca');
    end
    
end




               