clear all
close all

% This script will run an NBS analysis on the outputs of mediation
% modelling. It requires as input a series of critical values for the primary
% component-forming threshold, as obtained through
% getMediationCritValsEmpirical.m. 

%==========================================================================
% Settings
%==========================================================================

% add path to code
addpath /Users/alexfornito/Alex_Docs/Grants-Awards/Applications/2017/CIA/NHMRC-PGRS-psych-disorders/Analysis/code/

% path to BCT toolbox requires for NBS analysis
addpath(genpath('/Users/alexfornito/Alex_Docs/Matlab_functions/BCT/'))

% base directory
baseDir = '/Users/alexfornito/Alex_Docs/Grants-Awards/Applications/2017/CIA/NHMRC-PGRS-psych-disorders/Analysis/';

% directory where mediation modelling results are located
inDir = [baseDir,'results/mediation/PGRSoptimal/Prune03/'];

% list of mat files to analyse
inFiles = dir([inDir,'*bca0*']);

% primary component-forming threshold for NBS
primaryThr = [.05 .01];

% directory to write results
outDir = [baseDir,'results/NBS/PGRSoptimal/Prune03/'];

% prefix of output file name
outPref = 'mediationNbsExtent_';

% mediation parameters to analyse
params = {'a','b','c','c_prime','ab'};

% component-wide threshold for determining significance
compThr = .05;

% number of nulls run in mediation analysis
nNulls = 1000;

% intensity or extent NBS? 0 = extent; 1 = intensity.
intensity = 0;

% 1-sided or 2-sided hypothesis test? 1 = one-sided; 2 = two-sided. This
% will determine how the primary threshold is applied in the NBS.
% BAsically, if you are looking to threshold results in one direction
% (e.g., if edges in matrix represent correlations with some variable, 
% threshold such that all correlations are > r), then choose 1-sided. If you
% want to look at effects in both directions (e.g., threshold to retain
% strong positive and strong negative correlations), then choose 2-sided.
sides = 2;

%==========================================================================
% Run analysis
%==========================================================================

if sides == 2
   primaryThr = primaryThr./2;
end
   
for f = 1:length(inFiles)
    
    med = [];
    
    % load mediation results
    load([inDir,inFiles(f).name], 'med');

    if ~isempty(med)

        nbsMed = [];

        % for each latent variable
        for v = 1:length(med)

            tic;
            % format output and store basic parameters
            nbsMed(v).intensity = intensity;                            % analysis options
            nbsMed(v).sides = sides;

            % for each mediation parameter
            for par = 1:length(params)

                % for each primary component-forming threshold
                for t = 1:length(primaryThr)

                    % run NBS
                    nbsResults = mediationNBS(med(v), params{par}, primaryThr(t), compThr, nNulls, intensity);

                    % store results
                    results.primaryThr(t) = primaryThr(t);  % primary pval threhsold
                    results.nbs(t) = nbsResults;                  % NBS results

                end
                
                eval(['nbsMed(v).results.',params{par},'=results;']);

            end

            fprintf('Variable %d of %d for dataset %d of %d complete \n',...
                       v, length(med),f, length(inFiles)); toc;

        end

        % save output 
        outName = [outDir,outPref,inFiles(f).name(1:end-4)];
        save(outName,'nbsMed');
    end
  
end
