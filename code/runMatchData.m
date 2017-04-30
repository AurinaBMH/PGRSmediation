clear all 
close all

% base working directory
baseDir = '/Users/Aurina/GoogleDrive/Genetics_connectome/PGRSmediation/';

% directory where code is
addpath([baseDir,'code/']);

% directory where PRSice outputs are located
pgrsDir = [baseDir,'data/PRSiceOUTPUTfiles/'];

% where the SEM results are
latentFile = [baseDir,'data/Psychopathology-Factor-Scores-with-Bellgrove-IDs-DARIS-IDs.txt'];
latent = readtable(latentFile,'delimiter','\t');    % read in data

% where the covariate data are
covFile = [baseDir,'data/CovariateFile_AS.txt'];
cov = readtable(covFile, 'delimiter','\t');  
cov.Sex = num2cell(cov.Sex); 
% read in data
% column indices of covariates to use in analysis
covInds = [3 4];

% where the ID matching data are
idFile = [baseDir,'data/GenCogIDlist'];
load(idFile)                                        % load data

% tractography algorithm
tractAlgo = {'FACT'}%,'iFOD2'};
%tractAlgo = 'FACT';

%parcellation 
parc = {'aparcANDfirst'}% 'HCPMMP1ANDfslatlas20'};
%parc = 'HCPMMP1ANDfslatlas20';

%were data sifted?
sift = {'SIFT2'};
%sift = 'SIFT2';

%edge weight type
weight = {'standard'}% 'FA'};

%PGRS for which disorder? only choose 1
disorders = {'SCZ'}%,'ADHD','AUT','BIP','MDD'};

%imputation threshold for PGRS (choose 1)
pgrsProc = 'S179PRUNINGpheno';

%Latent variable ID used in PRSice output file names
latentNames = {'Pfact', 'Impuls', 'SCZpos', 'SCZneg','Instab'};

%GWAS imputation threshold
pgrsRThr = '03';

%output directory
outDir = [baseDir,'results/matchedData/'];

%prefix of output file name
outPref = 'matchedData_';

%==========================================================================
% run analysis
%==========================================================================

cnt = 1;
for d = 1:length(disorders)
    
        for a = 1:length(tractAlgo)
    
            for p = 1:length(parc)

                for s = 1:length(sift)

                    for w = 1:length(weight)
                        
                        % initialize output table for PGRS model stats
                         tic;
                         r2Max = {};
                         r2Max{1,1} = 'Phenotype';
                         r2Max{1,2} = 'Best fit model stats';
                         r2Max{1,3} = 'Scores at max';
                         r2Max{1,4} = 'Index at max in r2_out vector (add 1 for table of scores)';
                        
                         for lv = 1:length(latentNames)

                                % store latent variable names
                                pgrsModelFits{1,lv} = latentNames{lv};
                                
                                % read in PGRS model fits
                                pgrsModelFits{2,lv} = readtable([pgrsDir,pgrsProc,'/',disorders{d},pgrsRThr,'/PRSice_',latentNames{lv},'_RAW_RESULTS_DATA.txt'],'delimiter',' ');

                                % find best model fit
                                [maxVal maxInd] = max(pgrsModelFits{2,lv}.r2_out);
    
                                % store results for each latent variable
                                r2Max{lv+1,1} = latentNames{lv};
                                r2Max{lv+1,2} = pgrsModelFits{2,lv}(maxInd,:);

                                % read in PGRS scores and extract those for
                                % best fitting model
                                pgrsScores = readtable([pgrsDir,pgrsProc,'/',disorders{d},pgrsRThr,'/PRSice_',latentNames{lv},'_SCORES_AT_ALL_THRESHOLDS.txt'],'delimiter',' ');
                                r2Max{lv+1,3} = pgrsScores(:,[1 maxInd+1]);
                                r2Max{lv+1,4} = maxInd; 
                                
                                % load appropriate connectivity data
                                connName = [tractAlgo{a},'_',parc{p},'_',sift{s},'_',weight{w}];
                                connFile = [baseDir,'data/connectomes/',tractAlgo{a},'/',connName];
                                load(connFile); 

                                % match datasets for common subjects
                                data = matchData(r2Max{lv+1,3}, latent, conn, cov, IDlist, IDs, covInds);
                                
                                % separate outputs to store best PGRS scores
                                % separately for each latent variable
                                matchedData.missing = data.missing;
                                matchedData.latentNames = data.latentNames;
                                matchedData.latentNames(:,[1 2]) = [];                  % retain only variable names
                                matchedData.covNames = data.covNames;
                                
                                % stroe data for purely behavioural
                                % analysis
%                                 matchedData.behav.latent = data.behav{2,1};
%                                 matchedData.behav.pgrs{1,lv} = latentNames{lv};
%                                 matchedData.behav.pgrs{2,lv} = data.behav{2,2};
%                                 matchedData.behav.cov = data.behav{2,3};
%                                 matchedData.behav.conn = data.behav{2,4};
                                
                                % store data for subjects with MRI,
                                % genetics and behaviour
                                matchedData.mri.latent = data.mri{2,1};
                                matchedData.mri.pgrs{1,lv} = latentNames{lv};
                                matchedData.mri.pgrs{2,lv} = data.mri{2,2};
                                matchedData.mri.cov = data.mri{2,3};
                                matchedData.mri.conn = data.mri{2,4};
                                matchedData.pgrsR2Max = r2Max;
      
                         end
                         
                         % save output
                        save([outDir,outPref,disorders{d},'_',tractAlgo{a},'_',parc{p},'_',sift{s},'_',weight{w}], 'matchedData');
                         fprintf('Dataset %d of %d done \n',cnt,length(disorders)*length(tractAlgo)*length(parc)*length(sift)*length(weight)); toc;
                         cnt = cnt+1;
                         
                    end
                end
            end
        end
end

    
    
        
        
        