function matched = matchData(pgrs, latent, conn, cov, pgrsID, connID, covInds)

% matched = matchData(pgrs, latent, conn, cov, pgrsID, connID, covInds)
% 
% This function will identify the common subjects with SEM, MRI and PGRS
% data.
% 
% 
%
%==========================================================================

% get inds of subjects with DARIS ids
mriInds = find(latent.DARISID>0);
darisID = latent.DARISID(mriInds);

%==========================================================================
% Extract latent and PGRS data for subjects with DARIS IDs
%==========================================================================

cnt1 = 1;
cnt2 = 1;
for i = 1:length(mriInds)
    
    latDat = table2array(latent(mriInds(i),2:end));                     % convert latent variables to array
    idx = find(pgrsID.DARIS_ID==darisID(i));                          % get indices of daris subjects
    
    if ~isempty(idx)                                                    % if subject has both SEM and PGRS data

        latentMRI(cnt1,:) = latDat;                                     % add SEM data for subject
        
        pID = pgrsID.GWAS_ID(idx);                                   % get corresponding PGRS ID
        idx2 = find(pgrs.IID==pID);                                  % find row index in PGRS data
        pgrsMRI(cnt1,:) = table2array(pgrs(idx2,:));                    % extract PGRS data
        
        idx3 = find(cov.IID==pID);                                   % get corresponding covariate data
        covMRI(cnt1,:) = [pID table2array(cov(idx3, covInds))];
        cnt1 = cnt1+1;
        
    else                                                                % if both latent and PGRS data missing, note ID
        
        missingIDs(cnt2) = darisID(i);
        cnt2 = cnt2+1;

    end
    
end

matched.missing = missingIDs;
matched.latentNames = latent.Properties.VariableNames;

for i = 1:length(covInds)
    matched.covNames{i} = cov.Properties.VariableNames{covInds(i)};
end

matched.behav{1,1} = 'latent';
matched.behav{1,2} = 'pgrs';
matched.behav{1,3} = 'covariates';

matched.behav{2,1} = latentMRI;
matched.behav{2,2} = pgrsMRI;
matched.behav{2,3} = covMRI;

%==========================================================================
% Match latent and PGRS data to subjects with connectomes
%==========================================================================
      
% find common subjects with PGRS data based on DARIS IDs
commonInds = find(ismember(connID, latentMRI(:,1)))';
commonIDs = connID(commonInds)';

% find corresponding subjects in PGRS results
pgrsRowInds = find(ismember(latentMRI(:,1), commonIDs));

% extract common subjects from PGRS, latent and connectivty
% datasets
pgrsCommon = pgrsMRI(pgrsRowInds,:);
latentCommon = latentMRI(pgrsRowInds, :);
covCommon = covMRI(pgrsRowInds,:);
connCommon = conn(:,:,commonInds); 

matched.mri{1,1} = 'latent';
matched.mri{1,2} = 'pgrs';
matched.mri{1,3} = 'covariates';
matched.mri{1,4} = 'connectomes';

matched.mri{2,1} = latentCommon;
matched.mri{2,2} = pgrsCommon;
matched.mri{2,3} = covCommon;
matched.mri{2,4} = connCommon;
end



