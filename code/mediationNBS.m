function nbsMed = mediationNBS(med, param, primaryThr, compThr, nNulls, intensity)

% nbsMed = mediationNBS(med, param, primaryThr, compThr, nNulls, intensity, sides)
%
% This function will take the output of the mediation model and run an NBS
% analysis. It applies the primary threshold by setting a percentile
% cut-off on the null data, with the percentile determined by primaryThr.
%
% -------
% INPUTS:
% -------
% med -     a structured array that is the output of the mediation model 
%           (see runMediationFilter.m)
% 
% param -   string denoting the specific mediation parameter to analyse
%           ('a', 'b', 'ab', etc), as listed in med structure.
%
% primaryThr - the primary component forming threshold. This should be
%              specified as a p-value threshold for the parameter selected 
%              in param. Note that this threshold will be applied by taking
%              the absolute values in the matrix and looking for vals that
%              exceed the critical value at p = .05. If you are looking at
%              a component that includes edges significant in a two-tailed
%              sense (e.g., edges showing positive and negative
%              correlations with some variable), then primaryThr should be
%              havled. i.e., to obtain p = .05, p = .025 should be used.
%
% compThr -     p-val threshold for determining whether an NBS component is
%               significant
%
% nNulls -      Number of permutations forming the empirical null
%               distribution. This corresponds to the number of 'sims' run 
%               in the permutation model
%
% intensity -   1 = estimate NBS effect using intensity; 0 = extent.
%
% -------
% OUTPUTS:
% -------
% nbsMed -  a structured array containing the following fields:
%
%            .compP - pvals for each component
%
%            .compSigAdj - adjacency matrix for each significant component
%
%            .obsAdj - observed adjacency matrix
%
%            .nbsNullMaxDist - null distribution of max component sizes
%            used to estimat component-wide pval
%
%            .nullCompSz - the component sizes obtained for each
%            permutation
%
%            nulls - initial adjacency matrices for each permutation for
%            the chosen parameter 'param', as extracted from the input 'med'.
%
% ------------------
% UNRESOLVED ISSUES:
% ------------------
%
% In the inital dataset I examined, the distribution of observed values is
% substantially shifted to the right relative to the nulls, given by
% mediation bootstrapping, nearly always results in one large component
% that is highly significant (that is, p = 0 because no null gets even
% close). Need to verify if this behaviour is accurate or an artifact.
%
% Alex Fornito, Monash University, Jan 2017
%==========================================================================

% % testing
%  param = 'ab';
%  primaryThr = .025;
%  intensity = 0;
%  compThr = .05;
%  nNulls = 1000;
%  sides = 2;
%  med = medOrig;

%-------------------------------------------------------------------------
% Run analysis
%-------------------------------------------------------------------------

nbsMed = [];

for i = 1:length(med)
   
    % get observed data for chosen parameter
    eval(['obs = med(i).coeffs.',param,';']); 
    nNodes = size(obs,1);                               % number of nodes
    eval(['obsP = med(i).p.',param,';']);               % pvals
        
    % extract null data & compute critical value
    eval(['sims = [med(i).sims.',param,'];']);          % extract data
    nullData = horzcat(sims.data{:});                   % concatenate into matrix
    nullSrt = sort(abs(nullData),'descend');            % sort vals
    statThr = nullSrt(floor(primaryThr*nNulls),:)';     % get critical val
    
    %----------------------------------------------------------------------
    % Get component stats for observed data
    %----------------------------------------------------------------------
    
    % threshold and binarize
    obsMat = zeros(nNodes);                             % initialize
    edgeInds = find(~isnan(obs));                       % get indices of non-nan values
    edgeVals = abs(obs(edgeInds));                      % vectorize and take absolute
    edgeVals(edgeVals<=statThr) = 0;                    % threshold vector
    obsMat(edgeInds) = edgeVals;                        % put back in matrix
    obsAdj = double(logical(obsMat));                   % binarize
    obsAdj = obsAdj+obsAdj';                            % make symmetric

    % get components
    [obsComp obsCompSzNodes] = get_components(obsAdj);
        
    % find components with at least 1 edge (i.e., >1 node)
    indSz=find(obsCompSzNodes>1);
    obsCompSzNodesNontrivial = obsCompSzNodes(indSz);   % retain size of components with >1 node
    
    % initalize outputs
    szLinks=zeros(1,length(indSz));
    maxSz=0; 

    % for each component with >1 node
    for c = 1:length(indSz)

        % find nodes in that component
        nodes=find(indSz(c)==obsComp);

        % measure size of component
        if intensity == 1
            % Measure size as intensity
            % Andrew subtracts primary thresh from stat values. Not sure why. Have not done it here.
            szLinks(c)=sum(sum(logical(obsAdj(nodes,nodes)).*obsMat(nodes,nodes)))/2; 
        else
            %Measure size as extent
            szLinks(c)=sum(sum(logical(obsAdj(nodes,nodes))))/2;
        end
        
        % create output matrix where components are indicated by component
        % number
        nbsAdj(nodes,nodes)=obsAdj(nodes,nodes)*(c);

        % store max component size     
        if maxSz < szLinks(c)
            maxSz = szLinks(c);
        end
        
    end
    
    %----------------------------------------------------------------------
    % Get null distribution of max component sizes
    %----------------------------------------------------------------------
   
    % null = {};

    % initialize
    nullMat = zeros(nNodes);
  
    % get matrix indices for edges analysed
    indsNull = sub2ind([nNodes nNodes], sims.rowInd, sims.colInd); 
    
    % for each permutation
    for s = 1:nNulls
                   
        % generate thresholded adjacency matrix
        nullMat = zeros(nNodes);                                % initalize
        sigNull = abs(nullData(s,:))>statThr';                  % threshold
        nullMat(indsNull) = abs(nullData(s,:)).*sigNull;        % put 'significant' edges back in matrix 
        nullMat = nullMat+nullMat';                             % symmetrize
        nullMat(isnan(nullMat)) = 0;                            % set NaNs to zero (just in case)             
        nullAdj = double(logical(nullMat));                     % binarize  
 
        % for debugging
        % null{s} = nullMat;                                     % store null matrix 
        % nullAdjAll{s}=  nullAdj;                               % for debugging
        
        % get null components
        [nullComp nullCompSz] = get_components(nullAdj);

        % find components with at least 1 edge (i.e., >1 node)
        nullIndSz=find(nullCompSz>1);

        % initialize outputs
        nullSzLinks=zeros(1,length(nullIndSz));
        
        % initialize null max component size
        nullMaxSz = 0;

        % for each component with >1 node
        for c = 1:length(nullIndSz)

            % get nodes belonging to component
            nodes=find(nullIndSz(c)==nullComp);

            if intensity == 1
                % Measure size as intensity
                % Andrew subtracts primary thresh from stats. Not sure why. Have not done it here.
                nullSzLinks(c)=sum(sum(logical(nullAdj(nodes,nodes)).*nullMat(nodes,nodes)))/2; 
            else
                %Measure size as extent
                nullSzLinks(c)=sum(sum(logical(nullAdj(nodes,nodes))))/2;
            end
              
            % store max component size     
            if nullMaxSz < nullSzLinks(c)
                nullMaxSz = nullSzLinks(c);
            end

        end

        % store max size in null distribution
        nullMaxDist(s) = nullMaxSz;

        % store null component details 
        nullCompSzNodes{s} = nullCompSz;
        nullCompSzLinks{s} = nullSzLinks;
        
        % fprintf('Null %d of %d complete \n',s,nNulls);
    end

    %----------------------------------------------------------------------
    % Estimate p-val for each observed component
    %----------------------------------------------------------------------
  
    pComp = [];
    nbsComp = [];

    % for each observed component
    for c = 1:length(szLinks)    
        pComp(c) = sum(nullMaxDist>=szLinks(c))/nNulls;        % get pval
    end

    % find significant components
    sigInds = find(pComp<compThr);

    % if any components are significant
    if ~isempty(sigInds)
        
        nbsComp = zeros(nNodes, nNodes, length(sigInds));
        
        for c = 1:length(sigInds)   
            % define a matrix with edges comprising the component.
            m = zeros(nNodes);
            nodeInds = find(nbsAdj==sigInds(c));
            m(nodeInds) = 1;
            nbsComp(:,:,c) = m;
        end    
    end

    % store outputs
    nbsMed(i).compP = pComp;
    nbsMed(i).compSigAdj = nbsComp;
    nbsMed(i).compObsSzNodes = obsCompSzNodesNontrivial;
    nbsMed(i).compObsSzLinks = szLinks;
    nbsMed(i).compAdj = obsAdj;
    nbsMed(i).nbsNullMaxDist = nullMaxDist;
    nbsMed(i).nullCompSzNodes = nullCompSzNodes;
    nbsMed(i).nullCompSzLinks = nullCompSzLinks;
    %nbsMed(i).nulls = null;
    
end

