function [meanPVC,meanRSA,PVC,RSA] = popCoherence(dfof)
%% popCorr
% Calculates population activity coherence statisitcs.
%
% Population Vector Analysis (PVA): Provides a representation of how the
% cell population as a whole encodes each loacation.
%   10.1523/JNEUROSCI.4896-03.2004
%
% Representational Similarity Analysis (RSA): Similar to PVA
%   10.3389/neuro.06.004.2008
%
% Inputs:
%   dfof - activity matrix (laps x bins x cells). This dimension set is
%   that directly loaded from dfofM and related variables. Is converted
%   to bins x cells x laps.
%


%%

% clear; close all
% cd('D:\4m_optoPaper\TM221012\230115\loc1\pcaica')

% load('RunByRun_sig\dfofMInterpM_sig.mat')
% dfof = cat(3,dfofMInterpM_sig{:});

% load('RunByRun_dfof\dfofMInterpM.mat')
% dfof = cat(3,dfofMInterpM{:});


%% Compute PVA

% permute dfof
dfofPerm = permute(dfof,[2 3 1]);

% % normalize to peak of 1
mx = max(dfofPerm, [], 1, 'omitnan');
dfofNorm = dfofPerm./max(mx,eps);

% z-score normalize
% mu = mean(dfofPerm, 1, 'omitnan'); 
% sd = std(dfofPerm, 0, 1, 'omitnan');
% dfofNorm = (dfofPerm-mu)./max(sd,eps);

% define sizes
nBins = size(dfofNorm,1);
nLaps = size(dfofNorm,3);

% initialize output array
PVC = zeros(nBins,nBins,nLaps,nLaps);
RSA = zeros(nBins,nBins,nLaps,nLaps);

% calculate population vector correlations
for ii = 1:nLaps
    % get lap 1 activity
    dfofNorm1 = dfofNorm(:,:,ii)';
    dfof1 = dfofPerm(:,:,ii)';

    for jj = 1:nLaps
        if ii==jj
            PVC(:,:,ii,jj) = NaN;
            RSA(:,:,ii,jj) = NaN;
            continue
        end
        
        % get lap 2 activity
        dfofNorm2 = dfofNorm(:,:,jj)';
        dfof2 = dfofPerm(:,:,jj)';

        % calculate PVC correlation
        curPVC = corr(dfofNorm1,dfofNorm2,'rows','pairwise');
        PVC(:,:,ii,jj) = curPVC;

        % calculate RSA correlation
        curRSA = 1-corr(dfof1,dfof2,'rows','pairwise');
        RSA(:,:,ii,jj) = curRSA;
    end
end

% calculate mean
meanPVC = mean(PVC,[3 4],'omitnan');
meanRSA = mean(RSA,[3 4],'omitnan');


%%

% figure
% tiledlayout(1,2)
% 
% for kk = 1:2
%     if kk==1
%         plotData = meanPVC;
%     elseif kk==2
%         plotData = meanRSA;
%     end
% 
%     nexttile(kk)
%     imagesc(plotData)
%     axis('square')
%     clim([-0.2 1])
% end
% 


end

