function [PPS] = popSynchronyFull(dfof,winSz)
%% popCorr
% Calculates pairwise population activity synchrony statisitcs.
%
% Pairwise population synchrony is the pairwise mean correlation between
%   every pair of cells as a function of spatial location. 
%
% Inputs:
%   dfof - activity matrix (laps x bins x cells). This dimension set is
%   that directly loaded from dfofM and related variables. Is converted
%   to bins x cells x laps.
%   winSz - number of bins to use for each spatial window
%


%%

% clear; close all
% cd('D:\4m_optoPaper\TM221012\230115\loc1\pcaica')
% 
% % load('RunByRun_sig\dfofMInterpM_sig.mat')
% % dfof = cat(3,dfofMInterpM_sig{:});
% 
% load('RunByRun_dfof\dfofMInterpM.mat')
% dfof = cat(3,dfofMInterpM{:});
% 
% winSz = 5;
% 

%% Compute PVA

% permute dfof
dfofPerm = permute(dfof,[3 2 1]);

% define sizes
nCells = size(dfofPerm,1);
nWins = size(dfofPerm,2)-winSz+1;

% initialize output array
PPS = zeros(nWins,nWins);

% calculate pairwise synchrony
for ii = 1:nWins
    % get activity
    dfofWin1 = dfofPerm(:,ii:ii+winSz-1,:);
    dfofWinCat1 = reshape(dfofWin1,nCells,[])';

    for jj = 1:nWins
        % get activity
        dfofWin2 = dfofPerm(:,jj:jj+winSz-1,:);
        dfofWinCat2 = reshape(dfofWin2,nCells,[])';

        % calculate pairwise correlations
        dfofCorr = corr(dfofWinCat1,dfofWinCat2,'rows','pairwise');
        PPS(ii,jj) = mean(dfofCorr(~eye(size(dfofCorr))),'omitnan');
    end
end


%%

% figure; imagesc(PPS)
% ylim([0 1])

end

