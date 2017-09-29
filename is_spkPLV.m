function [sponSpkPLV,sacSpkPLV] = is_spkPLV(spkCell,sacTime,C,lfpFs,f)
% This function computes spike phase locking across the entire recording,
% as well as the time-resolved spike phase locking around saccades. 
% input:  spikeCell - spike times in seconds stored in a cell. Dimensions(1,numChans)
%         sacSamp   - saccades in samples
%         C         - The result of the convolution fo an lfp matrix with a
%                     complex Morlet wavelet. Dimensions(channels,samples)
%         lfpFs     - LFP sample frequency
%         f         - index of frequency
%
% output: sponSpkPLV - the estimated spike phase locking value
%                      computed using a predefined number of randomly 
%                      drawn spikes over various repetitions. Dimensions(numChans,1)
%         sacSpkPLV  - The time-resolved spike-PLV computed around the
%                      timing of saccades. Dimensions(channels,time bins) 
% I.S 2016

display(['processing spike PLV ' num2str(f)])

tres     = 2;   % time resolution (plus/minus)
winSize  = 0.5; % size of window to pool spikes
halfWin  = winSize/2;
stepSize = 0.01; % sliding window increment size
nSpks    = 100;  % number of spikes to compute spike PLV in each time bin
nReps    = 200;  % number of times to recompute spike PLV on randomly drawn spikes
sponSpks = 1000; % number of spikes to consider when computing spontaneous PLV
binC     = -tres:stepSize:tres; % vector of the center of spike PLV bins
numBins  = numel(binC);

numChans   = length(spkCell);
sponSpkPLV = nan(numChans,1);
sacSpkPLV  = nan(numChans,numBins);
for ichan = 1:numChans; 
    % display(['Spike PLV chan ' num2str(ichan)])
    % skip channel if there aren't enough spikes 
    if numel(spkCell{ichan}) < sponSpks; continue; end
    spks     = spkCell{ichan};
    spkSamp  = round(spks*lfpFs); % convert spike times to LFP phase
    spkPhase = angle(C(ichan,spkSamp)); % get LFP phase at time of spike
    
    % randomly draw spikes to compute spike PLV and take the mean
    tmpPLV = nan(nReps,1);
    for irep = 1:nReps
        rp = randperm(numel(spkPhase));
        tmpPLV(irep) = abs(nanmean(exp(1i*(spkPhase(rp(1:sponSpks)))))); % compute PLV on randomly drawn spikes
    end
    sponSpkPLV(ichan) = mean(tmpPLV); % mean PLV of all randomly drawn spikes
    
    spkMat   = nan(numel(sacTime),numBins); % Initialise matrix to count spikes per saccade and time bin
    sacPhase = cell(numel(sacTime),numBins); % initialise cell to fill with spike phases
    for isac = 1:numel(sacTime)
        for ibin = 1:numBins
            tmpWin  = [sacTime(isac)+binC(ibin)-halfWin sacTime(isac)+binC(ibin)+halfWin]; % time window for saccade and time bin
            binSpks = find(spks>tmpWin(1) & spks<tmpWin(2)); % get spikes that fit in timw window of interest
            if isempty(binSpks); spkMat(isac,ibin) = 0; else spkMat(isac,ibin) = numel(binSpks); end % count spikes in time window
            sacPhase{isac,ibin} = spkPhase(binSpks); % fill cell with spike phases
        end
    end
    
    nspksPerBin = sum(spkMat); % the number of spikes per time bin across all saccades
    if min(nspksPerBin) < nSpks;
        continue % if there aren't enough spikes to compute PLV then skip channel
    else
       
        for ibin = 1:numBins
            % display(num2str(ibin))
            % collect all spike phases for the time bin
            tmpPhase = [];
            for isac = 1:size(spkMat,1)
                tmpPhase = horzcat(tmpPhase,sacPhase{isac,ibin});
            end
            % loop through nReps times and compute PLV on nSpks number of
            % randomly drawn spikes
            tmpPLV = nan(1,nReps);
            for irep = 1:nReps
                rp           = randperm(numel(tmpPhase));
                tmpPLV(irep) = abs(nanmean(exp(1i*(tmpPhase(rp(1:nSpks))))));
            end
            sacSpkPLV(ichan,ibin) = mean(tmpPLV);
        end
    end
end
    

