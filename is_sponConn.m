function [PLV,AmpCorr] = is_sponConn(C)
% compute phase locking value and amplitude correlation for all
% channel pairs
% input :   C           - the result of the convolution of an LFP matrix with
%                         complex Morlet wavelets. Dimensions(channels,samples);
%
% ouput :   sponPLV     - Phase locking value computed across all samples.
%                         Dimensions(channels,channels)
%           sponAmpCorr - Linear correlation coefficient of the analytic
%                         amplitude of signal pairs. Dimensions(channels,channels)
% I.S 2016

numChans    = size(C,1);
PLV     = nan(numChans);
AmpCorr = nan(numChans);
for ichan = 1:(numChans-1)
    for jchan = (ichan+1):numChans
        ang                  = angle(C(ichan,:)) - angle(C(jchan,:)); % phase difference time series
        PLV(ichan,jchan)     = squeeze(abs(nanmean(exp(1i*(ang))))); % PLV formula
        [r,~]                = corrcoef(abs(C(ichan,:)),abs(C(jchan,:)),'rows','complete'); % linear correlation coefficient of analytic amplitudes
        AmpCorr(ichan,jchan) = r(1,2);
    end
end

