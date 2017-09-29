function funcCon = is_functionalConnectivity(xser,yser,fs,event,twin)
% This function computes an various metrics of functional and effective
% connectivity time-locked to certain events, including:
% 	* Spectral power
% 	* Phase locking value
% 	* Coherence
% 	* Imaginary coherence
% 	* Phase slope index
% 	* Granger causality
% 
% This code also generates a single plot where you can compare all of the
% above-mentioned metrics. 
% Input: - xser  (time series from one brain region)
%        - yser  (time series from another brain region)
%        - fs    (sample rate of time series)
%        - event (vector of event times in seconds)
%        - twin  (time window for spectral analysis, eg twin = [-2 5])
% Output - funcCon (structure containing all information on functional/effective connectivity analyses)
%
% For this script you will need to have the MVGC toolbox in your path. You
% can find the website for the toolbox here: http://users.sussex.ac.uk/~lionelb/MVGC/
% Or you can also just copy the code from the following folder in my
% codebase: C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeIain\mainCode\Code\mvgc_v1.0
% I.S. 2017

% initialize data structure
funcCon = struct;

% Define frequencies of interest. Linear spacing for Phase slope index, and
% logarithmic spacing for all other methods.
numFreqs = 100;
lowFreq  = 0.5;
highFreq = 128;
% foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing

% reject noise in recording (high amplitude noise can destroy coherence estimates in particular)
xsig = sub_rejectNoise(xser,fs,2,1); 
ysig = sub_rejectNoise(yser,fs,2,1); 

% Subtract event related potential to look at only induced oscillations
subtractERP = 1;
if subtractERP == 1
    % 
end

% Compute wavelets based on frequencies of interest
morWav = sub_makeWavelet(foi,fs);

% Convolve both signals with wavelets
xspec = nan(numFreqs,numel(xsig));
yspec = nan(numFreqs,numel(ysig));
for f = 1:numFreqs
    fprintf('convolving signals with wavelet %i/%i \n',f,numFreqs)
    xspec(f,:) = conv(xsig,morWav{f},'same');
    yspec(f,:) = conv(ysig,morWav{f},'same');
end

% Cut out spectral data around events
tsamps = round(twin*fs);
xmat = nan(numel(event),numFreqs,diff(tsamps)+1);
ymat = nan(numel(event),numFreqs,diff(tsamps)+1);
for iev = 1:numel(event)
    evSamp = round(event(iev)*fs);
    try % this may not work if the window looks for samples out of range
        xmat(iev,:,:) = xspec(:,evSamp+tsamps(1):evSamp+tsamps(2));
        ymat(iev,:,:) = yspec(:,evSamp+tsamps(1):evSamp+tsamps(2));
    catch
        % do nothing
    end
end
% clear spectral data from memory and compute event-triggered power spectrograms
clear xspec yspec
funcCon.xspec = squeeze(nanmean(abs(xmat).^2,1));
funcCon.yspec = squeeze(nanmean(abs(ymat).^2,1));
funcCon.tvec  = (tsamps(1):tsamps(2))/fs;

% % % % % % % % % % % % % % % %
% Compute phase locking value %
% % % % % % % % % % % % % % % %
fprintf('computing phase locking value \n')
phaseDiff   = angle(xmat) - angle(ymat); % compute phase lag between signals
plv         = squeeze(abs(nanmean(exp(1i*phaseDiff),1))); % PLV formula
funcCon.plv = plv; clear phaseDiff

% % % % % % % % % % %
% Compute coherence %
% % % % % % % % % % %
fprintf('computing coherence \n')
Sxy = squeeze(nanmean(xmat.*conj(ymat),1)); % multiply x by complex conjugate of y
Sxx = squeeze(nanmean(xmat.*conj(xmat),1)); % multiply x by complex conjugate of x
Syy = squeeze(nanmean(ymat.*conj(ymat),1)); % multiply y by complex conjugate of y
Cy  = Sxy./(sqrt(Sxx.*Syy)); % coherency formula
funcCon.coherency          = Cy; % might want to keep the complex part to look at phase lags 
funcCon.coherence          = abs(Cy);
funcCon.imaginaryCoherence = imag(Cy);

% calculate a z-transform of the imaginary part of coherence. Check Nolte
% et al 2004 for more details on this. 
Cxy           = Cy;
numTrials     = numel(event);
imagVariance  = ((1-abs(Cxy).^2)/(2*numTrials)).*((atanh(abs(Cxy)).^2)./abs(Cxy).^2);
imagSTD       = sqrt(imagVariance);
imagZ         = imag(Cxy)./imagSTD;
funcCon.imagZ = imagZ;

% % % % % % % % % % % % % % %
% Compute phase slope index %
% % % % % % % % % % % % % % %
fprintf('computing phase slope index \n')

% For more information on this method, please read Guido's paper:
% https://www.ncbi.nlm.nih.gov/pubmed/18643502
nFreqs = 5; % number of frequencies to compute PSI across
for f = 1:numel(foi) - nFreqs
    fb       = (1:nFreqs) + (f - 1); % frequency band to examine PSI 
    psi(f,:) = nansum(imag(conj(Cy(fb(1:end-1),:)).*Cy(fb(2:end),:))); % PSI formula
    fvec(f)  = mean(foi(fb)); % mean frequency of band of interest
end

% jacknife resampling to estimate the significance of PSI:
% subtract one trial n times and remeasure PSI. the standard deviation of
% PSI equals the standard deviation of the n-1 measured PSI values multiplied
% by the square root of the number of trials.
psiSTD   = nan(numTrials,size(psi,1),size(psi,2));
for n = 1:numTrials
    rt     = randperm(numTrials); % randomise trial indices
    nTrial = rt(1); % randomly select a trial to omit
    % subtract cross spectra from one trial (normalised by the number of trials)
    txy    = Sxy - squeeze(xmat(nTrial,:,:).*conj(ymat(nTrial,:,:)))/numTrials;
    txx    = Sxx - squeeze(xmat(nTrial,:,:).*conj(xmat(nTrial,:,:)))/numTrials;
    tyy    = Syy - squeeze(ymat(nTrial,:,:).*conj(ymat(nTrial,:,:)))/numTrials;
    % recompute coherency
    cy     = txy./sqrt(txx.*tyy);
    % recompute phase slope index
    for f = 1:numel(foi) - nFreqs
        fb            = (1:nFreqs) + (f - 1); % frequency band to examine PSI 
        psiSTD(n,f,:) = nansum(imag(conj(cy(fb(1:end-1),:)).*cy(fb(2:end),:))); % PSI formula
    end
end

% The standard deviation of psi equals the standard deviation of jacknifed
% psi values normalised by the square root of the number of trials:
ts = squeeze(nanstd(psiSTD,[],1)*sqrt(numTrials));

% save raw and std-normalised PSI values
funcCon.psi     = psi;
funcCon.psiNorm = psi./ts; % significant values are ± 2
funcCon.psiFreq = fvec;

funcCon = sub_grangerCausality(funcCon,xser,yser,isnan(xsig),isnan(ysig),event,fs,twin,foi);

doPlot = 1;
if doPlot == 1
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 screensize(3)-100 screensize(4)-150]);
    
    % Compute ticks for plotting
    fois = [0.5 1 2 4 8 16 32 64 128];
    for fi = 1:numel(fois)
        [bi,bb] = sort(abs(foi-fois(fi)));
        tickLoc(fi) = bb(1);
    end
    tickLabel = {'0.5','1','2','4','8','16','32','64','128'};
    
    % do the same for phase slope index frequencies
    fois = [0.6 1 2 4 8 16 32 64];
    for fi = 1:numel(fois)
        [bi,bb] = sort(abs(funcCon.psiFreq-fois(fi)));
        psitickLoc(fi) = bb(1);
    end
    psitickLabel = {'0.6','1','2','4','8','16','32','64'};
    
    % plot power spectrum for signal x
    subplot(2,4,1)
    imagesc(funcCon.tvec,1:numel(foi),pow2db(funcCon.xspec));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)'); % title('Signal X power')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    cl = colorbar('northoutside'); ylabel(cl,'Power (dB) signal X','FontSize',15)
    % plot power spectrum for signal y
    subplot(2,4,2)
    imagesc(funcCon.tvec,1:numel(foi),pow2db(funcCon.yspec));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Signal Y power')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    cl = colorbar('northoutside'); ylabel(cl,'Power (dB) signal Y','FontSize',15)
    % plot phase locking value
    subplot(2,4,3)
    imagesc(funcCon.tvec,1:numel(foi),funcCon.plv);
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('PLV')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',15)
    % plot coherence
    subplot(2,4,4)
    imagesc(funcCon.tvec,1:numel(foi),funcCon.coherence);
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Coherence')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',15)
    % plot imaginary coherence
    subplot(2,4,5)
    imagesc(funcCon.tvec,1:numel(foi),abs(funcCon.imagZ));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Imaginary coherence')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    cl = colorbar('northoutside'); ylabel(cl,'Imaginary coherence (z-score)','FontSize',15)
    % plot phase slope index
    subplot(2,4,6)
    imagesc(funcCon.tvec,1:numel(funcCon.psiFreq),funcCon.psiNorm);
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Phase slope index')
    set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel)
    cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z-score)','FontSize',15)
    %caxis([-5 5])
    % plot granger causality X to Y
    subplot(2,4,7)
    imagesc(funcCon.grangerCausality.tvec,1:numel(foi),funcCon.grangerCausality.X_to_Y);
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
    caxis([0 0.3]); ylim([1 90])
    % plot granger causality Y to X
    subplot(2,4,8)
    imagesc(funcCon.grangerCausality.tvec,1:numel(foi),funcCon.grangerCausality.Y_to_X);
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: Y to X')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: Y to X','FontSize',15)
    caxis([0 0.3]); ylim([1 90])
    colormap(jet)
end

return

function wav = sub_makeWavelet(foi,Fs)
% This function generates complex morlet wavelets that are Gaussian shaped
% in both the time and frequency domains. The equtions used to generate
% wavelets are taken from Tallon-Baundry et al (1996) J Neuroscience.
% 
% Inputs:  foi - vector of center frequencies to generate wavelets
%          Fs  - sample frequency of signal that is to be convolved with wavelets
% Outputs: wav - a cell array that contains the complex morlet wavelets 
% I.S. 2016 

q  = 7; % Wavelet width in cycles
w  = 3; % Width of wavelet taper in terms of standard deviation

wav = cell(numel(foi),1);
for f = 1:numel(foi)
    sf     = foi(f)/q;                   % standard deviation in the frequency domain
    st     = 1/(2*pi*sf);                % standard deviation in the time domain
    t      = -w*st:1/Fs:w*st;            % time vector for wavelet calculation
    A      = 1/sqrt(st*sqrt(pi));        % normalisation factor
    tap    = (A*exp(-t.^2/(2*st^2)));    % Gaussian window
    wav{f} = tap.*exp(2*1i*pi*foi(f)*t); % wavelet formula
    E      = sum(abs(wav{f}).^2);        % energy of wavelet
    wav{f} = wav{f}./sqrt(E);            % normalise wavelet energy to 1
end
return

function sigOut = sub_rejectNoise(sigIn,fs,winLen,rejThresh)
sigOut      = sigIn;
% delWin      = ones(1,round(fs*winLen));                   % window to cut out
% delInd      = (abs(zscore(sigIn)) > rejThresh);           % Samples that surpass threshold
% delVec      = (conv(double(delInd),delWin,'same') > 0);   % Convolve to smooth window and detect non-zero samples
% sigOut(delVec) = nan;                                     % Noisy samples change to NaN's

rejThresh = 200; % 500uV threshold for rejection
[b,a] = butter(4,40/(fs/2),'high'); % define highpass filter at 40Hz
hfSig = filtfilt(b,a,sigIn); % filtyer signal
delWin      = ones(1,round(fs*winLen));                   % window to cut out
delInd      = (abs((hfSig)) > rejThresh);           % Samples that surpass threshold
delVec      = (conv(double(delInd),delWin,'same') > 0);   % Convolve to smooth window and detect non-zero samples
sigOut(delVec) = nan;

return

function funcCon = sub_grangerCausality(funcCon,xsig,ysig,xnan,ynan,event,fs,twin,foi)

% set up priors
regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation
acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)
fres      = [];     % frequency resolution (empty for automatic calculation)
segLength = 1;      % window length for computation of Granger Causality
newFs     = 200; % This is very important for GC! 

% define low pass filter at 100Hz
nyqFreq = fs/2;
[b,a]   = butter(2,100/nyqFreq,'low');
xfilt   = filtfilt(b,a,xsig);
yfilt   = filtfilt(b,a,ysig);

% resample data to sample rate of newFs (defined above)
idat = resample(xfilt,newFs,fs);
jdat = resample(yfilt,newFs,fs);
inan = round(resample(double(xnan),newFs,fs));
jnan = round(resample(double(ynan),newFs,fs));

stepSize  = 0.1; % sliding window increments
stepCen   = twin(1):stepSize:twin(2); % all window centers
recLength = numel(idat)/newFs;
halfWinSamp = (segLength*newFs)/2;

X2Y = nan(numel(foi),numel(stepCen));
Y2X = nan(numel(foi),numel(stepCen));
for istep = 1:numel(stepCen)
    c = 0; % counting variable
    clear imat jmat
    % fill up matrices with data 
    for iev = 1:numel(event)
        % skip if we have window range issues
        if event(iev) < abs(twin(1))+segLength; continue; end
        if event(iev) > recLength-twin(2)-segLength; continue; end
        samp      = round((stepCen(istep)+event(iev))*newFs);
        
        itmp = inan(samp-halfWinSamp:samp+halfWinSamp);
        jtmp = jnan(samp-halfWinSamp:samp+halfWinSamp);
        if sum(itmp) + sum(jtmp) == 0 % only use data that have no noise (ie, no nan's)
            c = c + 1;
            imat(:,c) = idat(samp-halfWinSamp:samp+halfWinSamp);
            jmat(:,c) = jdat(samp-halfWinSamp:samp+halfWinSamp);
        else
            continue
        end
    end
    clear X
    X(1,:,:)  = imat;
    X(2,:,:)  = jmat;
    numSeg    = c;
    
    % compute information criterion
    nvars = size(X,1);
    [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
    amo = 10; % actual model order
    
    % Select model order.
    if strcmpi(morder,'actual')
        morder = amo;
        fprintf('\nusing actual model order = %d\n',morder);
    elseif strcmpi(morder,'AIC')
        morder = moAIC;
        fprintf('\nusing AIC best model order = %d\n',morder);
    elseif strcmpi(morder,'BIC')
        morder = moBIC;
        fprintf('\nusing BIC best model order = %d\n',morder);
    else
        fprintf('\nusing specified model order = %d\n',morder);
    end
    
    % Fit autoregressive model to data
    [A,SIG] = tsdata_to_var(X,morder,regmode);
    assert(~isbad(A),'VAR estimation failed');
    
    % return autocovariance sequence
    [G,info] = var_to_autocov(A,SIG,acmaxlags);
    
    var_info(info,true); % report results (and bail out on error)
    
    % compute Granger Causality based on autocovariance
    f = autocov_to_spwcgc(G,fres);
    assert(~isbad(f,false),'spectral GC calculation failed');
    freqRes = size(f,3)-1;
    freqs   = linspace(0,newFs/2,freqRes+1)'; % compute frequencies
    
    % interpolate to frequencies of interest
    X2Y(:,istep) = interp1(freqs,squeeze(f(2,1,:)),foi,'spline');
    Y2X(:,istep) = interp1(freqs,squeeze(f(1,2,:)),foi,'spline');    
end
funcCon.grangerCausality.X_to_Y = X2Y;
funcCon.grangerCausality.Y_to_X = Y2X;
funcCon.grangerCausality.tvec   = stepCen;
return


% Let's generate some surrogate data for testing this code
numEvs = 100; % define number of events
fs     = 1e3; % sample rate
fband  = [30 60]; % frequency band of surrogate interaction
evDur  = 2; % surrogate stimulus duration
evs    = (1:numEvs)*10+rand(1,numEvs); % Define event times
reclength = evs(end)+evDur+5; % length of vector to make (in seconds)
recSamps  = round(reclength*fs); % number of samples in surrogate data vector
[c,d]     = butter(2,0.5/(fs/2),'high'); % highpass filter to remove low freq components
x         = filter(c,d,pinknoise(recSamps)); % highpass filter pink noise (1/f)
y         = filter(c,d,pinknoise(recSamps)); % highpass filter pink noise (1/f)
[b,a]     = butter(2,fband/(fs/2),'bandpass'); % bandpass filter for adding band limited signals
s         = filter(b,a,randn(1,recSamps))*2; % surrogate band limited signal
timeLag   = -round(0.005*fs); % time lag between two surrogate signals (in seconds)
randSamps = round((rand(1,numEvs)*(reclength-10))*fs); % samples to take surrogate oscillatory signal

% Loop through events and add the band-limited surrogate data
for iev = 1:numEvs
    samp  = round(evs(iev)*fs);
    rsamp = randSamps(iev);
    x(samp:samp+(evDur*fs)) = x(samp:samp+(evDur*fs)) + s(rsamp:rsamp+(evDur*fs)); % add band limited data
    shiftSamp = rsamp + timeLag;
    y(samp:samp+(evDur*fs)) = y(samp:samp+(evDur*fs)) + s(shiftSamp:shiftSamp+(evDur*fs)) + rand(1,evDur*fs+1); % add band limited data with offset plus some noise
end

funcCon = is_functionalConnectivity(x,y,fs,evs,[-2 4])

