This function computes various metrics of functional and effective
connectivity time-locked to certain behavioral/stimulation events, including:
	* Spectral power
	* Phase locking value
 	* Coherence
 	* Imaginary coherence
 	* Phase slope index
 	* Granger causality
 
This code also generates a single plot where you can compare all of the
above-mentioned metrics. 
Input: - xser  (time series from one brain region)
       - yser  (time series from another brain region)
       - fs    (sample rate of time series)
       - event (vector of event times in seconds)
       - twin  (time window for spectral analysis, eg twin = [-2 5])
Output - funcCon (structure containing all information on functional/effective connectivity analyses)
For this script you will need to have the MVGC toolbox in your path. You
can find the website for the toolbox here: http://users.sussex.ac.uk/~lionelb/MVGC/
Or you can also just copy the code from the following folder in my
codebase: C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeIain\mainCode\Code\mvgc_v1.0
I.S. 2017
