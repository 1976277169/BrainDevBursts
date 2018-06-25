function cluster_sweep(jobid)
% cluster_sweep wrapper function for use with the cluster
% Input parameter: 
% jobid: integer from 1 to N where N = length(Hs)*length(Thrs)*length(taus)*nb_repeats 

% NB: In producing paper, this script was duplicated as follows. 
% In the first, lines 42 and 43 were commented out so as to produce time
% series as fast as possible.
% In the second, lines 40 and 41 were commented out and network statistics
% were calculated over the pre-calculated files. 

    % Variables over which a range was explored
    Hs = [0.5,0.6,0.7,0.8];  % Hurst exponents
    Thrs = 2:0.2:3; % Losslikelihood thresholds (gain/loss thresholds) -- set range
    tauLs = 100:10:150; % L decay -- set range
    
    % Variables with fix value (could be a range if needed)
    tau_STDP = 10; % STDP time constants
    nb_repeats = 20; % Number of repeats (used as seed)
    SL = 5*10^6; % Length of time series
    N = 200; % Size of network
    Ap = 0.5; % LTP
    Ad = 0.55; % LTD
   
    % Code below assumes three variables with range: Hs, Thrs, tauLs
    hind = floor((jobid-1)/(length(Thrs)*length(tauLs)*nb_repeats)); % index of Hs (-1)
    H = Hs(hind+1);
    
    jobid = mod(jobid-1, length(Thrs)*length(tauLs)*nb_repeats); % Take remainder of jobid
    thrsind = floor(jobid/(length(tauLs)*nb_repeats)); % index of Thrs (-1)
    thr = Thrs(thrsind+1);

    jobid = mod(jobid-1, length(tauLs)*nb_repeats); % Take remainder of jobid
    tauLsind = floor(jobid/nb_repeats); % index of tauLs (-1)
    tau_L = tauLs(tauLsind+1);
    
    iter = mod(jobid, nb_repeats) + 1; % repeat index
   
    %fprintf('%d %d %d\n', hind+1, thrsind+1, iter);
    runSim(H, SL, N, thr, Ap, Ad, tau_STDP, tau_L, iter, 0); % First run the LRTC case (which needs to be run before shuffled case)
    runSim(H, SL, N, thr, Ap, Ad, tau_STDP, tau_L, iter, 1); % Needed to get the required data 
    no_rand = 50; % number of null networks to generate at each time point
    calculateNormalisedNetworkStatistics(H, N, Ap, Ad, thr, tau_STDP, tau_L, iter, no_rand); % Calculate normalised network statistics
end
