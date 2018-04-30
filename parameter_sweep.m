function parameter_sweep()
% parameter_sweep performs a sweep over the Hurst exponents
% Input parameter: 
% shuffled: optional parameter -- 0 for 
    
    Hs = [0.5, 0.6, 0.7, 0.8];  % Hurst exponents
    Thrs = 1.4:0.2:3; % Losslikelihood thresholds (gain/loss thresholds)
    SL = 1*10^6; % Length of time series
    N = 200; % Size of network
    nb_repeats = 10; % Number of repeats (used as seed)

    for hind = 1:length(Hs)
        H = Hs(hind);
        for thrsind = 1:length(Thrs)
            thr = Thrs(thrsind); 
            for iter=1:nb_repeats % Number of repeats
                runSim(H, SL, N, thr, iter, 0); % First run the LRTC case (which needs to be run before shuffled case)
                runSim(H, SL, N, thr, iter, 1); % Needed to get the required data     
            end
        end
    end    
end