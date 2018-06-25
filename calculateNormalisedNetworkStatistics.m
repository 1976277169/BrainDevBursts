
function calculateNormalisedNetworkStatistics(H, N, Ap, Ad, thr, tau_STDP, tau_L, iter, no_rand)
% generates null netw and calc normalised clustering and mean path length
% Input parameters:
% H: Hurst exponent
% N: Network size
% Ap: Long-Term Potentiation parameter
% Ad: Long-Term Depression parameter
% thr: Losslikelihood threshold 
% tau_stdp: STDP time constant
% tau_L: L decay
% iter: Seed
% no_rand: Number of null networks to generate for each case

    % File containing IBI sequences + key network metrics in LRTC case
    filename1 = strcat('./data/sims_N', num2str(N), '_H', num2str(H*100), '_Ap', num2str(Ap), '_Ad', num2str(Ad), '_thr', num2str(thr), '_tauSTDP', num2str(tau_STDP), '_tauL', num2str(tau_L), '_s', num2str(iter),'.mat');
    
    if ~exist(filename1,'file')
        fprintf('compare_to_random: File %s could not be found\n', filename1);
        return
    else
        load(filename1); % This will load all content in workspace. Variables of interest are IBIseq and results.nocount, results.clust, results.mpl
        nocount = results.nocount; 
        clust = results.clust; 
        mpl = results.mpl; 
    end
    
    % File containing key network metrics in shuffled case
    filename2 = strcat('./data/sims_shuff_N', num2str(N), '_H', num2str(H*100), '_Ap', num2str(params.Ap), '_Ad', num2str(params.Ad), '_thr', num2str(params.gain_thres), '_tauSTDP', num2str(tau_STDP), '_tauL', num2str(tau_L), '_s', num2str(iter),'.mat');

    if ~exist(filename2,'file')
        fprintf('compare_to_random: File %s could not be found\n', filename2);
        return
    else
        load(filename2); % This will load all content in workspace. Variables of interest are results_shuff.nocount, results_shuff.clust, results_shuff.mpl
        nocount_shuff = results_shuff.nocount; 
        nocount_shuff = results_shuff.nocount; 
        clust_shuff = results_shuff.clust; 
        mpl_shuff = results_shuff.mpl; 

    end
        
    % Storage
    clustrand_store = zeros(length(nocount), 1); % global clustering for null networks
    mpl_store = zeros(length(nocount), 1); % mean path length
    clustrand_shuff_store = zeros(length(nocount), 1); % global clustering
    mpl_shuff_store = zeros(length(nocount), 1); % mean path length

    
    for i = 1:length(nocount) % for each of the recorded time points
        disp(i/length(nocount)*100); % progress update
        
        % Storage (not really needed, could calculate mean incrementally)
        clustrand = zeros(no_rand, 1);
        mplrand = zeros(no_rand,1);
        clustrand_shuff = zeros(no_rand, 1);
        mplrand_shuff = zeros(no_rand, 1);

        for j = 1:no_rand % for each of the null models

            % Generate a random network with proportion of links identical 
            % to network with LRTC input
            % and calculate global clustering coefficient and mean path
            % length. This follows the same process as used in runSim
            CIJ = createRandomNetwork(N, floor(nocount(i)*(N*(N-1))));
            cc = clustering_coef_bd(CIJ);
            clustrand(j) = mean(cc);
            [~,D] = reachdist(CIJ);
            [lambda, ~, ~, ~, ~] = charpath(D);
            mplrand(j) = lambda;

            % As above but for the network with shuffled input
            CIJ_shuff = createRandomNetwork(N, floor(nocount_shuff(i)*(N*(N-1))));
            cc = clustering_coef_bd(CIJ_shuff);
            clustrand_shuff(j) = mean(cc);
            [~,D] = reachdist(CIJ_shuff);
            [lambda, ~, ~, ~, ~] = charpath(D);
            mplrand_shuff(j) = lambda;

        end

        % Store averages
        clustrand_store(i) = mean(clustrand);
        mplrand_store(i) = mean(mplrand);
        clustrand_shuff_store(i) = mean(clustrand_shuff);
        mplrand_shuff_store(i) = mean(mplrand_shuff);

    end

    % Calculate normalised clustering and mean path length
    norm_clust = clust./clustrand_store;
    norm_mpl = mpl./mplrand_store';

    norm_clust_shuff = clust_shuff./clustrand_shuff_store;
    norm_mpl_shuff = mpl_shuff./mplrand_shuff_store';

    % Save to file
    filename3 = strcat('./data/sims_norm_clust_N', num2str(N), '_H', num2str(H*100), '_Ap', num2str(Ap), '_Ad', num2str(Ad), '_thr', num2str(thr), '_tauSTDP', num2str(tau_STDP), '_tauL', num2str(tau_L), '_s', num2str(iter),'.mat');
    save(filename3, 'norm_clust', 'norm_mpl', 'norm_clust_shuff', 'norm_mpl_shuff');  
end