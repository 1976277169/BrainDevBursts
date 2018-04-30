function runSim(H, SL, N, thr, Ap, Ad, iter, shuffled)
% runSim runs one simulation for a given set of parameters
% Input parameters: 
% - Hexp: Hurst exponent
% - SL: Length of record
% - N: Size of network
% - thr: Losslikelihood threshold 
% - iter: Seed for random generators
% - shuffled: Whether IBI sequence has LRTC or is to be shuffled
% This code require functions from BCT toolbox (not included here)
% https://sites.google.com/site/bctnet/ 
% Version used: 2017_01_15_BCT

    %=====================================================================
    % Initialisation
    %=====================================================================

    % Check directory ./data exists if not create it
    if ~exist('data', 'dir')
        mkdir('./data');
    end
        
    % Add path to progressbar routine
    addpath('./progressbar/');
    
    % Parameters shared across all simulations
    % Could be loaded from file... 
    params.H = H; % Hurst exponent
    params.SL = SL; % Length of record
    params.N = N; % Size of network
    params.iter = iter; % Seed
    params.burst_amp = 0.8; % Amplitude of the bursts
    params.burst_len = 5; % Duration (in number of spikes) of the bursts
    params.p = 0.4; % Probability of connection in the network
    params.mean_dur = 5; % Target mean duration for IBIs
    params.Vinit = -60; % Initial voltage
    params.Rm = -70; % Resting potential (V_r in paper)
    params.reset = -60; % Potential after spike
    params.fire_thres = -54; %Firing threshold (V_{thres} in paper)
    params.tau_stdp = 5; % time constant for STDP
    params.tau_L = 100; % decay of likelihood (\tau_L in paper)
    params.gl_mu = 0.025; % mean conductance
    params.gl_std = 0.005; % std conductance
    params.gain_thres = thr; % Likelihood threshold over which connection is gained
    params.loss_thres = -thr; % Likelihood threshold under which connection is lost
    params.Ap = Ap; % Long-Term Potentiation parameter
    params.Ad = Ad; % Long-Term Depression parameter
    params.ccfind = 100; % Period with which to collect network information (= number of iterations)
    params.win = 200; % Number of samples over which to calculate input rate
    
    %=====================================================================
    % LRTC input + Random network
    %=====================================================================

    % About random number generation: either the data is being loaded
    % (e.g., shuffled case OR when the data have already been generated and
    % can be loaded from file -- code checks automatically) or the data
    % need generating and the random number generators will be seeded by
    % the generateFDN. This seeding will apply to both the generation of
    % the LRTC input and the adjacency matrix of the starting network. 
    % The random number generator will then be reset to the seed when
    % initialising the neuronal simulations. This will guarantee the same
    % sequence of events for both LRTC and shuffled cases. 
        
    if shuffled == 0 % most data need to be generated
        % Load (or generate) IBIs
        filename = strcat('./data/H_', num2str(round(H*100)), '_', num2str(iter), '.mat');
        if exist(filename, 'file') 
            load(filename, 'IBIseq'); % This will load variable 'IBIseq'
        else % File has not been saved already
            IBIseq = generateIBI_from_FDN(H, params.mean_dur, SL, iter); 
        end

        % Input
        [h, IBIseq] = calc_external_input(SL, IBIseq, params.burst_amp, params.burst_len, shuffled); % Generate LRTC input and its IBI sequence
        
        % Network structure
        CIJ = createRandomNetwork(N, params.p*(N*(N-1))); % Create random network with N nodes and p*(N*(N-1)) links -- do not allows self-loops so p=1 means all but diagonal links exist 
        CIJstart = CIJ; % Store initial matrix (which will be used in the shuffled case)
    else
        % Load newIEI and CIJ from file
        filename = strcat('./data/sims_N', num2str(N), '_H', num2str(H*100), '_Ap', num2str(params.Ap), '_Ad', num2str(params.Ad), '_thr', num2str(params.gain_thres), '_s', num2str(iter),'.mat');
        if ~exist(filename, 'file')
            fprintf('runSim: Could not find %s. Exiting...', filename); 
            return
        end
        load(filename, 'IBIseq', 'CIJstart'); 
                
        % Need to generate input corresponding to this IBI sequence
        rng(iter, 'twister'); % Reset generator
        shindIBI = randperm(length(IBIseq)); % Random permutation of the IEI
        IBIseq = IBIseq(shindIBI);
        h = calc_external_input(SL, IBIseq, params.burst_amp, params.burst_len, shuffled); % No need to return IBIseq since it is unchanged
        CIJ = CIJstart; % Load initial matrix from LRTC case
    end

    
    %=====================================================================
    % Neuronal simulation initialisation
    %=====================================================================

    % Reset random generators to seed
    rng(iter, 'twister'); % Reset generator
    
    % Initialisation of the neuronal simulations' variables + storage
    v = params.Vinit*ones(N, 1); % Initial values of v
    gl = abs(params.gl_mu+params.gl_std*randn(N, 1)); % Conductances (make sure they are positive)
    %firings = []; % storage for spike timings -- see code for lines to comment out if interested in firing times
    x = zeros(N, 1); % decay variables for each from spike
    losslikelihood = zeros(N); % Variable L in paper
    tlast = zeros(N, 1); %time of last spike

    % Prepare storage for network parameters and collect stats on initial network        
    if SL >= params.ccfind % Only needed if SL > 100 (which should always be the case)
        clust = zeros(SL/params.ccfind, 1); % storage for clustering coefficient
        mpl = zeros(SL/params.ccfind, 1); % storage for mean path length
        nocount = zeros(SL/params.ccfind+1, 1); % storage for proportion of connections (with respect to N*(N-1)) -- referred to as ratio in paper
        modularity = zeros(SL/params.ccfind, 1); % storage for modularity (via Louvain algorithm)
        noconn_gained = []; % list containing number of connections gained over the last ccfind iterations
        noconn_lost = []; % list containing number of connections lost over the last ccfind iterations
        nocomponents = zeros(SL/params.ccfind+1, 1); % storage for number of connected components
        count = 1; % Index (in units of SL/params.ccfind)
        nocount(count) = params.p; % initial proportion of connections
        cc = clustering_coef_bd(CIJ); % local clustering coefficient for binary / directed networks
        clust(count) = mean(cc); % initial global clustering coefficient 
        [~,D] = reachdist(CIJ);
        [lambda,~,~,~,~] = charpath(D);
        mpl(count) = lambda; % initial mean path length
        [~,Q] = community_louvain(CIJ); 
        modularity(count) = Q; % initial modularity
        [~,comp_sizes] = get_components(CIJ);
        nocomponents(1) = length(comp_sizes); % initial number of connected components
    end
    %firingrate = zeros(SL, 1); % storage for number of spikes (not actually a rate) -- see code for lines to comment out if interested in firing rate
        
    %=====================================================================
    % Simulation main loop
    %=====================================================================

    progressbar(0);
    for t = 1:SL % loop over time
        progressbar(t/SL);

        % Collect statistics about the network at this point
        if mod(t,params.ccfind) == 0 % collect info every ccfind
            % disp(t/SL*100) % unless progressbar not available
            count = count+1; % current index (in units of SL/params.ccfind)
            cc = clustering_coef_bd(CIJ);
            clust(count) = mean(cc); % current global clustering coefficient
            [~,D] = reachdist(CIJ);
            [lambda,~,~,~,~] = charpath(D);
            mpl(count) = lambda; % current mean path length
            nocount(count) = sum(sum(CIJ))/(N*(N-1)); % current proportion of connections
            [~,Q] = community_louvain(CIJ);
            modularity(count) = Q; % current modularity
            [~,comp_sizes] = get_components(CIJ);
            nocomponents(count) = length(comp_sizes); % current number of connected components
        end

        % Update x and losslikelihood (decays)
        xupdate = exp(-(t-tlast)/params.tau_stdp).*x;
        losslikelihood = losslikelihood*exp(-1/params.tau_L);

        % Current input
        input=h(t);

        % Find neurons that have fired
        fired = find(v > params.fire_thres);
        % firingrate(t) = length(fired); % a count rather than a rate
        % firings =[firings; t+0*fired, fired]; 
        
        if ~isempty(fired) % There were spikes
            v(fired) = params.reset; % Reset firing neuron to reset voltage

            cur_p = sum(sum(CIJ))/(N*(N-1)); % current proportion of connections
            weight = (2/(cur_p*N));
            input = input+(sum(weight*CIJ(fired,:), 1))'; % synaptic input

            tlast(fired) = t; % update last firing time
            x(fired) = xupdate(fired)+1; % update x variable
            losslikelihood(:, fired) = losslikelihood(:, fired)+params.Ap*repmat(xupdate, 1, length(fired)); % LTP if post>pre
            losslikelihood(fired, :) = losslikelihood(fired, :)-params.Ad*repmat(xupdate', length(fired), 1); % LTD if pre>post

            [radding,cadding] = find((losslikelihood(:, fired)-repmat(xupdate, 1, length(fired))) < params.gain_thres & losslikelihood(:, fired) >= params.gain_thres); % Candidates for link addition
            [rlosing,closing] = find((losslikelihood(fired, :)+repmat(xupdate', length(fired), 1)) > params.loss_thres & losslikelihood(fired, :) <= params.loss_thres); % Candidates for link removal

            if ~isempty(radding) % If there are candidates for addition
                nb_added = 0; % counter of how many links get added
                for cidx = 1:length(radding) % cycle through the list
                    if (radding(cidx) ~= fired(cadding(cidx))) && (CIJ(radding(cidx), fired(cadding(cidx))) == 0) % link is not self-loop and link exists
                        CIJ(radding(cidx), fired(cadding(cidx))) = 1; % Add link to adjacency matrix
                        nb_added = nb_added + 1;
                    end
                end
                noconn_gained = [noconn_gained, nb_added]; % Keep track of number of links added
            else
                noconn_gained = [noconn_gained, 0]; % No links added
            end

            if ~isempty(rlosing) % If there are candidates for removal
                nb_removed = 0; % counter of how many links get removed
                for cidx = 1:length(rlosing) % cycle through the list
                    if CIJ(fired(rlosing(cidx)),closing(cidx)) == 1 % link exists (no need to check if it's self-loop since they shouldn't be created)
                        CIJ(fired(rlosing(cidx)),closing(cidx)) = 0; % Remove link from adjacency matrix
                        nb_removed = nb_removed + 1;
                    end
                end
                noconn_lost = [noconn_lost, nb_removed]; % Keep track of number of links removed
            else
                noconn_lost = [noconn_lost,0]; % No links removed
            end
        else % No firings this time so no change to adjacency matrix
            noconn_gained = [noconn_gained,0]; 
            noconn_lost = [noconn_lost,0];
        end

        % Voltage update using small steps (a la Izhikevich)
        v = v+0.25*(-gl.*(v-params.Rm)+input);
        v = v+0.25*(-gl.*(v-params.Rm)+input);
        v = v+0.25*(-gl.*(v-params.Rm)+input);
        v = v+0.25*(-gl.*(v-params.Rm)+input);
    end
    progressbar(1);

    %=====================================================================
    % Plotting
    %=====================================================================

    figure; plot(nocount); title('Number of connections')
    
    figure; plot(modularity); title('Modularity')
    % [M,Q]=community_louvain(CIJ);
    % [On,CIJr] = reorder_mod(CIJ,M);
    %  figure; imagesc(CIJr);
    %figure(); plot(firings(:,1),firings(:,2),'.')

    figure; plot(nocomponents); title('Number of components')
    
    figure; 
    p1=subplot(3,1,1); plot(noconn_gained); ylabel('No. connections gained'); ylim([0,50])
    p2=subplot(3,1,2); plot(noconn_lost); ylabel('No. connections lost'); ylim([0,150])

    % Calculate rate of input 
    rate_input = NaN*ones(length(h), 1); 
    for i=params.win:length(h)
        rate_input(i)=sum(h(i-params.win+1:i))/params.win;
    end
    p3=subplot(3,1,3); plot(rate_input); ylabel('Rate of input'); ylim([0.2,0.8])
    
    linkaxes([p1,p2,p3],'x')

   
    %=====================================================================
    % Saving data
    %=====================================================================

    if shuffled == 0
        filename = strcat('./data/sims_N', num2str(N), '_H', num2str(H*100), '_Ap', num2str(params.Ap), '_Ad', num2str(params.Ad), '_thr', num2str(params.gain_thres), '_s', num2str(iter),'.mat');

        % Put simulation results into structure 'results'
        results.nocount = nocount;
        results.clust = clust;
        results.mpl = mpl;
        results.CIJ = CIJ;
        results.noconn_gained = noconn_gained;
        results.noconn_lost = noconn_lost;
        results.rate_input = rate_input;
        results.nocomponents = nocomponents;
        results.modularity = modularity; 
        % Save all relevant data
        save(filename,'SL','N','IBIseq','CIJstart','params','results');
    else
        % Save variables under different names for comparison purposes
        results_shuff.nocount = nocount; 
        results_shuff.clust = clust;
        results_shuff.mpl = mpl;
        results_shuff.CIJ = CIJ;
        results_shuff.noconn_gained = noconn_gained;
        results_shuff.noconn_lost = noconn_lost; 
        results_shuff.rate_input = rate_input;
        results_shuff.nocomponents = nocomponents; 
        IBIseq_shuff = IBIseq; 
        filename = strcat('./data/sims_shuff_N', num2str(N), '_H', num2str(H*100), '_Ap', num2str(params.Ap), '_Ad', num2str(params.Ad), '_thr', num2str(params.gain_thres), '_s', num2str(iter),'.mat');
        save(filename, 'IBIseq_shuff', 'results_shuff');
    end
    
end

function [h, IBIseq] = calc_external_input(len, IBIseq, burst_amp, burst_len, shuffled) 
% calc_external_input calculates the input provided a sequence of IBIs
% Input parameters
% len: max length of record
% burst_amp: amplitude of input
% burst_len: burst len (= number of consecutive spikes)
% Output: 
% h: the external input
% IBIseq: the exact sequence of IBIs used to construct the input. This is 
%   important so that in the shuffled case, we are guaranteed to use the very
%   same IBIs (albeit in a different order)

    h = burst_amp*ones(len,1); % h initialised to constant input of amplitude burst_amp
    pos = burst_len+1; % index of first IBI sample, i.e., there is a burst at the start
    IBIidx = 1; % IBI index
    while pos <= len % This could result in an error if sum(IBIseq)+burst_len*len(IBIseq)<SL -- currently unchecked
        h(pos:pos+IBIseq(IBIidx)-1) = 0; % set h to 0 during the IBI 
        pos = pos+IBIseq(IBIidx)+burst_len; % index of next IBI sample
        IBIidx = IBIidx+1; 
    end

    if shuffled == 0 % In the shuffled case, we have the exact sequence to use but in th LRTC case, we may have too many IBIs
        IBIseq = IBIseq(1:IBIidx-1); % We exited loop because last IEI + burst took us past SL. 
        % NB: Improvement would be to differentiate between the IBI or the burst taking us past SL
        % In the case of it being the burst, we should admit the IBI.
        % Effect likely negligible. 
    end
end


