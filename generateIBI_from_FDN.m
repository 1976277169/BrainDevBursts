function IBIseq = generateIBI_from_FDN(Hexp, mean_dur, SL, aseed)
% generateIBI_from_FDN returns sequence of IBIs from FDN with given Hurst
% This function checks whether the data has already been generated/saved.
% Input parameters: 
% Hexp: Hurst exponent
% mean_dur: Mean duration of the IBIs (integer = number of spikes)
% SL: Length of data, no default provided but used 1*10^6 in paper.  
% aseed: seed for random generators
% Output: 
% IBIseq: IBI sequence

    filename = strcat('./data/FDN_',num2str(Hexp*100),'_',num2str(aseed),'.mat'); 
    if exist(filename,'file') % Data have already been generated
        load(filename, 'xpn'); % This will load 'xpn' which contains FDN noise
    else
        xpn = generateFDN(Hexp, SL/mean_dur, aseed); % generate the FDN. 
        % NB: Given mean_dur, we only generate sequences with SL/mean_dur elements. 
        % Code could be made more robust by checking that sum of all IBIs >= SL
    end

    % Generate exponentially distributed integer waiting times with mean_dur mean duration
    IBIh = ceil(abs(mean_dur*randn(length(xpn),1)));  
    if sum(IBIh)<SL 
        fprintf('generateIBI_from_FDN: Warning: Sequence is short %d vs %d but not expected to be a problem.\n', sum(IBIh), SL); 
    end

    % Principle: We want to order the IBIh so they have the same ordering as the FDN. 
    Is=sort(IBIh); % IBIs sorted in increasing length
    [~,Ix]=sort(xpn); % Ix is index of FDN noise (xpn) sorted by increasing value, so xpn(Ix) would be FDN sorted in increasing length
    [~,IIx]=sort(Ix); % IIx is telling us in which order this sorted FDN was taken to form xpn. NB: Ix(IIx) = 1..length(IBIh). 
    IBIseq=Is(IIx); % Use that order to sequence the sorted IBIh. 
    
    % Save IBI sequence (variable 'IBIseq')
    save(strcat('./data/H_',num2str(round(Hexp*100)),'_',num2str(aseed),'.mat'),'IBIseq')

end