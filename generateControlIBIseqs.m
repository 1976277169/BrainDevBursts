function generateControlIBIseqs(refH, refseed, Hs, itermax, newdir)
% generateControlIBIseqs generate IBIseqs with same values but different order
% Input parameters: 
% - refH: Hurst exponent for the reference IBIseq
% - refseed: Seed for the reference IBIseq
% - Hs: Range of Hurst exponents for which to generate IBIseqs, e.g., [0.5, 0.6, 0.7, 0.8]
% - itermax: Max number of seeds for each Hurst exponent
% - newdir: Directory in which to store files -- user to create before!
% NB: Reference file + all FDN files assumed to be in ./data

    % Load reference IBIseq
    fpath = './data/';
    load(strcat(fpath, 'H_', num2str(refH*100), '_', num2str(refseed)))

    % For each Hs
    for h = 1:length(Hs)
        Hexp = Hs(h);
        
        % for each iter
        for iter = 1:itermax

            % Load noise values for each configuration
            load(strcat(fpath, 'FDN_', num2str(Hexp*100), '_', num2str(iter)));
            xpn = xpn(1:length(IBIseq));

            % Re-order
            Is = sort(IBIseq); % Index from the reference sequence
            [~, Ix] = sort(xpn); % Index from the noise
            [~, IIx] = sort(Ix); % Order from the noise

            IBIseq = Is(IIx); % Re-ordered sequence (same values as original but different order)

            % Check directory newdir exists if not create it
            if ~exist(newdir, 'dir')
                mkdir(newdir);
            end
            
            % save files
            save(strcat(newdir, 'H_', num2str(Hexp*100), '_', num2str(iter)), 'IBIseq');

        end
    end
end