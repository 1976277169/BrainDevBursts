function xpn = generateFDN( H, len, aseed )
% generateFDN returns Fractional Differencing Noise with given Hurst exponent
% This is a wrapping function for FDN. 
% The function will only save the data to file is a seed is provided. 
% Input parameters: 
% H: Hurst exponent
% len: Length of sequence
% aseed (optional): seed for random generator (required for file saving)
% Output: 
% xpn: FDN sequence

    % If seed is provided
    if exist('aseed', 'var')
        rng(aseed, 'twister');
    end

    xpn = FDN(len, H-0.5); % Fractional differencing parameter d = H-0.5. 
    
    % Save only if seed is provided
    if exist('aseed', 'var')
        save(strcat('./data/FDN_', num2str(H*100), '_', num2str(aseed), '.mat'), 'xpn'); % Save time series (to avoid unnecessary repeat calculations)
    end
end
