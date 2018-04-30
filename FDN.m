function X = FDN( N, d )
% FDN Returns fractional differencing noise with parameter d
% It is calculated using Hosking 1981's method (equation 3.1),
% actually the so-called direct method. 
% Input parameters: 
% N: number of samples (integer)
% d = H-0.5 where H is the Hurst exponent (as measured by DFA)

    X = zeros(N, 1); % the FDN
    a = randn(N, 1); % source of gaussian noise 
    psicoeff = FDNcoeff(N, d)';

    progressbar(0);
    for n = 2:N
        progressbar((n/N));
        % Progress update
        % if mod(n, 10000)==0
        %    fprintf('FDN: percentage complete: %3d\n', (n/N)*100)
        % end
        X(n) = a(n) + psicoeff(1:n-1) * X(n-1:-1:1); 
    end
    progressbar(1);
end

function psicoeff = FDNcoeff( n, d )
% psicoeff as defined in the so-called direct method which really is a 
% rearrangement of Equation 3.1 from Hosking 1981.  

    psicoeff = ones(n, 1);
    psicoeff(1) = d;
    for k=2:n
        psicoeff(k) = (k-1-d)/k*psicoeff(k-1);
    end
end
