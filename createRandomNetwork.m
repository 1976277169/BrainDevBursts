function CIJ=createRandomNetwork( N, kc)
% createRandomNetwork returns a random (directed) network, no self-loop
% Input parameters: 
% N: number of neurons
% kc: total number of connections in the network

    CIJ=zeros(N); % empty adjacency matrix
    ind=~eye(N); % matrix of 1 for all non-autapse
    i = find(ind); % index of all non-autapse
    rp = randperm(length(i)); % random permutation of indices
    irp = i(rp); % reordered sequence 
    CIJ(irp(1:kc)) = 1; % set links for first kc element of reordered sequence
end

