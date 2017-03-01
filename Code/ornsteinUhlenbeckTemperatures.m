function temperatures = ornsteinUhlenbeckTemperatures(X,dimX) 
% Computes the diffusion coefficients of a TxN time-series up to a scaling
% factor

if dimX == 1
temperatures = squeeze(mean((X(2:end,:)- X(1:end-1,:)).^2,1)'/2);
end

if dimX == 2
    temperatures = zeros(size(X,2));
    for k=1:size(X,2)
        for l=k:size(X,2)
            temperatures(k,l) = mean((X(2:end,k)- X(1:end-1,k)).*(X(2:end,l)- X(1:end-1,l)),1)'/2;
        end
    end
    temperatures = temperatures + temperatures' - diag(diag(temperatures));
end
