function couplings = ornsteinUhlenbeckInverseMethod(covariance, temperatures)
% Solves the Lyapunov Equation [1] and finds the couplings of an heterogenous Ornstein-Uhlenbeck process 
% covariance   = covariance matrix of the zero-mean signals 
% temperatures = diffusion coefficient of each of the signals, see ornsteinUhlenbeckTemperatures()

% [1] Localization in covariance matrices of coupled heterogenous
% Ornstein-Uhlenbeck processes -http://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.062129

if size(temperatures,2)==1
    T = diag(temperatures);
else if size(temperatures,1)==size(temperatures,2)
        T = temperatures;
    end
end

nVars = length(T);

%Spectral method 
[eVec eVal] = eig(-covariance);

eMat = zeros(nVars);

% for k=1:nVars
%     for l=k:nVars      
%         eMat(k,l) = sum(eVec(:,k).*T.*eVec(:,l));       
%     end
% end
% eMat = eMat + eMat'-diag(diag(eMat));

eMat = eVec'*T*eVec;
eVal = ones(nVars)*eVal;
eVal = (eVal+eVal').^-1;
couplings = -eVec*(2*eVal.*eMat)*eVec';
