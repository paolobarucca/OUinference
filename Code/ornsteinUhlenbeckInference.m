function [couplings temperatures covariance indeces]= ornsteinUhlenbeckInference(X)
% P. Barucca 11.10.16 
% Solves the Lyapunov Equation [1] and finds the couplings of an heterogenous Ornstein-Uhlenbeck process 
% X = is the TxN time-series (warning: for long time-series it is crucial not to insert the transpose NxT time-series, 
% that would result in computing a TxT matrix)
% [1] Localization in covariance matrices of coupled heterogenous
% Ornstein-Uhlenbeck processes - http://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.062129

temperatures = mean((X(2:end,:)- X(1:end-1,:)).^2,1)'/2;

indeces = temperatures>0;

yTemperatures = temperatures(indeces);
Y = X(:,indeces);
yCovariance = cov(Y);
yCouplings = ornsteinUhlenbeckInverseMethod(-yCovariance, yTemperatures);

covariance = zeros(size(X,2));
covariance(indeces,indeces) = yCovariance;

couplings = zeros(size(X,2));
couplings(indeces,indeces) = yCouplings;
