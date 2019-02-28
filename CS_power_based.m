% @Author: Francesco Picariello email: fpicariello@unisannio.it
% Updated: 28/02/2019
% MATLAB version R2018a
%
% [xest,xcut,s,y,Phi,x_th,Psi] = CS_power_based(x,N,M,res,type_algorithm) 
% is the fuction that implements the CS power-based algorithm reported in 
% the paper: "A Novel Method for Compressed Sensing based Sampling of ECG 
% Signals in Medical-IoT era".
%
% Input parameters:
% x is the input signal to be compressed
% N is the number of samples to be compressed for each ECG data frame
% M is the number of compressed samples, the compression ratio (CR) is N/M
% res is the expected residual of the OMP algorithm
% type_algorithm is the type of algorithm that wil be used for estimating
% the alpha coefficients according to the Phi and Psi matrices (OMP-default
% or cvx)
%
% Output parameters:
% xest is the estimate of the input signal x according to the compressed
% xcut is the input vector x cut according to an integer number of frames
% of N samples
% samples y, the sensing matrix Phi, and the orthonormal matrix Psi that
% defines the domain where the input signal will be reconstructed
% alpha is the vector of the reconstructed coefficients in the domain
% defined by Psi
% Phi is the sensing matrix
% x_th vector of threshold values estimated for each frame of samples
% Psi is the matrix that define the domain where the signal is sparse

function [xest,xcut,alpha,y,Phi,x_th,Psi] = CS_power_based(x,N,M,res,type_algorithm)

    xcut = x(1:(floor(length(x)/N)*N)); % cutting of the input vector accordiing to an entire number of frames according to N
    xmat = zeros(floor(length(xcut)/N),N); % initialization of the matrix that will contain on each row a frame of the input signal
% definition of the input matrix ccording to the number of samples for
% frame
    for ii = 1:floor(length(xcut)/N)
        xmat(ii,:) = xcut(1+(ii-1)*N:ii*N); 
    end
% definition of the Mexican Hat matrix
    Psi = mexican_hat_matrix(N);

% CS processing of all the frames
    ymat = zeros(floor(length(xcut)/N),M); % initialization of the matrix containing the compressed M samples, y
    x_th = zeros(1,floor(length(xcut)/N)); % initialization of the threshold vector
    
    f = waitbar(0,'Please wait...');
    for ii = 1:floor(length(xcut)/N)
        x_th(ii) = threshold_value(xmat(ii,:)); % definition of the threshold value used for genereting the sampling vector p
        Phi = sensing_matrix_power_based(xmat(ii,:),x_th(ii),M); % definition of the MxN sensing matrix Phi according to the threshold value x_th
        ymat(ii,:) = Phi*xmat(ii,:)'; % compression of the ii-th frame of x
        if (strcmp(type_algorithm,'cvx'))
            % cvx estimation of the dictionary coefficients
            s = size(Psi); % size of the matrix Psi
            yd = ymat(ii,:)';
            cvx_begin
            variable alpha_est(s(2),1);
            minimize (norm(alpha_est,1)); % function to be minimized according to alpha_est    
            subject to
            yd == Phi*Psi*alpha_est % minimization constraint
            cvx_end
            alpha(ii,:) = alpha_est;
        else
            alpha(ii,:) = OMP(Phi*Psi,ymat(ii,:)',{res},[],[]); % OMP estimation of the dictionary coefficients
        end
        xest((ii-1)*N+1:ii*N) = Psi*alpha(ii,:)'; % reconstruction of the signal in the discrete-time domain according to the estimated Mexican Hat coefficients
        y((ii-1)*M+1:ii*M) = ymat(ii,:); % definition of the compressed vector y;
        waitbar(ii/floor(length(x)/N),f,'Processing ECG data');
    end
    delete(f)
end