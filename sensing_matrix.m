% @Author: Francesco Picariello email: fpicariello@unisannio.it
% Updated: 15/02/2019
% MATLAB version R2018a
%
% Phi = sensing_matrix(x,x_th,M), this function returns the sensing matrix
% Phi according to the threshold value x_th as reported in the paper titled
% "A Novel Method for Compressed Sensing based Sampling of ECG Signals in
% Medical-IoT era".
%
% Input parameters:
% x is the ii-th frame of the input signal
% x_th is the threshold value
% M is the number of compressed samples
%
% Output parameters:
% Phi is the sensing matrix defined according to the threshold value x_th

function Phi = sensing_matrix(x,x_th,M)
    I = find(abs(x-mean(x))>=x_th); % identification of the indices of the samples of |x-x_avg| that overcomes the threshold values x_th
    Phi(1,:) = zeros(1,length(x)); % inizialization of the first row of the matrix Phi
    Phi(1,I) = 1; % positioning of the ones in the first row of the sensing matrix according to the threshold value
    for ii = 2:M
        Phi(ii,:) = circshift(Phi(ii-1,:),round(length(x)/M)); % definition of the other rows of the sensing matrix
    end
end