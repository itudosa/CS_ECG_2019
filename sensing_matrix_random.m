% @Author: Francesco Picariello email: fpicariello@unisannio.it
% Updated: 25/02/2019
% MATLAB version R2018a
%
% Phi = sensing_matrix_random(M,N,type), this function returns the sensing matrix
% Phi (MxN) acording to the Bernoulli (+1, -1) distribution and the
% Gaussian distribution.
%
% Input parameters:
% N is the number of samples acquired for each frame
% M is the number of compressed samples
% type is the type of used distribution ('Bernoulli' or 'Gaussian')
%
% Output parameters:
% Phi is the sensing matrix

function Phi = sensing_matrix_random(M,N,type)
    if strcmp(type,'Bernoulli')
        Phi = randi(2,M,N)-1;
    elseif strcmp(type,'Gaussian')
        Phi = randn(M,N);
    end
end