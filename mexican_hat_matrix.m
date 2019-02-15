% @Author: Francesco Picariello email: fpicariello@unisannio.it
% Updated: 15/02/2019
% MATLAB version R2018a
%
% Psi = mexican_hat_matrix(N) is the function that returns the NxN Mexican
% Hat wavelet matrix (Psi) defined according to the paper titled "ECG analysis
% using the Mexican-Hat wavelet"
% 
% Input parameter:
% N number of samples for frame
%
% Output parameter:
% Psi, the NxN Mexican_Hat wavelet matrix

function Psi = mexican_hat_matrix(N)
    n = 0:N-1;
    rr = 1;
    for mm = 1:log2(N)
        a = 2^mm; % definition of the Mexican Hat scaling factor
        for b = 0:a:N-1 % definition of the delay factor
            Psi(:,rr) = (1/sqrt(a)) * (2/sqrt(3)) * pi^(-1/4) * (1 - ((n-b)/a).^2 ) .* exp( -((n-b)/a).^2/2 ); % definition of the Mexican Hat vector psi(a,b)
            rr = rr + 1;
        end
    end
    Psi(:,N) = ones(N,1); % last column vector of the matrix Psi which consists of ones
end