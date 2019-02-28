% @Author: Francesco Picariello email: fpicariello@unisannio.it
% Updated: 25/02/2019
% MATLAB version R2018a
%
% Psi = dictionary_matrix(N,type)is the function that generated the NxK
% matrix Psi according to the choesn dictionary (e.g. Mexican Hat, Symlet4,
% etc.)
%
% Input parameters:
% N is the number of samples to be compressed for each ECG data frame
% type is the type of dictionary, Mexican Hat or Symlet 4 or others
%
% Output parameters:
% Psi is the matrix that define the domain where the signal is sparse

function Psi = dictionary_matrix(N,type)
    if(strcmp(type, 'Mexican Hat'))
    % Mexican Hat defined according to the paper "Compressed Sensing for
    % bioelectric signals: a review, Mexican Hat with two intervals: (i) N
    % shifted with {-5,5}, and (ii) N shifted with {-1,1}.
        Psi(:,1) = mexihat(-5,5,N)';
        for ii = 2:N
            Psi(:,ii) = circshift(Psi(:,ii-1),1);
        end
        Psi(:,N+1) = mexihat(-1,1,N)';
        for ii = N+2:2*N
            Psi(:,ii) = circshift(Psi(:,ii-1),1);
        end
    elseif(strcmp(type, 'Mexican Hat Scaled'))
    % Hat wavelet matrix (Psi) defined according to the paper titled "ECG analysis
    % using the Mexican-Hat wavelet"
        Psi = mexican_hat_matrix(N);
    else
    % Definition of the matrix Psi according to the wavelet kernel
    % available in MATLAB
        dwtmode('per'); % dwt with periodized extension
        IN = eye(N);
        for ii = 1:N
            [Psib(:,ii),l] = wavedec(IN(ii,:),fix(log2(N)),type);
        end
        Psi = Psib';
    end
end