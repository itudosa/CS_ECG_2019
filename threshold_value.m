% @Author: Francesco Picariello email: fpicariello@unisannio.it
% Updated: 15/02/2019
% MATLAB version R2018a
%
% x_th = threshold_value(x), this function returns the threshold value on
% each frame according to the amplitude value of the x signal that exhibits
% the highest probability (i.e. it corresponds to the noise floor).
%
% Input parameters:
% x the frame of the acquired samples
%
% Output parameters:
% x_th is the estimated threshold value

function x_th = threshold_value(x)
    [N_hist,X_hist] = hist(abs(x-mean(x)),sqrt(length(x))); % histogram of the input signal x
    [~,I] = max(N_hist); % index of the class which exhibits the highest probability
    x_th = X_hist(I); % identification of the threshold value
end