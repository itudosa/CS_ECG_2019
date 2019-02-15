% @Author: Francesco Picariello email: fpicariello@unisannio.it
% Updated: 15/02/2019
% MATLAB version R2018a
%
% Compression and Reconstruction of the signal No. 106 from the MIT-BIH
% arrhythmia database for N = 720, M = 180 => CR = 4

clear all, close all, clc

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% Database reading
file_name = 'mitdb/106'; % file name related to the ECG signal from the MIT-BIH arrhythmia database
[x, Fs, tm] = rdsamp(['/database/' file_name],1); % reading of the ECG data from the database

%% CS processing of the ECG signal
[xest,xcut,alpha,y,Phi,x_th,Psi] = CS_power_based(x,720,180,0.2); % compression and reconstruction of the entire ECG signal

%% Results displaying
% plotting of the last sensing matrix Phi and of the Mexican Hat wavelet
% matrix Psi
figure
subplot(2,1,1)
imagesc(Phi)
title('The last sensing matrix \Phi')
set(gca,'FontSize',16)
subplot(2,1,2)
imagesc(Psi)
title('The Mexican Hat wavelet matrix \Psi')
set(gca,'FontSize',16)
% plotting of the compressed vector y
figure
plot(y,'LineWidth',2)
xlabel('Sample')
ylabel('Amplitude [mV]')
set(gca,'FontSize',16)
grid on
% plotting of the input signal xcut and of the estimated signal xest
figure
t = 0:1/Fs:(length(xcut)-1)*1/Fs;
plot(t,xcut,'LineWidth',2)
hold on
plot(t,xest,'LineWidth',2)
xlabel('Time [s]')
ylabel('Amplitude [mV]')
grid on
set(gca,'FontSize',16)
% evaluation of the Percentage of Root-mean-squared Difference (PRD)
disp('The Percentage of Root-mean-squared Difference:')
PRD = norm(xcut-xest')/norm(xcut)*100

