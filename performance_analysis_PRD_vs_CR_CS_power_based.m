% @Author: Francesco Picariello email: fpicariello@unisannio.it
% Updated: 28/02/2019
% MATLAB version R2018a
%
% Performance analysis in terms of PRD and PRD per frame (PRDf) for the
% compression and the reconstruction of the signal No. 100 from the
% MIT-BIH arrhythmia database for N = 720, CR = 2,3,4,5,6,7,8,9,10 by 
% adopting the CS algorithm presented in the paper: "A Novel Method for 
% Compressed Sensing based Sampling of ECG Signals in Medical-IoT era".

clear all, close all, clc

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
% Create the folder "Results", where the performance results will be stored
mkdir Results

%% Database reading
file_name = 'mitdb/100'; % file name related to the ECG signal from the MIT-BIH arrhythmia database
wfdbdownload(file_name); % downloading of the ECG data from the database
[xa, Fs, tm] = rdsamp(['database/',file_name],1); % loading of the ECG signal in MATLAB
x = xa(1:360*60); % selection of 1 min of ECG data

% Parameters of the analysis
N = 720; % number of samples for acquired frame to be compressed
CR = 2:5; % compression ratios CRs
M = round(N./CR);
CR = N./M;

for ii=1:length(CR)
% CS processing of the ECG signal
    [xest,xcut,alpha,y,Phi,x_th,Psi] = CS_power_based(x,N,M(ii),0.2,'cvx'); % compression and reconstruction of the entire ECG signal
    PRD(ii) = norm(xcut-xest')/norm(xcut)*100; % evaluation of the PRD value for the specific CR
    % evaluation of the PRDf value
    mat_xcut = vec2mat(xcut',length(xcut)/N)';
    mat_xest = vec2mat(xest,length(xest)/N)';
    for jj = 1:length(xest)/N
        PRDf(ii,jj)= norm(mat_xcut(jj,:)-mat_xest(jj,:))/norm(mat_xcut(jj,:))*100;
    end
end

save('./Results/PRD_vs_CR_CS_power_based.mat')

for ii=1:length(CR)
    figure
    hist(PRDf(ii,:))
    xlabel('PRD per frame')
    ylabel('Counts')
    grid on
    set(gca,'FontSize',16)
    % percentile evaluation of the dataset
    PRD_95(ii) = prctile(PRDf(ii,:),95);
    % mean evaluation of the dataset PRDf
    PRD_avg(ii) = mean(PRDf(ii,:));
end

% PRD at 95th percentile vs CR
figure
plot(CR,PRD_95,'*-','LineWidth',2)
grid on
xlabel('CR')
ylabel('PRD at 95^{th} percentile [%]')
set(gca,'FontSize',16)

% Mean of PRDf vs CR
figure
plot(CR,PRD_avg,'*-','LineWidth',2)
grid on
xlabel('CR')
ylabel('PRD_{avg} [%]')
set(gca,'FontSize',16)

% Total PRD vs CR
figure
plot(CR,PRD,'*-','LineWidth',2)
grid on
xlabel('CR')
ylabel('PRD [%]')
set(gca,'FontSize',16)