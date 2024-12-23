clear 
close all
clc

data_Path="C:\Users\rober\OneDrive\Desktop\Uni\Bioinformatics\Neurosciences\project\V6a_mef24\V6A_mef24\";
cells_in_Directory = dir(data_Path);
cells_in_Directory ([1,2],:) = [];

event_Name= "Saccade-Off";    % the event marker for the first PSTH 
time_Window = [500,1000];                  % time before and after the event


%SPIKE DENSITY FUNCTION ---------------------------------------------------
sDF_bin_Size = 100;              % the width of window for spike density function

% SIMOULTANEOUS RECORDINGS ------------------------------------------------
ifSimultaneousRecording = true;  % change this to simulate simultaneous recording
                                 % (imply same number of trials for each neuron)
dataname="V6A_mef24 - Saccade-Off";


% CALCULATE FIRING RATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[firingRates, trialNum] = A_general_calculate_firing_rates_dpca(data_Path, cells_in_Directory, time_Window, sDF_bin_Size, event_Name);

% Calculate FiringRates considering the information about hemisphere
%[firingRates, trialNum] = A_hem_calculate_firing_rates_dpca(data_Path, cells_in_Directory, time_Window, sDF_bin_Size, event_Name);

N = size(firingRates, 1);    % number of neurons
S = size(firingRates, 2);    % number of stimuli conditions       (F1 frequencies in Romo's task, TARGET POSITION in this task)
D = size(firingRates, 3);    % number of decisions                (D=2)
T = size(firingRates, 4);    % number of time points              (all trials should have the same length in time!)
E = size(firingRates, 5);    % maximal number of trial repetitions
time = (1:T); %/ 10;

% For Romo-like task, data are joined in three arrays of the following sizes:
% trialNum: N x S x D
% firingRates: N x S x D x T x maxTrialNum
% firingRatesAverage: N x S x D x T

% trialNum -----  n# of trials for each neuron in each S,D condition (usually
%                 different for different conditions and different sessions)
% firingRates --- all single-trial data together, massive array. 
% maxTrialNum --- here is the maximum value in trialNum.
%                 E.g. if the number of trials per condition varied between
%                 1 and 20, then maxTrialNum = 20.
% For the neurons and conditions with less trials, fill remaining entries 
% in firingRates with zeros or nans.

disp(['Calculate_firing_rates_dpca, Cell #: ' num2str(size(trialNum,1))])

% Computing PSTHs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% firingRatesAverage -- average of firingRates over trials (5th dimension).
firingRatesAverage = mean(firingRates, 5); %nanmean

% If the firingRates is filled up with nans, then it's simply
%    firingRatesAverage = nanmean(firingRates,5)
% If it's filled up with zeros (as is convenient if it's stored on hard drive
% as a sparse matrix), then 
%    firingRatesAverage = bsxfun(@times, mean(firingRates,5), ...
%        size(firingRates,5)./trialNum)

%'combinedParams'- cell array of cell arrays specifying 
%                  which marginalizations should be added up together
combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Target', 'Hand', 'Condition-independent', 'T/H Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = time_Window(1);

disp(['dPCA with regularization'])

optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, dataname, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 2, ...  % (2) increase this number to ~10 for better accuracy ***
    'filename', 'tmp_optimalLambdas.mat');

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

numComp = 15; % numComp=20 default

[W,V,whichMarg] = dpca(firingRatesAverage, numComp, ... % numComp=20
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

save("results", "W", "V", "whichMarg", "firingRatesAverage");