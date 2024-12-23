%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   dPCA    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
clc
tic

% PSTH properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATASET ------------------------------------------------------------------
list_dataset= {'PE_mef24', 'PE_mef25', 'PEc_mef24', 'PEc_mef25', 'V6A_mef22', 'V6A_mef24', 'prova_PE_mef24', 'tot_PE', 'tot_PEc', 'tot_V6A', 'tot_mef24', 'tot_mef24_PE_PEc', 'tot_mef25'};
[indx_ds] = listdlg('ListString',list_dataset, 'SelectionMode', 'single');
% The folder with cells is present in the same path as this function
data_Path=strcat(string(cd),'\Cells_data\', string(list_dataset(indx_ds)), '\');
data_Path="C:\Users\rober\OneDrive\Desktop\Uni\Bioinformatics\Neurosciences\project\V6A_mef24\V6A_mef24\"
cells_in_Directory = dir(data_Path);
cells_in_Directory ([1,2],:) = [];

%EVENT DESCRIPTION AND TIME WINDOW ----------------------------------------
% ['-', 'KeyDown', 'LED1', 'Saccade-Off', '-' , '-', 'GO', 'KeyUp', '-' , 'TOUCH1', '-' , 'RedOff', 'TOUCH2' , '-', 'keydown', '-'];
list_marker= {'KeyDown', 'LED1', 'Saccade-Off', 'GO', 'KeyUp', 'TOUCH1', 'RedOff', 'TOUCH2', 'keydown'};
[indx_am] = listdlg('ListString',list_marker, 'SelectionMode', 'single');

event_Name= string(list_marker(indx_am));    % the event marker for the first PSTH 
if event_Name == 'TOUCH1'
    time_Window = [500, 1800];               % time before and after the event
else
    time_Window = [500,1000];                  % time before and after the event
end

%SPIKE DENSITY FUNCTION ---------------------------------------------------
sDF_bin_Size = 100;              % the width of window for spike density function

% SIMOULTANEOUS RECORDINGS ------------------------------------------------
ifSimultaneousRecording = true;  % change this to simulate simultaneous recording
                                 % (imply same number of trials for each neuron)

dataname= strcat(string(list_dataset(indx_ds)), {' - '}, string(list_marker(indx_am)));
disp([strcat('Dataset and event marker used:', {' '}, dataname)])


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


%% Define parameter grouping
% K=E;S;Q=D;T;
% firingRates array has [N S D T E] size; here we ignore the 1st dimension
% (neurons), i.e. we have the following parameters:
%    1 - stimulus 
%    2 - decision
%    3 - time
% There are three pairwise interactions:
%    [1 3] - stimulus/time interaction
%    [2 3] - decision/time interaction
%    [1 2] - stimulus/decision interaction
% And one three-way interaction:
%    [1 2 3] - rest
% As explained in the eLife paper (Kobak 5et al.,2016), we group stimulus with
% stimulus/time interaction etc.: {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}}.

%'combinedParams'- cell array of cell arrays specifying 
%                  which marginalizations should be added up together
combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Target', 'Hand', 'Condition-independent', 'T/H Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = time_Window(1);

% Check consistency between trialNum and firingRates
for n = 1:size(firingRates,1)
    for s = 1:size(firingRates,2)
        for d = 1:size(firingRates,3)
            assert(isempty(find(isnan(firingRates(n,s,d,:,1:trialNum(n,s,d))), 1)), ...
                'Something is wrong!')
        end
    end
end
clear n s d 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: PCA of the dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Step 1: PCA of the dataset'])
fig_count=0;

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

% SVD ---------------------------------------------------------------------
%svd returns the singular values of matrix X in descending order.
%[___] = svd(A,"econ") produces an economy-size decomposition of A 

% W = matrix of size #neuron x #neuron
[W,~,~] = svd(X, 'econ'); 
% W = matrix of size #neuron x 1:20 rows = component
W = W(:,1:20);

% Computing explained variance --------------------------------------------
% explVar = dpca_explainedVariance(X, W, V)      X is the data matrix
%                                                W is the decoder matrix
%                                                V is the encoder matrix
% computes various measures of % explained variance and returns them in a
% structure explVar.Returned values:
%
%  * explVar.totalVar             - total variance
%  * explVar.totalMarginalizedVar - total variance in each marginalization
%  * explVar.componentVar         - variance of each component (%)
%  * explVar.margVar              - variance of each component in each marginalization (%)
%  * explVar.cumulativePCA        - cumulative variance of the PCA components (%)
%  * explVar.cumulativeDPCA       - cumulative variance of the dPCA components (%)

explVar1 = dpca_explainedVariance(firingRatesAverage, W, W, ...
    'combinedParams', combinedParams);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2: PCA in each marginalization separately %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Step 2: PCA in each marginalization separately'])

% dpca_perMarginalization --> dpca_marginalize (ENTER)
% dpca_perMarginalization = PLOT FIGURE (PCA for each marginaliztion)
fig_count= fig_count + 1;
dpca_perMarginalization(firingRatesAverage, dataname, fig_count, ...
    @dpca_plot_default, ...
    'combinedParams', combinedParams, ...
    'timeEvents', timeEvents);

%saveas(1,'C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\dPCA_img\Fig1_PCAxMarginalization_SaccadeOff_tot_mef24_PEPEc')  %+++

% dpca_perMarginalization(X, plotFunction, ...) performs PCA in each
% marginalization of X and plots the components using plotFunction, a
% pointer to the function that plots one component
%
% --> YY = dpca_marginalize(X) computes data marginalized over all combinations
%   of parameters. 
%   X = multi-dimensional array of dimensionality D+1, where the first dimension
%       corresponds to neurons and the rest D dimensions to various parameters.
%   YY = cell array of marginalized datasets,containing 2^D-1 arrays,
%       marginalized over all combinations of D parameters, excluding empty set.
%   For each i size(YY{i}) equals size(X).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3: dPCA without regularization and ignoring noise covariance %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Step 3: dPCA without regularization and ignoring noise covariance'])

% CORE FUNCTION -----------------------------------------------------------
% [W, V, whichMarg] = dpca(X, numComp, ...) performs dPCA on the data in X
% and returns decoder matrix W and encoder matrix V.
numComp = 15; % numComp=20 default

%tic
[W,V,whichMarg] = dpca(firingRatesAverage, numComp, ... 
    'combinedParams', combinedParams);
%toc

%INPUT
%    X       multi-dimensional array of dimensionality D+1, where the first
%            dimension corresponds to neurons and the rest D dimensions to
%            various parameters.
% numComp    specifies the number of dPCA components to be extracted (can be
%            either one number or a list of numbers for each marginalization)
% whichMarg  is an array of integers providing the 'type' of each component
%            (which marginalization it describes).
% If the total number of required components is S=sum(numComp),then W and V
% are of NxS size, and whichMarg has length S.

%OUTPUT
%    W      is the decoder
%    V      is the encoder (ordered by explained variance),

% dpca_explainedVariance --> dpca_marginalize  (ENTER)
explVar2 = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

fig_count= fig_count + 1;
fig_name="dPCA without regularization and Cnoise";
myFig = figure('Render', "painters", 'Position', [0 0 1800 1000], ...
    'Name', 'Demixed-Principal Component Analysis');

dpca_plot(firingRatesAverage, W, V, dataname,fig_name, fig_count, ...
    @dpca_plot_default, ...
    'explainedVar', explVar2, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg, ...
    'time', time, ...
    'timeEvents', timeEvents, ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);

%saveas(2,'C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\dPCA_img\Fig2_dpca_NoReg_SaccadeOff_tot_mef24_PEPEc')  %+++

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4: dPCA with regularization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Step 4: dPCA with regularization'])

% This function takes some minutes to run. It will save the computations 
% in a .mat file with a given name. Once computed, you can simply load 
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself
% (even with lambda set to zero).

fig_count= fig_count + 1;
optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, dataname, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 2, ...  % (2) increase this number to ~10 for better accuracy ***
    'filename', 'tmp_optimalLambdas.mat');


% optimalLambda = dpca_optimizeLambda(X, Xtrial, numOfTrials, ...)
% computes optimal regularization parameter. X is the data array. Xtrial 
% is an array storing single trials. It has one extra dimension as compared 
% with X and stores individual single trial firing rates, as opposed to the 
% trial average. numOfTrials has one dimension fewer than X and for each 
% neuron and combination of parameters (without time) specifies the number 
% of available trials in X_trial. All entries have to be larger than 1.
%
% This code assumes that time parameter is stored in the last dimension of
% X. For datasets without time, some other cross-validation needs to be
% used.
%
% [optimalLambda, optimalLambdas] = dpca_optimizeLambda(...) additionally
% returns a list of optimal lambdas found separately for each
% marginalization

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W,V,whichMarg] = dpca(firingRatesAverage, numComp, ... % numComp=20
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar3 = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

%saveas(3,'C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\dPCA_img\Fig3_relative_crossval_errors_SaccadeOff_tot_mef24_PEPEc')  %+++

%--------------------------------------------------------------------------
fig_count= fig_count + 1;
fig_name="dPCA with Regularization";
myFig = figure('Render', "painters", 'Position', [0 0 1800 1000], ...
    'Name', 'Demixed-Principal Component Analysis with Regularization');

dpca_plot(firingRatesAverage, W, V, dataname, fig_name, fig_count, ...
    @dpca_plot_default, ...
    'explainedVar', explVar3, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);

%saveas(4,'C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\dPCA_img\Fig4_dpca_regularization_SaccadeOff_tot_mef24_PEPEc')  %+++

%% Optional: estimating "signal variance"

disp(['Step 4.1 Optional: estimating "signal variance"'])

explVar4 = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams, ...
    'Cnoise', Cnoise, 'numOfTrials', trialNum);

% Note how the pie chart changes relative to the previous figure.
% That is because it is displaying percentages of (estimated) signal PSTH
% variances, not total PSTH variances. See paper for more details.

fig_count= fig_count + 1;
fig_name = "dPCA with Regularization - Estimating signal variance";
myFig = figure('Render', "painters", 'Position', [0 0 1800 1000], ...
    'Name', 'Demixed-Principal Component Analysis with Regularization - Estimating signal variance');

dpca_plot(firingRatesAverage, W, V, dataname, fig_name, fig_count, ...
    @dpca_plot_default, ...
    'explainedVar', explVar4, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);

%saveas(5,'C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\dPCA_img\Fig5_dpca_regularization_signalVar_SaccadeOff_tot_mef24_PEPEc')  %+++


%% Optional: decoding

decodingClasses = {[(1:S)' (1:S)'], repmat([1:2], [S 1]), [], [(1:S)' (S+(1:S))']};

accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 2, ...        % (5) increase to 100 ***
    'filename', 'tmp_classification_accuracy.mat');

%dpca_classificationPlot(accuracy, [], [], [], decodingClasses)

accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 2, ...        % (5) increase to 100 ***
    'numShuffles', 2, ...  % (20) increase to 100 *** (time consuming)
    'filename', 'tmp_classification_accuracy.mat');

fig_count= fig_count + 1;
dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses)

componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);

%saveas(6,'C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\dPCA_img\Fig6_classification_acuracy_xMarg_SaccadeOff_tot_mef24_PEPEc')  %+++

%--------------------------------------------------------------------------
fig_count= fig_count + 1;
fig_name="dPCA with Regularization: Decoding";
myFig = figure('Render', "painters", 'Position', [0 0 1800 1000], ...
    'Name', 'Demixed-Principal Component Analysis with Regularization - Decoding');

dpca_plot(firingRatesAverage, W, V, dataname, fig_name, fig_count, ...
    @dpca_plot_default, ...
    'explainedVar', explVar4, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16,                ...
    'componentsSignif', componentsSignif);

%saveas(7,'C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\dPCA_img\Fig7_dpca_regularization_decoding_SaccadeOff_tot_mef24_PEPEc')        %+++
%saveas(7,'C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\dPCA_img\Fig7_dpca_regularization_decoding_SaccadeOff_tot_mef24_PEPEc','epsc') %+++

toc