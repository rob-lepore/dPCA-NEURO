%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Calculate Firing Rates from data %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [firingRates, trialNum] = A_general_calculate_firing_rates_dpca(data_Path, cells_in_Directory, time_Window, sDF_bin_Size, event_Name, n_decisions);

% Event Marker --> 16 events
%   1)-,     2)KeyDown,      3)LED1,     4) Saccade-Off,  5)-,         6)-,        7)GO,         8) KeyUp,
%   9)-,     10) TOUCH1,     11)-,      12)RedOff,        13)TOUCH2,   14)-,       15)keydown,   16)-
Event_Marker.n16 = {'-', 'KeyDown', 'LED1', 'Saccade-Off', '-' , '-', 'GO', 'KeyUp', '-' , 'TOUCH1', '-' , 'RedOff', 'TOUCH2' , '-' , 'keydown', '-'};
%                           2         3          4                      7       8              9              10         11                12                            

N = length(cells_in_Directory);  % number of neurons   --> neurons from MEF animals 
T = sum(time_Window);            % number of time points
S = 9;                           % number of stimuli   --> 9-led panel with fixed position
D = n_decisions;                           % number of decisions --> animal starting point (near-far)
E = 10;                          % maximal number of trial repetitions (all trials should have the same length in time)

%--------------------------------------------------------------------------
%  1) trialNum: N x S x D                      ----------------------------
%  2) firingRates: N x S x D x T x maxTrialNum ----------------------------
%  3) firingRatesAverage: N x S x D x T        ----------------------------
%--------------------------------------------------------------------------

firingRates = zeros(N, S, D, T, E);
trialNum = zeros(N, S, D);
all_Trial = [];
count_Cell = 1;


for cell = 1 :length(cells_in_Directory) 
    data_Path_cell=strcat(string(data_Path), string(cells_in_Directory(cell).name));
    load([data_Path_cell])
    event_Marker = find(strcmp(Event_Marker.n16, event_Name));

        % fprintf(['Cell #:' num2str(cell)])
        for iHand_Position = 1 : D  % we have two hand positions
           % fprintf(['   Hand #: ' num2str(iHand_Position)])

            for iTarget_Position = 1 : 9  % we have nine target positions
                %fprintf(['   target #: ' num2str(iTarget_Position) ', '])

                num_Trial = sum(~cellfun(@isempty,Data.SpkTimeCS{iHand_Position, iTarget_Position}));  % how many trials we have for this condition
                index_Trial = ~cellfun(@isempty,Data.SpkTimeCS{iHand_Position, iTarget_Position});  % trial ind: 1: we had trial 0: we didn't have trial
                this_Trial = zeros(sum(index_Trial), round(max(cell2mat(cellfun(@max,Data.SpkTimeCS{iHand_Position, iTarget_Position}, 'UniformOutput', false))))+5000);  
                % to be filled by spikes 1ms resolution. 5000: just add some zeros at the end
                %fprintf(num2str(length(index_Trial), '%0.2d'))
                
                all_Trial = [all_Trial length(index_Trial)];
                
                trialNum(count_Cell, iTarget_Position, iHand_Position) = num_Trial;
                trial_Counter = 1;

                for iTrial = 1 : length(index_Trial)  % loop over all trials for one condition
                    if index_Trial(iTrial)
                        
                        % align spikes to this event
                        spike_Times = Data.SpkTimeCS{iHand_Position, iTarget_Position}{iTrial} - Data.EventsTimeMS{iHand_Position, iTarget_Position}{iTrial}(event_Marker);  
                        % then add spikes times with a time_Window(1), time before the event, as MATALB doesn't accept negative indexes
                        spike_Times = round(spike_Times + time_Window(1)); 
                        % only consider this_Trial = conv(sDF, ones(1, sDF_bin_Size), 'same')*(1/(sDF_bin_Size/1000));spikes after 0 ms (Wind(1))
                        spike_Times = spike_Times(spike_Times > 0);
                        % spike train 1 ms resolution for one trial
                        sDF = zeros(1, size(this_Trial, 2)); 
                        sDF(spike_Times) = 1;  % spike train
                        
                        this_Trial = conv(sDF, ones(1, sDF_bin_Size), 'same')*(1/(sDF_bin_Size/1000));
                        
                        firingRates(count_Cell, iTarget_Position, iHand_Position, :, trial_Counter) = this_Trial(1 : T);
                        trial_Counter = trial_Counter + 1;
                        
                    end
                end
            end
            %fprintf('\n')
        end
        count_Cell = count_Cell + 1;
end