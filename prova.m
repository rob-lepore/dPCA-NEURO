% clc;
% clear;
% close all;
% 
% data_Path = 'C:\Users\rober\OneDrive\Desktop\Uni\Bioinformatics\Neurosciences\project\V6a_mef24\V6A_mef24\';
% cells_in_Directory = dir(data_Path);
% cells_in_Directory ([1,2],:) = [];
% 
% markers = {'Saccade-off', 'GO'};
% events = [200 800; 700 0];
% sDF_bin_Size = 100;              
% ifSimultaneousRecording = true;  
% 
% for i=1:length(markers)
%     eventName = markers(i);
%     timeWindow = events(i,:);
%     [firingRates, trialNum] = A_general_calculate_firing_rates_dpca(data_Path, cells_in_Directory, timeWindow, sDF_bin_Size, eventName);
%     firingRatesAverage = mean(firingRates, 5); %nanmean
%     size(firingRatesAverage)
% end
% 

eventPositions = [4,7,9];
eventsTime = [];

for iEvent=1:length(eventPositions)
    eventPosition = eventPositions(iEvent);
    
    % Initialize a variable to store the sum of the fourth values
    sumFourthValues = 0;
    
    % Initialize a counter to track the number of elements processed
    numElements = 0;
    
    % Loop over the 2x9 cell array
    for i = 1:size(Data.EventsTimeMS, 1)
        for j = 1:size(Data.EventsTimeMS, 2)
            % Access the 1xN cell array (N can vary)
            cellArray1xN = Data.EventsTimeMS{i, j};
            
            % Loop over the 1xN cell array
            for k = 1:length(cellArray1xN)
                % Access the array of doubles
                arrayDoubles = cellArray1xN{k};
                
                % Ensure the array has at least 4 elements
                if length(arrayDoubles) >= 4
                    % Extract the fourth value and add it to the sum
                    sumFourthValues = sumFourthValues + arrayDoubles(eventPosition);
                    
                    % Increment the counter
                    numElements = numElements + 1;
                end
            end
        end
    end
    
    % Check if numElements is greater than 0 to avoid division by zero
    if numElements > 0
        % Calculate the average
        averageFourthValue = sumFourthValues / numElements;
        
        % Display the result
        fprintf('Average value of the %d number: %f\n',eventPosition, averageFourthValue);
        eventsTime(iEvent) =  averageFourthValue;
    else
        disp('No valid fourth values found.');
    end
end

fprintf("Fixation: %f %f\n",eventsTime(1)-200,eventsTime(1)+800)
fprintf("Plan: %f %f\n",eventsTime(2)-700,eventsTime(2))
fprintf("Reach: %f %f\n",eventsTime(2),eventsTime(3))
fprintf("Hold: %f %f\n",eventsTime(3),eventsTime(3)+800)

