%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Selection Mean Firing Rates from data %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

%--------------------------------------------------------------------------
% Dataset
list_dataset= {'PE_mef24', 'PE_mef25', 'PEc_mef24', 'PEc_mef25', 'V6A_mef22', 'V6A_mef24', 'tot_PE', 'tot_PEc', 'tot_V6A', 'tot_mef24', 'tot_mef24_PE_PEc', 'tot_mef25'};
%                   1           2            3          4              5           6          7          8          9           10             11                 12 

list_marker={'Saccade-Off', 'GO', 'KeyUp', 'TOUCH1', 'RedOff', 'TOUCH2'};
list_marker_tot= {'-', 'KeyDown', 'LED1', 'Saccade-Off', '-' , '-', 'GO', 'KeyUp', '-' , 'TOUCH1', '-' , 'RedOff', 'TOUCH2' , '-' , 'keydown', '-'};
%                  1       2         3          4         5     6    7      8      9       10       11     12        13       14     15        16

marker_alignment_tot = table;

for ds = 1:length(list_dataset)
    data_Path=strcat(string(cd),'\Cells_data\', string(list_dataset(ds)), '\');
    cells_in_Directory = dir(data_Path);
    cells_in_Directory ([1,2],:) = [];
    disp([strcat({'Dataset: '}, string(list_dataset(ds)))])
    disp([strcat({'# cell: '}, string(length(cells_in_Directory)))])
    cell=0;

    sz=[length(cells_in_Directory), 6];
    varTypes = ["double","double","double","double","double","double"];
    varNames = ["Saccade-Off", "GO", "KeyUp", "TOUCH1", "RedOff", "TOUCH2"];
    marker_alignment_cell=table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

    for cell = 1:length(cells_in_Directory)
        data_Path_cell=strcat(string(data_Path), string(cells_in_Directory(cell).name));
        load(data_Path_cell)
        %marker_align_all{cell} = [];

        for hand_pos = 1 : 2  % we have two hand positions
            for target_pos = 1 : 9  % we have nine target positions
                index_trial = ~cellfun(@isempty,Data.SpkTimeCS{hand_pos, target_pos});

                for trial = 1 : length(index_trial)  % loop over all trials for one condition
                    if index_trial(trial) == 1

                        for am = 1:length(list_marker) %alignment marker
                            %{'Saccade-Off', 'GO', 'KeyUp', 'TOUCH1', 'RedOff', 'TOUCH2'}
                            %        1         2      3        4          5         6
                            % {'-', 'KeyDown', 'LED1', 'Saccade-Off', '-' , '-', 'GO', 'KeyUp', '-' , 'TOUCH1', '-' , 'RedOff', 'TOUCH2' , '-' , 'keydown', '-'};
                            %   1       2         3          4         5     6     7      8      9       10       11     12        13       14     15        16
                            indx_marker= find(strcmp((list_marker(am)), list_marker_tot));

                            % dati della V6A allineati su Saccade-Off
                            if ~startsWith(string(list_dataset(ds)), "V6A")
                                marker_alignment_cell{cell, am} = Data.EventsTimeMS{hand_pos, target_pos}{1, trial}(indx_marker);
                            else
                                % Allineamento sul KeyUp per poter calcolare media
                                marker_alignment_cell{cell, am}= Data.EventsTimeMS{hand_pos, target_pos}{1, trial}(indx_marker) - Data.EventsTimeMS{hand_pos, target_pos}{1, trial}(8);
                            end
                        end
                    end
                end
            end
        end
    end

    marker_alignment_tot= vertcat(marker_alignment_tot,marker_alignment_cell);

    %Compute Mean and STDV for each cell
    for ii = 1:size(marker_alignment_cell,2)
        mean_am=mean(table2array(marker_alignment_cell(:,ii),1));
        stdev_am=std(table2array(marker_alignment_cell(:,ii),1));
        my_field = strcat(string(list_dataset(ds)));
        marker_alignment_mean.(my_field)(1,ii) = mean_am;
        marker_alignment_mean.(my_field)(2,ii) = stdev_am;
    end
end


