%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Selection Mean Firing Rates from data %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

%--------------------------------------------------------------------------
%dei neuroni rimasti elimino quelli che non fanno almeno un numero minimo di spike (minimum_spike)
list_dataset= {'PE_mef24', 'PE_mef25', 'PEc_mef24', 'PEc_mef25', 'V6A_mef22', 'V6A_mef24', 'prova_PE_mef24', 'tot_PE', 'tot_PEc', 'tot_V6A', 'tot_mef24', 'tot_mef24_PE_PEc', 'tot_mef25'};
%                   1           2            3          4              5           6              7              8         9          10         11               12               13

list_epoch={'Saccade-Off -- KeyUp', 'KeyUp -- TOUCH1', 'TOUCH1 -- TOUCH1 + 1800'};
list_marker= {'-', 'KeyDown', 'LED1', 'Saccade-Off', '-' , '-', 'GO', 'KeyUp', '-' , 'TOUCH1', '-' , 'RedOff', 'TOUCH2' , '-' , 'keydown', '-'};
%               1       2         3          4         5     6     7      8      9       10       11     12        13       14     15        16
neu2remove=[];
max_spike= 50;

for ds = 1:length(list_dataset)
    data_Path=strcat(string(cd),'\Cells_data\', string(list_dataset(ds)), '\');
    cells_in_Directory = dir(data_Path);
    cells_in_Directory ([1,2],:) = [];
    disp([strcat({'Dataset: '}, string(list_dataset(ds)))])
    disp([strcat({'# cell: '}, string(length(cells_in_Directory)))])
    cell=0;

    for cell = 1:length(cells_in_Directory)
        data_Path_cell=strcat(string(data_Path), string(cells_in_Directory(cell).name));
        load(data_Path_cell)
        frequenza_all{cell} = []; %imposto il vettore frequenza che si riempirà con le frequenze di ogni ciclo

        for hand_pos = 1 : 2  % we have two hand positions
            for target_pos = 1 : 9  % we have nine target positions
                index_trial = ~cellfun(@isempty,Data.SpkTimeCS{hand_pos, target_pos});

                for trial = 1 : length(index_trial)  % loop over all trials for one condition
                    if index_trial(trial) == 1

                        for ep = 1:length(list_epoch)
                            % {'-', 'KeyDown', 'LED1', 'Saccade-Off', '-' , '-', 'GO', 'KeyUp', '-' , 'TOUCH1', '-' , 'RedOff', 'TOUCH2' , '-' , 'keydown', '-'};
                            %   1       2         3          4         5     6     7      8      9       10       11     12        13       14     15        16

                            if isequal(string(list_epoch(ep)),'Saccade-Off -- KeyUp')
                                m_i=4; %Saccade-Off
                                m_f=8; % KeyUp
                                marker_i = find(strcmp(list_marker(m_i), list_marker)); %imposto un marker iniziale per l'epoch
                                marker_f = find(strcmp(list_marker(m_f), list_marker)); %imposto un marker finale per l'epoch
                            elseif isequal(string(list_epoch(ep)),'KeyUp -- TOUCH1')
                                m_i=8; %KeyUp
                                m_f=10; %TOUCH1
                                marker_i = find(strcmp(list_marker(m_i), list_marker)); %imposto un marker iniziale per l'epoch
                                marker_f = find(strcmp(list_marker(m_f), list_marker)); %imposto un marker finale per l'epoch
                            elseif isequal(string(list_epoch(ep)),'TOUCH1 -- TOUCH1 + 1800')
                                m_i=10; %TOUCH1
                                marker_i = find(strcmp(list_marker(m_i), list_marker)); %imposto un marker iniziale per l'epoch
                                marker_f = find(strcmp(list_marker(m_i), list_marker)); %imposto un marker finale per l'epoch
                            end

                            i_ep = Data.EventsTimeMS{hand_pos, target_pos}{1, trial}(marker_i) ; %imposto l'inizio dell'epoch
                            f_ep = Data.EventsTimeMS{hand_pos, target_pos}{1, trial}(marker_f); %imposto la fine dell'epoch

                            if isequal(string(list_epoch(ep)),'TOUCH1 -- TOUCH1 + 1800')
                                f_ep = f_ep + 1800;
                            end

                            %trovo le spikes comprese nelle epoch considerate
                            idx = find(Data.SpkTimeCS{hand_pos, target_pos}{1, trial} > i_ep & Data.SpkTimeCS{hand_pos, target_pos}{1,trial} < f_ep);
                            num_spikes = length(idx); %calcolo la lunghezza del vettore idx, cioè il numero di scariche nell'epoch
                            durata_ep = f_ep-i_ep; %calcolo la durata dell'epoch in s

                            frequenza = num_spikes*1000/durata_ep; %calcolo la frequenza convertendo in ms
                            frequenza_all{cell} = [frequenza_all{cell} frequenza]; %inserisco i valori della frequenza nel vettore vuoto
                            mean_firing_rate(cell) = mean(frequenza_all{cell},2); %calcolo la media per ogni neurone della frequenza di scarica
                        end
                    end
                end
            end
        end
    end
    neu2remove=unique([neu2remove find(mean_firing_rate > max_spike)]);
    neu_toremove.(string(list_dataset(ds)))= neu2remove;
    mean_fr_areas.(string(list_dataset(ds)))=mean_firing_rate; %***
end
save('mean_fr_areas.mat', 'mean_fr_areas')
save('neu_toremove.mat', 'neu_toremove')