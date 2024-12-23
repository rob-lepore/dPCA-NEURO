%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Hemisphere information %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%list_dataset= {'PE_mef24', 'PE_mef25', 'PEc_mef24', 'PEc_mef25', 'prova_PE_mef24', 'tot_PE', 'tot_PEc', 'tot_mef24_PE+PEc', 'tot_mef25'};
list_dataset= {'PE_mef24', 'PE_mef25', 'PEc_mef24', 'PEc_mef25', 'V6A_mef22', 'V6A_mef24', 'tot_PE', 'tot_PEc', 'tot_V6A', 'tot_mef24', 'tot_mef24_PE_PEc', 'tot_mef25'};
%                   1           2            3          4              5           6          8         9          10         11               12               13
[indx_ds] = listdlg('ListString',list_dataset, 'SelectionMode', 'single');
data_Path=strcat(string(cd),'\Cells_data\', string(list_dataset(indx_ds)), '\');
cells_in_Directory = dir(data_Path);
count_cells=0;
cells_in_Directory ([1,2],:) = [];

disp(string(list_dataset(indx_ds)))
dataset=string(list_dataset(indx_ds));

    if dataset == 'PE_mef24';
    T=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF24\MEF24_PE+HEM.xlsx');
    T([4,7,9],:) = [];     % neurons to be deleted
elseif dataset == 'PEc_mef24';
    T=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF24\MEF24_PEc+HEM.xlsx');
    T([4,9,10],:) = [];    % neurons to be deleted
elseif dataset == 'PE_mef25';
    T=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF25\MEF25_PE+HEM.xlsx');
    T([2,27],:) = [];      % neurons to be deleted
elseif dataset =='PEc_mef25';
    T=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF25\MEF25_PEc+HEM.xlsx');
    T([11,25,32],:) = [];   % neurons to be deleted
elseif dataset =='tot_PE';  %Folder order: PE_mef24 + PE_mef25 
    T1=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF24\MEF24_PE+HEM.xlsx');
    T1([4,7,9],:) = [];     % neurons to be deleted
    T2=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF25\MEF25_PE+HEM.xlsx');
    T2([2,27],:) = [];      % neurons to be deleted
    T= vertcat(T1,T2);
elseif dataset =='tot_PEc';  %Folder order: PEc_mef24 + PEc_mef25 
    T1=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF24\MEF24_PEc+HEM.xlsx');
    T1([4,9,10],:) = [];    % neurons to be deleted
    T2=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF25\MEF25_PEc+HEM.xlsx');
    T2([11,25,32],:) = [];  % neurons to be deleted
    T= vertcat(T1,T2);
elseif dataset =='tot_mef24_PE_PEc'; %Folder order: PE_mef24 + PEc_mef24
    T1=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF24\MEF24_PE+HEM.xlsx');
    T1([4,7,9],:) = [];     % neurons to be deleted
    T2=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF24\MEF24_PEc+HEM.xlsx');
    T2([4,9,10],:) = [];    % neurons to be deleted
    T= vertcat(T1,T2);
elseif dataset == 'tot_mef25'; %Folder order: PE_mef25 + PEc_mef25
    T1=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF25\MEF25_PE+HEM.xlsx');
    T1([2,27],:) = [];      % neurons to be deleted
    T2=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF25\MEF25_PEc+HEM.xlsx');
    T2([11,25,32],:) = [];  % neurons to be deleted
    T= vertcat(T1,T2);
elseif dataset == 'V6A_mef22';
    T=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\MEF22_onlyV6A+HEM.xls');
    T([8],:) = [];          % neurons to be deleted
elseif dataset == 'V6A_mef24';
    T=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\MEF22_onlyV6A+HEM.xls');
elseif dataset == 'tot_V6A'; %Folder order: V6A_mef22 + V6A_mef24
    T1=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\MEF22_onlyV6A+HEM.xls');
    T1([8],:) = [];          % neurons to be deleted
    T2=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\MEF24_onlyV6A+HEM.xls');
    T= vertcat(T1,T2);
elseif dataset == 'tot_mef24'; %Folder order: PE_mef24 + PEc_mef24 + V6A_mef24
    T1=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF24\MEF24_PE+HEM.xlsx');
    T1([4,7,9],:) = [];     % neurons to be deleted
    T1=table2cell(T1)
    T2=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF24\MEF24_PEc+HEM.xlsx');
    T2([4,9,10],:) = [];    % neurons to be deleted
    T2=table2cell(T2)
    T3=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\MEF24_onlyV6A+HEM.xls');
    T3 = renamevars(T3,'UNIT_2','Var3');
    T3 = renamevars(T3,'UNIT_1','UNIT');
    T3=table2cell(T3)
    T= vertcat(T1,T2,T3);
    T=cell2table(T)
else dataset == 'prova_PE_mef24'; 
    T=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_MEF24\MEF24_PE+HEM.xlsx');
    T([4,7,9],:) = []; % neurons to be deleted

end

path_tosave=strcat(string('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_hem\'))
T= T(:,6);
T=table2array(T);

for cell = 1 :length(cells_in_Directory)
    data_Path_cell=strcat(string(data_Path), string(cells_in_Directory(cell).name));
    load(data_Path_cell);
    Data.Hem=T(cell,1);

    filename=strcat(path_tosave, list_dataset(indx_ds),'_Cell_', num2str(count_cells, '%0.3d'), '.mat');
    save(filename, 'Data');
    disp(strcat('Saved: ', filename))
 
    hem_info_dataset(cell)=T(cell,1);
    count_cells= count_cells+1;
end