%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Hemisphere information %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

% General import table  ---------------------------------------------------
%exported_data= uigetfile('', 'Select the Excel for the Exported MEF ** data', 'MultiSelect','off')
%db_hem_info= uigetfile('', 'Select the Database File Maker Pro', 'MultiSelect','off')
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
% path = C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\
db_hem_info=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\DATA BASE file maker pro 22+24.xlsx','ReadVariableNames',false);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT TABLE x mef22 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEF         exported_data 1st column 
% UNIT_1      exported_data 2nd column
% UNIT_1      exported_data 3rd column
% CITO        db_hem_info 2th column
% PROG        db_hem_info 4th column
% HEM         db_hem_info 3st column

% % IMPORT V6A mef22 == esportazioneMEF22tot.xls -------------------------------
exported_data=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\esportazioneMEF22tot.xls','ReadVariableNames',false);
sz=[(size(exported_data, 1)), 6];
varTypes = ["string","double","double","string","string","string"];
varNames = ["MEF","UNIT_1","UNIT_2","CITO", "PROG","HEM"];
MEF22_V6A=table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

% MEF ---------------------------------------------------------------------
MEF22_V6A(:,1) = exported_data(:,1);       %mef22
mef_name= table2array(exported_data(:,1)); %mef22
% UNIT_1 lettera ----------------------------------------------------------
MEF22_V6A(:,2) = exported_data(:,2);     %mef22
unit1 = table2array(exported_data(:,2)); %mef22
unit1= alphabet(unit1);
% UNIT_2 ------------------------------------------------------------------
MEF22_V6A(:,3) = exported_data(:,3);     %mef22
unit2 = table2array(exported_data(:,3)); %mef22

%name as in db_hem_info, i.e. 22,B084 -------------------------------------
for ii = 1:length(unit2)
    unit_value = strcat(num2str(mef_name(ii)), ',', unit1(ii), num2str(unit2(ii),'%03.f'));
    MEF22_V6A{ii,1} = {unit_value};
end
unit_value = table2array(MEF22_V6A(:,1));
unit_data= string(table2array(db_hem_info(:,1)));

%Find the number of row in db_hem_info I want to extract ------------------
for ii = 1:length(unit_value)
    indx_name(ii) = [find(strcmp(unit_value(ii), unit_data))];
end
%Table of only extracted cells, in the correct order
db_hem_extracted=db_hem_info(indx_name,:);

% CITO --------------------------------------------------------------------
MEF22_V6A(:,4) = db_hem_extracted(:,2);

% PROG --------------------------------------------------------------------
MEF22_V6A(:,5) = db_hem_extracted(:,4);

% HEM ---------------------------------------------------------------------
MEF22_V6A(:,6) = db_hem_extracted(:,3);

% writetable(MEF22_V6A, 'C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\MEF22_V6A+HEM.xls','WriteVariableNames', true); 


% % Delete the cell that are not labelled as APO* = V6A -------------------
cito= table2array(MEF22_V6A(:,"CITO"));
c=0;
for ii = 1:length(cito)
    if contains(string(cito(ii)), "PE")
        c=c+1;
        indx_cell2del(c,1)= [ii];
    end
end

MEF22_onlyV6A = MEF22_V6A
MEF22_onlyV6A(indx_cell2del,:) = [];
%writetable(MEF22_onlyV6A, 'C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\MEF22_onlyV6A+HEM.xls','WriteVariableNames', true); 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT TABLE x mef24 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IMPORT V6A mef24 = esportazioneMEF24.xls ------------------------------
% exported_data=readtable('C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\esportazioneMEF24.xls','ReadVariableNames',false);
% exported_data([169:173],:) = []; %solo per mef24, ci sono delle note sul fondo dell'Excel
% sz=[(size(exported_data, 1)), 6];
% varTypes = ["string","double","double","string","string","string"];
% varNames = ["MEF","UNIT_1","UNIT_2","CITO", "PROG","HEM"];
% MEF24_V6A=table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
% 
% % MEF ---------------------------------------------------------------------
% MEF24_V6A(:,1) = exported_data(:,2);  %mef24, scambiati le righe nell'Excel
% mef_name= table2array(exported_data(:,2));
% % UNIT_1 lettera ----------------------------------------------------------
% MEF24_V6A(:,2) = exported_data(:,3);  %mef24, scambiati le righe nell'Excel
% unit1 = table2array(exported_data(:,3));
% unit1= alphabet(unit1);
% % UNIT_2 ------------------------------------------------------------------
% MEF24_V6A(:,3) = exported_data(:,4);  %mef24, scambiati le righe nell'Excel
% unit2 = table2array(exported_data(:,4));
% 
% % name as in db_hem_info, i.e. 22,B084 -------------------------------------
% for ii = 1:length(unit2)
%     unit_value = strcat(num2str(mef_name(ii)), ',', unit1(ii), num2str(unit2(ii),'%03.f'));
%     MEF24_V6A{ii,1} = {unit_value};
% end
% unit_value = table2array(MEF24_V6A(:,1));
% unit_data= string(table2array(db_hem_info(:,1)));
% 
% %Find the number of row in db_hem_info I want to extract ------------------
% for ii = 1:length(unit_value)
%     indx_name(ii) = [find(strcmp(unit_value(ii), unit_data))];
% end
% %Table of only extracted cells, in the correct order
% db_hem_extracted=db_hem_info(indx_name,:);
% 
% % CITO --------------------------------------------------------------------
% MEF24_V6A(:,4) = db_hem_extracted(:,2);
% 
% % PROG --------------------------------------------------------------------
% MEF24_V6A(:,5) = db_hem_extracted(:,4);
% 
% % HEM ---------------------------------------------------------------------
% MEF24_V6A(:,6) = db_hem_extracted(:,3);
% 
% %writetable(MEF24_V6A, 'C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\MEF24_V6A+HEM.xls','WriteVariableNames', true); 

% % Delete the cell that are not labelled as APO* = V6A -------------------
% cito= table2array(MEF24_V6A(:,"CITO"));
% c=0;
% for ii = 1:length(cito)
%     if contains(string(cito(ii)), "PE")
%         c=c+1;
%         indx_cell2del(c,1)= [ii];
%     end
% end
% 
% MEF24_onlyV6A = MEF24_V6A
% MEF24_onlyV6A(indx_cell2del,:) = [];
% %writetable(MEF24_onlyV6A, 'C:\Users\Martina\Desktop\Tesi\dpca_matlab\Cells_data\DATI_V6A\MEF24_onlyV6A+HEM.xls','WriteVariableNames', true); 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATABASE V6A mef22-mef24 RD task %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Extract table rows from db_hem_info relative to V6A and RD task
% find_APO = table2cell(db_hem_info(:,"Var2"));
% find_task = table2cell(db_hem_info(:,"Var4"));
% c=0; db_V6A= table;
% 
% for ii = 1:length(find_APO)
%     if contains(string(find_APO(ii)), "APO") & contains(string(find_task(ii)), "RD")
%         c=c+1
%         db_V6A = [db_V6A; db_hem_info(ii, :)]
%     end
% end