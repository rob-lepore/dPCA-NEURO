clc;
clear;
close all;

% path = 'C:\Users\rober\OneDrive\Documents\Uni\Bioinformatics\Neurosciences\project\mef25_V6A\mef25_V6A\';
path = 'C:\Users\rober\OneDrive\Documents\Uni\Bioinformatics\Neurosciences\project\V6A_mef24\V6A_mef24\';

% % DATASET WITH ONLY ONE HAND
% principal_angles = full_dpca_analysis(path);

% % DATASET WITH TWO HANDS, CONSIDER ONLY ONE HAND
% principal_angles = full_dpca_analysis(path, 'hand_position', 1);
% principal_angles = full_dpca_analysis(path, 'hand_position', 2);

% DATASET WITH TWO HANDS, CONSIDER BOTH HANDS
principal_angles = full_dpca_analysis_both_hands(path);

disp(principal_angles)
