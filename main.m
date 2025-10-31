clc;
clear;
close all;

path = 'C:\Users\rober\OneDrive\Documents\Uni\Bioinformatics\Neurosciences\project\V6A_mef24\V6A_mef24\';
hand_position = 1;

principal_angles = full_dpca_analysis(path, hand_position);
disp(principal_angles)