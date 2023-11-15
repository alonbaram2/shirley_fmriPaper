clear all
close all
clc

root = 'C:\Users\User\Documents\fMRI_EXP\Alon\';

cells_data = fullfile(root,'from_git','cells','grid_and_place_cells','smthRm_grid_all_animals.mat');

load(cells_data)