% main.m
%clear;clc

global Nx dimPk NumGLP NumEq flux_type
%Nx = 601;
k = 2;
NumGLP = 5;
dimPk = k + 1;
NumEq = 3;
CFL = 0.15;
%flux_type = 1;
%Limit_type = 2;

get_GLP

init_data

get_basis
 
L2_Pro
 
RK3

calculate_L2_Error

%draw_solution
