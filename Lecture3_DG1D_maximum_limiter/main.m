% main.m
clear;clc

global Nx dimPk NumGLP
Nx = 160;
k = 2;
NumGLP = 5;
dimPk = k + 1;
CFL = 0.05;

get_GLP

init_data

get_basis

L2_Pro

RK3

calculate_L2_Error

plot(Xc,uh(:,1),'b-','linewidth',1.3); axis([-1,1,-0.1,1.1])
