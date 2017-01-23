% Causal Inference Measurement
% This example program estimates the directed informations between two
% processes linked by a forward channel and a backward channel. The results
% show that the directed information correctly reflects the directions of
% causal inferences between processes.

close all;
clear all;
clc

%% Setting the parameters
D=3;
Nx=2;
start_show_sample=100;
s_plot='b';
s_plot2='k';
jump=10;
n_data=10^5;
X=-1+zeros(n_data,1)';
Y=X;

%----------------------------------------

alpha=1-0.1;  % transition probability of the clean data
beta=1-0.2;   % transition probability from clean data to observation;

true_value_DI=ctwentropy1D(alpha*(1-beta)+beta*(1-alpha))-ctwentropy1D(alpha);
true_value_inv_DI=ctwentropy1D(alpha*(1-beta)+beta*(1-alpha))-ctwentropy1D(beta);
true_value_MI=2*ctwentropy1D(alpha*(1-beta)+beta*(1-alpha))-ctwentropy1D(alpha)-ctwentropy1D(beta);

axis_cor1=[start_show_sample, n_data, 0, true_value_MI+0.15];

figure;

for i_plots=1:1,

X(1)=(rand(1,1)>0.5);
Y(1)=xor(X(1),rand(1,1)>alpha);
for i=2:n_data
    X(i)=xor(Y(i-1),rand(1,1)>beta);
    Y(i)=xor(X(i),rand(1,1)>alpha);
end;

disp('---simulated data generated.')

[B_MI, B_DI, B_inv_DI]=compute_DI_MI(X,Y,Nx,D,'E1',0,0,0);

disp('---directed information estimator 1 calculated.')

subplot(2,2,1)

semilogx([start_show_sample:jump:length(B_DI)],B_DI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_DI)],s_plot)
hold on
semilogx([start_show_sample,length(B_DI)],true_value_DI*ones(1,2),s_plot2)
xlabel('n')
axis(axis_cor1)
semilogx([start_show_sample:jump:length(B_DI)],B_MI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_MI)],s_plot)
semilogx([start_show_sample,length(B_MI)],true_value_MI*ones(1,2),s_plot2)
semilogx([start_show_sample:jump:length(B_inv_DI)],B_inv_DI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_inv_DI)],s_plot)
semilogx([start_show_sample,length(B_inv_DI)],true_value_inv_DI*ones(1,2),s_plot2)
text(length(B_DI)*1.05,true_value_DI,'DI')
text(length(B_MI)*1.05,true_value_MI,'MI')
text(length(B_inv_DI)*1.05,true_value_inv_DI,'revDI')

title('Estimator 1')

[B_MI, B_DI, B_inv_DI]=compute_DI_MI(X,Y,Nx,D,'E2',0,0,0);

disp('---directed information estimator 2 calculated.')
subplot(2,2,2)

semilogx([start_show_sample:jump:length(B_DI)],B_DI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_DI)],s_plot)
hold on
semilogx([start_show_sample,length(B_DI)],true_value_DI*ones(1,2),s_plot2)
xlabel('n')
axis(axis_cor1)
semilogx([start_show_sample:jump:length(B_DI)],B_MI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_MI)],s_plot)
semilogx([start_show_sample,length(B_MI)],true_value_MI*ones(1,2),s_plot2)
semilogx([start_show_sample:jump:length(B_inv_DI)],B_inv_DI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_inv_DI)],s_plot)
semilogx([start_show_sample,length(B_inv_DI)],true_value_inv_DI*ones(1,2),s_plot2)
text(length(B_DI)*1.05,true_value_DI,'DI')
text(length(B_MI)*1.05,true_value_MI,'MI')
text(length(B_inv_DI)*1.05,true_value_inv_DI,'revDI')

title('Estimator 2')

[B_MI, B_DI, B_inv_DI]=compute_DI_MI(X,Y,Nx,D,'E3',0,0,0);

disp('---directed information estimator 3 calculated.')
subplot(2,2,3)
semilogx([start_show_sample:jump:length(B_DI)],B_DI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_DI)],s_plot)
hold on
semilogx([start_show_sample,length(B_DI)],true_value_DI*ones(1,2),s_plot2)
xlabel('n')
axis(axis_cor1)
semilogx([start_show_sample:jump:length(B_DI)],B_MI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_MI)],s_plot)
semilogx([start_show_sample,length(B_MI)],true_value_MI*ones(1,2),s_plot2)
semilogx([start_show_sample:jump:length(B_inv_DI)],B_inv_DI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_inv_DI)],s_plot)
semilogx([start_show_sample,length(B_inv_DI)],true_value_inv_DI*ones(1,2),s_plot2)
text(length(B_DI)*1.05,true_value_DI,'DI')
text(length(B_MI)*1.05,true_value_MI,'MI')
text(length(B_inv_DI)*1.05,true_value_inv_DI,'invDI')
title('Estimator 3')

[B_MI, B_DI, B_inv_DI]=compute_DI_MI(X,Y,Nx,D,'E4',0,0,0);

disp('---directed information estimator 4 calculated.')

subplot(2,2,4)

semilogx([start_show_sample:jump:length(B_DI)],B_DI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_DI)],s_plot)
hold on
semilogx([start_show_sample,length(B_DI)],true_value_DI*ones(1,2),s_plot2)
xlabel('n')
axis(axis_cor1)
semilogx([start_show_sample:jump:length(B_DI)],B_MI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_MI)],s_plot)
semilogx([start_show_sample,length(B_MI)],true_value_MI*ones(1,2),s_plot2)
semilogx([start_show_sample:jump:length(B_inv_DI)],B_inv_DI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_inv_DI)],s_plot)
semilogx([start_show_sample,length(B_inv_DI)],true_value_inv_DI*ones(1,2),s_plot2)
text(length(B_DI)*1.05,true_value_DI,'DI')
text(length(B_MI)*1.05,true_value_MI,'MI')
text(length(B_inv_DI)*1.05,true_value_inv_DI,'revDI')
title('Estimator 4')
end; 



