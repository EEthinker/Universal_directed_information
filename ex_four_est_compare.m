% Stationary Hidden Markov Processes
% For a statoinary hidden Markov process, this example program estimates the directed information 
% from the process observations to the process states. This program shows that the estimators
% converge to the analytical value asymptotically.


close all;
clear all;
clc

%% Setting parameters
D=3;
Nx=2;
start_show_sample=100;
s_plot='b';
s_plot2='k';
n_data=10^5;
X=-1+zeros(n_data,1)';
Y=X;

%---------------------------------------------
p_t=0.3;  % cross-over probability of the clean data
o_t=0.2;  % cross-over probability from clean data to observation

true_value_DI=ctwentropy1D(p_t)-(((1-p_t)*(1-o_t)+p_t*o_t)*ctwentropy1D(p_t*o_t/((1-p_t)*(1-o_t)+p_t*o_t))+...
    ((p_t)*(1-o_t)+(1-p_t)*o_t)*ctwentropy1D((1-p_t)*o_t/((p_t)*(1-o_t)+(1-p_t)*o_t)));
axis_cor=[start_show_sample, n_data,  true_value_DI-0.2, true_value_DI+0.2];

n_plots=3;

for i_plots=1:n_plots

data1=(rand(1,n_data)<p_t);
data(1)=(rand(1,1)>0.5);
for i=2:n_data
    data(i)=xor(data(i-1),data1(i));
end;
sum(data)/n_data;

Y=data;
X=xor(Y,(rand(1,n_data)<o_t));

figure(1)
jump=10;
disp('---simulated data generated.')


[B_MI, B_DI, B_rev_DI]=compute_DI_MI(X,Y,Nx,D,'E1',0,0,0);
disp('---directed information estimator 1 calculated.')
subplot(2,2,1)
semilogx([start_show_sample:jump:length(B_DI)],B_DI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_DI)],s_plot)
hold on
if i_plots==1
semilogx([start_show_sample,length(B_DI)],true_value_DI*ones(1,2), s_plot2);
end;
title('Estimator 1')
xlabel('n')
axis(axis_cor)

[B_MI, B_DI, B_rev_DI]=compute_DI_MI(X,Y,Nx,D,'E2',0,0,0);
disp('---directed information estimator 2 calculated.')
subplot(2,2,2)
semilogx([start_show_sample:jump:length(B_DI)],B_DI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_DI)],s_plot)
hold on
if i_plots==1
semilogx([start_show_sample,length(B_DI)],true_value_DI*ones(1,2), s_plot2)
end;
title('Estimator 2')
xlabel('n')
axis(axis_cor)

[B_MI, B_DI, B_rev_DI]=compute_DI_MI(X,Y,Nx,D,'E3',0,0,0);
disp('---directed information estimator 3 calculated.')
subplot(2,2,3)
semilogx([start_show_sample:jump:length(B_DI)],B_DI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_DI)],s_plot)
hold on
if i_plots==1
semilogx([start_show_sample,length(B_DI)],true_value_DI*ones(1,2), s_plot2)
end;
title('Estimator 3')
xlabel('n')
axis(axis_cor)
hold on


[B_MI, B_DI, B_rev_DI]=compute_DI_MI(X,Y,Nx,D,'E4',0,0,0);
disp('---directed information estimator 4 calculated.')
subplot(2,2,4)
semilogx([start_show_sample:jump:length(B_DI)],B_DI(start_show_sample:jump:end)./[start_show_sample:jump:length(B_DI)],s_plot)
hold on
if i_plots==1
semilogx([start_show_sample,length(B_DI)],true_value_DI*ones(1,2), s_plot2)
end;
title('Estimator 4')
xlabel('n')
axis(axis_cor)
end; 


