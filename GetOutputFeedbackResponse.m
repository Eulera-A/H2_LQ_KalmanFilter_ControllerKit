function [y,t,z] = GetOutputFeedbackResponse(A,B1,B2,C1,C2,dn,D21,KH2,FH2,time)
%% this function computes the closed-loop feedback control plant response with controller K and estimator F
%% October 21 2021
%% By Owner of Eulera-A github
%% input:
%% A: NxN system dynamic matrix
%% B1: N x num_disturbances disturbance matrix
%% B2: N x num_controllers control B
%% C1: 1 X M selected controlled states
%% C2: 1 x q measurement/sensor states
%% dn: scalar or 1 x M controll noise
%% D21 : 1 x q sensor/measurement noise
%% KH2 : num_controllers x N control gain 
%% FH2 : q x 1 Kalman/estimator gain 
%% time: vector of discrete simulated time stamps, i.e. 0:0.1:100


close all
clc
clear all

%% obtaining dimensions:
num_dist = dim(B1,2);
dim_zp = dim(A,1);
row_B1 = dim(B1,1);

%% setting up the closed-loop state-space:


sys_A = [A           -B2*KH2;
         FH2*C2 (A-B2*KH2-FH2*C2)];
 

sys_B = [B1 zeros(row_B1,dim(dn,2));
         zeros(row_B1,num_dist) FH2*dn];
  

sys_C2 = [C2 zeros(row_C2,col_C2);
          zeros(row_C2,col_C2) C2];
      
sys_D = [D21;
         zeros(size(D21))]; 

%% obtaining state space:
sys_ofb = ss(sys_A,sys_B,sys_C2,sys_D);   

      
nt=length(time);

% assuming each noise dim 1 here
v1=randn(nt,1)*1;
v2=randn(nt,1)*sqrt(1);
//v3 = randn(nt,1)*sqrt(1);


%% specifiy dimensionalities of the plant to be controlled:
%% dim_weights = dim_W1+dim_W_dist+dim_eta;
dim_plant = dim_zp;

%% setting up your initial condition. default is at rest:
z0=[zeros(dim_plant,1);zeros(dim_plant,1)];   

% solving ode using lsim
[y,t,z]=lsim(sys_ofb,[v1 v2],time,z0);



end



