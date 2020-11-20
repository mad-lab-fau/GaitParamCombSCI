% This code apply LDA on the chosen gait parameters
%
% This code is related to the publication Timotius et.al, 
% "Combination of defined CatWalk gait parameters for predictive locomotion 
% alterations in experimental thoracic spinal cord injured rat models", 2020.
%
% Ivanna K. Timotius (2020)

clear all
close all
clc

% Assign where to save the result:
file_hasil = [pwd,'\Result_parameter_combination_2n.xlsx'];

% Read Study 2 run data:
    [file_number,file_text] = xlsread([pwd,'\CatWalk_Study2_rundata.xlsx'],1,'A2:AE51');  

    group_D = file_text(:,1);
    animal_D = file_text(:,2);
    data_D = file_number;

    print_D = 0.5*(data_D(:,20) + data_D(:,21));
    speed_D = 0.5*(data_D(:,11)./(data_D(:,3)+data_D(:,5))+data_D(:,12)./(data_D(:,4)+data_D(:,6)));
    data_D_1 = [data_D(:,2:19),print_D,speed_D,data_D(:,22:28)];
    
    no_parameter_D = [4,10,12,17,1,20,18,8,11]; %This is the chosen parameter for the LDA
    baca_D = data_D_1(:,no_parameter_D);
    baca_D_I400 = baca_D(1:18,:);    
    baca_D_C60 = baca_D(19:50,:);
    
% Read Study 1 run data:
    File_Br = [pwd,'\CatWalk_Study1_rundata.xlsx'];
    data_Br = xlsread(File_Br,1,'H2:LN58');
    [n1,group_Br] = xlsread(File_Br,1,'C2:C58');

    data_Br_kanan_kiri = 0.5*(data_Br(:,[11:54,59:102]) + data_Br(:,[107:150,155:198]));
    body_speed_Br = 0.25*(data_Br(:,[55:58])+data_Br(:,[103:106])+data_Br(:,[151:154])+data_Br(:,[199:202]));
    print_Br = 0.5*(data_Br(:,221) + data_Br(:,222));
    
    baca_Br_1 = [data_Br(:,[2:4]),body_speed_Br,data_Br_kanan_kiri,data_Br(:,[203:220]),print_Br,data_Br(:,[223:319])];
    no_parameter_B = [36,40,44,105,103,4,100,12,84]; %This is the chosen parameter for the LDA
    baca_Br = baca_Br_1(:,no_parameter_B);
    baca_Br_I = baca_Br(1:27,:);   
    baca_Br_C2 = baca_Br(28:57,:);    

% Combine data    
data_1 = [baca_Br_I;baca_D_I400];
data_2 = [baca_Br_C2;baca_D_C60];
  
% LDA:
[beta,lambda]= LDA2d(data_1,data_2);

xlswrite(file_hasil,{'Study1 param.no', 'Study2 param.no', 'Weight'},1,'A1');
xlswrite(file_hasil,[no_parameter_B',no_parameter_D',beta(:,1)],1,'A2');