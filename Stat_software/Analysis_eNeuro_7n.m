% This code combines gait parameter based on the calculated weight.
% Then, the results are analyzed.
%
% The parameter included in the combination:
% F_Swing_Time_(s)
% F_Stride_Length_(cm)
% F_Duty_Cycle_(%)
% BOS_Hind_(cm)
% Regularity_Index
% Body_speed
% Seq_AB_(%)
% F_MaxContactAt_(%)
% H_Stride_Length_(cm)
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
file_hasil = [pwd,'\Result_eNeuro_10n.xlsx'];

% Load the parameter combination
file_beta = [pwd,'\Result_parameter_combination_2n.xlsx'];
beta = xlsread(file_beta,1,'C2:C10');

% Study 1:
    File_B = [pwd,'\Study1_average.xlsx'];
    data_B = xlsread(File_B,1,'H2:LN39');
    [n1,group_B] = xlsread(File_B,1,'C2:C39');
    
    data_B_kanan_kiri = 0.5*(data_B(:,[11:54,59:102]) + data_B(:,[107:150,155:198]));
    body_speed_B = 0.25*(data_B(:,[55:58])+data_B(:,[103:106])+data_B(:,[151:154])+data_B(:,[199:202]));
    print_B = 0.5*(data_B(:,221) + data_B(:,222));
    baca_B_1 = [data_B(:,[2:4]),body_speed_B,data_B_kanan_kiri,data_B(:,[203:220]),print_B,data_B(:,[223:319])];
    baca_B_bbb = data_B(:,1);
    
    no_parameter = [36,40,44,105,103,4,100,12,84]; %Single parameters
    baca_B_2 = baca_B_1(:,no_parameter);
    
    LDA_P_B = baca_B_2*beta(:,1); % Combine gait parameters
    
    % Stat Test:
    baca_B = [LDA_P_B,baca_B_2,baca_B_bbb];
    baca_B_I = baca_B(1:6,:);
    baca_B_E2 = baca_B(7:15,:);    % Epo Severe
    baca_B_E1 = baca_B(16:23,:);   % Epo Mild
    baca_B_C2 = baca_B(24:30,:);   % Control_Severe
    baca_B_C1 = baca_B(31:38,:);   % Control_Mild
    group_B_I = group_B(1:6,:);
    group_B_E2 = group_B(7:15,:);    % Epo Severe
    group_B_E1 = group_B(16:23,:);   % Epo Mild
    group_B_C2 = group_B(24:30,:);   % Control_Severe
    group_B_C1 = group_B(31:38,:);   % Control_Mild
    
    % Study 1: 3 class: I, C1, C2
    baca_B_2 = [baca_B_I;baca_B_C1;baca_B_C2];
    group_B_2 = [group_B_I;group_B_C1;group_B_C2];
    % Study 2: C1 E1
    baca_B_3 = [baca_B_C1;baca_B_E1];
    group_B_3 = [group_B_C1;group_B_E1];
    % Study 2: C2 E2
    baca_B_4 = [baca_B_C2;baca_B_E2];
    group_B_4 = [group_B_C2;group_B_E2];    
    
for k = 1:11, % parameter
    stat_B_I(k,:) = [nanmean(baca_B_I(:,k));nanstd(baca_B_I(:,k));nanstd(baca_B_I(:,k))./sqrt(6)];
    stat_B_C1(k,:) = [nanmean(baca_B_C1(:,k));nanstd(baca_B_C1(:,k));nanstd(baca_B_C1(:,k))./sqrt(8)];
    stat_B_C2(k,:) = [nanmean(baca_B_C2(:,k));nanstd(baca_B_C2(:,k));nanstd(baca_B_C2(:,k))./sqrt(7)];
    stat_B_E1(k,:) = [nanmean(baca_B_E1(:,k));nanstd(baca_B_E1(:,k));nanstd(baca_B_E1(:,k))./sqrt(8)];
    stat_B_E2(k,:) = [nanmean(baca_B_E2(:,k));nanstd(baca_B_E2(:,k));nanstd(baca_B_E2(:,k))./sqrt(9)];    
   
    [p2,t2,s2] = anova1(baca_B_2(:,k),group_B_2,'off');
    [c2,m2,h2,nms2] = multcompare(s2,'CType','bonferroni');
    simpan_B_2(k,:) = [p2;c2(:,6)];
   
    [p3,t3,s3] = anova1(baca_B_3(:,k),group_B_3,'off');
    [c3,m3,h3,nms3] = multcompare(s3,'CType','bonferroni');
    simpan_B_3(k,:) = [p3;c3(:,6)];
    [h,pTTest_3,ci,stats] = ttest2(baca_B_C1,baca_B_E1,'Vartype','unequal');
    
    [p4,t4,s4] = anova1(baca_B_4(:,k),group_B_4,'off');
    [c4,m4,h4,nms4] = multcompare(s4,'CType','bonferroni');
    simpan_B_4(k,:) = [p4;c4(:,6)];
    [h,pTTest_4(:,k),ci,stats] = ttest2(baca_B_C2(:,k),baca_B_E2(:,k),'Vartype','unequal');
end
    % Save the results:    
    judul1 = {'LDA','F_swing_time_(s)','F_stride_length_(cm)','F_duty_cycle_(%)',...
    'H_base_of_support_(cm)','Regularity_index_(%)','Body_speed_(cm/s)',...
    'AB_sequence_(%)','F_max_contact_at_(%)','H_stride_length_(cm)','BBB'};
    
    % Statistic
    xlswrite(file_hasil,{'Study 1';'ave_I';'std_I';'sem_I';...
        '';'ave_C_Mild';'std_C_Mild';'sem_C_Mild';'';'ave_C_Severe';'std_C_Severe';'sem_C_Severe';...
        '';'ave_E_Mild';'std_E_Mild';'sem_E_Mild';'';'ave_E_Severe';'std_E_Severe';'sem_E_Severe'},1,'A1');
    xlswrite(file_hasil,judul1,1,'B1');
    xlswrite(file_hasil,stat_B_I',1,'B2');
    xlswrite(file_hasil,stat_B_C1',1,'B6');
    xlswrite(file_hasil,stat_B_C2',1,'B10');
    xlswrite(file_hasil,stat_B_E1',1,'B14');
    xlswrite(file_hasil,stat_B_E2',1,'B18');

    % 3 group: I, C1, C2
    xlswrite(file_hasil,{'Study 1';'p_3groups';'Intact_Control-Mild';...
        'Intact_Control-Severe';...
        'Control-Mild_Control-Severe';'';'ave_I';'std_I';'sem_I';...
        '';'ave_C_Mild';'std_C_Mild';'sem_C_Mild';...
        '';'ave_C_Severe';'std_C_Severe';'sem_C_Severe'},2,'A1');
    xlswrite(file_hasil,judul1,2,'B1');
    xlswrite(file_hasil,simpan_B_2',2,'B2');
    xlswrite(file_hasil,stat_B_I',2,'B7');
    xlswrite(file_hasil,stat_B_C1',2,'B11');
    xlswrite(file_hasil,stat_B_C2',2,'B15');

    % 2 group: C1 E1
    xlswrite(file_hasil,{'Study 1';'p_2groups';'Epo-Mild_Control-Mild';...
        'pTTest';'';'ave_C_Mild';'std_C_Mild';'sem_C_Mild';...
        '';'ave_E_Mild';'std_E_Mild';'sem_E_Mild'},3,'A1');
    xlswrite(file_hasil,judul1,3,'B1');
    xlswrite(file_hasil,simpan_B_3',3,'B2');
    xlswrite(file_hasil,pTTest_3,3,'B4');
    xlswrite(file_hasil,stat_B_C1',3,'B6');
    xlswrite(file_hasil,stat_B_E1',3,'B10');    
    
   % 2 group: C2 E2
    xlswrite(file_hasil,{'Study 1';'p_2groups';'Epo-Severe_Control-Severe';...
        'pTTest';'';'ave_C_Severe';'std_C_Severe';'sem_C_Severe';...
        '';'ave_E_Severe';'std_E_Severe';'sem_E_Severe'},4,'A1');
    xlswrite(file_hasil,judul1,4,'B1');
    xlswrite(file_hasil,simpan_B_4',4,'B2');
    xlswrite(file_hasil,pTTest_4,4,'B4');
    xlswrite(file_hasil,stat_B_C2',4,'B6');
    xlswrite(file_hasil,stat_B_E2',4,'B10'); 
    
% Study 2:
    File_D = [pwd,'\Study2_average_4.xlsx'];
    [data_D_1,data_D_2] = xlsread(File_D,1,'A2:AD30');
    group_D = data_D_2(:,1);
    animal_D = data_D_2(:,2);
    bbb_D = data_D_1(:,1);

    File_Dp = [pwd,'\Study2_average_Prog_2.xlsx'];
    [data_Dp_1,data_Dp_2] = xlsread(File_Dp,1,'A2:AD29');
    group_Dp = data_Dp_2(:,1);
    animal_Dp = data_Dp_2(:,2); 
    bbb_Dp = data_Dp_1(:,1);
    
    print_D = 0.5*(data_D_1(:,20) + data_D_1(:,21));
    speed_D = 0.5*(data_D_1(:,11)./(data_D_1(:,3)+data_D_1(:,5))+data_D_1(:,12)./(data_D_1(:,4)+data_D_1(:,6)));
    baca_D_1 = [data_D_1(:,2:19),print_D,speed_D,data_D_1(:,22:28)];

    print_Dp = 0.5*(data_Dp_1(:,20) + data_Dp_1(:,21));
    speed_Dp = 0.5*(data_Dp_1(:,11)./(data_Dp_1(:,3)+data_Dp_1(:,5))+data_Dp_1(:,12)./(data_Dp_1(:,4)+data_Dp_1(:,6)));
    baca_Dp_1 = [data_Dp_1(:,2:19),print_Dp,speed_Dp,data_Dp_1(:,22:28)];
    
    no_parameter_D = [4,10,12,17,1,20,18,8,11]; %Single parameters
    baca_LDA_D = baca_D_1(:,no_parameter_D);
    baca_LDA_Dp = baca_Dp_1(:,no_parameter_D);
    
    LDA_P_D = baca_LDA_D*beta(:,1); % Combine gait parameters    
    LDA_P_Dp = baca_LDA_Dp*beta(:,1); 
    
    % Stat Test:
    baca_D = [LDA_P_D,baca_LDA_D,bbb_D];
    baca_D_I3 = baca_D(1:6,:);
    baca_D_I4 = baca_D(7:11,:);
    baca_D_C7 = baca_D(12:15,:);
    baca_D_C30 = baca_D(16:22,:);
    baca_D_C60 = baca_D(23:29,:);

    baca_Dp = [LDA_P_Dp,baca_LDA_Dp,bbb_Dp];
    baca_Dp_7d = baca_Dp(1:7,:);
    baca_Dp_14d = baca_Dp(8:13,:);
    baca_Dp_30d = baca_Dp(14:21,:);
    baca_Dp_60d = baca_Dp(22:28,:);
    
for k = 1:11,
    % stat:
    stat_D_I3(k,:) = [nanmean(baca_D_I3(:,k));nanstd(baca_D_I3(:,k));nanstd(baca_D_I3(:,k))./sqrt(6)];
    stat_D_I4(k,:) = [nanmean(baca_D_I4(:,k));nanstd(baca_D_I4(:,k));nanstd(baca_D_I4(:,k))./sqrt(5)];
    stat_D_C7(k,:) = [nanmean(baca_D_C7(:,k));nanstd(baca_D_C7(:,k));nanstd(baca_D_C7(:,k))./sqrt(4)];
    stat_D_C30(k,:) = [nanmean(baca_D_C30(:,k));nanstd(baca_D_C30(:,k));nanstd(baca_D_C30(:,k))./sqrt(7)];
    stat_D_C60(k,:) = [nanmean(baca_D_C60(:,k));nanstd(baca_D_C60(:,k));nanstd(baca_D_C60(:,k))./sqrt(7)];
    
    stat_Dp_7d(k,:) = [nanmean(baca_Dp_7d(:,k));nanstd(baca_Dp_7d(:,k));nanstd(baca_Dp_7d(:,k))./sqrt(7)];
    stat_Dp_14d(k,:) = [nanmean(baca_Dp_14d(:,k));nanstd(baca_Dp_14d(:,k));nanstd(baca_Dp_14d(:,k))./sqrt(6)];
    stat_Dp_30d(k,:) = [nanmean(baca_Dp_30d(:,k));nanstd(baca_Dp_30d(:,k));nanstd(baca_Dp_30d(:,k))./sqrt(8)];
    stat_Dp_60d(k,:) = [nanmean(baca_Dp_60d(:,k));nanstd(baca_Dp_60d(:,k));nanstd(baca_Dp_60d(:,k))./sqrt(7)];
    
    % T-test:
    [h_7,pTTest_7(:,k),ci_7,stats_7] = ttest2(baca_D_C7(:,k),baca_Dp_7d(:,k),'Vartype','unequal');
    [h_30,pTTest_30(:,k),ci_30,stats_30] = ttest2(baca_D_C30(:,k),baca_Dp_30d(:,k),'Vartype','unequal');
    [h_60,pTTest_60(:,k),ci_60,stats_60] = ttest2(baca_D_C60(:,k),baca_Dp_60d(:,k),'Vartype','unequal');
    
    %Two-way ANOVA:
    data_DD = [baca_D_I4(:,k);baca_D_C7(:,k);baca_D_C30(:,k);baca_D_C60(:,k);baca_D_I4(:,k);baca_Dp_7d(:,k);baca_Dp_30d(:,k);baca_Dp_60d(:,k)];
    Study_DD = [zeros(5+4+7+7,1);ones(5+7+8+7,1)];
    Time_DD = [zeros(5,1);7*ones(4,1);30*ones(7,1);60*ones(7,1);zeros(5,1);7*ones(7,1);30*ones(8,1);60*ones(7,1)];   
    [p_DD(:,k),t_DD,s_DD] = anovan(data_DD,{Study_DD Time_DD},'alpha',0.05,'display','off','model',2,'varnames',{'Study','Time'});
    
    % Anova1/Mult compare per time point:
    [pD_a_7(:,k),tD,statsD_a_7] = anova1([baca_D_C7(:,k);baca_Dp_7d(:,k)],[zeros(4,1);ones(7,1)],'off');
    [pD_a_30(:,k),tD,statsD_a_30] = anova1([baca_D_C30(:,k);baca_Dp_30d(:,k)],[zeros(7,1);ones(8,1)],'off');
    [pD_a_60(:,k),tD,statsD_a_60] = anova1([baca_D_C60(:,k);baca_Dp_60d(:,k)],[zeros(7,1);ones(7,1)],'off');
end
    % Save the results:           
    % Between group
    xlswrite(file_hasil,{'Study 2';'D7-P7';'D30-P30';'D60-P60';...
        '';'ave_I300';'std_I300';'sem_I300';'';'ave_I400';'std_I400';'sem_I400';...
        '';'ave_C7';'std_C7';'sem_C7';'';'ave_C30';'std_C30';'sem_C30';'';...
        'ave_C60';'std_C60';'sem_C60';'';'ave_P7';'std_P7';'sem_P7';'';...
        'ave_P14';'std_P14';'sem_P14';'';'ave_P30';'std_P30';'sem_P30';...
        '';'ave_P60';'std_P60';'sem_P60'},5,'A1');
    xlswrite(file_hasil,judul1,5,'B1');
    xlswrite(file_hasil,pTTest_7,5,'B2');
    xlswrite(file_hasil,pTTest_30,5,'B3');
    xlswrite(file_hasil,pTTest_60,5,'B4');
    xlswrite(file_hasil,stat_D_I3',5,'B6');
    xlswrite(file_hasil,stat_D_I4',5,'B10');
    xlswrite(file_hasil,stat_D_C7',5,'B14');
    xlswrite(file_hasil,stat_D_C30',5,'B18');
    xlswrite(file_hasil,stat_D_C60',5,'B22'); 
    xlswrite(file_hasil,stat_Dp_7d',5,'B26');
    xlswrite(file_hasil,stat_Dp_14d',5,'B30');
    xlswrite(file_hasil,stat_Dp_30d',5,'B34');
    xlswrite(file_hasil,stat_Dp_60d',5,'B38');   
    
    xlswrite(file_hasil,{'Two-way ANOVA'},5,'A42');
    xlswrite(file_hasil,{'Study';'Time';'Interaction';'p_7';'p_30';'p_60'},5,'A43');
    xlswrite(file_hasil,p_DD,5,'B43');
    xlswrite(file_hasil,pD_a_7,5,'B46');
    xlswrite(file_hasil,pD_a_30,5,'B47');
    xlswrite(file_hasil,pD_a_60,5,'B48');
  
% Study 3
    [file_number,file_text] = xlsread([pwd,'\Study3_average_1.xlsx'],1,'A2:AD96');  

    group = file_text(:,1);
    animal = file_number(:,1);
    dpi = file_number(:,3);
    bbb_L = file_number(:,5);
    data = file_number(:,[5:27]);

    nama_group = intersect(group,group);
    nama_dpi = intersect(dpi,dpi);

    jumlah_group = length(nama_group);
    jumlah_dpi = length(nama_dpi);
    
    F_Swing = 0.5*(data(:,6)+data(:,13));
    F_Stride = 0.5*(data(:,7)+data(:,14));
    F_Duty = 0.5*(data(:,8)+data(:,15));
    Speed = 0.25*(data(:,9)+data(:,11)+data(:,16)+data(:,18));
    F_Contact = 0.5*(data(:,5)+data(:,12));
    H_Stride = 0.5*(data(:,10)+data(:,17));

    baca_L_T8 = [F_Swing,F_Stride,F_Duty,data(:,[21,20]),Speed,data(:,19),F_Contact,H_Stride];

    LDA_P_L = baca_L_T8(:,1:9)*beta(:,1); % Combine gait parameters
    
    % Stat Test:
    baca_L = [LDA_P_L,baca_L_T8,bbb_L];
    baca_L_YI = baca_L(56:74,:);
    baca_L_Y15 = baca_L(28:29,:);
    baca_L_Y22 = baca_L(30:35,:);
    baca_L_Y29 = baca_L(36:41,:);
    baca_L_Y36 = baca_L(42:47,:);
    baca_L_Y43 = baca_L(48:54,:);
    baca_L_OI = baca_L(16:27,:);
    baca_L_O15 = baca_L(1:3,:);
    baca_L_O22 = baca_L(4:9,:);
    baca_L_O29 = baca_L(10:15,:);

for k = 1:11,
    stat_L_YI(k,:) = [nanmean(baca_L_YI(:,k));nanstd(baca_L_YI(:,k));nanstd(baca_L_YI(:,k))./sqrt(19)];
    stat_L_Y15(k,:) = [nanmean(baca_L_Y15(:,k));nanstd(baca_L_Y15(:,k));nanstd(baca_L_Y15(:,k))./sqrt(2)];
    stat_L_Y22(k,:) = [nanmean(baca_L_Y22(:,k));nanstd(baca_L_Y22(:,k));nanstd(baca_L_Y22(:,k))./sqrt(6)];
    stat_L_Y29(k,:) = [nanmean(baca_L_Y29(:,k));nanstd(baca_L_Y29(:,k));nanstd(baca_L_Y29(:,k))./sqrt(6)];
    stat_L_Y36(k,:) = [nanmean(baca_L_Y36(:,k));nanstd(baca_L_Y36(:,k));nanstd(baca_L_Y36(:,k))./sqrt(6)];
    stat_L_Y43(k,:) = [nanmean(baca_L_Y43(:,k));nanstd(baca_L_Y43(:,k));nanstd(baca_L_Y43(:,k))./sqrt(7)];
    stat_L_OI(k,:) = [nanmean(baca_L_OI(:,k));nanstd(baca_L_OI(:,k));nanstd(baca_L_OI(:,k))./sqrt(12)];
    stat_L_O15(k,:) = [nanmean(baca_L_O15(:,k));nanstd(baca_L_O15(:,k));nanstd(baca_L_O15(:,k))./sqrt(3)];
    stat_L_O22(k,:) = [nanmean(baca_L_O22(:,k));nanstd(baca_L_O22(:,k));nanstd(baca_L_O22(:,k))./sqrt(6)];
    stat_L_O29(k,:) = [nanmean(baca_L_O29(:,k));nanstd(baca_L_O29(:,k));nanstd(baca_L_O29(:,k))./sqrt(6)];
    
    % T-test:
    [h_LI,pTTest_LI(:,k),ci_LI,stats_LI] = ttest2(baca_L_YI(:,k),baca_L_OI(:,k),'Vartype','unequal');
    [h_L15,pTTest_L15(:,k),ci_L15,stats_L15] = ttest2(baca_L_Y15(:,k),baca_L_O15(:,k),'Vartype','unequal');
    [h_L22,pTTest_L22(:,k),ci_L22,stats_L22] = ttest2(baca_L_Y22(:,k),baca_L_O22(:,k),'Vartype','unequal');
    [h_L29,pTTest_L29(:,k),ci_L29,stats_L29] = ttest2(baca_L_Y29(:,k),baca_L_O29(:,k),'Vartype','unequal');
    
    %Two-way ANOVA:
    data_LL = [baca_L_YI(:,k);baca_L_Y15(:,k);baca_L_Y22(:,k);baca_L_Y29(:,k);baca_L_OI(:,k);baca_L_O15(:,k);baca_L_O22(:,k);baca_L_O29(:,k)];
    Study_LL = [zeros(19+2+6+6,1);ones(12+3+6+6,1)];
    Time_LL = [zeros(19,1);15*ones(2,1);30*ones(6,1);60*ones(6,1);zeros(12,1);15*ones(3,1);30*ones(6,1);60*ones(6,1)];    
    [p_LL(:,k),t_LL,s_LL] = anovan(data_LL,{Study_LL Time_LL},'alpha',0.05,'display','off','model',2,'varnames',{'Study','Time'});
    
    % Anova1/Mult compare per time point:
    [pL_a_0(:,k),tL,statsL_a_0] = anova1([baca_L_YI(:,k);baca_L_OI(:,k)],[zeros(19,1);ones(12,1)],'off');
    [pL_a_15(:,k),tL,statsL_a_15] = anova1([baca_L_Y15(:,k);baca_L_O15(:,k)],[zeros(2,1);ones(3,1)],'off');  
    [pL_a_22(:,k),tL,statsL_a_22] = anova1([baca_L_Y22(:,k);baca_L_O22(:,k)],[zeros(6,1);ones(6,1)],'off');
    [pL_a_29(:,k),tL,statsL_a_29] = anova1([baca_L_Y29(:,k);baca_L_O29(:,k)],[zeros(6,1);ones(6,1)],'off');
end
    % Save the results: 
    % Between group
    xlswrite(file_hasil,{'Study 3';'p Uninjured';'p 15 dpi';'p 22 dpi';'p 29 dpi';...
        '';'ave_young0';'std_young0';'sem_young0';'';'ave_young15';'std_young15';'sem_young15';...
        '';'ave_young22';'std_young22';'sem_young22';'';'ave_young29';'std_young29';'sem_young29';'';...
        'ave_young36';'std_young36';'sem_young36';'';'ave_young43';'std_young43';'sem_young43';'';...
        'ave_old0';'std_old0';'sem_old0';'';'ave_old15';'std_old15';'sem_old15';...
        '';'ave_old22';'std_old22';'sem_old22';'';'ave_old29';'std_old29';'sem_old29'},6,'A1');
    xlswrite(file_hasil,judul1,6,'B1');
    xlswrite(file_hasil,pTTest_LI,6,'B2');
    xlswrite(file_hasil,pTTest_L15,6,'B3');
    xlswrite(file_hasil,pTTest_L22,6,'B4');
    xlswrite(file_hasil,pTTest_L29,6,'B5');
    xlswrite(file_hasil,stat_L_YI',6,'B7');
    xlswrite(file_hasil,stat_L_Y15',6,'B11');
    xlswrite(file_hasil,stat_L_Y22',6,'B15');
    xlswrite(file_hasil,stat_L_Y29',6,'B19');
    xlswrite(file_hasil,stat_L_Y36',6,'B23'); 
    xlswrite(file_hasil,stat_L_Y43',6,'B27');
    xlswrite(file_hasil,stat_L_OI',6,'B31');
    xlswrite(file_hasil,stat_L_O15',6,'B35');
    xlswrite(file_hasil,stat_L_O22',6,'B39'); 
    xlswrite(file_hasil,stat_L_O29',6,'B43'); 

    xlswrite(file_hasil,{'Two-way ANOVA'},6,'A47');
    xlswrite(file_hasil,{'Study';'Time';'Interaction';'p_0';'p_15';'p_22';'p_29'},6,'A48');
    xlswrite(file_hasil,p_LL,6,'B48');   
    xlswrite(file_hasil,pL_a_0,6,'B51');
    xlswrite(file_hasil,pL_a_15,6,'B52');
    xlswrite(file_hasil,pL_a_22,6,'B53');
    xlswrite(file_hasil,pL_a_29,6,'B54');
    
% Study 4    
    [file_number2,file_text2] = xlsread([pwd,'\Study4_average_week_0_4_8_12.xlsx'],1,'A2:LN60');  

    animal2 = file_text2(:,2);
    dpi2 = file_number2(:,1);
    data2 = file_number2(:,[5:end]);
    BBB2 = file_number2(:,5);

    F_Swing2 = 0.5*(data2(:,39)+data2(:,135));
    F_Stride2 = 0.5*(data2(:,43)+data2(:,139));
    F_Duty2 = 0.5*(data2(:,47)+data2(:,143));
    Speed2 = 0.25*(data2(:,55)+data2(:,103)+data2(:,151)+data2(:,199));
    F_Contact2 = 0.5*(data2(:,15)+data2(:,111));
    H_Stride2 = 0.5*(data2(:,91)+data2(:,187));

    baca2 = [F_Swing2,F_Stride2,F_Duty2,data2(:,[212,210]),Speed2,data2(:,207),F_Contact2,H_Stride2];

    nama_dpi = intersect(dpi2,dpi2);
    jumlah_dpi = length(nama_dpi);
    
    LDA_P_V = baca2(:,1:9)*beta(:,1); % Combine gait parameters
    
    % Stat Test:    
    baca_V = [LDA_P_V,baca2,BBB2];
    baca_V_0 = baca_V(1:19,:);
    baca_V_4 = baca_V(20:35,:);
    baca_V_8 = baca_V(36:47,:);
    baca_V_12 = baca_V(48:59,:);
    
for k = 1:11,
    stat_V_0(k,:) = [nanmean(baca_V_0(:,k));nanstd(baca_V_0(:,k));nanstd(baca_V_0(:,k))./sqrt(19)];
    stat_V_4(k,:) = [nanmean(baca_V_4(:,k));nanstd(baca_V_4(:,k));nanstd(baca_V_4(:,k))./sqrt(16)];
    stat_V_8(k,:) = [nanmean(baca_V_8(:,k));nanstd(baca_V_8(:,k));nanstd(baca_V_8(:,k))./sqrt(12)];
    stat_V_12(k,:) = [nanmean(baca_V_12(:,k));nanstd(baca_V_12(:,k));nanstd(baca_V_12(:,k))./sqrt(12)];
end
    % Save the results:     
    xlswrite(file_hasil,{'Study 4';'ave_0';'std_0';'sem_0';'';'ave_4';'std_4';'sem_4';...
        '';'ave_8';'std_8';'sem_8';'';'ave_12';'std_12';'sem_12'},7,'A1');
    xlswrite(file_hasil,judul1,7,'B1');
    xlswrite(file_hasil,stat_V_0',7,'B2');
    xlswrite(file_hasil,stat_V_4',7,'B6');
    xlswrite(file_hasil,stat_V_8',7,'B10');
    xlswrite(file_hasil,stat_V_12',7,'B14');
    
% Study 2 & 4: 0, 4 and 8 wpi
for k = 1:11,
    % T-test:
    [h_C0,pTTest_C0(:,k),ci_C0,stats_C0] = ttest2(baca_D_I4(:,k),baca_V_0(:,k),'Vartype','unequal');
    [h_C30,pTTest_C30(:,k),ci_C30,stats_C30] = ttest2(baca_D_C30(:,k),baca_V_4(:,k),'Vartype','unequal');
    [h_C60,pTTest_C60(:,k),ci_C60,stats_C60] = ttest2(baca_D_C60(:,k),baca_V_8(:,k),'Vartype','unequal');
    
    %Two-way ANOVA:
    data_DV = [baca_D_I4(:,k);baca_D_C30(:,k);baca_D_C60(:,k);baca_V_0(:,k);baca_V_4(:,k);baca_V_8(:,k)];
    Study_DV = [zeros(5+7+7,1);ones(19+16+12,1)];
    Time_DV = [zeros(5,1);30*ones(7,1);60*ones(7,1);zeros(19,1);30*ones(16,1);60*ones(12,1)];    
    p_DV(:,k) = anovan(data_DV,{Study_DV Time_DV},'alpha',0.05,'display','off','model',2,'varnames',{'Study','Time'});
    
    %Multiple comparison/ANOVA1:
    [pDV_a_0(:,k),tDV,statsDV_a_0] = anova1([baca_D_I4(:,k);baca_V_0(:,k)],[zeros(5,1);ones(19,1)],'off');
    [pDV_a_30(:,k),tDV,statsDV_a_30] = anova1([baca_D_C30(:,k);baca_V_4(:,k)],[zeros(7,1);ones(16,1)],'off');
    [pDV_a_60(:,k),tDV,statsDV_a_60] = anova1([baca_D_C60(:,k);baca_V_8(:,k)],[zeros(7,1);ones(12,1)],'off');
end  
    % Save the results: 
    xlswrite(file_hasil,{'Study 2 & 4';'p_I';'p_30dpi';'p_60dpi';...
        '';'ave_D0';'std_D0';'sem_D0';'';'ave_D30';'std_D30';'sem_D30';'';'ave_D60';'std_D60';'sem_D60';...
        '';'ave_V0';'std_V0';'sem_V0';'';'ave_V30';'std_V30';'sem_V30';'';'ave_V60';'std_V60';'sem_V60'},8,'A1');
    xlswrite(file_hasil,judul1,8,'B1');
    xlswrite(file_hasil,pTTest_C0,8,'B2');
    xlswrite(file_hasil,pTTest_C30,8,'B3');
    xlswrite(file_hasil,pTTest_C60,8,'B4');
    xlswrite(file_hasil,stat_D_I4',8,'B6');
    xlswrite(file_hasil,stat_D_C30',8,'B10');
    xlswrite(file_hasil,stat_D_C60',8,'B14');
    xlswrite(file_hasil,stat_V_0',8,'B18');
    xlswrite(file_hasil,stat_V_4',8,'B22');
    xlswrite(file_hasil,stat_V_8',8,'B26');

    xlswrite(file_hasil,{'Two-way ANOVA'},8,'A30');
    xlswrite(file_hasil,{'Study';'Time';'Interaction';'0';'30';'60'},8,'A31');
    xlswrite(file_hasil,p_DV,8,'B31');
    xlswrite(file_hasil,pDV_a_0,8,'B34');
    xlswrite(file_hasil,pDV_a_30,8,'B35');
    xlswrite(file_hasil,pDV_a_60,8,'B36');
    
    xlswrite(file_hasil,{'Note: BBB can not be compared directly'},8,'A38');
    
% Study 5:
    [file_number,file_text] = xlsread([pwd,'\Study5_average_0_28_49.xlsx'],1,'A2:Z248');  

    group = file_text(:,1);
    animal = file_text(:,2);
    dpi = file_text(:,3);
    data = file_number;

    nama_group = intersect(group,group);
    nama_dpi = intersect(dpi,dpi);

    jumlah_group = length(nama_group);
    jumlah_dpi = length(nama_dpi);
    
    F_Swing = 0.5*(data(:,2)+data(:,9));
    F_Stride = 0.5*(data(:,3)+data(:,10));
    F_Duty = 0.5*(data(:,4)+data(:,11));
    Speed = 0.25*(data(:,5)+data(:,7)+data(:,12)+data(:,14));
    F_Contact = 0.5*(data(:,1)+data(:,8));
    H_Stride = 0.5*(data(:,6)+data(:,13));

    baca_L_C4 = [F_Swing,F_Stride,F_Duty,data(:,[17,16]),Speed,data(:,15),F_Contact,H_Stride];
    
    LDA_P_L_C4 = baca_L_C4(:,1:9)*beta(:,1); % Combine gait parameters
    
    % Stat Test:  
    baca_L_C4 = [LDA_P_L_C4 ,baca_L_C4];
    baca_L_C4_0 = baca_L_C4(1:6,:);
    baca_L_C4_28 = baca_L_C4(7:12,:);
    baca_L_C4_49 = baca_L_C4(13:18,:);
    baca_L_Sham_0 = baca_L_C4(19:27,:);
    baca_L_Sham_28 = baca_L_C4(28:36,:);
    baca_L_Sham_49 = baca_L_C4(37:45,:);
    
    for k = 1:10,
    stat_L_C4_0(k,:) = [nanmean(baca_L_C4_0(:,k));nanstd(baca_L_C4_0(:,k));nanstd(baca_L_C4_0(:,k))./sqrt(6)];
    stat_L_C4_28(k,:) = [nanmean(baca_L_C4_28(:,k));nanstd(baca_L_C4_28(:,k));nanstd(baca_L_C4_28(:,k))./sqrt(6)];
    stat_L_C4_49(k,:) = [nanmean(baca_L_C4_49(:,k));nanstd(baca_L_C4_49(:,k));nanstd(baca_L_C4_49(:,k))./sqrt(6)];
       
    [h_L_C4_0_28,pTTest_L_C4_0_28(:,k),ci_L_C4_0_28,stats_L_C4_0_28] = ttest2(baca_L_C4_0(:,k),baca_L_C4_28(:,k),'Vartype','unequal');
    
    %Study 2,4,5 dpi 30:
    [p245,t245,s245] = anova1([baca_D_C30(:,k);baca_V_4(:,k);baca_L_C4_28(:,k)],...
        [ones(length(baca_D_C30(:,k)),1);2*ones(length(baca_V_4(:,k)),1);3*ones(length(baca_L_C4_28(:,k)),1)],'off');
    [c245,m245,h245,nms245] = multcompare(s245,'CType','bonferroni');
    simpan_L_C245(k,:) = [p245;c245(:,6)]; 
    
    %Study 2,4,5 dpi 0:
    [p0245,t0245,s0245] = anova1([baca_D_I4(:,k);baca_V_0(:,k);baca_L_C4_0(:,k)],...
        [ones(length(baca_D_I4(:,k)),1);2*ones(length(baca_V_0(:,k)),1);3*ones(length(baca_L_C4_0(:,k)),1)],'off');
    [c0245,m0245,h0245,nms0245] = multcompare(s0245,'CType','bonferroni');
    simpan_L_C0245(k,:) = [p0245;c0245(:,6)];   
    
    % Two Way Anova:
    data_DVL = [baca_D_I4(:,k);baca_D_C30(:,k);baca_V_0(:,k);baca_V_4(:,k);baca_L_C4_0(:,k);baca_L_C4_28(:,k)];
    Study_DVL = [zeros(5+7,1);ones(19+16,1);2*ones(6+6,1)];
    Time_DVL = [zeros(5,1);30*ones(7,1);zeros(19,1);30*ones(16,1);zeros(6,1);30*ones(6,1)];    
    p_DVL(:,k) = anovan(data_DV,{Study_DV Time_DV},'alpha',0.05,'display','off','model',2,'varnames',{'Study','Time'});
end
    % Save the results:    
    judul2 = {'LDA','F_swing_time_(s)','F_stride_length_(cm)','F_duty_cycle_(%)',...
    'H_base_of_support_(cm)','Regularity_index_(%)','Body_speed_(cm/s)',...
    'AB_sequence_(%)','F_max_contact_at_(%)','H_stride_length_(cm)'};
    
    xlswrite(file_hasil,{'Study 5: C4';'pTTest_0_28';'';'ave_0';'std_0';'sem_0';'';'ave_28';'std_28';'sem_28';...
        '';'ave_49';'std_49';'sem_49'},9,'A1');
    xlswrite(file_hasil,judul2,9,'B1');
    xlswrite(file_hasil,pTTest_L_C4_0_28,9,'B2');
    xlswrite(file_hasil,stat_L_C4_0',9,'B4');
    xlswrite(file_hasil,stat_L_C4_28',9,'B8');
    xlswrite(file_hasil,stat_L_C4_49',9,'B12');
    
    xlswrite(file_hasil,{'Study 2,4,5, 30dpi';'p_30dpi';'D_V';'D_C4';'V_C4';...
        '';'ave_D0';'std_D0';'sem_D0';'';'ave_D30';'std_D30';'sem_D30';...
        '';'ave_V0';'std_V0';'sem_V0';'';'ave_V30';'std_V30';'sem_V30';...
        '';'ave_L0';'std_L0';'sem_L0';'';'ave_L30';'std_L30';'sem_L30'},10,'A1');
    xlswrite(file_hasil,judul2,10,'B1');
    xlswrite(file_hasil,simpan_L_C245',10,'B2');
    xlswrite(file_hasil,stat_D_I4',10,'B7');
    xlswrite(file_hasil,stat_D_C30',10,'B11');
    xlswrite(file_hasil,stat_V_0',10,'B15');
    xlswrite(file_hasil,stat_V_4',10,'B19');
    xlswrite(file_hasil,stat_L_C4_0',10,'B23');
    xlswrite(file_hasil,stat_L_C4_28',10,'B27');
    
    xlswrite(file_hasil,{'Study 2,4,5, 0dpi';'p_0dpi';'D_V';'D_C4';'V_C4';...
        '';'ave_D0';'std_D0';'sem_D0';'';'ave_D30';'std_D30';'sem_D30';...
        '';'ave_V0';'std_V0';'sem_V0';'';'ave_V30';'std_V30';'sem_V30';...
        '';'ave_L0';'std_L0';'sem_L0';'';'ave_L30';'std_L30';'sem_L30'},11,'A1');
    xlswrite(file_hasil,judul2,11,'B1');
    xlswrite(file_hasil,simpan_L_C0245',11,'B2');
    xlswrite(file_hasil,stat_D_I4',11,'B7');
    xlswrite(file_hasil,stat_D_C30',11,'B11');
    xlswrite(file_hasil,stat_V_0',11,'B15');
    xlswrite(file_hasil,stat_V_4',11,'B19');
    xlswrite(file_hasil,stat_L_C4_0',11,'B23');
    xlswrite(file_hasil,stat_L_C4_28',11,'B27');
    
    xlswrite(file_hasil,{'Two-way ANOVA Study 2,4,5'},12,'A1');
    xlswrite(file_hasil,{'Study';'Time';'Interaction'},12,'A2');
    xlswrite(file_hasil,judul2,12,'B1');
    xlswrite(file_hasil,p_DVL,12,'B2');
    
    xlswrite(file_hasil,{'Study 4 have NaN data, therefore the SEM are:'},7,'A20');
    xlswrite(file_hasil,{'30dpi';'pLDA';'BBB'},7,'A21');
    V_SEM_LDA = nanstd(baca_V_4(:,1))./15;
    V_SEM_BBB = nanstd(baca_V_4(:,11))./13;
    xlswrite(file_hasil,V_SEM_LDA,7,'B22');
    xlswrite(file_hasil,V_SEM_BBB,7,'B23');