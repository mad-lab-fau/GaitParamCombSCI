function [beta,lambda]= LDA2d(data1,data2)
% LDA between two classes
% input: the two classes
%
% example:
% load('IRIS.DAT');
% input = IRIS/10;
% data1 = input(1:50,1:2);
% data2 = input(51:100,1:2);
% [beta,lambda]= LDA2c(data1,data2)

data = [data1;data2];
n1 = size(data1,1);
n2 = size(data2,1);
n=n1+n2;

% Within-class scatter:
data1mean = data1 - repmat(mean(data1),n1,1);
data2mean = data2 - repmat(mean(data2),n2,1);
S1 = data1mean'*data1mean;
S2 = data2mean'*data2mean;
Sw = S1 + S2;

% Between-class scatter:
Sb1 = n1*(mean(data1) - mean(data))'*(mean(data1) - mean(data));
Sb2 = n2*(mean(data2) - mean(data))'*(mean(data2) - mean(data));
Sb = Sb1 + Sb2;

%Resolving the eigenvalue problem
S = inv(Sw)*Sb;
[beta,lambda] = eig(S);
