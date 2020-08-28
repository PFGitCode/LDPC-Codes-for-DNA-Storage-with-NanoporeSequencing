addpath('/home/peng/research/LDPC2/LDPC');
clear
addpath distribution
% load dis_rate075.mat;
vard = [0,0,1];
chkd = [0,0,0,0,0,0,0,0,0,1];
sumvard = 0;
sumchkd = 0;
for i = 1: numel(vard)
    sumvard = sumvard+ vard(i)/i;
end
for i = 1: numel(chkd)
    sumchkd = sumchkd+ chkd(i)/i;
end
Rate = 1 - sumchkd/sumvard;

% [vard,chkd] = edge2Node(vard,chkd); % transfer edge degree to node degree
% rate = ;
% L = [1];
% R = [1];
% Rate = 0.5;
% k = 500;
% n = round(k/(1-Rate));
n = 1000;
L = vard(vard~=0);
R = chkd(chkd~=0);
% L = vard;
% R = chkd;
Deg_Lambda = find(vard ~= 0);
Deg_Rho = find(chkd ~= 0);
DEBUG_FLAG = 0;
FileName = [];

[H1,exactRate,nk,ERROR_FLAG] = hgen(Rate,n,L,Deg_Lambda,R,Deg_Rho,0,DEBUG_FLAG,FileName);
[b,rearranged_cols]=rearrange_cols(H1);
k = size(H1,1);
strlen = n-k;
H1(:,1:strlen) = b(:,k+1:n);
H1(:,strlen+1:n) = b(:,1:k);
hEnc = comm.LDPCEncoder(H1);

[H2,exactRate,nk,ERROR_FLAG] = hgen(Rate,n,L,Deg_Lambda,R,Deg_Rho,1,DEBUG_FLAG,FileName);
[b,rearranged_cols]=rearrange_cols(H2);
k = size(H2,1);
strlen = n-k;
H2(:,1:strlen) = b(:,k+1:n);
H2(:,strlen+1:n) = b(:,1:k);
hEnc = comm.LDPCEncoder(H2);
save 310Sep H1 H2