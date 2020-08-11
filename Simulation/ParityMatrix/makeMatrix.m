addpath('/home/peng/research/LDPC2/LDPC');
clear
addpath distribution
load dis_rate0875.mat;

sumvard = 0;
sumchkd = 0;
for i = 1: numel(vard)
    sumvard = sumvard+ vard(i)/i;
end
for i = 1: numel(chkd)
    sumchkd = sumchkd+ chkd(i)/i;
end
Rate = 1 - sumchkd/sumvard;
% vard = [0,0,1];
% chkd = [0,0,0,0,0,1];
% [vard,chkd] = edge2Node(vard,chkd); % transfer edge degree to node degree
rate = 0.8754;
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
seed = 1;
% rand('state',seed);
% for i = 1 : 1: 1
[H1,exactRate,nk,ERROR_FLAG] = hgen(Rate,n,L,Deg_Lambda,R,Deg_Rho,0,DEBUG_FLAG,FileName);
% [H1, rate]= createIrH(L,R,k);
% H1 = removeCircle(H1);
[b,rearranged_cols]=rearrange_cols(H1);
H1(:,1:k) = b(:,k+1:n);
H1(:,k+1:n) = b(:,1:k);
hEnc = comm.LDPCEncoder(H1);

[H2, rate]= createIrH(L,R,k);
% [H2,exactRate,nk,ERROR_FLAG] = hgen(Rate,n,L,Deg_Lambda,R,Deg_Rho,0,DEBUG_FLAG,FileName);
% H2 = removeCircle(H2);
% [b,rearranged_cols]=rearrange_cols(H2);
% H2(:,1:k) = b(:,k+1:n);
% H2(:,k+1:n) = b(:,1:k);
hEnc = comm.LDPCEncoder(H2);
%         H = removeCircle(H);
%         for i = 1:n
%             if sum(H(:,i)) == 0
%                 H(randperm(k,1),i) = 1;
%             end
%         end
%     [G,H,succ] = ((systematic(full(H))));
%     if(succ == 1)

%         break;
%     end
%     seed = seed + 1;
% end
% H = sparse(H);
save 36matrix2k1 H1 H2