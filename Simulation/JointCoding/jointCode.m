clear
matrix_file = "rate36Joint.mat";
error = 0.8;
bitDis = [0.25,0.25,0.25,0.25]; % distribution for A,T,C,G
testNum = 10000;
decodeMethod = 'method4';
% hardOrSoft = 'hard';
hardOrSoft = 'soft';
modelDimen = '2D';

avgError = simulation_joint(matrix_file, error, bitDis, testNum, decodeMethod, hardOrSoft,modelDimen);

fileName = "/home/peng/research/DNA_LDPC/simulation/result/"+decodeMethod+"_"+hardOrSoft + modelDimen + "/"+num2str(error*100)+".mat";
save(fileName, 'avgError');