clear
addpath(genpath('Function'),genpath('ChannelData'),genpath('ParityMatrix'));
matrix_file = "rate36Sep2.mat";
channelError = 0.65;
bitDis = [0.25,0.25,0.25,0.25]; % distribution for A,T,C,G
testNum = 1;
decodeMethod = 'method4'; % method1, method2, method3, method4
hardOrSoft = 'hard'; % hard, soft

fprintf("%s, %s, channelError: %d, testNum: %d\n",decodeMethod,hardOrSoft, channelError, testNum);
 
avgError = simulation(matrix_file, channelError, bitDis, testNum, decodeMethod, hardOrSoft);


% fileName = "/home/peng/research/DNA_LDPC/simulation/result/"+decodeMethod+"_"+hardOrSoft + modelDimen + "/"+num2str(error*100)+".mat";
% save(fileName, 'avgError');