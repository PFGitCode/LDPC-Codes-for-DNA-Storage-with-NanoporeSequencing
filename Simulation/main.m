clear
addpath(genpath('Function'),genpath('ChannelData'),genpath('ParityMatrix'));
matrix_file = "quater36rate.mat"; %quater36rate, rate36Sep
channelError = 0.65;
bitDis = [0.25,0.25,0.25,0.25]; % distribution for A,T,C,G
testNum = 10000;
decodeMethod = 'quater'; % method1, method2, method3, method4, quater
hardOrSoft = 'hard'; % hard, soft

fprintf("%s, %s, channelError: %d, testNum: %d\n",decodeMethod,hardOrSoft, channelError, testNum);
 
avgError = simulation(matrix_file, channelError, bitDis, testNum, decodeMethod, hardOrSoft);

fileName = "Results/"+decodeMethod+"_"+hardOrSoft + "_"+num2str(channelError*100)+".mat";
save(fileName, 'avgError');
