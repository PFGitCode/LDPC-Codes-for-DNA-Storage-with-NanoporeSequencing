clear
addpath(genpath('Function'),genpath('ChannelData'),genpath('ParityMatrix'));
matrix_file = "rate36Sep.mat"; %quater36rate, rate36Sep
channelError = 0.70;
alpha = -1; %for average error model parameter, -1 means using the original channel
bitDis = [0.25,0.25,0.25,0.25]; % distribution for A,T,C,G
testNum = 10000;
decodeMethod = 'method2'; %baseline, method1, method2, method3, method4, quater
hardOrSoft = 'hard'; % hard, soft

fprintf("%s, %s, channelError: %d, testNum: %d\n",decodeMethod,hardOrSoft, channelError, testNum);
 
[avgError,avgTime] = simulation(matrix_file, channelError, bitDis, testNum, decodeMethod, hardOrSoft,alpha);
if alpha > 0
    fileName = "Results/"+decodeMethod+"_"+hardOrSoft + "_"+ "alpha" + num2str(alpha) + "_"+num2str(channelError*100)+".mat";
else
    fileName = "Results/"+decodeMethod+"_"+hardOrSoft + "_"+num2str(channelError*100)+".mat";
end
save(fileName, 'avgError', 'avgTime');
