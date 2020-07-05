function avgError = simulation(matrix_file, channelError, bitDis, testNum, decodeMethod, hardOrSoft);
load(matrix_file, 'H1', 'H2'); % get parity matrix H1, H2
errorModel = "data/testData" + num2str(channelError)+".mat";
load(errorModel,'mu','theta','matrix');
if hardOrSoft == "hard"
    P = hardError(mu,matrix); % P = [Pac,Pat,Pag,Pct,Pcg,Ptg] 
end
n = size(H1,2);
k = size(H1,1);

strlen = n-k;
hEnc1 = comm.LDPCEncoder(H1);
hEnc2 = comm.LDPCEncoder(H2);
H1 = full(H1);
H2 = full(H2);

runningTime = 0;
sumError = 0;
for iter = 1:testNum
    tic
    rawData = (randsrc(n-k,1,[0,1,2,3;bitDis]))';
    [biData1, biData2] = quaternary2bi(rawData);
    encodedData1 = step(hEnc1, biData1');
    encodedData2 = step(hEnc2, biData2');
    encodedData = bi2quaternary(encodedData1,encodedData2);
    
    if hardOrSoft == "hard"
        receivedSignal = hardChannel(double(encodedData), P);
        llr = hardLLR(receivedSignal, P);
    elseif hardOrSoft == "soft" && modelDimen == "2D"
        [receivedSignalTime, receivedSignalCurrent] = GaussianChannel2D(double(encodedData),mu,matrix);
        receivedSignal = [receivedSignalTime; receivedSignalCurrent];
        llr = GaussianLLR2D(receivedSignalTime,receivedSignalCurrent,mu,matrix);
    end
    if  decodeMethod == "method4"
        [bireceivedBits1, bireceivedBits2]= Method4Hard(llr(1:n),llr(n+1:end), H1,H2,receivedSignal,P);
    end
    
    runningTime = runningTime +toc;
    receivedBits = bi2quaternary(bireceivedBits1,bireceivedBits2);
    error = sum(double(rawData) - double(receivedBits(1:strlen))~=0)/strlen;
    sumError = sumError+error;
    
    if error > 0
        fprintf("%s, %s: %d\n",decodeMethod,hardOrSoft,iter);
        avgError = sumError/iter;
        avgTime = runningTime/iter;
        fprintf("error rate: %f,%f, %f\n",channelError,avgError, avgTime);
    end
end
avgError = sumError/testNum;
end


