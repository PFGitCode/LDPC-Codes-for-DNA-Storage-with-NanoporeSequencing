function avgError = simulation_joint(matrix_file,channelError,bitDis,testNum,decodeMethod,hardOrSoft, modelDimen)
load(matrix_file, 'H');  % get genate matrix G and parity matrix H
H = sparse(H);
hEnc = comm.LDPCEncoder(H);
H = full(H);
errorModel = "data/testData" + num2str(channelError)+".mat";
load(errorModel,'mu','matrix');
Pac = 0;
Pat = 0;
Pag = 0;
Pct = 0;
Pcg = 0;
Ptg = 0;

n = size(H,2);
k = size(H,1);
strlen = (n-k)/2;

sumError = 0;
% Date = datestr(today('datetime'));
% fid = fopen("result/"+decodeMethod+"_"+hardOrSoft+modelDimen+num2str(errorRate)+"_"+Date+".txt",'w');
runningTime = 0;
if hardOrSoft == "hard"
    [Pac,Pat,Pag,Pct,Pcg,Ptg] = hardErrorFun(mu,matrix);
    P = [Pac,Pat,Pag,Pcg,Pct,Ptg];
end
% q1 = 0.37;
% q2 = 0.19;
% q3 = 0.19;
% q4 = 0.25;
for iter = 1:testNum
    rawData = (randsrc((n-k)/2,1,[0,1,2,3;bitDis]))';% 0-A,1-C,2-T,3-G
    binaryData = quaternary2bi_joint(rawData);
    biEncodedData = step(hEnc, binaryData');
    encodedData = bi2quaternary_joint(biEncodedData');
    if hardOrSoft == "hard"
        receivedSignal = hardChannel_Gau(double(encodedData), Pac, Pat,Pag,Pct,Pcg,Ptg);
        llr = hardLLR_Gau(receivedSignal, Pac, Pat,Pag,Pct,Pcg,Ptg);
    elseif hardOrSoft == "soft"
        [receivedSignalTime, receivedSignalCurrent] = GaussianChannel2D(double(encodedData),mu,matrix);
        receivedSignal = [receivedSignalTime; receivedSignalCurrent];
        llr = GaussianLLR2D(receivedSignalTime,receivedSignalCurrent,mu,matrix);
    else
        errorRate("Soft or Hard coding?");
    end
    tic
    if decodeMethod == "baseLine"
        hDec = comm.LDPCDecoder(sparse(H));
        bireceivedBits   = step(hDec, llr);
    elseif decodeMethod == "Method1"&& hardOrSoft == "hard"
        bireceivedBits   = method1DecoderHard(llr,H,receivedSignal);
    elseif decodeMethod == "Method1"&& hardOrSoft == "soft"
        bireceivedBits   = method1DecoderSoft2D(llr,H,receivedSignal,mu,matrix);
        
    elseif decodeMethod == "Method2" && hardOrSoft == "hard"
        bireceivedBits = method2DecoderHard(llr,H,receivedSignal,Pac,Pat,Pag,Pct,Pcg,Ptg);
        
    elseif decodeMethod == "Method2" && hardOrSoft == "soft" 
        bireceivedBits = method2DecoderSoft2D(llr, H,receivedSignal,mu,matrix);
        
    elseif decodeMethod == "Method3" && hardOrSoft == "hard"
        bireceivedBits = method3DecoderHard(llr,H,receivedSignal,Pac,Pat,Pag,Pct,Pcg,Ptg);
        
    elseif decodeMethod == "Method3" && hardOrSoft == "soft"
        bireceivedBits = method3DecoderSoft2D(llr, H,receivedSignal,mu,matrix);
        
    elseif decodeMethod == "method4" && hardOrSoft == "hard"
        bireceivedBits = JointHard(llr, H,receivedSignal,P);     
    elseif decodeMethod == "method4" && hardOrSoft == "soft"
        bireceivedBits = method4_soft_joint(llr, H,receivedSignal,mu,matrix);
        
    end
    runningTime = runningTime +toc;
    receivedBits = bi2quaternary_joint(bireceivedBits(1:2*strlen));
    error = sum(double(rawData) - double(receivedBits)~=0)/strlen;
    sumError = sumError+error;
    %     if error > 0
    fprintf("%s, %s: %d\n",decodeMethod,hardOrSoft,iter);
    avgError = sumError/iter;
    avgTime = runningTime/iter;
    %         fprintf(fid, '%f %f %d\n', [avgError avgTime iter]');
    fprintf("error rate: %f,%f, %f\n ",error,avgError, avgTime);
    %     end
end

avgError = sumError/testNum;
fprintf("%s, %s\n",decodeMethod,hardOrSoft);
fprintf("Avrage error rate:%f\n",avgError);

fprintf("%s,%s: %d\n",decodeMethod,hardOrSoft,testNum);
avgError = sumError/testNum;
avgTime = runningTime/testNum;
fprintf(fid, '%f %f %d\n', [avgError avgTime testNum]');
fprintf("%f\n",avgError);
end

