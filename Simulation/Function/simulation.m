function [avgError,avgTime] = simulation(matrix_file, channelError, bitDis, testNum, decodeMethod, hardOrSoft, alpha)
errorModel = "data/testData" + num2str(channelError)+".mat";
load(errorModel,'mu','matrix');
if hardOrSoft == "hard"
    P = hardError(mu,matrix);% P = [Pac,Pat,Pag,Pct,Pcg,Ptg]
    if alpha > 0
        lamda = sum(P)/(5+alpha);
        P = [lamda,lamda,lamda,alpha*lamda,lamda,lamda];
    end
end


if  decodeMethod == "quater"
    load(matrix_file, 'H', 'g'); % get parity matrix H1, H2
    n = size(H,2);
    k = size(H,1);
    H = full(H);
    
    maxlen = 0;
    H = H.x;
    for i = 1:k
        hi = H(i,:);
        non = find(hi);
        if length(non) > maxlen
            maxlen = length(non);
        end
    end
    for i  = 1:n-k
        hi = H(:,i);
        nonZerosElement = find(hi~=0);
        nonZerosElements(i,1:length(nonZerosElement)) = nonZerosElement;
    end
    
    [nonZerosElements,s]=decInit(H,maxlen);
    
else
    load(matrix_file, 'H1', 'H2'); % get parity matrix H1, H2
    n = size(H1,2);
    k = size(H1,1);
    hEnc1 = comm.LDPCEncoder(H1);
    hEnc2 = comm.LDPCEncoder(H2);
    H1 = full(H1);
    H2 = full(H2);
    
end

strlen = n-k;
runningTime = 0;
sumError = 0;

for iter = 1:testNum
    tic
    rawData = (randsrc(n-k,1,[0,1,2,3;bitDis]))';
    if  decodeMethod == "quater"
        encodedData = LDPC_Encoder(rawData,g);
    else
        [biData1, biData2] = quaternary2bi(rawData);
        encodedData1 = step(hEnc1, biData1');
        encodedData2 = step(hEnc2, biData2');
        encodedData = bi2quaternary(encodedData1,encodedData2);
    end
    
    
    if hardOrSoft == "hard"
        receivedSignal = hardChannel(double(encodedData), P);
        llr = hardLLR(receivedSignal, P);
        
        if decodeMethod == "baseline"
            %         hDec1 = comm.LDPCDecoder(sparse(H1));
            %         hDec2 = comm.LDPCDecoder(sparse(H2));
            %         bireceivedBits1 = step(hDec1, llr(1:n));
            %         bireceivedBits2 = step(hDec2, llr(n+1:end));
            [bireceivedBits1, bireceivedBits2]= baseline(llr(1:n),llr(n+1:end), H1,H2,receivedSignal);
        elseif decodeMethod == "method1"
            [bireceivedBits1, bireceivedBits2]= method1Hard(llr(1:n),llr(n+1:end), H1,H2,receivedSignal);
        elseif decodeMethod == "method2"
            [bireceivedBits1, bireceivedBits2]= method2Hard(llr(1:n),llr(n+1:end), H1,H2,receivedSignal,P);
        elseif decodeMethod == "method3"
            [bireceivedBits1, bireceivedBits2]= method3Hard(llr(1:n),llr(n+1:end), H1,H2,receivedSignal,P);
        elseif  decodeMethod == "method4"
            [bireceivedBits1, bireceivedBits2]= method4Hard(llr(1:n),llr(n+1:end), H1,H2,receivedSignal,P);
        elseif decodeMethod == "quater"
            receivedBits = quaternaryDecoder(receivedSignal,H,P,nonZerosElements,s);
        end
        
    elseif hardOrSoft == "soft"
        [receivedSignalTime, receivedSignalCurrent] = GaussianChannel2D(double(encodedData),mu,matrix);
        receivedSignal = [receivedSignalTime; receivedSignalCurrent];
        llr = GaussianLLR2D(receivedSignalTime,receivedSignalCurrent,mu,matrix);
        if decodeMethod == "baseline"
            %         hDec1 = comm.LDPCDecoder(sparse(H1));
            %         hDec2 = comm.LDPCDecoder(sparse(H2));
            %         bireceivedBits1 = step(hDec1, llr(1:n));
            %         bireceivedBits2 = step(hDec2, llr(n+1:end));
            [bireceivedBits1, bireceivedBits2]= baseline(llr(1:n),llr(n+1:end), H1,H2,receivedSignal);
        elseif decodeMethod == "method1"
            [bireceivedBits1, bireceivedBits2]= method1Hard(llr(1:n),llr(n+1:end), H1,H2,receivedSignal);
        elseif decodeMethod == "method2"
            [bireceivedBits1, bireceivedBits2]= method2Hard(llr(1:n),llr(n+1:end), H1,H2,receivedSignal,P);
        elseif decodeMethod == "method3"
            [bireceivedBits1, bireceivedBits2]= method3Soft(llr(1:n),llr(n+1:end), H1,H2,receivedSignal,mu,matrix);
        elseif  decodeMethod == "method4"
            [bireceivedBits1, bireceivedBits2]= method4Soft(llr(1:n),llr(n+1:end), H1,H2,receivedSignal,mu,matrix);
        elseif decodeMethod == "quater"
            receivedBits = quaternaryDecoderSoft(receivedSignal,H,mu,matrix,nonZerosElements,s);
        end
    end

    
    if decodeMethod ~= "quater"
        receivedBits = bi2quaternary(bireceivedBits1,bireceivedBits2);
    end
    
    runningTime = runningTime +toc;
    error = sum(double(rawData) - double(receivedBits(1:strlen))~=0)/strlen;
    sumError = sumError+error;
    %     if error > 0
    fprintf("%s, %s: %d\n",decodeMethod,hardOrSoft,iter);
    avgError = sumError/iter;
    avgTime = runningTime/iter;
    fprintf("error rate: %f,%f, %f\n",channelError,avgError, avgTime);
    %     end
end
avgError = sumError/testNum;
avgTime = runningTime/testNum;
end


