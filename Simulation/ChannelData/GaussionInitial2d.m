clear
for r = 1.01:0.01:1.1
    muA = [0.5,0.62];
    muC = [0.06,0.31];
    muT = [0.09,0.49];
    muG = [0.15,0.83];
    mu = [muA,muC,muT,muG];
    
    thetaA = [0.1,0.1]*r;
    thetaT = [0.2,0.2]*r;
    thetaC = [0.15,0.14]*r;
    thetaG = [0.1,0.15]*r;
    theta = [thetaA, thetaC, thetaT, thetaG];
    
    angleA = -pi/6;
    angleT = 0;
    angleC = 0;
    angleG = -pi/3;
    
    mA = generateMatrix(angleA,thetaA);
    mT = generateMatrix(angleT,thetaT);
    mC = generateMatrix(angleC,thetaC);
    mG = generateMatrix(angleG,thetaG);
    matrix = [mA,mC,mT,mG];
    
    filename = "data/testData" + num2str(r)+".mat";
    save(filename,'mu','theta','matrix');
end

function m = generateMatrix(sigma,theta)
theta_X = theta(1);
theta_Y = theta(2);
s = [theta_X^2, 0; 0,theta_Y^2];
B = [cos(sigma), sin(sigma);-sin(sigma), cos(sigma)];
m = B'*s*B;
end
