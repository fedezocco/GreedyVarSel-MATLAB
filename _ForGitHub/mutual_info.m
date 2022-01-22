%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

function MI = mutual_info(X,Y)
    X = X - mean(X);
    Y = Y - mean(Y);
    Xtot = [X Y];
    n = size(X,1);
    sigmaX = X'*X/n;
    sigmaY = Y'*Y/n;
    sigma = Xtot'*Xtot/n;
    MI = 1/2*log2(det(sigmaX)*det(sigmaY)/det(sigma));
end

