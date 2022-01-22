%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

function [X_dataStore] = DataGenerator(TestSelector, NumOfSimulations)

switch TestSelector
    
    case 1
        %% Constants:
        NumOfObservations = 1000;
        mean0 = 0;
        sigma0 = 1;
        sigmaNoise = 0.1;
        sigmaLargeNoise = 0.4;
        
        for i = 1:NumOfSimulations
            %% Base variables (4 in total):
            w0 = normrnd(mean0*ones(NumOfObservations,1), sigma0*ones(NumOfObservations,1));
            x0 = normrnd(mean0*ones(NumOfObservations,1), sigma0*ones(NumOfObservations,1));
            y0 = normrnd(mean0*ones(NumOfObservations,1), sigma0*ones(NumOfObservations,1));
            z0 = normrnd(mean0*ones(NumOfObservations,1), sigma0*ones(NumOfObservations,1));
            
            %% Noise variables (20 in total):
            e1 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e2 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e3 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e4 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e5 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e6 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e7 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e8 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e9 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e10 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e11 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e12 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e13 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e14 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e15 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e16 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e17 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e18 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e19 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            e20 = normrnd(mean0*ones(NumOfObservations,1), sigmaNoise*ones(NumOfObservations,1));
            
            %% Larger noise variables (2 in total):
            e21 = normrnd(mean0*ones(NumOfObservations,1), sigmaLargeNoise*ones(NumOfObservations,1));
            e22 = normrnd(mean0*ones(NumOfObservations,1), sigmaLargeNoise*ones(NumOfObservations,1));
            
            
            
            %% Elements of X_data (26 in total):
            
            w1 = w0 + e1;
            w2 = w0 + e2;
            w3 = w0 + e3;
            w4 = w0 + e4;
            w5 = w0 + e5;
            
            x1 = x0 + e6;
            x2 = x0 + e7;
            x3 = x0 + e8;
            x4 = x0 + e9;
            x5 = x0 + e10;
            
            y1 = y0 + e11;
            y2 = y0 + e12;
            y3 = y0 + e13;
            y4 = y0 + e14;
            y5 = y0 + e15;
            
            z1 = z0 + e16;
            z2 = z0 + e17;
            z3 = z0 + e18;
            z4 = z0 + e19;
            z5 = z0 + e20;
            
            h1 = w0 + x0 + e21;
            h2 = y0 + z0 + e22;
            
            
            X_dataStore(:,:,i) = [w0 w1 w2 w3 w4 w5 x0 x1 x2 x3 x4 x5 y0 y1 y2 y3 y4 y5 z0 z1 z2 z3 z4 z5 h1 h2];
        end
        
    case 2
        %% Constants:
        NumOfObservations = 1000; % n
        NumOfColumnsOfX0 = 25; % u
        v = 50; % NOTE: Must be > NumOfColumnsOfX0
        mean = 0;
        sigma = 1;
        sigmaNoise = 0.1;
        
        for i = 1:NumOfSimulations
            X0 = normrnd(mean*ones(NumOfObservations,NumOfColumnsOfX0), sigma*ones(NumOfObservations,NumOfColumnsOfX0));
            fi = normrnd(mean*ones(NumOfColumnsOfX0, v - NumOfColumnsOfX0), sigma*ones(NumOfColumnsOfX0, v - NumOfColumnsOfX0));
            e = normrnd(mean*ones(NumOfObservations, v - NumOfColumnsOfX0), sigmaNoise*ones(NumOfObservations, v - NumOfColumnsOfX0));
            X1 = X0 * fi + e;
            
            X_dataStore(:,:,i) = [X0 X1]; % It has NumOfColumnsOfX0+(v-NumOfColumnsOfX0) colums, NumOfObservations rows
        end
        
    case 3
        load('medical.mat');
        X_dataStore = arrhythmia;
        
    case 4
        load('X50sites.mat'); 
        X_dataStore = X50sites;

    case 5
        load('Xpitprops.mat'); 
        X_dataStore = Xpitprops;

    case 6
        load('Xmusic.mat'); 
        X_dataStore = Xmusic; 
  
    case 7
        load('XgasBatch3.mat');
        X_dataStore = XgasBatch3; 
        
    case 8
        load('financial.mat');
        X_dataStore = sales;
        
    case 9 
        load('MofERdata.mat');
        X_dataStore = M;
        
   case 10
        load('wafer.mat');
        X_dataStore = X;
        
   case 11
        load('YaleB.mat');
        X_dataStore = X;
        
   case 12
        load('USPS.mat');
        X_dataStore = X;
        
end
