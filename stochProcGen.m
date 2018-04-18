function [stochProc] = stochProcGen(pdf, meanNull, meanPref, stdNull, stdPref,...
    iPrefCh, noObs, noCh)
% generator of parallel stochastic processes (stochProc), dim 1: steps,
% dim 2: streams.
%
% syntax:
% [stochProc] = stochProcGen(pdf, meanNull, meanPref, 
%     stdNull, stdPref, iPrefCh, noObs, noCh)
%
% where:
% % pdf: string for probability density function among... 
% %......n, Gaussian
% %......w, inverse Gaussian (Wald)
% %......l, lognormal
% %......g, gamma
% %......i, inverse-gamma
% %......e, exponential... only valid when the mean is equal to the
%           standard deviation!
% 
% meanNull: mean of null channels
% meanPref: mean of preferred channel
% stdNull: standard deviation of null channels
% stdPref: standard deviation of preferred channel
% iPrefCh: index of the preferred channel
% noObs: number of observations
% noCh: number of channels

switch pdf
    case 'n'% Gaussian
        stochProc = random('norm', meanNull, stdNull, noObs, noCh);
        stochProc(:, iPrefCh) = random('norm', meanPref, stdPref, noObs, 1);
        
    case 'w'% Inverse Gaussian
        lambdaNull = meanNull^3/stdNull^2;
        lambdaPref = meanPref^3/stdPref^2;
        stochProc = randInvGauss(meanNull, lambdaNull, noObs, noCh);
        stochProc(:, iPrefCh)= randInvGauss(meanPref, lambdaPref, noObs, 1);
        
    case 'l'% Lognormal
        stochProc = lognrnd( (log(meanNull) - ((log((stdNull^2/meanNull^2) + 1))/2)),...
            sqrt(log((stdNull^2/meanNull^2) + 1)), noObs, noCh);
        stochProc(:, iPrefCh) = lognrnd((log(meanPref) -...
            (log((stdPref^2/meanPref^2) + 1)/2)),...
            sqrt(log((stdPref^2/meanPref^2) + 1)), noObs, 1);
        
    case 'g'% Gamma
        stochProc = random('gam', meanNull^2/stdNull^2, stdNull^2/meanNull, noObs,...
            noCh);
        stochProc(:, iPrefCh) = random('gam', meanPref^2/stdPref^2, stdPref^2/meanPref,...
            noObs, 1);
        
    case 'i'% Inverse-Gamma
        stochProc = random('gam', (meanNull/stdNull)^2 + 2,...
            1/((meanNull^3/stdNull^2) + meanNull), noObs, noCh);
        stochProc(:, iPrefCh) = random('gam', (meanPref/stdPref)^2 + 2,...
            1/((meanPref^3/stdPref^2) + meanPref), noObs, 1);        
        stochProc = 1./stochProc;
        
    case 'e'% Exponential, using only the mean!
        stochProc = random('exp', meanNull, noObs, noCh);
        stochProc(:, iPrefCh) = random('exp', meanPref, noObs, 1);
        
    otherwise% Lognormal by default...
        stochProc = lognrnd( (log(meanNull) - ((log((stdNull^2/meanNull^2) + 1))/2)),...
            sqrt(log((stdNull^2/meanNull^2) + 1)), noObs, noCh);
        stochProc(:, iPrefCh) = lognrnd((log(meanPref) -...
            (log((stdPref^2/meanPref^2) + 1)/2)),...
            sqrt(log((stdPref^2/meanPref^2) + 1)), noObs, 1);
        warning(['Density "' pdf '" not recognized, output lognormally'...
            ' distributed by default'])
end