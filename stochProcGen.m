% @fileName: stochProcGen.m
% @purpose: generator of parallel stochastic processes (stochProc), 
%                dim 1: observations, dim 2: channels.
%
% @syntax:
%           [stochProc] = stochProcGen(pdf, meanNull, meanPref, 
%                          stdNull, stdPref, iPrefCh, noObs, noCh)
%
% where:
%
% pdf: string for probability density function among... 
%............n: Gaussian
%............w: inverse Gaussian (Wald)
%............l:   lognormal
%............g:  gamma
%............i:   inverse-gamma/scaled-inverse-chi-square
%............e:  exponential... only the mean will be used!
% 
% meanNull: mean of null channels
% meanPref: mean of preferred channel
% stdNull: standard deviation of null channels
% stdPref: standard deviation of preferred channel
% iPrefCh: index of the preferred channel
% noObs: number of observations
% noCh: number of channels
%  
% @version: 0.1
% @created:  12-Jan-2016       
% @lastModified: 19-Apr-2018   
% @bugs: no known bugs  
% @author: Javier A. Caballero (jcaballeroo)
% @lastEditor: Javier A. Caballero (jcaballeroo)
% 
% @license: MIT license
%  
%  Copyright (c) 2018, Javier A. Caballero
%  
%  Permission is hereby granted, free of charge, to any person obtaining a copy
%  of this software and associated documentation files (the "Software"), to deal
%  in the Software without restriction, including without limitation the rights
%  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%  copies of the Software, and to permit persons to whom the Software is
%  furnished to do so, subject to the following conditions:
%  
%  - The above copyright notice and this permission notice shall be included in all
%  copies or substantial portions of the Software.
% 
%  - Reference to the original source should be given.
%  
%  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%  SOFTWARE.

function [stochProc] = stochProcGen(pdf, meanNull, meanPref, stdNull, stdPref,...
    iPrefCh, noObs, noCh)

switch pdf
    case 'n'% Gaussian
        stochProc = random('norm', meanNull, stdNull, noObs, noCh);
        stochProc(:, iPrefCh) = random('norm', meanPref, stdPref, noObs, 1);
        
    case 'w'% inverse gaussian
        lambdaNull = meanNull^3/stdNull^2;
        lambdaPref = meanPref^3/stdPref^2;
        stochProc = randInvGauss(meanNull, lambdaNull, noObs, noCh);
        stochProc(:, iPrefCh)= randInvGauss(meanPref, lambdaPref, noObs, 1);
        
    case 'l'% lognormal
        stochProc = lognrnd( (log(meanNull) - ((log((stdNull^2/meanNull^2) + 1))/2)),...
            sqrt(log((stdNull^2/meanNull^2) + 1)), noObs, noCh);
        stochProc(:, iPrefCh) = lognrnd((log(meanPref) -...
            (log((stdPref^2/meanPref^2) + 1)/2)),...
            sqrt(log((stdPref^2/meanPref^2) + 1)), noObs, 1);
        
    case 'g'% gamma
        stochProc = random('gam', meanNull^2/stdNull^2, stdNull^2/meanNull, noObs,...
            noCh);
        stochProc(:, iPrefCh) = random('gam', meanPref^2/stdPref^2, stdPref^2/meanPref,...
            noObs, 1);
        
    case 'i'% inverse-gamma
        stochProc = random('gam', (meanNull/stdNull)^2 + 2,...
            1/((meanNull^3/stdNull^2) + meanNull), noObs, noCh);
        stochProc(:, iPrefCh) = random('gam', (meanPref/stdPref)^2 + 2,...
            1/((meanPref^3/stdPref^2) + meanPref), noObs, 1);        
        stochProc = 1./stochProc;
        
    case 'e'% exponential, using only the mean!
        stochProc = random('exp', meanNull, noObs, noCh);
        stochProc(:, iPrefCh) = random('exp', meanPref, noObs, 1);
        
    otherwise% lognormal by default...
        stochProc = lognrnd( (log(meanNull) - ((log((stdNull^2/meanNull^2) + 1))/2)),...
            sqrt(log((stdNull^2/meanNull^2) + 1)), noObs, noCh);
        stochProc(:, iPrefCh) = lognrnd((log(meanPref) -...
            (log((stdPref^2/meanPref^2) + 1)/2)),...
            sqrt(log((stdPref^2/meanPref^2) + 1)), noObs, 1);
        warning(['Density "' pdf '" not recognized, output lognormally'...
            ' distributed by default'])
end