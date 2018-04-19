% @fileName: paramRepository.m
% @purpose: repository of parameters for simulation of the recursive 
%                multi-hypothesis sequential probability ratio test (rMSPRT), 
%                as introduced in Caballero JA, Humphries MD, 
%                Gurney KN (2018) A probabilistic, distributed, recursive mechanism for 
%                decision-making in the brain. PLoS Comput Biol 14(4): e1006033. 
%                https://doi.org/10.1371/journal.pcbi.1006033. 
%
% @syntax:
%       [meansNullCh, meansPrefCh, stdsNullCh, stdsPrefCh] = paramRepository(paramCase)
%
% where paramCase will retrieve paramters...
% %......'mt': ... computed from MT data
% %......'mtInfoDeplN2': ... after information depletion for N = 2
% %......'mtInfoDeplN4': ... after information depletion for N = 4
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

function [meansNullCh, meansPrefCh, stdsNullCh, stdsPrefCh] =...
paramRepository(paramCase)
% parameters for preferred channel (same in all cases)
% means of preferred channel
meansPrefCh = [54.1067812950675;...
  51.9555167200400;46.1073617021512;...
  37.7431198991798;29.8554021786382];
% standard deviations of preferred channel
stdsPrefCh = [33.1325606188147;...
  32.1756309152138;30.4920600270835;...
  28.0096721079346;25.9790487446720];

% choose
switch paramCase
    case 'mt'
    % means of null channel
    meansNullCh = [59.3736631745858;...
        62.9023641601734;65.5267578611166;...
        70.1465464403720;83.4519960352983];
    % standard deviations of null channel
    stdsNullCh = [34.5386065496385;...
        35.3233041337630;36.1434994400871;...
        37.1454736402515;40.6329345670797];     

    case 'mtInfoDeplN2'
    % means of null channel
    meansNullCh = [58.98362;...
        60.59836;60.30201;...
        62.00465;75.51389];
    % standard deviations of null channel
    stdsNullCh = [34.43448;...
        34.66081;34.62299;...
        34.84995;38.46257];
        
    case 'mtInfoDeplN4'
    % means of null channel
    meansNullCh = [58.53778;...
        59.83197;58.31644;...
        59.87656;71.80696];
    % standard deviations of null channel
    stdsNullCh = [34.31546;...
        34.44044;34.04515;...
        34.24996;37.44906];
        
    otherwise% return 'mt' parameters by default
    % means of null channel
    meansNullCh = [59.3736631745858;...
        62.9023641601734;65.5267578611166;...
        70.1465464403720;83.4519960352983];
    % standard deviations of null channel
    stdsNullCh = [34.5386065496385;...
        35.3233041337630;36.1434994400871;...
        37.1454736402515;40.6329345670797];        
        
        % give warning of default
        warning(strcat('paramCase =  ', paramCase, ' not recognized, parameters for paramCase = ''mt'' returned instead'))
end