% @fileName: rMSPRT.m
% @purpose: implementation of the recursive multi-hypothesis sequential probability 
%                ratio test (rMSPRT), as introduced in Caballero JA, Humphries MD, 
%                Gurney KN (2018) A probabilistic, distributed, recursive mechanism for 
%                decision-making in the brain. PLoS Comput Biol 14(4): e1006033. 
%                https://doi.org/10.1371/journal.pcbi.1006033. Also enabling to effectively
%                simulate the (non-recursive) traditional MSPRT for comparison.
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



clc
clear



%% PARAMETERS
% general 
thCtxOn = 1;% thalamocortical connection 1 = on (rMSPRT) or 0 = off (effectively MSPRT)
iCoherence = 4;% index of coherence % within {3.2, 6.4, 12.8, 25.6, 51.2}
noHyp = 4;% No. hypotheses, N (here equal to No. channels, C)
iPrefCh = 1;% index of preferred channel
mtPdf = 'l';% probability density function of ISIs (keep as 'l', which stands for lognormal)
decThresh = 1/100;% decision threshold for basal ganglia output, \theta

% free parameters
constBaseline = 15;% constant baseline, l    
ctxThWeight = 0.4;% corticothalamic weitght, w_{yu}
dataScalingFactor = 40;% data scaling factor, n

% pick MT data parameter set among:
%......'mt': computed from MT data (set \Omega)
%......'mtInfoDeplN2': after information depletion (\Omega_d) for N = 2 
%......'mtInfoDeplN4': after information depletion (\Omega_d) for N = 4
paramCase = 'mtInfoDeplN2';       

% retrieve parameter set
[meanNullCh, meanPrefCh, stdNullCh, stdPrefCh] =  paramRepository(paramCase);

% keep only parameters corresponding to chosen 
% coherence % and scale them by n 
meanNullCh = meanNullCh(iCoherence)/dataScalingFactor;
meanPrefCh = meanPrefCh(iCoherence)/dataScalingFactor;
stdNullCh = stdNullCh(iCoherence)/dataScalingFactor;
stdPrefCh = stdPrefCh(iCoherence)/dataScalingFactor;  
    
% latencies
delayCtxTh = 1;% cortex to thalamus, \delta_{yu}
delayThCtx = 1;% thalamus to cortex, \delta_{uy}
delayCtxBg = 1;% cortex to basal ganglia output, \delta_{yb}
delayBgTh = 1;% basal ganglia output to thalamus, \delta_{bu}
loopDelay = delayCtxBg + delayBgTh +...
    delayThCtx;% overall loop delay, \Delta
fbExtraSteps = delayBgTh +...
    delayThCtx;% extra steps to allow feedback to reach cortex post-termination

    

%% SET-UP 
% pre-decision start-up values
bg = (-log(1/noHyp)) * ones(1,noHyp);% basal ganglia
th = (log(1/noHyp) +...% thalamus
    ctxThWeight*(...
    (ctxThWeight*constBaseline +...
    log(1/noHyp)) +...
    constBaseline)) .*...
    ones(1,noHyp);
ctx = ((1 - thCtxOn)*log(1/noHyp) +...% sensorimotor cortex
    thCtxOn *...
    (ctxThWeight*constBaseline +...
    log(1/noHyp)) +...
    constBaseline) .* ones(1,noHyp);

% indicator of decision termination, 0 = not terminated, 1 = terminated
indTerm = 0;

% indicator of basal ganglia output at termination reaching cortex
indFedback = 0;

% default/error outputs 
obsTerm = 0;% observation at which decision is terminated
iSelHyp  = 0;% index of hypothesis selected

% start observation counter
t = 0;



%% RUN OBSERVATION-BY-OBSERVATION 
while indTerm == 0 || indFedback == 0 
    % update observation index, t
    t = t + 1;
    
    % sensory cortex, x(t)
    data(t, :) = stochProcGen(mtPdf,...
        meanNullCh, meanPrefCh,...
        stdNullCh, stdPrefCh, ...
        iPrefCh, 1, noHyp);
    
    % lognormal simplified log-likelihood
    % gains
    g1 = (1/(2*log(((stdPrefCh^2)/(meanPrefCh^2))+1))) -...% g_1
        (1/(2*log(((stdNullCh^2)/(meanNullCh^2 ))+1)));
    g2 = (log(((meanPrefCh^2))/sqrt((stdPrefCh^2) +...% g_2
        (meanPrefCh^2)))/log((((stdPrefCh^2))/...
        ((meanPrefCh^2)))+1)) - (log(((meanNullCh^2))/...
        sqrt((stdNullCh^2  )+(meanNullCh^2  )))/...
        log((((stdNullCh^2  ))/((meanNullCh^2  )))+1));
    % cumulative factors
    if t == 1% first step with priors only
        y1 = zeros(1, noHyp);
        y2 = zeros(1, noHyp);
    elseif  t <= (loopDelay + 1) || thCtxOn == 0% to execute MSPRT                
        y1(2:t,:) = cumsum((log(data(2:t,:))).^2, 1);
        y2(2:t,:) = cumsum(log(data(2:t,:)), 1);
    else
        dummy = cumsum((log(data((t - loopDelay + 1):t, :))).^2, 1);
        y1(t,:) = dummy(loopDelay, :);
        dummy = cumsum(log(data((t - loopDelay + 1):t, :)), 1);
        y2(t,:) = dummy(loopDelay, :);
    end
    # simplified log-likelihood, y_i(t)
    simpLL(t,:) = g2*y2(t,:) - g1*y1(t,:);
    
    
    % sensorimotor cortex...
    % ... not yet influenced by thalamus?
    if t <= delayThCtx
        ctx(t,:) = simpLL(t,:) + (1 - thCtxOn)*log(1/noHyp) +...
            thCtxOn * (ctxThWeight*constBaseline +...
            log(1/noHyp)) + constBaseline;
    else% thalamo-cortical input onset
        ctx(t,:) = simpLL(t,:) + (1 - thCtxOn)*log(1/noHyp) +...
            thCtxOn*th(t - delayThCtx,:) + constBaseline;
    end
    
    % thalamus...
    if t <= delayCtxTh% ... not yet influenced by cortex?
        % ... not yet influenced by basal ganglia?
        if t <= delayBgTh
            th(t,:) = log(1/noHyp) + ctxThWeight*(...
                (ctxThWeight*constBaseline + log(1/noHyp)) +...
                constBaseline);
        else% gangliothalamic feedback onset
            th(t,:) = ctxThWeight*((ctxThWeight*constBaseline +...
                log(1/noHyp)) + constBaseline) - bg(t - delayBgTh,:);
        end
        
    else% corticothalamic feedback onset
        % thalamus not yet influenced by basal ganglia?
        if t <= delayBgTh
            th(t,:) = log(1/noHyp) + ctxThWeight *...
                mean(ctx(t - delayCtxTh, :));            
        else% gangliothalamic feedback onset
            th(t,:) = ctxThWeight * mean(ctx(t - delayCtxTh,:)) -...
                bg(t - delayBgTh,:);               
        end
        
    end
    
    % basal ganglia
    if t <= delayCtxBg % cortical signal has not arrived
        bg(t,:) = -log(1/noHyp);
    else% cortical signal onset in the basal ganglia
        bg(t,:) =  - ctx(t - delayCtxBg,:) +...
            log(sum(exp(ctx(t - delayCtxBg,:))));
    end
   
    % if feedback is complete after termination
    if t == (obsTerm + delayCtxBg + 1 + fbExtraSteps)...
            && indTerm == 1
        indFedback = 1;% flag it
    end
       
    
    % decision detection 
    if min(bg(t,:)) <= decThresh && indTerm == 0
        % termination observation (at the basal ganglia), T
        obsTerm = t - delayCtxBg - 1;
        % flag termination
        indTerm = 1;
        % index of selected hypothesis
        iSelHyp  = find(bg(t,:) == min(bg(t,:)));
        
        % in the rare case that several alternatives are 
        % selected, pick one at random
        if length(iSelHyp ) > 1
            iSelHyp  =  iSelHyp (unidrnd(length(iSelHyp )));
            warning('Multiple hypotheses were selected and one of them was reported randomly')
        end       
  
    end  
    
end



%% PLOT
% plot parameters
tLine = 1:t;% timeline
fontSize = 10;% font size 
lineW = 1.5;% line width 
lMargDec = 0.03; % left margin decrement
widthScale = 0.98; % width scaling factor 

figure(1)
% sensory cortex
h = subplot(2,2,1);
p = get(h, 'pos');
p(1) = p(1) - lMargDec;
p(3) = p(3)*widthScale;
set(h, 'pos', p);
set(gca, 'FontSize', fontSize)
plot(tLine, data*dataScalingFactor, 'k', 'lineWidth', lineW);
hold on
plot(tLine, data(:, iPrefCh)*dataScalingFactor, 'r', 'lineWidth', lineW);
ylabel('sensory cortex', 'FontName', 'Arial')
axis tight
xlim([0 size(data, 1)])
hold off
% sensorimotor cortex
h = subplot(2, 2, 2);
p = get(h, 'pos');
p(1) = p(1) - lMargDec;
p(3) = p(3)*widthScale;
set(h, 'pos', p);
set(gca, 'FontSize', fontSize)
plot(tLine, ctx, 'k', 'lineWidth', lineW);
hold on
plot(tLine, ctx(:,iPrefCh), 'r', 'lineWidth', lineW);
ylabel('sensorimotor cortex', 'FontName', 'Arial')
axis tight
xlim([0 size(data, 1)])
hold off
% basal ganglia output
h = subplot(2, 2, 3);
p = get(h, 'pos');
p(1) = p(1) - lMargDec;
p(3) = p(3)*widthScale;
set(h, 'pos', p);
set(gca, 'FontSize', fontSize)
plot(tLine, bg, 'k', 'lineWidth', lineW);
hold on
plot(tLine, bg(:,iPrefCh), 'r', 'lineWidth', lineW);
ylabel('basal ganglia output', 'FontName', 'Arial')
xlabel('time step')
axis tight
xlim([0 size(data, 1)])
hold off
% thalamus
h = subplot(2, 2, 4);
p = get(h, 'pos');
p(1) = p(1) - lMargDec;
p(3) = p(3)*widthScale;
set(h, 'pos', p);
set(gca, 'FontSize', fontSize)
plot(tLine, th, 'k', 'lineWidth', lineW);
hold on
plot(tLine, th(:, iPrefCh), 'r', 'lineWidth', lineW);
ylabel('thalamus')
axis tight
xlim([0 size(data, 1)])
hold off

 