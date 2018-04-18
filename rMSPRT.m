% rMSPRT

clc
clear



%% PARAMETERS
% general 
iCoherence = 4;% index of coherence % within {3.2, 6.4, 12.8, 25.6, 51.2}
out.params.noHyp = 2;% No. hypotheses, N (here equal to No. channels, C)
out.params.iPrefCh = 1;% index of preferred channel
out.params.pdf = 'l';% probability density function of inputs (keep as 'l', which stands for lognormal)
out.params.thresh = 1/100;% decision threshold for basal ganglia output

% free parameters
out.params.constBaseline = 15;% constant baseline, l    
out.intControls.ctxThWeight = 0.4;% corticothalamic weitght, w_yt 
dataScalingFactor = 40;% data scaling factor, n

% MT data parameter set among:
%......'mt': computed from MT data (set \Omega)
%......'mtInfoDeplN2': after information depletion (\Omega_d) for N = 2 
%......'mtInfoDeplN4': after information depletion (\Omega_d) for N = 4
out.params.pcase = 'mtInfoDeplN2';       

% retrieve parameter set
[meansNullCh, meansPrefCh, stdsNullCh, stdsPrefCh] =  paramRepository(out.params.pcase);

% pick parameters corresponding to chosen coherence %
% and scale them by n 
out.params.meanNullCh = meansNullCh(iCoherence)/dataScalingFactor;
out.params.meanPrefCh = meansPrefCh(iCoherence)/dataScalingFactor;
out.params.stdNullCh = stdsNullCh(iCoherence)/dataScalingFactor;
out.params.stdPrefCh = stdsPrefCh(iCoherence)/dataScalingFactor;  
    
    


% Path of the m-file producing these results
out.generatingMFile = mfilename('fullpath');

%% INTERNAL CONTROLS
% Decision detection, on = 1, off = 0
out.intControls.decdelayBgTh = 1;

% Thalamocortical input on. rMSPRT = 1, MSPRT = 0
out.intControls.thCtxSwitch = 1;

%% Intermediate MT statistics for lateral choices when noHyp=4?, 1: yes, 0: no
%indInterLat = 0;
%if indInterLat == 1 && noHyp == 4
%    iLat = [3 4];% Index of lateral hypotheses
%end



% Delays
out.intControls.delayCtxTh = 1;% cortex to thalamus
out.intControls.delayThCtx = 1;% thalamus to cortex
out.intControls.delayCtxBg = 1;% cortex to BG output
out.intControls.delayBgTh = 1;% BG output to thalamus
loopDelay = out.intControls.delayCtxBg + out.intControls.delayBgTh +...
    out.intControls.delayThCtx;% Overall loop delay in observations, 4
fbExtraSteps = out.intControls.delayBgTh +...
    out.intControls.delayThCtx;% Extra steps to complete feedback to cortex post-termination

%% Scaling fraction of marginal distribution if Levy processes are assumed
%% 1: no Levy scaling
%out.intControls.dt = 1;

% % Observations after termination
% extraObs = 0;


% pre-decision values
out.comps.bg = (-log(1/out.params.noHyp)) * ones(1,out.params.noHyp);
out.comps.th = (log(1/out.params.noHyp) +...
    out.intControls.ctxThWeight*(...
    (out.intControls.ctxThWeight*out.params.constBaseline +...
    log(1/out.params.noHyp)) +...
    out.params.constBaseline)) .*...
    ones(1,out.params.noHyp);
out.comps.ctx = ((1 - out.intControls.thCtxSwitch)*log(1/out.params.noHyp) +...
    out.intControls.thCtxSwitch *...
    (out.intControls.ctxThWeight*out.params.constBaseline +...
    log(1/out.params.noHyp)) +...
    out.params.constBaseline) .* ones(1,out.params.noHyp);

% Indicator of decision termination, 0 = not terminated (starting value), 1 = terminated
indTerm = 0;

% Indicator of BG signal at termination being fed back to cortex
indFedback = 0;

% % Indicator of extra steps taken after termination
% indExtraObs = 0;

% Default/error output 
out.results.obsTerm = 0;% Observation at which decision is terminated
out.results.iSelHyp = 0;% Index of hypothesis selected

% Start observation counter
t = 0;

% % Start extra observation counter
% tExtra = 0;

%% RUN OBSERVATION-BY-OBSERVATION 
while indTerm == 0 || indFedback == 0 
    % Index of new observation
    t = t + 1;
    
    % INPUT STOCHASTIC PROCESSES
    % Observing the stochastic processes at t
%     if indTerm == 0% If decision has not been terminated
        % Make random observations for all channels
        out.comps.data(t, :) = stochProcGen(out.params.pdf,...
            out.params.meanNullCh, out.params.meanPrefCh,...
            out.params.stdNullCh, out.params.stdPrefCh, ...
            out.params.iPrefCh, 1, out.params.noHyp);
%        % Correct observations to simulate lateral choices
%        if indInterLat == 1 && noHyp == 4
%            fracRedStatDif = 0.4;% Fraction that statistics are reduced ~0.4 from STA92 data
%            out.comps.data(t, iLat) = input_sp_gen(out.params.pdf,...
%                out.params.meanNullCh -...
%                (fracRedStatDif*(out.params.meanNullCh - out.params.meanPrefCh)),...
%                out.params.meanNullCh -...
%                (fracRedStatDif*(out.params.meanNullCh - out.params.meanPrefCh)),...
%                out.params.stdNullCh -...
%                (fracRedStatDif*(out.params.stdNullCh - out.params.stdPrefCh)),...
%                out.params.stdNullCh -...
%                (fracRedStatDif*(out.params.stdNullCh - out.params.stdPrefCh)),...
%                out.intControls.dt, 1, 1, 2);
%        end
%     else% After termination...
% %         % ... feed a value close to 0 to all channels
% %         out.comps.data(t, :) = 1e-10 * ones(1, out.params.noHyp);
% %         % ... make the random observations equal to out.those of the null
% %         % channel
% %         out.comps.data(t, :) = input_sp_gen(out.params.pdf, out.params.meanNullCh, out.params.meanNullCh, out.params.stdNullCh, out.params.stdNullCh, out.intControls.dt, out.params.iPrefCh, 1, out.params.noHyp);
%         % ... observe normally
%         out.comps.data(t, :) = input_sp_gen(out.params.pdf,...
%             out.params.meanNullCh, out.params.meanPrefCh,...
%             out.params.stdNullCh, out.params.stdPrefCh,...
%             out.intControls.dt, out.params.iPrefCh, 1,...
%             out.params.noHyp);
%         % Correct observations to simulate lateral choices
%         if indInterLat == 1 && noHyp == 4
%             fracRedStatDif = 0.4;% Fraction that statistics are reduced ~0.4 from STA92 data
%             out.comps.data(t, iLat) = input_sp_gen(out.params.pdf,...
%                 out.params.meanNullCh -...
%                 (fracRedStatDif*(out.params.meanNullCh - out.params.meanPrefCh)),...
%                 out.params.meanNullCh -...
%                 (fracRedStatDif*(out.params.meanNullCh - out.params.meanPrefCh)),...
%                 out.params.stdNullCh -...
%                 (fracRedStatDif*(out.params.stdNullCh - out.params.stdPrefCh)),...
%                 out.params.stdNullCh -...
%                 (fracRedStatDif*(out.params.stdNullCh - out.params.stdPrefCh)),...
%                 out.intControls.dt, 1, 1, 2);
%         end
%     end
    
      
% lognormal integration
% Gains
                g4 = (log(((out.params.meanPrefCh^2))/...
                    sqrt((out.params.stdPrefCh^2) +...
                    (out.params.meanPrefCh^2)))/...
                    log((((out.params.stdPrefCh^2))/...
                    ((out.params.meanPrefCh^2)))+1)) -...
                    (log(((out.params.meanNullCh^2))/...
                    sqrt((out.params.stdNullCh^2  )+...
                    (out.params.meanNullCh^2  )))/...
                    log((((out.params.stdNullCh^2  ))/...
                    ((out.params.meanNullCh^2  )))+1));
                g5 = (1/(2*log(((out.params.stdPrefCh^2)/...
                    (out.params.meanPrefCh^2))+1))) -...
                    (1/(2*log(((out.params.stdNullCh^2)/...
                    (out.params.meanNullCh^2 ))+1)));
                % Cumulative factors
                if t == 1% First step with priors only
                    yy1 = zeros(1, out.params.noHyp);
                    yy2 = zeros(1, out.params.noHyp);
                elseif  t <= (loopDelay + 1) || out.intControls.thCtxSwitch == 0                
                    yy1(2:t,:) = cumsum(log(out.comps.data(2:t,:)), 1);
                    yy2(2:t,:) = cumsum((log(out.comps.data(2:t,:))).^2, 1);
                else
                    dummy = cumsum(log(out.comps.data((t - loopDelay + 1):t, :)), 1);
                    yy1(t,:) = dummy(loopDelay, :);
                    dummy = cumsum((log(out.comps.data((t - loopDelay + 1):t, :))).^2, 1);
                    yy2(t,:) = dummy(loopDelay, :);
                end
                % Integration
                int(t,:) = g4*yy1(t,:) - g5*yy2(t,:);
    
   
    
    % EXECUTION
%     % If the decision has not been terminated yet, simulate all sites
%     if indTerm == 0
        % Integrator cortex
        % Cortex not yet influenced by thalamus?
        if t <= out.intControls.delayThCtx
            out.comps.ctx(t,:) = int(t,:) +...
                (1 - out.intControls.thCtxSwitch)*log(1/out.params.noHyp) +...
                out.intControls.thCtxSwitch *...
                (out.intControls.ctxThWeight*out.params.constBaseline +...
                log(1/out.params.noHyp)) +...
                out.params.constBaseline;
        else% Thalamo-cortical input onset
            out.comps.ctx(t,:) = int(t,:) +...
                (1 - out.intControls.thCtxSwitch)*log(1/out.params.noHyp) +...
                out.intControls.thCtxSwitch *...
                out.comps.th(t - out.intControls.delayThCtx,:) +...
                out.params.constBaseline;
        end
        
        % Thalamus
        if t <= out.intControls.delayCtxTh% Thalamus not yet influenced by cortex?
%             out.comps.ctx(t,:) = int(t,:) + log(1/out.params.noHyp) +...
%                 out.intControls.ctxThWeight*out.params.constBaseline +...
%                 out.params.constBaseline;
            % Thalamus not yet influenced by basal ganglia?
            if t <= out.intControls.delayBgTh
                out.comps.th(t,:) = log(1/out.params.noHyp) +...
                    out.intControls.ctxThWeight*(...
                    (out.intControls.ctxThWeight*out.params.constBaseline +...
                    log(1/out.params.noHyp)) +...
                    out.params.constBaseline);
            else% Gangliothalamic feedback onset
                out.comps.th(t,:) =...
                    out.intControls.ctxThWeight*(...
                    (out.intControls.ctxThWeight*out.params.constBaseline +...
                    log(1/out.params.noHyp)) +...
                    out.params.constBaseline) -...
                    out.comps.bg(t - out.intControls.delayBgTh,:);
            end
            
        else% Corticothalamic input onset
            % Thalamus not yet influenced by basal ganglia?
            if t <= out.intControls.delayBgTh
                out.comps.th(t,:) = log(1/out.params.noHyp) +...
                    out.intControls.ctxThWeight *...
                    mean(out.comps.ctx(t - out.intControls.delayCtxTh, :));
                
            else% Gangliothalamic feedback onset
                out.comps.th(t,:) = out.intControls.ctxThWeight *...
                    mean(out.comps.ctx(t - out.intControls.delayCtxTh,:)) -...
                    out.comps.bg(t - out.intControls.delayBgTh,:);

            end
        end
        
        %Basal ganglia
        if t <= out.intControls.delayCtxBg %Cortical signal has not arrived
            out.comps.bg(t,:) = -log(1/out.params.noHyp);
        else%Cortical onset in the BG
            out.comps.bg(t,:) =...
                - out.comps.ctx(t - out.intControls.delayCtxBg,:) +...
                log(sum(exp(out.comps.ctx(t - out.intControls.delayCtxBg,:))));
        end
       

        % If feedback is complete after termination
        if t == (out.results.obsTerm + out.intControls.delayCtxBg + 1 + fbExtraSteps)...
                && indTerm == 1
            indFedback = 1;
        end
       
    
    %Decision detection for MSPRT
    if out.intControls.decdelayBgTh == 1 &&...
            min(out.comps.bg(t,:)) <= (1/out.params.thresh) &&...
            indTerm == 0
        % Observation at which decision is made, from the basal ganglia POV
        out.results.obsTerm = t - out.intControls.delayCtxBg - 1;
        % Change indicator of termination
        indTerm = 1;
        % Index of chosen hypothesis
        out.results.iSelHyp = find(out.comps.bg(t,:) ==...
            min(out.comps.bg(t,:)));
        
        % If the minimum bg is equal for multiple alternatives, choose 
        % one at random
        if length(out.results.iSelHyp) > 1
            out.results.iSelHyp =...
                out.results.iSelHyp(unidrnd(length(out.results.iSelHyp)));
            warning('More than one hypotheses was selected as equally good and among them one was reported randomly')
        end
        
  
    end
    

    
    
end



%% PLOT
% plot paramters
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
plot(tLine, out.comps.data*dataScalingFactor, 'k', 'lineWidth', lineW);
hold on
plot(tLine, out.comps.data(:, out.params.iPrefCh)*dataScalingFactor, 'r', 'lineWidth', lineW);
ylabel('sensory cortex', 'FontName', 'Arial')
axis tight
xlim([0 size(out.comps.data, 1)])
hold off
% sensorimotor cortex
h = subplot(2, 2, 2);
p = get(h, 'pos');
p(1) = p(1) - lMargDec;
p(3) = p(3)*widthScale;
set(h, 'pos', p);
set(gca, 'FontSize', fontSize)
plot(tLine, out.comps.ctx, 'k', 'lineWidth', lineW);
hold on
plot(tLine, out.comps.ctx(:,out.params.iPrefCh), 'r', 'lineWidth', lineW);
ylabel('sensorimotor cortex', 'FontName', 'Arial')
axis tight
xlim([0 size(out.comps.data, 1)])
hold off
% basal ganglia output
h = subplot(2, 2, 3);
p = get(h, 'pos');
p(1) = p(1) - lMargDec;
p(3) = p(3)*widthScale;
set(h, 'pos', p);
set(gca, 'FontSize', fontSize)
plot(tLine, out.comps.bg, 'k', 'lineWidth', lineW);
hold on
plot(tLine, out.comps.bg(:,out.params.iPrefCh), 'r', 'lineWidth', lineW);
ylabel('basal ganglia output', 'FontName', 'Arial')
axis tight
xlim([0 size(out.comps.data, 1)])
hold off
% thalamus
h = subplot(2, 2, 4);
p = get(h, 'pos');
p(1) = p(1) - lMargDec;
p(3) = p(3)*widthScale;
set(h, 'pos', p);
set(gca, 'FontSize', fontSize)
plot(tLine, out.comps.th, 'k', 'lineWidth', lineW);
hold on
plot(tLine, out.comps.th(:, out.params.iPrefCh), 'r', 'lineWidth', lineW);
ylabel('thalamus')
axis tight
xlim([0 size(out.comps.data, 1)])
hold off

 