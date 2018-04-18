%function [out] =...
%    rMSPRT_online(pdf, noHyp, iPrefCh, stdNullCh, stdPrefCh, meanNullCh, meanPrefCh, thresh, constBaseline)
% rMSPRT

clc
clear


%% PARAMETERS
    out.params.pdf = 'l';% probability density function for the inputs
    out.params.noHyp = 2;% No. hypotheses, N (equal to No. channels, C)
    out.params.iPrefCh = 1;% index of preferred channel
    out.params.constBaseline = 15;% constant baseline, l
    %Input process parameters
    
%    out.params.pcase = 'mtNorm';%Parameterization case
%    % NOTE: non-normalized parameterizations, like 'mt' make the model ring a
%    %       lot, disrupting the inference.
%    %......bg, Bogacz & Gurneout.y, 2007
%    %......ev, "Neuro-biological" equi-variant cannels
%    %......sp, Single parameter, assumed to be the mean
%    %......nb, Neuro-biological, different means and variances
%    %......ev144, ISI equi-variant cannels (based on Britten1992 w144, 12.8%)
%    %......sp144, ISI single-degree (based on Britten1992 w144, 12.8%), assumed to be the mean (it would not change wiout.th the population sampling correction)
%    %......nb144, ISI neuro-biological (based on Britten1992 w144, 12.8%), different means and variances
%    %......nb144f, IFR neuro-biological (based on Britten1992 w144, 12.8%)
%    %......nb144l, ISI neuro-biological (based on Britten1992 w144, 12.8%), and Levy scaling for the preferred channel
%    %......nb144ps,  " " " " " " ", correcting as in population sampling out.params.stdPrefCh/sqrt(No. ISIs in bin)
%    %......nb144psc,  " " " " " " " " " " " ",wiout.th Gaussian correction for all  cN*out.params.stdPrefCh/sqrt(No. ISIs in bin)
%    %......nb144pss, " " " " " " " " " " " ", swapping (for the Race model) the statistics of the preferred and background channels
%    %......nb144n, ISI neuro-biological (based on Britten1992 w144, 12.8%), statistics not normalized by StepMeanings
%    %......nb144nps, " " " " " " " " " " " ", correcting as in population sampling out.params.stdPrefCh/sqrt(No. ISIs in bin)
%    %......nb144npsc, " " " " " " " " " " " " " " " " ", wiout.th Gaussian correction for all cN*out.params.stdPrefCh/sqrt(No. ISIs in bin)
%    %......sb, Sandbox...
%    %Retrieving the parameters from the "repositorout.y" script
%    [StepMeanings, duration, dt, meansNullCh, meansPrefCh, stdsNullCh, stdsPrefCh,...
%        desERs] =...
%        inputparamsDB(out.params.pcase);
        
    % NOTE: this parameter was estimated from the whole sample of
        %       Britten et al. 1992. The estimates are made out of between
        %       ~70 and 213 MT neurons (less recordings for high coherences
        %       in the null/background direction). The parameters are
        %       vectors/matrices where the rows are the coherence % in
        %       order (0, 3.2, 6.4, 12.8, 25.6, 51.2). To estimate them,
        %       trials are classified according to the direction of the
        %       dots. In this variation of the parameterizations all the
        %       statistics are normalized by the StepMeanings(normCohNo) in
        %       order for the integrations in the rMSPRT to be positive
        
        dataScalingFactor = 40;% data scaling factor, n
        

                % NOTE: Both correct and error trials were used to estimate 
                %       these statistics 
                
                % Error rates (fraction)
                desERs = [0.5; 0.373293137;...
                    0.311515152; 0.189561878;...
                    0.079708826; 0.020252118];

%                % How much time a sampling step is equivalent to (mean ISI in
%                % preferred channel)? (ms)
%                StepMeanings = [60.2704696282276;54.1067812950675;...% sensory
%                    51.9555167200400;46.1073617021512;...
%                    37.7431198991798;29.8554021786382]...
%                    ./dataScalingFactor;
%                StepMeanings = dataScalingFactor*ones(6,1);
                % Mean of NULL channel
                meansNullCh = [59.8388569087918;59.3736631745858;...% sensory
                    62.9023641601734;65.5267578611166;...
                    70.1465464403720;83.4519960352983]...
                    ./dataScalingFactor;
                % Standard deviation of NULL channel
                stdsNullCh = [34.6210540304759;34.5386065496385;...% sensory
                    35.3233041337630;36.1434994400871;...
                    37.1454736402515;40.6329345670797]...
                    ./dataScalingFactor;
                % Mean of PREFERRED channel
                meansPrefCh = [60.2704696282276;54.1067812950675;...% sensory
                    51.9555167200400;46.1073617021512;...
                    37.7431198991798;29.8554021786382]...
                    ./dataScalingFactor;
                % Standard deviation of PREFERRED channel
                stdsPrefCh = [35.0144580577522;33.1325606188147;...% sensory
                    32.1756309152138;30.4920600270835;...
                    28.0096721079346;25.9790487446720]...
                    ./dataScalingFactor;

        duration = 3000;% Maximum number of observations
%        dt = 1;% Ratio of observation stats vs. marginal ones (Levy processes)  
    
    
    % Index of the coherence level to be used
    iCoherence = 4;
    
    % Is out.this the mt parameterization (composed of vectors)?
%    if strcmp(out.params.pcase(1:2), 'mt')
        % Select the appropriate parameters
%        StepMeanings = StepMeanings(iCoherence);
        out.params.meanNullCh = meansNullCh(iCoherence);
        out.params.meanPrefCh = meansPrefCh(iCoherence);
        out.params.stdNullCh = stdsNullCh(iCoherence);
        out.params.stdPrefCh = stdsPrefCh(iCoherence);
%         desiER = round(desERs(iCoherence)*100)/100;
%    else% Take parameters as such
%        StepMeanings = StepMeanings;
%        out.params.meanNullCh = meansNullCh;
%        out.params.meanPrefCh = meansPrefCh;
%        out.params.stdNullCh = stdsNullCh;
%        out.params.stdPrefCh = stdsPrefCh;
%%         desiER = round(desERs*100)/100;
%    end
    
    % Diagnostics plot (input densities)
    figure(1)
    clf
    xx = linspace(1e-10,out.params.meanNullCh*5, 200);
    plot(xx, pdfAtX(out.params.pdf, out.params.meanNullCh, out.params.stdNullCh, xx),'k')
    hold on
    plot(xx, pdfAtX(out.params.pdf, out.params.meanPrefCh,...
        out.params.stdPrefCh, xx),'r')
    legend('Null','Pref')
    hold off
    
    out.params.thresh = 100;%Decision threshold


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

% Weight from cortex to thalamus, 0.7 good (CGD14) 0.25 (CHG 15), 0.1 good without 'm' and 1 
out.intControls.ctxThWeight = 0.4;

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
                g4 = (log(((out.params.meanPrefCh^2)*(out.intControls.dt^2))/...
                    sqrt((out.params.stdPrefCh^2)*out.intControls.dt +...
                    (out.params.meanPrefCh^2)*(out.intControls.dt^2)))/...
                    log((((out.params.stdPrefCh^2)*out.intControls.dt)/...
                    ((out.params.meanPrefCh^2)*(out.intControls.dt^2)))+1)) -...
                    (log(((out.params.meanNullCh^2) * (out.intControls.dt^2))/...
                    sqrt((out.params.stdNullCh^2  )*out.intControls.dt +...
                    (out.params.meanNullCh^2  )*(out.intControls.dt^2)))/...
                    log((((out.params.stdNullCh^2  )*out.intControls.dt)/...
                    ((out.params.meanNullCh^2  )*(out.intControls.dt^2)))+1));
                g5 = (1/(2*log(((out.params.stdPrefCh^2*out.intControls.dt)/...
                    (out.params.meanPrefCh^2*out.intControls.dt^2))+1))) -...
                    (1/(2*log(((out.params.stdNullCh^2*out.intControls.dt)/...
                    (out.params.meanNullCh^2  *out.intControls.dt^2))+1)));
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

spms = mean(sum(exp(-out.comps.bg')));%Mean of the sum of all the probabilities



    %Timeline for graphs
    tLine = 1:t;

    fsize = 10;%Font size for plots
    linew = 1.5;%Line widout.th for plots
    ldecr = 0.03; %left margin decrement for subplots
    wfact = 0.98; %widout.th scaling factor for subplots
    
    
    figure(1)
    h=subplot(2,3,1);%%%% Cumulative x_i
    p = get(h, 'pos');
    p(1) = p(1)-ldecr;
    p(3) = p(3)*wfact;
    set(h, 'pos', p);
    set(gca,'FontSize',fsize)
    plot(tLine, yy1,'k','LineWidth',linew);
    hold on
    plot(tLine, yy1(:,out.params.iPrefCh),'r','LineWidth',linew);
    ylabel('Cumulative $x_i$','Interpreter','latex')
    axis tight
    xlim([0 size(out.comps.data, 1)])
    hold off
    h=subplot(2,3,2);%%%% Cumulative 1/x_i
    p = get(h, 'pos');
    p(1) = p(1)-ldecr;
    p(3) = p(3)*wfact;
    set(h, 'pos', p);
    set(gca,'FontSize',fsize)
    plot(tLine, yy2,'k','LineWidth',linew);
    hold on
    plot(tLine, yy2(:,out.params.iPrefCh),'r','LineWidth',linew);
    ylabel('Cumulative $\frac{1}{x_i}$','Interpreter','latex')
    axis tight
    xlim([0 size(out.comps.data, 1)])
    hold off
    h=subplot(2,3,3);%%%% Raw out.y_i
    p = get(h, 'pos');
    p(1) = p(1)-ldecr;
    p(3) = p(3)*wfact;
    set(h, 'pos', p);
    set(gca,'FontSize',fsize)
    switch out.params.pdf
        case 'w'
            rw=-g1*yy1-g3*yy2;
        case 'l'
            rw=g4*yy1-g5*yy2;
        otherwise
            rw=-g1*yy1-g3*yy2;
    end
    plot(tLine, rw,'k','LineWidth',linew);
    hold on
    plot(tLine, rw(:,out.params.iPrefCh),'r','LineWidth',linew);
    ylabel('Raw $out.y_i$','Interpreter','latex')
    axis tight
    xlim([0 size(out.comps.data, 1)])
    hold off
    h=subplot(2,3,4);%%%% out.y_i+prior+initialFR
    p = get(h, 'pos');
    p(1) = p(1)-ldecr;
    p(3) = p(3)*wfact;
    set(h, 'pos', p);
    set(gca,'FontSize',fsize)
    plot(tLine, out.comps.ctx,'k','LineWidth',linew);
    hold on
    plot(tLine, out.comps.ctx(:,out.params.iPrefCh),'r','LineWidth',linew);
    ylabel('$out.y_i+lnP_i(T-\loopDelay t)+out.params.constBaseline$','Interpreter','latex')
    axis tight
    xlim([0 size(out.comps.data, 1)])
    hold off
    h=subplot(2,3,5);%%%% GPi/EP/SNr
    p = get(h, 'pos');
    p(1) = p(1)-ldecr;
    p(3) = p(3)*wfact;
    set(h, 'pos', p);
    set(gca,'FontSize',fsize)
    plot(tLine, out.comps.bg,'k','LineWidth',linew);
    hold on
    plot(tLine, out.comps.bg(:,out.params.iPrefCh),'r','LineWidth',linew);
    ylabel('GPi/EP/SNr','Interpreter','latex')
    axis tight
    xlim([0 size(out.comps.data, 1)])
    hold off
    h=subplot(2,3,6);%%%% Choice probabilitout.y
    p = get(h, 'pos');
    p(1) = p(1)-ldecr;
    p(3) = p(3)*wfact;
    set(h, 'pos', p);
    set(gca,'FontSize',fsize)
    es=exp(-out.comps.bg);
    plot(tLine,sum(exp(-out.comps.bg')),'b--','LineWidth',linew)
    hold on
    plot(tLine, es,'k','LineWidth',linew);
    plot(tLine, es(:,out.params.iPrefCh),'r','LineWidth',linew);
    ylabel('Choice Probabilitout.y','Interpreter','latex')
    axis tight
    xlim([0 size(out.comps.data, 1)])
    ylim([0 1.1])
    legend('Sum of all probabilities','Location','NorthOutside')
    hold off
    
    %StepMeanings=1;
    tLine = tLine*StepMeanings;
    
    figure(2)
    h=subplot(2,2,1);%%%% Sensory cortex
    p = get(h, 'pos');
    p(1) = p(1)-ldecr;
    p(3) = p(3)*wfact;
    set(h, 'pos', p);
    set(gca,'FontSize',fsize)
    plot(tLine, out.comps.data*StepMeanings,'k','LineWidth',linew);
    hold on
    plot(tLine, out.comps.data(:,out.params.iPrefCh)*StepMeanings,'r','LineWidth',linew);
    ylabel('Sensorout.y afferent', 'FontName', 'Arial')
    axis tight
    xlim([0 size(out.comps.data, 1)*StepMeanings])
    hold off
    h=subplot(2,2,2);%%%% Integrator cortex
    p = get(h, 'pos');
    p(1) = p(1)-ldecr;
    p(3) = p(3)*wfact;
    set(h, 'pos', p);
    set(gca,'FontSize',fsize)
    plot(tLine, out.comps.ctx,'k','LineWidth',linew);
    hold on
    plot(tLine, out.comps.ctx(:,out.params.iPrefCh),'r','LineWidth',linew);
    ylabel('Integrator cortex', 'FontName', 'Arial')
    axis tight
    xlim([0 size(out.comps.data, 1)*StepMeanings])
    hold off
    h=subplot(2,2,3);%%%% BG output
    p = get(h, 'pos');
    p(1) = p(1)-ldecr;
    p(3) = p(3)*wfact;
    set(h, 'pos', p);
    set(gca,'FontSize',fsize)
    plot(tLine, out.comps.bg,'k','LineWidth',linew);
    hold on
    plot(tLine, out.comps.bg(:,out.params.iPrefCh),'r','LineWidth',linew);
    ylabel('BG output nuclei', 'FontName', 'Arial')
    axis tight
    xlim([0 size(out.comps.data, 1)*StepMeanings])
    %xlabel('($t$)', 'FontName', 'Arial','Interpreter','latex')
    hold off
    h=subplot(2,2,4);%%%% Thalamus
    p = get(h, 'pos');
    p(1) = p(1)-ldecr;
    p(3) = p(3)*wfact;
    set(h, 'pos', p);
    set(gca,'FontSize',fsize)
    switch out.params.pdf
        case 'w'
            rw=-g1*yy1-g3*yy2;
        case 'l'
            rw=g4*yy1-g5*yy2;
        otherwise
            rw=-g1*yy1-g3*yy2;
    end
    plot(tLine, out.comps.th,'k','LineWidth',linew);
    hold on
    plot(tLine, out.comps.th(:,out.params.iPrefCh),'r','LineWidth',linew);
    ylabel('Thalamus')
    axis tight
    xlim([0 size(out.comps.data, 1)*StepMeanings])
    hold off
    
    % plot(tLine,yy1)
    % title('Cumulative out.comps.data')
    % grid
    %
    % figure(2)
    % plot(tLine, yy2);
    % title('Cumulative 1/out.comps.data')
    % grid
    %
    % figure(3)
    % plot(tLine, -g1*yy1-g2*yy2);
    % title('out.comps.ctx')
    % grid
    %
    % figure(4)
    % plot(tLine, out.comps.ctx);
    % title('out.comps.ctx+lnPi')
    % grid
    %
    % figure (5)
    % plot(tLine,out.comps.bg)
    % title('EP/SNr')
    % grid
    %
    % out.comps.ctx=-g1*yy1-g2*yy2;
    % out.yp=out.comps.ctx;
    % out.yp(:,out.params.iPrefCh)=[];
    % P=exp(-out.comps.bg);
    %
    % figure(6)
    % plot(tLine,sum(exp(-out.comps.bg')),'b--')
    % hold on
    % plot(tLine,exp(-out.comps.bg))
    % title('exp(-EP/SNr)')
    % legend('Sum of all probabilities')
    % grid
    % hold off
    %
    % figure(7)
    % hist(out.comps.ctx(:,out.params.iPrefCh), 25)
    % h = findobj(gca,'Tout.ype','patch');
    % set(h,'FaceColor','g')
    % hold on
    % hist(out.yp, 25)
    % title('out.comps.ctx distribution (green = winner)')
    % hold off
    
    % figure(7)
    % hist(out.comps.data(:,out.params.iPrefCh), 25)
    % h = findobj(gca,'Tout.ype','patch');
    % set(h,'FaceColor','g')
    % hold on
    % hist(xp, 25)
    % title('out.comps.data distribution (green = winner)')
    % hold off
    
    % tile
