% Copyright (C) 2015 by Mark A. Davenport, Andrew K. Massimino, 
%     Deanna Needell, and Tina Woolf

% This file is part of ``Constrained adaptive sensing.''

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Setup
% optAlg1 : 1 - Run CoSaMP to locate support
%           2 - Run OMP to locate support
%           3 - Run L1-Minimization to locate support
% opt :  1 - Run algorithm 1 to locate support, then use pseudo-inverse to
%            reconstruct
%        2 - Run algorithm 1 to locate support, then algorithm 2 to
%            reconstruct
% optAlg2 : 1 - If opt=2, run CoSaMP to reconstruct signal
%           2 - If opt=2, run OMP to reconstruct signal
%           3 - If opt=2, run L1-Minimization to reconstruct signal
% optRound: 1 - Simple rounding
%           2 - S Probability (typically always use this)
%           3 - Preconditioning S Probability (not implemented)
%           4 - Water emptying rounding
%
% NOTES: 
% 1) optAlg2 will be used for nonadaptive reconstruction
% 2) Pseudoinverse is always used for Oracle adaptive sensing

% Clean up and Setup
clc;
close all;
clearvars;
addpath TFOCS/
addpath spgl1/
global F Hstar;

USINGPARFOR = 0;        % must also change the parfor to for

QUIET     = 0;
numTrials = 200;
nlist     = 2.^[8 9 10 11 12 13 14];
mfactor   = 0.6; % mfactor*n = number of measurement to use
noisevar  = 0.01;
optAlg1   = 1; % locate support with CoSaMP (1), OMP (2), L1-Minimization (3) 
optAlg2   = 1; % reconstruct signal with CoSaMP (1), OMP (2), L1-Minimiztion (3) (assuming opt = 2)
opt       = 1; % use pseudoinverse with reconstruction (1) or Alg2 (2)
optRound  = 2;
variableDensity = 1; % use variable density sampling for nonadaptive measurements
preconditioning = 0; % use preconditioning for nonadaptive measurements
k         = 10;
suppType  = 'genSig'; % top, random, bottom, genSig
pNonadapt = 0.5;   % percent of measurements to use in finding support in adaptive sensing
resThresh = 10^10;  % residual threshold in support estimation. If norm(residual) > resThresh, run again. 
                    % Setting this value very large turns off the loop, will just run once
saveFigs  = 0;
saveMat   = 0;
seed      = 1; % random seed

rng('default')

%% Run through mlist and numTrials
% Initialize
errorsAdapVDS = zeros(length(nlist),numTrials);
errorsAdapUnif = zeros(length(nlist),numTrials);
correctSuppVDS = zeros(length(nlist),numTrials);
correctSuppUnif = zeros(length(nlist),numTrials);
optErrors = zeros(length(nlist),numTrials);
errorsNonVDS = zeros(length(nlist),numTrials);
errorsNonUnif = zeros(length(nlist),numTrials);
avgErrGoodTrialsVDS = zeros(length(nlist),1);
avgErrGoodTrialsMSEVDS = zeros(length(nlist),1);
avgErrGoodTrialsUnif = zeros(length(nlist),1);
avgErrGoodTrialsMSEUnif = zeros(length(nlist),1);

for nn = 1:length(nlist)
    n = nlist(nn);
    mu        = sqrt(n);
    m = round(mfactor*n);
    if ~QUIET
        display(['*** n = ', num2str(n)]);
        display(['*** m = ', num2str(m)]);
    end
    % F = 1/sqrt(n)*dftmtx(n);
    F = 1/sqrt(n)*fft(eye(n));
    H = haarmtx(n);
    Hstar = H';
    
    if USINGPARFOR
        globaldata = struct;
        globaldata.F = F;
        globaldata.Hstar = Hstar;
    else
        globaldata = [];
    end
        
    % Define Variable Density Sampling Probabilities once and pass v and d into
    % functions to save time
    if (variableDensity)
        kappa = max(abs(F*Hstar),[],2);
        v = kappa.^2/norm(kappa,2)^2;
        d = norm (kappa,2)./kappa;
    else
        v = [];
        d = [];
    end
    
    % for trial=1:numTrials
    for trial=1:numTrials
        if USINGPARFOR && isempty(getCurrentTask())
            error(['USINGPARFOR option specified but parfor is not ' ...
                'used!  You may need to open a pool.']);
        end

       successVDS  = 0; % success of TFOCS for adaptive sensing: 0 (no), 1 (yes). Continue until success.
       successUnif = 0;
       while (~(successVDS && successUnif))
           if ~QUIET
               display(['*** trial = ', num2str(trial)]);
           end
           switch suppType      
                case 'top'
                    supp = 1:k;          % support starting from top of signal
                    x = zeros(n,1);
                    % x(supp) = ones(k, 1)*10;
                    x(supp) = randn(k, 1)+mu;
                    %x = x/norm(x);
                case  'random'
                    supp = randperm(n);
                    supp = supp(1:k);    % random support of size k
                    x = zeros(n,1);
                    % x(supp) = ones(k, 1)*10;
                    x(supp) = randn(k, 1)+mu;
                    %x = x/norm(x);
                case 'bottom'
                    supp = n:-1:n-k+1;   % support starting from bottom of signal
                    x = zeros(n,1);
                    % x(supp) = ones(k, 1)*10;
                    x(supp) = randn(k, 1)+mu;
                    %x = x/norm(x);
               case 'genSig'
                   x = genSignal(n, k, randn(k, 1)+mu);
                   x = x';
                   %x = x'/norm(x);
                   supp = find(abs(x)>0);
           end
           
           % Run Adaptive: with VDS
           [xhatAdapVDS, supp1VDS, optErr, str1, successVDS] = runAdaptive(x, round(m*pNonadapt), m-round(m*pNonadapt), noisevar, optAlg1, opt, optRound, optAlg2, 1, resThresh, preconditioning, v, d, globaldata);
           % Run Adaptive: with Uniform Sampling
           [xhatAdapUnif, supp1Unif, optErr, str1, successUnif] = runAdaptive(x, round(m*pNonadapt), m-round(m*pNonadapt), noisevar, optAlg1, opt, optRound, optAlg2, 0, resThresh, preconditioning, v, d, globaldata);
           if (~(successVDS && successUnif))
               display('In tester: TFOCS failed somewhere - need to try new signal');
           end

       end
       errorsAdapVDS(nn, trial)  = norm(x-xhatAdapVDS);
       correctSuppVDS(nn, trial) = all(ismember(supp, supp1VDS));
       errorsAdapUnif(nn, trial) = norm(x-xhatAdapUnif);
       correctSuppUnif(nn, trial)= all(ismember(supp, supp1Unif));
       optErrors(nn, trial)   = optErr; % this could be from either adaptive run, with or without VDS
       
       % Run NonAdaptive: with VDS
       disp('Running VDS Nonadaptive');
       [xhatNon, str2] = runNonAdaptive(x, m, noisevar, optAlg2, 1, preconditioning, v, d, globaldata);
       errorsNonVDS(nn, trial) = norm(x-xhatNon);
       % Run NonAdaptive: with Uniform Sampling
       disp('Running Uniform Nonadaptive');
       [xhatNon, str2] = runNonAdaptive(x, m, noisevar, optAlg2, 0, preconditioning, v, d, globaldata);
       errorsNonUnif(nn, trial) = norm(x-xhatNon);
       
    end
    goodTrialsVDS = find(correctSuppVDS(nn,:)==1);
    avgErrGoodTrialsVDS(nn) = mean(errorsAdapVDS(nn, goodTrialsVDS));
    avgErrGoodTrialsMSEVDS(nn) = mean(errorsAdapVDS(nn, goodTrialsVDS).^2);
    goodTrialsUnif = find(correctSuppUnif(nn,:)==1);
    avgErrGoodTrialsUnif(nn) = mean(errorsAdapUnif(nn, goodTrialsUnif));
    avgErrGoodTrialsMSEVDS(nn) = mean(errorsAdapUnif(nn, goodTrialsUnif).^2);
    
    filename = sprintf('SignalRecoveryVaryDimension2_%dk_%1.1fmfactor_%dtrials_%sSupport_%1.2fsigma_%doptAlg1_%doptAlg2_%dopt_%dprecond',k,mfactor,numTrials,suppType,noisevar,optAlg1,optAlg2,opt,preconditioning)
    if (saveMat)
       save([filename '.mat']); 
    end
end

str = [str1 ', ' str2];

%% Plot

figure
set(gcf,'WindowStyle','normal')
semilogy(nlist, median(errorsNonVDS.^2, 2), '-ro','linewidth',2);hold on
semilogy(nlist, median(errorsNonUnif.^2, 2), '--ro','linewidth',2);
semilogy(nlist, median(errorsAdapVDS.^2, 2), '-b^','linewidth',2);
semilogy(nlist, median(errorsAdapUnif.^2, 2), '--b^','linewidth',2);
semilogy(nlist, median(optErrors.^2, 2), '-ks','linewidth',2);
hl = legend('VDS Nonadaptive','Uniform Nonadaptive', 'VDS Adaptive','Uniform Adaptive', 'Oracle Adaptive');
set(hl,'fontsize',20);
ylabel('Median Squared Error','fontsize',22,'FontName','Helvetica');
xlabel('Signal Dimension','fontsize',22,'FontName','Helvetica');
set(gca,'FontSize',20);
grid on;
% title({str;['n=' num2str (n) ', k=' num2str (k) ', Support type=' suppType ' sigma= ' num2str (noisevar) ', Trials=', num2str (numTrials) ', p=' num2str (pNonadapt) ', VDS = ' num2str (variableDensity) ', precond = ' num2str(preconditioning)]});
hold off
set(gcf,'PaperPositionMode','Auto');
filename3 = sprintf('SignalRecoveryVaryDimensionMedian_%dk_%1.1fmfactor_%dtrials_%sSupport_%1.2fsigma_%doptAlg1_%doptAlg2_%dopt_%dprecond',k,mfactor,numTrials,suppType,noisevar,optAlg1,optAlg2,opt,preconditioning)
if (saveFigs)  
    print('-depsc','-r0',[filename3 '.eps'])
    saveas(gcf,[filename3 '.png'],'png');  
    saveas(gcf,[filename3 '.fig'],'fig');
end

figure
set(gcf,'WindowStyle','normal')
semilogy(nlist, median(errorsNonVDS.^2, 2)./median(errorsAdapVDS.^2, 2), '-go','linewidth',2);hold on
semilogy(nlist, median(errorsNonUnif.^2, 2)./median(errorsAdapUnif.^2, 2), '--go','linewidth',2);hold on
semilogy(nlist, log(nlist), '-m^','linewidth',2);
semilogy(nlist, nlist./k, '--m^','linewidth',2);
hl = legend('VDS Nonadaptive/Adaptive','Uniform Nonadaptive/Adaptive', 'Log(n)','n/s');
set(hl,'fontsize',20,'location','southeast');
ylabel('Median Squared Error Ratio','fontsize',22,'FontName','Helvetica');
xlabel('Signal Dimension','fontsize',22,'FontName','Helvetica');
set(gca,'FontSize',20);
xlim([0 nlist(end)])
grid on;
% title({str;['n=' num2str (n) ', k=' num2str (k) ', Support type=' suppType ' sigma= ' num2str (noisevar) ', Trials=', num2str (numTrials) ', p=' num2str (pNonadapt) ', VDS = ' num2str (variableDensity) ', precond = ' num2str(preconditioning)]});
hold off
set(gcf,'PaperPositionMode','Auto');
filename8 = sprintf('SignalRecoveryVaryDimensionRatioMedianVDSNonVSAdapt_%dk_%1.1fmfactor_%dtrials_%sSupport_%1.2fsigma_%doptAlg1_%doptAlg2_%dopt_%dprecond',k,mfactor,numTrials,suppType,noisevar,optAlg1,optAlg2,opt,preconditioning)
if (saveFigs)  
    print('-depsc','-r0',[filename8 '.eps']) 
    saveas(gcf,[filename8 '.png'],'png');  
    saveas(gcf,[filename8 '.fig'],'fig');
end
