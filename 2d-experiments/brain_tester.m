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

%% Brain_tester.m
% Test adaptive sensing with MRI image

% Clean up
clearvars;
close all;

% Add paths
addpath ../../RWT/bin;
addpath ../../TFOCS;
addpath ../../spgl1;

%% Setup Parameters
n1 = 64;            % image dimension is n1 x n2
n2 = n1;            %
opt_method = 'l1';  % Reconstruction method: l1 or cosamp (currenly only working for l1)
mlist = [4000:-500:500]; % List of number of measurements to loop through (should be even values)
noisevar = 1e-2;    % Gaussian noise standard deviation
klist  = [500 1000 1500];       % estimated sparsity level
ntrials_Nonadapt = 50;   % number of trials for nonadaptive sensing
ntrials_Adapt1 = 50; % number of trials in stage 1 of adaptive sensing (i.e., nonadaptive measurements)
ntrials_Adapt2 = 50; % number of trials in stage 2 of adaptive sensing (i.e., adaptive measurements)
run_Nonadapt = 1;   % Flag to run VDS Nonadaptive simulations (1=yes,0=no)
run_Adapt    = 1;   % Flag to run VDS Adaptive simulations (1=yes,0=no)
run_TFOCS    = 1;   % Flag to run TFOCS in adaptive simulations (1=yes,0=no and load from existing mat file)
saveMat      = 1;   % Flag to save data to mat files
saveFigs     = 0;   % Flag to save figures
plotFinalImages.flag = 0; % Flag to plot final "average" images
 plotFinalImages.kk = 1;  % mlist index for which measurement level to plot
 plotFinalImages.mm = 1;  % klist index for which sparsity level to plot

%% Signal loading and wavelet setup

% Load image
load ../brain

use_rwt = 1; % set to 1 if using the Rice Wavelet Transform package

if use_rwt
    hh = daubcqf(6, 'min');
else
    wname = 'haar';
end
nlevels = log2(n1);

% Original image = X, wavelet coefficients = a
X = imresize(abs(im), n1/512);
if use_rwt
    [a, ~] = mdwt(X, hh, nlevels);
    a = a(:);
else
    dwtmode per;
    [a, s] = wavedec2(X, nlevels, wname);
    a = a(:);
end

% Plot wavelet coefficients
% figure; clf;
% stem(real(a));
% title('Coefficients');

% Plot sorted absolute value coefficients
% figure; clf;
[sortedabscoef, Imax] = sort(abs(a), 'descend');
% semilogy(sortedabscoef);
% title('Sorted absolute value coefficients');

ncoef = length(a);

%% Generate Master Sensing Matrix
% Full n1*n2 by n1*n2 sensing matrix is FHstar
if ~exist('FHstar','var')
    if use_rwt
        % computes F*Hstar
        n = n1*n2;
        FHstar = zeros(n);

        fprintf('generating FHstar...\n');

        tic
        for ii = 1:n
            aa = zeros(n1, n2);
            aa(ii) = 1;
            HH = midwt(aa, hh, nlevels);        % computes 2d inverse wavelet
            XX = fft2( HH ) / sqrt(n1*n2);      % computes 2d FFT
            FHstar(:,ii) = XX(:);
        end
        toc
    else
        error('not supported');
    end
end

%% closest k-sparse approximation
for kk = 1:length(klist)
    k   = klist(kk);
    ind = Imax(1:k);
    accuracy = norm(a(ind)) ./ norm(a);
    fprintf('%g%% energy contained in %d (%g%%) largest support\n', ...
        accuracy*100, k, k/ncoef*100);
    FWinv = FHstar(:,ind);
    ahat_bestk(kk).ahat = zeros(ncoef, 1);
    ahat_bestk(kk).ahat(ind) = a(ind);
    Xhat_bestk(kk).Xhat = ifft2(reshape(FWinv * a(ind), n1, n2)) * sqrt(n1*n2);
    
    figure; clf;
    imagesc(abs(Xhat_bestk(kk).Xhat));
    PSNRdB = psnr(Xhat_bestk(kk).Xhat, X);
    title({['Best ' num2str(k) '-term (' num2str(k/ncoef*100) '%) Approximation'];...
        ['PSNR = ' num2str(PSNRdB) ' dB, MSE = ' num2str(norm(X-Xhat_bestk(kk).Xhat))]},'fontsize',16)
    if (saveFigs)
       filename = sprintf('Best_%d-term_approximation',k);
       %print('-depsc','-r0',[filename '.eps']) 
       saveas(gcf,[filename '.png'],'png');  
       %saveas(gcf,[filename '.fig'],'fig');
    end
end

%% Define VDS Probabilities
% compute VDS parameters
% v is the pmf for sampling
% d is the weighting (if being used)
% i.e., norm(diag(sqrt(v))*FHstar,'fro') = 1 <-- it's fair!
kappa = max(abs(FHstar), [], 2);
v = kappa.^2/norm(kappa,2)^2;
d = norm(kappa,2)./kappa;

%% Nonadaptive Sensing

if (run_Nonadapt)
    
    % allocate result variables
    Xhat_Nonadapt   = cell(length(klist),length(mlist),ntrials_Nonadapt);
    psnr_Nonadapt   = zeros(length(klist),length(mlist),ntrials_Nonadapt);
    error2_Nonadapt = zeros(length(klist),length(mlist),ntrials_Nonadapt);
    nnz_Nonadapt    = zeros(length(klist),length(mlist),ntrials_Nonadapt);
    
    for kk = 1:length(klist)
        k = klist(kk);
        for mm = 1:length(mlist)
            m = mlist(mm);
            for nn = 1:ntrials_Nonadapt

                % Generate VDS samples
                OmegaN = zeros(m,1);
                for ii = 1:m
                   OmegaN(ii) = GenerateRandomIndex(v); 
                end
                PhiHstarN  = FHstar(OmegaN,:);     % VDS nonadaptive sensing matrix
                noiseTermN = randn(m,1)*noisevar;  % Gaussian noise
                yN = PhiHstarN*a + noiseTermN;     % a is not exactly sparse

                % Perform Reconstruction
                switch opt_method
                    case 'cosamp'
                        tic
                        ahat = cosamp(PhiHstarN, yN, k, a); 
                        toc
                    case 'l1'
                        D = speye(n1*n2);
                        D = D(OmegaN,:);
                        fwd = @(a) D*reshape(fft2(midwt(...
                            reshape(a,n1,n2), hh, nlevels)),n1*n2,1) ./ sqrt(n1*n2);
                        adj = @(y) reshape(mdwt(ifft2(...
                            reshape(D'*y,n1,n2)), hh, nlevels),n1*n2,1) .* sqrt(n1*n2);
                        Aop = Aopfast(fwd,adj);

                        if (noisevar > 0)
                            bpdn_sigma = norm(noiseTermN) + norm(PhiHstarN*a - PhiHstarN*ahat_bestk(kk).ahat);
                        else
                            error('Need to decide how to define bpdn_sigma when there is no noise');
                            bpdn_sigma = 1;
                        end
                        spgl1_opts = spgSetParms('verbosity', 0);

                        fprintf('starting spgl1\n');
                        tic
                        ahat = spg_bpdn(Aop, yN, bpdn_sigma, spgl1_opts);
                        toc

                    otherwise 
                        error('invalid method');
                end

                % Estimated image, PSNR, and squared error
                Xhat_Nonadapt{kk,mm,nn}   = ifft2(reshape(FHstar*ahat,n1,n2))*sqrt(n1*n2);
                psnr_Nonadapt(kk,mm,nn)   = psnr(Xhat_Nonadapt{kk,mm,nn}, X);
                error2_Nonadapt(kk,mm,nn) = norm(X-Xhat_Nonadapt{kk,mm,nn})^2;
                nnz_Nonadapt(kk,mm,nn)    = nnz(ahat);

                disp(['Nonadaptive m ' num2str(m) ' trial ' num2str(nn) ' sparsity k ' num2str(k)]);
            end
        end
    end
    
    if (saveMat)
       filenameMat = sprintf('NonadaptiveImages_%dn_%dtrials_%1.4fsigma_%sopt_%dk1_%dkn_%dm1_%dmn',n1,ntrials_Nonadapt,noisevar,opt_method,klist(1),klist(end),mlist(1),mlist(end));
       save([filenameMat '.mat'],'Xhat_Nonadapt','psnr_Nonadapt','error2_Nonadapt');
       clear Xhat_Nonadapt;
    end

end

%% Adaptive Sensing

if (run_Adapt)
    
    % allocate result variables for stage 2
    Xhat_Adapt   = cell(length(klist),length(mlist),ntrials_Adapt2);
    psnr_Adapt   = zeros(length(klist),length(mlist),ntrials_Adapt2);
    error2_Adapt = zeros(length(klist),length(mlist),ntrials_Adapt2);
    
    for kk = 1:length(klist)
        k = klist(kk);
        for mm = 1:length(mlist)
            m = mlist(mm);
            m1 = m/2; % m1 = number of measurements used in stage 1
            m2 = m1;  % m2 = number of measurements used in stage 2
            
            %parfor trial=1:numTrials % make sure to open a pool
            for nn = 1:ntrials_Adapt1

                % First half of measurements are VDS nonadaptive
                Omega1 = zeros(m1,1);
                for ii = 1:m1
                    Omega1(ii) = GenerateRandomIndex(v); 
                end
                PhiHstar1 = FHstar(Omega1, :);     % m/2 VDS nonadaptive measurements
                noiseTerm1 = randn(m1,1)*noisevar; % Gaussian noise
                y1 = PhiHstar1*a + noiseTerm1;     % a is not exactly sparse

                % Perform Support Estimation
                switch opt_method
                    case 'cosamp'
                        tic
                        [ahat, supp_est, residual] = cosamp(PhiHstar1, y1, k, a); 
                        toc
                    case 'l1'
                        D = speye(n1*n2);
                        D = D(Omega1,:);
                        fwd = @(a) D*reshape(fft2(midwt(...
                            reshape(a,n1,n2), hh, nlevels)),n1*n2,1) ./ sqrt(n1*n2);
                        adj = @(y) reshape(mdwt(ifft2(...
                            reshape(D'*y,n1,n2)), hh, nlevels),n1*n2,1) .* sqrt(n1*n2);
                        Aop = Aopfast(fwd,adj);

                        if (noisevar > 0)
                            bpdn_sigma = norm(noiseTerm1) + norm(PhiHstar1*a - PhiHstar1*ahat_bestk(kk).ahat);
                        else
                            error('Need to decide how to define bpdn_sigma when there is no noise');
                            bpdn_sigma = 1;
                        end
                        spgl1_opts = spgSetParms('verbosity', 0);

                        fprintf('starting spgl1\n');            
                        tic
                        ahat = spg_bpdn(Aop, y1, bpdn_sigma, spgl1_opts);
                        toc

                    otherwise 
                        error('invalid method');
                end
                
                % Get indices of largest min(k,nnz(ahat)) coefficients
                nnz_ahat = nnz(ahat); % number of nonzeros in ahat
                [sortedabscoef_ahat, Imax_ahat] = sort(abs(ahat), 'descend');
                Imax_ahat_all(kk,mm).data(nn,1:min(k,nnz_ahat)) = Imax_ahat(1 : min(k,nnz_ahat));
                Imax_ahat_all(kk,mm).nnz(nn) = min(k,nnz_ahat);
                Imax_ahat_all(kk,mm).accuracy(nn) = length(intersect(Imax_ahat_all(kk,mm).data(nn,:),Imax(1:klist(kk)))) ...
                    / Imax_ahat_all(kk,mm).nnz(nn); % Accuracy = number in common with true top k support / number in this support estimate
                
                adapt1_data_all(kk,mm).y1(1:m1,nn) = y1;
                adapt1_data_all(kk,mm).noiseTerm1(1:m1,nn) = noiseTerm1;
                adapt1_data_all(kk,mm).PhiHstar1{nn} = PhiHstar1;
                
                % Store VDS measurement indices
                Omega1_all(kk,mm).data(nn,:) = Omega1;
                
            end
            
            % Idea: Want to select a "representative" support estimate
            % Choose support estimate with average accuracy over the trials
            [junk stage1_idx] = min(abs(Imax_ahat_all(kk,mm).accuracy - mean(Imax_ahat_all(kk,mm).accuracy)));
            supp_est_stage1 = Imax_ahat_all(kk,mm).data(stage1_idx,1:Imax_ahat_all(kk,mm).nnz(stage1_idx));
            % Get the corresponding VDS measurements
            Omega1_stage1   = Omega1_all(kk,mm).data(stage1_idx,:);
            y1_stage1       = adapt1_data_all(kk,mm).y1(:,stage1_idx);
            noiseTerm1_stage1 = adapt1_data_all(kk,mm).noiseTerm1(:,stage1_idx);
            PhiHstar1_stage1  = adapt1_data_all(kk,mm).PhiHstar1{stage1_idx};
            
            % Run TFOCS only once per m and k selection!
            % Second half of measurements are selected adaptively
            % select adaptive samples by solving relaxation with TFOCS
            if (run_TFOCS)
                G0 = FHstar(:,supp_est_stage1);
                disp('Running TFOCS');
                tic
                [omega, S, out, success] = bestrows_tfocs2(G0, 1, 0, 1, []);
                toc

                S = real(S);

                if (saveMat)
                   filenameMat = sprintf('TFOCS_Data_%dn_%dk_%dm_%dtrials_%1.4fsigma',n1,k,m,ntrials_Adapt1,noisevar); 
                   save([filenameMat '.mat'],'S','supp_est_stage1','Omega1_stage1','y1_stage1','noiseTerm1_stage1','PhiHstar1_stage1');
                end
            else
                filenameMat = sprintf('TFOCS_Data_%dn_%dk_%dm_%dtrials_%1.4fsigma',n1,k,m,ntrials_Adapt1,noisevar);
                load([filenameMat '.mat']);
            end
                
            %for nn = 1:ntrials_Adapt2
            nn = 1;
            while (nn <= ntrials_Adapt2)
                Omega2 = zeros(m2,1);
                for ii = 1:m2
                   Omega2(ii) = GenerateRandomIndex(S);
                end

                PhiHstar = FHstar([Omega1_stage1'; Omega2], :);   % Complete adaptive sensing matrix
                noiseTerm2 = randn(m2,1)*noisevar;                % Gaussian noise
                y2 = FHstar(Omega2,:)*a + noiseTerm2;             % a is not exactly sparse
                y = [y1_stage1; y2];                              % Complete observation vector

                % Perform Reconstruction: Algorithm + Pseudoinverse
                switch opt_method
                    case 'cosamp'
                        tic
                        [ahat_cosamp_Adapt, suppRec, residual] = cosamp(PhiHstar, y, k, a);
                        %Xhat_opt_Adapt = ifft2(reshape(FHstar*ahat_cosamp_Full,n1,n2))*sqrt(n1*n2);
                        toc
                    case 'l1'
                        D = speye(n1*n2);
                        D = D([Omega1_stage1'; Omega2],:);
                        fwd = @(a) D*reshape(fft2(midwt(...
                            reshape(a,n1,n2), hh, nlevels)),n1*n2,1) ./ sqrt(n1*n2);
                        adj = @(yy) reshape(mdwt(ifft2(...
                            reshape(D'*yy,n1,n2)), hh, nlevels),n1*n2,1) .* sqrt(n1*n2);
                        Aop = Aopfast(fwd,adj);

                        if (noisevar > 0)
                            bpdn_sigma = norm([noiseTerm1_stage1; noiseTerm2]) + norm(PhiHstar*a - PhiHstar*ahat_bestk(kk).ahat);
                        else
                            error('Need to decide how to define bpdn_sigma when there is no noise');
                            bpdn_sigma = 1;
                        end
                        spgl1_opts = spgSetParms('verbosity', 0);

                        fprintf('starting spgl1\n');
                        tic
                        ahat_l1_Adapt = spg_bpdn(Aop, y, bpdn_sigma, spgl1_opts);
                        toc

                        suppRec = find(abs(ahat_l1_Adapt)>0); 

                    otherwise 
                        error('invalid method');
                end

                % Psuedoinverse for adaptive reconstruction
                % Observed that sometimes pinv fails with the error "SVD
                % failed to converge." Repeating the trial if this happens
                try
                    atemp = pinv(PhiHstar(:,suppRec))*y;
                catch
                    disp('pinv failed - Maybe SVD failed to converge. Repeating this trial.');
                    nn = nn - 1;
                    continue;
                end
                ahat_Final = zeros(length(a),1);
                ahat_Final(suppRec) = atemp;

                % Estimated image, PSNR, and squared error
                Xhat_Adapt{kk,mm,nn}   = ifft2(reshape(FHstar*ahat_Final,n1,n2))*sqrt(n1*n2);
                psnr_Adapt(kk,mm,nn)   = psnr(Xhat_Adapt{kk,mm,nn}, X);
                error2_Adapt(kk,mm,nn) = norm(X-Xhat_Adapt{kk,mm,nn})^2;
                
                disp(['Adaptive m ' num2str(m) ' trial ' num2str(nn) ' sparsity k ' num2str(k)]);              
                nn = nn + 1;
            end
            if (saveMat)
               filenameMat = sprintf('AdaptiveImages_%dn_%dtrials_%1.4fsigma_%sopt_%dk1_%dkn_%dm1_%dmn',n1,ntrials_Adapt2,noisevar,opt_method,klist(1),klist(end),mlist(1),mlist(end));
               save([filenameMat '.mat'],'Xhat_Adapt','psnr_Adapt','error2_Adapt');
               %clear Xhat_Adapt;
            end
        end
    end

end

%% Save and Plot Results

if (saveMat)
   filenameMat = sprintf('Brain_Tester_%dn_%dtrialsNonadapt_%dtrialsAdapt1_%dtrialsAdapt2_%1.4fsigma_%sopt_%dk1_%dkn_%dm1_%dmn',n1,ntrials_Nonadapt,ntrials_Adapt1,ntrials_Adapt2,noisevar,opt_method,klist(1),klist(end),mlist(1),mlist(end));
   save([filenameMat '.mat']);
end

% Plot Original Image
figure; clf;
imagesc(abs(X));
title(['Original ' num2str(n1) ' by ' num2str(n2) ' Brain Image'],'fontsize',20);

% Plot Nonadaptive Results
if (run_Nonadapt)
    
    for kk = 1:length(klist)
    
        % Plot Average number of nonzeros
        figure;
        set(gcf,'WindowStyle','normal')
        plot(mlist, mean(squeeze(nnz_Nonadapt(kk,:,:)),2), '-ro','linewidth',2)
        ylabel('Mean Number of Nonzeros','fontsize',22,'FontName','Helvetica');
        xlabel('Number of Measurements','fontsize',22,'FontName','Helvetica');
        set(gca,'FontSize',20);
        grid on;
        title(['Nonadaptive, k=' num2str(klist(kk))],'fontsize',20);
        hold off
        set(gcf,'PaperPositionMode','Auto');

        % Plot PSNR
        figure;
        set(gcf,'WindowStyle','normal')
        plot(mlist, mean(squeeze(psnr_Nonadapt(kk,:,:)),2), '-ro','linewidth',2); hold on
        plot(mlist, median(squeeze(psnr_Nonadapt(kk,:,:)),2), '--r+','linewidth',2);
        hl = legend('Mean PSNR','Median PSNR');
        set(hl,'fontsize',20);
        ylabel('Mean or Median PSNR','fontsize',22,'FontName','Helvetica');
        xlabel('Number of Measurements','fontsize',22,'FontName','Helvetica');
        set(gca,'FontSize',20);
        grid on;
        title(['Nonadaptive, k=' num2str(klist(kk))],'fontsize',20);
        hold off
        set(gcf,'PaperPositionMode','Auto'); 

        % Plot MSE
        figure;
        set(gcf,'WindowStyle','normal')
        semilogy(mlist, mean(squeeze(error2_Nonadapt(kk,:,:)),2), '-ro','linewidth',2); hold on
        semilogy(mlist, median(squeeze(error2_Nonadapt(kk,:,:)),2), '--r+','linewidth',2);
        hl = legend('Mean Squared Error','Median Squared Error');
        set(hl,'fontsize',20);
        ylabel('Mean or Median Squared Error','fontsize',22,'FontName','Helvetica');
        xlabel('Number of Measurements','fontsize',22,'FontName','Helvetica');
        set(gca,'FontSize',20);
        grid on;
        title(['Nonadaptive, k=' num2str(klist(kk))],'fontsize',20);
        hold off
        set(gcf,'PaperPositionMode','Auto'); 
        
    end
    
end

% Plot Adaptive Results
if (run_Adapt)
    
    for kk = 1:length(klist)
    
        % Plot PSNR
        figure;
        set(gcf,'WindowStyle','normal')
        plot(mlist, mean(squeeze(psnr_Adapt(kk,:,:)),2), '-bo','linewidth',2); hold on
        plot(mlist, median(squeeze(psnr_Adapt(kk,:,:)),2), '--b+','linewidth',2);
        hl = legend('Mean PSNR','Median PSNR');
        set(hl,'fontsize',20);
        ylabel('Mean or Median PSNR','fontsize',22,'FontName','Helvetica');
        xlabel('Number of Measurements','fontsize',22,'FontName','Helvetica');
        set(gca,'FontSize',20);
        grid on;
        title(['Adaptive, k=' num2str(klist(kk))],'fontsize',20);
        hold off
        set(gcf,'PaperPositionMode','Auto'); 

        % Plot MSE
        figure;
        set(gcf,'WindowStyle','normal')
        semilogy(mlist, mean(squeeze(error2_Adapt(kk,:,:)),2), '-bo','linewidth',2); hold on
        semilogy(mlist, median(squeeze(error2_Adapt(kk,:,:)),2), '--b+','linewidth',2);
        hl = legend('Mean Squared Error','Median Squared Error');
        set(hl,'fontsize',20);
        ylabel('Mean or Median Squared Error','fontsize',22,'FontName','Helvetica');
        xlabel('Number of Measurements','fontsize',22,'FontName','Helvetica');
        set(gca,'FontSize',20);
        grid on;
        title(['Adaptive, k=' num2str(klist(kk))],'fontsize',20);
        hold off
        set(gcf,'PaperPositionMode','Auto'); 
    
    end

end

if (plotFinalImages.flag)
    % To plot the nonadaptive and adaptive images corresponding to selected
    % levels of m and k
    mm = plotFinalImages.mm;
    kk = plotFinalImages.kk;

    % Plot Nonadaptive
    mean_error2_Nonadapt = mean(squeeze(error2_Nonadapt(kk,mm,:)));
    [junk imageIdx_Nonadapt] = min(squeeze(error2_Nonadapt(kk,mm,:)) - mean_error2_Nonadapt);
    figure; 
    set(gcf,'WindowStyle','normal')
    imagesc(abs(Xhat_Nonadapt{kk,mm,imageIdx_Nonadapt}));
    set(gcf,'visible','on');
    set(gca,'FontSize',20);
    set(gcf,'PaperPositionMode','Auto'); 
    title({[sprintf('Nonadaptive VDS Recovery.  PSNR = %2.2f dB, MSE = %1.4f', psnr_Nonadapt{kk,mm,imageIdx_Nonadapt}, error2_Nonadapt{kk,mm,imageIdx_Nonadapt})];...
        ['n1 = ' num2str(n1) ', n2 = ' num2str(n2) ', k = ' num2str(klist(kk)) ',  m = ' num2str(mlist(mm)) ', sigma = ' num2str(noisevar)]});

    % Plot Adaptive
    mean_error2_Adapt = mean(squeeze(error2_Adapt(kk,mm,:)));
    [junk imageIdx_Adapt] = min(squeeze(error2_Adapt(kk,mm,:)) - mean_error2_Adapt);
    figure;
    set(gcf,'WindowStyle','normal')
    imagesc(abs(Xhat_opt_Adapt));
    set(h,'visible','on');
    set(gca,'FontSize',20);
    set(gcf,'PaperPositionMode','Auto');
    title({[sprintf('Adaptive Recovery.  PSNR = %2.2f dB, MSE = %1.4f', psnr_Adapt{kk,mm,imageIdx_Adapt}, error2_Adapt{kk,mm,imageIdx_Adapt})];...
        ['n1 = ' num2str(n1) ', n2 = ' num2str(n2) ', k = ' num2str(klist(kk)) ', m = ' num2str(mlist(mm)) ', sigma = ' num2str(noisevar) ', k=' num2str(k)]});

end
