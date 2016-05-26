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

function [xhat, supp1, optErr, str, success] = runAdaptive(x, m1, m2, noisevar, optAlg1, opt, optRound, optAlg2, variableDensity, resThresh, preconditioning, v, d, globaldata)

% x : n-dimensional SPARSE signal (often called alpha)
% m1 : how many measurements to use in support identification
% m2 : how many additional (adaptive) measurements to use for reconstruction
% noisevar : variance of noise to add to measurements
% optAlg1 : 1 - Run CoSaMP to locate support
%           2 - Run OMP to locate support
% optAlg2 : 1 - If opt=2, run CoSaMP to reconstruct signal
%           2 - If opt=2, run OMP to reconstruct signal
% opt     : 1 - Run algorithm 1 to locate support, then use pseudo-inverse to
%               reconstruct
%           2 - Run algorithm 1 to locate support, then algorithm 2 to
%               reconstruct
% optRound: 1 - Simple rounding
%           2 - S Probability
%           3 - Preconditioning S Probability
%           4 - Water emptying rounding
% variableDensity - variable density sampling in nonadaptive measurements flag (1 = yes, 0 = no)
% resThresh: In support estimation step, if norm(residual) is greater than
% resThresh, try again with new nonadaptive measurements
% v : variable density sampling probabilities array
% d : variable desnity smapling conditioning parameter array
%
% xhat : estimate to x
% supp1 : support found by initial support identification
% success: if 0, need to try new signal - TFOCS failed somewhere

global F Hstar;

% sometimes parfor ignores globals; you can set them explicitely here
if nargin > 13 && ~isempty(globaldata)   
    F = globaldata.F;
    Hstar = globaldata.Hstar;
else
    globaldata = [];
end

QUIET = 0;

totargs = 8;
if nargin < totargs
    optAlg2 = 1;
end
if nargin < totargs - 1
    optRound = 1;
end
if nargin < totargs - 2
    opt = 1;
end
if nargin < totargs - 3
    optAlg1 = 1;
end
if nargin < totargs - 4
    noisevar = 0;
end
if nargin < totargs - 5
    display('ERROR: Must define parameters m and x');
    return;
end

k = nnz(x);
n = length(x);

%str = ['Adaptive: '];
str = [];

% To test against the oracle case, when we know support
truesupp = find(abs(x)>0);
if k ~= length(truesupp)
    display('Warning vector not actually k-sparse');
end
if ~QUIET
    display('Running adaptive rows for Oracle adaptive case');
end
%[truePhi, S, success] = adaptiveRows(m2+m1, truesupp, n, optRound);
[truePhi, S, success] = adaptiveRows(m2+m1, truesupp, n, optRound, globaldata);
if (~success)
    display('adaptiveRows did not succeed in oracle case - need to try new signal');
    xhat = []; supp1 = []; optErr = [];
    return;
end
S = real(S);
H_truesupp = Hstar(:, truesupp);
ytrue = truePhi*H_truesupp*x(truesupp) + randn(m1+m2, 1)*noisevar; % a little faster

residual = resThresh + 1;
count = 0;
while (norm(residual) > resThresh && count < 100)
    % Choose first m1 rows nonadaptively
    if (variableDensity)
        Omega = zeros(m1,1);
        for ii = 1:m1
           Omega(ii) = GenerateRandomIndex(v); 
        end
        if (preconditioning)
            D    = diag(d(Omega));
            Phi1 = 1/sqrt(m1)*D*F(Omega, :);
            %y1 = Phi1*Hstar*x + 1/sqrt(m1)*D*randn(m1,1)*noisevar;
            noiseTerm1 = 1/sqrt(m1)*D*randn(m1,1)*noisevar;
            y1 = Phi1*Hstar(:,truesupp)*x(truesupp) + noiseTerm1; % a little faster
        else
            Phi1 = F(Omega,:);
            %y1 = Phi1*Hstar*x + randn(m1,1)*noisevar;
            noiseTerm1 = randn(m1,1)*noisevar;
            y1 = Phi1*Hstar(:,truesupp)*x(truesupp) + noiseTerm1; % a little faster
        end
    else
        rows1 = randperm(n);
        rows1 = rows1(1:m1);
        Phi1  = F(rows1, :);
        %y1    = Phi1*Hstar*x + randn(m1,1)*noisevar;
        noiseTerm1 = randn(m1,1)*noisevar;
        y1    = Phi1*Hstar(:,truesupp)*x(truesupp) + noiseTerm1; % a little faster
    end

    switch optAlg1
        case 1 % CoSaMP
            if ~QUIET
                display(['Running CoSaMP to find support using m=', num2str(m1)]);
            end
            [xhat1, supp1, residual] = cosamp(Phi1*Hstar, y1, k, x);
            if (count == 0)
                str = [str, 'Using CoSaMP for support, '];
            end
            if ~QUIET
                display('Done.');   
            end
        case 2 % OMP
            if ~QUIET
                display(['Running OMP to find support using m=', num2str(m1)]);
            end
            [xhat1, supp1, residual] = singleMatrixOMP(y1, Phi1*Hstar, k, []);
            if (count == 0)
               str = [str, 'Using OMP for support, '];
            end
             if ~QUIET
                display('Done.');   
             end     
        case 3 % L1-Minimization
            if ~QUIET
                display(['Running L1-Minimization (SPGL1) to find support using m=', num2str(m1)]);
            end
            spgl1_opts = spgSetParms('verbosity',0,'optTol',1e-10);
            xhat = spg_bpdn(Phi1*Hstar, y1, norm(noiseTerm1), spgl1_opts);
            [tmp idx] = sort(abs(xhat),'descend');
            supp1 = idx(1:k); % support of largest k entries (in magnitude)
            residual = resThresh-1; % Force loop on residual to end after 1 iteration. 
            if (count == 0)
               str = [str, 'Using L1-Min for support, '];
            end
            if ~QUIET
                display('Done.');
            end
    end
    count = count+1;
    if ~QUIET  
       if (norm(residual) <= resThresh)
           if (all(ismember(truesupp, supp1)))
               display('residual less than resThresh and correct support');
           else
               display('residual less than resThresh and NOT correct support');
           end
       else
           if (all(ismember(truesupp, supp1)))
               display('residual greater than resThresh and correct support, repeating support estimation step');
           else
               display('residual greater than resThresh and NOT correct support, repeating support estimation step');
           end
       end
       
       if (count >= 100)
          display('NOTE: After 100 tries, residual still not below resThresh!'); 
       end
    end
end

if ~QUIET
    display(['Running adaptive rows for adaptive case using m2=' num2str(m2)]);
end
[Phi2, S, success] = adaptiveRows(m2, supp1, n, optRound);
if (~success)
    display('adaptiveRows did not succeed in adaptive case - need to try new signal');
    xhat = []; optErr = [];
    return;
end
S = real(S);

if k==1  % Do a sanity check
    U = unique(Phi2, 'rows');
    if size(U, 1) > 1
        display('k=1 but picking more than one row');
    end
end
Phi = [Phi1; Phi2];
%y2 = Phi2*Hstar*x + randn(m2,1)*noisevar;
noiseTerm2 = randn(m2,1)*noisevar;
y2 = Phi2*Hstar(:,truesupp)*x(truesupp) + noiseTerm2; % a little faster
y = [y1; y2];
%H_supp = Hstar(:, supp1);
if length(supp1) > m1 + m2
    display('WARNING : Support size found is larger than m!');
end

% Reconstruction 
% For Oracle case, always use Pseudo-inverse
xtemptrue = pinv(truePhi*H_truesupp)*ytrue;
xhattrue = zeros(n,1);
xhattrue(truesupp) = xtemptrue;
optErr = norm(x - xhattrue);
if ~QUIET
    display(['Oracle Recovery error is : ', num2str(optErr)]);
end
% For purely adaptive, use psuedo-inverse, or run OMP or CoSaMP
switch opt
    case 1 % Use Pseudo-inverse
        % get new estimate of support, using all m1+m2 measurements
        switch optAlg1
            case 1 % CoSaMP
                if ~QUIET
                    display(['Running CoSaMP to find support for reconstruction using m=', num2str(m1+m2)]);
                end
                [xhatRec, suppRec, residual] = cosamp(Phi*Hstar, y, k, x);
                if ~QUIET
                    display('Done.');   
                end
            case 2 % OMP
                if ~QUIET
                    display(['Running OMP to find support for reconstruction using m=', num2str(m1+m2)]);
                end
                [xhatRec, suppRec, residual] = singleMatrixOMP(y, Phi*Hstar, k, []);
                 if ~QUIET
                    display('Done.');   
                 end      
            case 3 % L1-Minimization
                if ~QUIET
                    display(['Running L1-Minimization (SPGL1) to find support for reconstruction using m=', num2str(m1+m2)]);
                end
                spgl1_opts = spgSetParms('verbosity',0,'optTol',1e-10);
                xhatRec = spg_bpdn(Phi*Hstar, y, sqrt(norm(noiseTerm1)^2+norm(noiseTerm2)^2), spgl1_opts);
                [tmp idx] = sort(abs(xhatRec),'descend');
                suppRec = idx(1:k); % support of largest k entries (in magnitude)
                if ~QUIET
                    display('Done.');
                end
        end
        if ~QUIET
            if ( all(ismember(truesupp, suppRec)) )
               display('Correct support estimate in reconstruction step!');
            else
               display('Incorrect support estimate in reconstruction step!');
            end
        end
        suppRec = sort(suppRec); 
        H_supp = Hstar(:,suppRec);
        xtemp = pinv(Phi*H_supp)*y;
        xhat = zeros(n,1);
        xhat(suppRec) = xtemp;
        str = [str, 'Pseudo-inv for Adapt rec'];
    case 2 % Use algorithm 2
        switch optAlg2
            case 1 % CoSaMP
                if ~QUIET
                    display(['Running CoSaMP for reconstruction using m=', num2str(m1+m2)]);
                end
                xhat = cosamp(Phi*Hstar, y, k, x);
                str = [str, 'CoSaMP for Adapt rec'];
                if ~QUIET
                    display('Done.');
                end
            case 2 % OMP
                if ~QUIET
                    display(['Running OMP for reconstruction using m=', num2str(m1+m2)]);
                end
                xhat = singleMatrixOMP(y, Phi*Hstar, k, []);
                str = [str, 'OMP for Adapt rec'];  
                if ~QUIET
                    display('Done.');
                end
            case 3 % L1-Minimization
                if ~QUIET
                    display(['Running L1-Minimization (SPGL1) for reconstruction using m=', num2str(m1+m2)]);
                end
                spgl1_opts = spgSetParms('verbosity',0,'optTol',1e-10);
                xhat = spg_bpdn(Phi*Hstar, y, sqrt(norm(noiseTerm1)^2+norm(noiseTerm2)^2), spgl1_opts);
                str = [str, 'L1-Minimization for Adapt rec'];
                if ~QUIET
                    display('Done.');
                end
        end   
end

if ~QUIET
    display(['Adaptive Recovery error is : ', num2str(norm(x-xhat))]);
end

