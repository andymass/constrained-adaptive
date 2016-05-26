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

function [xhat str] = runNonAdaptive(x, m, noisevar, optAlg, variableDensity, preconditioning, v, d, globaldata)

% x : n-dimensional SPARSE signal (often called alpha)
% m : how many measurements to use in total
% noisevar : variance of noise to add to measurements
% optAlg : 1 - Use CoSaMP to reconstruct
% v : variable density sampling probabilities array
% d : variable desnity smapling conditioning parameter array

global F Hstar;

% sometimes parfor ignores globals; you can set them explicitely here
if nargin > 8 && ~isempty(globaldata)   
    F = globaldata.F;
    Hstar = globaldata.Hstar;
else
    globaldata = [];
end

QUIET = 0;

totargs = 4;
if nargin < totargs
    optAlg = 1;
end
if nargin < totargs - 1
    noisevar = 0;
end

k = nnz(x);
n = length(x);
truesupp = find(abs(x)>0);

if (variableDensity)
    Omega = zeros(m,1);
    for ii = 1:m
       Omega(ii) = GenerateRandomIndex(v); 
    end
    if (preconditioning)
        D = diag(d(Omega));
        Phi = 1/sqrt(m)*D*F(Omega, :);
        %y = Phi*Hstar*x + 1/sqrt(m)*D*randn(m,1)*noisevar;
        noiseTerm = 1/sqrt(m)*D*randn(m,1)*noisevar
        y = Phi*Hstar(:,truesupp)*x(truesupp) + noiseTerm; % a little faster
    else
        Phi = F(Omega,:);
        %y = Phi*Hstar*x + randn(m,1)*noisevar;
        noiseTerm = randn(m,1)*noisevar;
        y = Phi*Hstar(:,truesupp)*x(truesupp) + noiseTerm; % a little faster
    end
else
    % uniformly random rows
    rows = randperm(n);
    rows = rows(1:m);
    Phi = F(rows, :);
    %y = Phi*Hstar*x + randn(m,1)*noisevar;
    noiseTerm = randn(m,1)*noisevar;
    y = Phi*Hstar(:,truesupp)*x(truesupp) + noiseTerm; % a little faster
end

switch optAlg
    case 1 % CoSaMP
        if ~QUIET
            display(['Running CoSaMP for NonAdaptive reconstruction']);
        end
        xhat = cosamp(Phi*Hstar, y, k, x); 
        if ~QUIET
            display('Done.');
        end
        str = ['CoSaMP for Nonadapt rec'];
    case 2 % OMP
        if ~QUIET
            display(['Running OMP for NonAdaptive reconstruction']);
        end
        xhat = singleMatrixOMP(y, Phi*Hstar, k, []);
        if ~QUIET
            display('Done.');
        end
        str = ['OMP for Nonadapt rec'];
    case 3 % L1-Minimization
        if ~QUIET
            display(['Running L1-Minimization (SPGL1) for NonAdaptive reconstruction']);
        end
        %spgl1_opts = spgSetParms('verbosity',0,'optTol',1e-10);
        spgl1_opts = spgSetParms('verbosity',0);
        xhat = spg_bpdn(Phi*Hstar, y, norm(noiseTerm), spgl1_opts);
        if ~QUIET
            display('Done.');
        end
        str = ['L1-Minimization for Nonadapt rec'];
        
end

if ~QUIET
    display(['NonAdaptive Recovery error is : ', num2str(norm(x-xhat))]);
end



