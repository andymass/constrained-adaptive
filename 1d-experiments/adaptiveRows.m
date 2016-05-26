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

function [Phi, S, success] = adaptiveRows(m, supp, n, opt, globaldata)

% m : Number of rows to return (integer)
% n : signal dimension (integer)
% supp : indices of the (estimated) signal support (k dimensional vector)
% Phi : m times n matrix containing adaptively chosen rows
% opt : 1 - Simple rounding
% opt : 2 - S Probability
% opt : 3 - Preconditioning S Probability (not implemented)
% opt : 4 - Water emptying rounding

global F Hstar;

% parfor basically ignores globals; set them explicitly    
if nargin > 4 && ~isempty(globaldata)   
    F = globaldata.F;
    Hstar = globaldata.Hstar;
else
    globaldata = [];
end

QUIET = 0;

G0 = F*Hstar(:,supp);

k = length(supp);
K = zeros(k*k, n);
for kk = 1:k
    for ell = 1:k
        K(kk+k*(ell-1), :) = G0(:,kk) .* conj(G0(:,ell));
        %K(kk+k*(ell-1), :) = G0(:,kk) .* (G0(:,ell));
    end
end

if ~QUIET
    display('Solving optimization problem...');
end

success = 0;

while (~success)
    [omega2, S, out, success] = bestrows_tfocs2(G0, 1, 0, 1, K);
    %[omega2, S] = bestrows_cvx(G0);
    if (~success)
       fprintf('\n ERROR: TFOCS didnt finish with 100 attempts\n');
       Phi = []; S = [];
       return;
    else
        success = 1;
        S = real(S);
        if ~QUIET
            display('TFOCS done');
        end
    end
end

switch opt
    case 1  % Simple rounding 
        if ~QUIET
            display('Using simple rounding');
        end
        Sorig= S*n;
        [Ssort_new I_new] = sort(Sorig,'descend');
        Omega = [];
        numLeft = m;
        while numLeft > 0
            num = round(Ssort_new(1));
            num = min(num, numLeft);
            Omega = [Omega, I_new(1)*ones(1,num)]; % Add these rows
            numLeft = numLeft - num;
            Sorig(I_new(1)) = 0;
            [Ssort_new I_new] = sort(Sorig,'descend');

        end
        Phi = F(Omega,:);
    case 2 % S Probabilities
        if ~QUIET
            display('Using S Probabilities');
        end
        % continue until k unique rows are selected
        kuniqueRows = 0;
        while(~kuniqueRows)
            Omega = zeros(m,1);
            for ii = 1:m
               Omega(ii) = GenerateRandomIndex(S); 
            end
            Phi = F(Omega,:);
            if (rank(Phi) < k)
                disp('Rank of Phi less than k, selecting rows again');
            else
                kuniqueRows = 1;
            end
        end
        
    case 3 % S Probabilities with Preconditioning
        if ~QUIET
            display('Using S Probabilities with Preconditioning');
        end
    case 4 % Water emptying
        if ~QUIET
            display('Using water emptying method');
        end
        [Ssort_new I_new] = sort(S*n,'descend');
        clear Omega;
        for ind = 1:m
           Omega(ind) = I_new(1);
           Ssort_new(1) = Ssort_new(1)-1;
           [Ssort_new I_new] = sort(Ssort_new(I_new),'descend');
        end
        Phi = F(Omega,:);
end
    
if (length(Omega) ~= m)
    error('Length of Omega is not m!');
end

end
