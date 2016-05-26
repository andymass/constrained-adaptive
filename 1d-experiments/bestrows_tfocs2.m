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

function [omega, S, out, success] = bestrows_tfocs2(G0, m, u, tau, K)

[n, k] = size(G0);

if nargin == 5      % If K is provided, use that
    herm_part = @(X) (X+X')/2;
    Af = @(S) herm_part( reshape( K*S, k, k) );
    At = @(Q) reshape( K'*Q(:), n, 1);
else                % No K: ``memory saving'' formulation
    Af = @(S) G0'*diag(S)*G0;
    At = @(Q) real( diag(G0*Q'*G0') );
%    Af = @(S) G0'*diag(abs(S) + 1e-12)*G0;
end

h = linop_handles({[n 1], [k k]}, Af, At, 'C2C'); 


% positive definite starting point
s0 = rand(n,1);
% s0 = ones(n, 1);
% s0 = sqrt(diag(G0*G0'));
s0 = s0 / norm(s0,1) * tau;

opts = [];
opts.maxIts = 600;
%opts.printEvery = 0;
opts.tol = 1e-15; % 1e-8;
opts.maxIts = 1000;

% opts.alg = 'GRA';

projectorF = proj_l1(tau);


 go = 1;
 countBad = 0;
 while (go && countBad < 100)
    [S, out] = tfocs(smooth_traceinv, { h, 0 }, ...
        projectorF, s0, opts);
%     if (strcmp(out.status,'Iteration limit reached'))
      if (out.niter < 5)
         disp('Running tfocs again because less than 5 iterations!');
        % get new positive definite starting point
         s0 = rand(n,1);
         s0 = s0 / norm(s0,1) * tau;
         countBad = countBad + 1;
         continue;
     else
         go = 0;
     end
 end
 
 if (countBad == 100)
     success = 0;
 else
     success = 1;
 end

[~, omega] = sort(abs(S), 'descend');
omega = omega(1:m);
