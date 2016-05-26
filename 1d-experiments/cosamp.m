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

function[xhat, Last_Set, v] = cosamp(Phi, u, s, x)

% Phi : measurement matrix
% u : measurement vector
% s : sparsity level
% x : real solution (for experimental purposes)

                PRINTOUT = 0;
                
                d = size(Phi,2);
                Last_Set = [];
                %Initialize
                a = zeros(d,1);
                v = u;
                it=0;
                stop = 0;
                while ~stop
                    
                    %Signal Proxy
                    y = Phi'*v;
                    [tmp, ix] = sort(abs(y), 'descend');
                    Omega = ix(1:2*s); 
                    [tmp, ix] = sort(abs(a), 'descend');
                    T = union(Omega, Last_Set);
                    if PRINTOUT
                        display(['Current index set is:']);
                        T
                        pause
                    end
                    
                    %Signal Estimation
                    b = zeros(d, 1);
                    b(T) = Phi(:, T) \ u;
                    
                    %Prune
                    [tmp, ix] = sort(abs(b), 'descend');
                    a = zeros(d, 1);
                    a(ix(1:s), 1) = b(ix(1:s), 1);
                    Last_Set = ix(1:s);
                    if PRINTOUT
                        display('Pruned index set is:');
                        ix(1:s)
                        display('Current estimation on this set:');
                        a(ix(1:s))
                        pause
                    end
                    
                    %Sample Update
                    %Edit: This is the residual, edit to return this variable 
                    v = u - Phi*a;
                    
                    %Iteration counter
                    it = it + 1;
                    
                    %Check Halting Condition (change to residual check if
                    %desired)
                    if norm(a-x) <= 10^(-4) || it > max(8*s, 100)
                        stop = 1;
                    end
                    
                end %End CoSaMP iteration
                xhat = a;
