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

function a = genSignal(n, k, values)

[f,~] = log2(abs(n));
if f ~= 0.5
    error('n must be a power of two');
end

if nargin < 3
    values = randn(k, 1);
end

a = zeros(1,n);
%parent = [0 1 2.^(ceil(log2(3:n))-2)+1]; % original parent , not binary tree
parent = [0 1 2 2 floor((2:n-3)/2)+2];    % Update: binary tree

a(1) = values(1);
for ii = 2:k
    I = find( a(2:end) == 0 & a( parent(2:end) ) ~= 0 ) + 1;
    jj = I(randi(end));
    a(jj) = values(ii);
    % fprintf('ii: %d, from %s, choose %d\n', ...
    %     ii, num2str(I), jj);
end

