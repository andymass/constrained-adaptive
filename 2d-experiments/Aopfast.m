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

function op = Aopfast(fwd, adj)

op = @(varargin) Aop_impl(fwd, adj, varargin{:});

end

function v = Aop_impl(fwd, adj, w, mode)
%   v = A(w,mode)   which returns  v = A *w  if mode == 1;
%                                  v = A'*w  if mode == 2. 
%
if mode == 1
    v = fwd(w);
elseif mode == 2
    v = adj(w);
else
    error('invalid mode');
end

end

