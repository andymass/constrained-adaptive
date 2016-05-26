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

function op = smooth_traceinv(q, C)

if nargin < 1, q = 1; end
if nargin < 2, C = []; end
if ~isreal(q) || q <= 0
    error('smooth_traceinv: First argument must be real and positive');
end

if isempty(C)
    op = @(varargin)smooth_traceinv_impl( q, varargin{:} );
else
    error('smooth_traceinv: C argument not yet supported');
end

function [ v, g ] = smooth_traceinv_impl( q, x, t )
if size(x,1) ~= size(x,2)
    error('smooth_traceinv_impl: input must be a square matrix');
end
switch nargin
    case 2
        if any(isinf(x(:)))
            error('inf');
        end
        if any(isnan(x(:)))
            error('nan');
        end
        d = eig(x);
        if any(abs(imag(d)) > 1e-15) || any(d <= 0)
            %warning('x must be positive definite');
            if nargout > 1
                error('can''t handle this!');
            end
            v = inf;
            return;
        end
        
        % the function is being used in a "smooth" fashion        
        if nargout > 1 
            % g = -x^(-2)
            invx = inv(x);
            g = -q*invx^2;    % Boyd, p120
            v = q*trace(invx);
            
            0;
        else
            % v = tr(inv(x));
           % d = eig(x);
            v = q*sum(1./d);
        
            0;
            
            % v = q*sum(1./diag(chol(x)).^2);
        end
        
    case 3
        error('smooth_traceinv_impl: non smooth not yet supported');

        % the function is being used in a "nonsmooth" fashion
        % i.e. return g = argmin_g  -q*log(det(g)) + 1/(2t)||g-x||^2
        [V,D]   = eig(x);
        d       = diag(D);
        if any(d<=0),
            v   = Inf;
            g   = nan(size(x));
            return;
%             error('log_det requires a positive definite point'); 
        end
        l       = ( d + sqrt( d.^2 + 4*t*q ) )/2;
        g       = V*diag(l)*V';
        v       = -q*sum(log(l));
    otherwise
        error('Wrong number of arguments');
end
