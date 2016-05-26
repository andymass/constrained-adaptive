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

%% GenerateRandomIndex.m
function [r] = GenerateRandomIndex(p)
% Description: Given an array of probabilities, return a random index with
% the probabilities specified in the array.
%
% Usage: r = GenerateRandomIndex(p)
%  Inputs:
%       p = vector of probabilities
%  Outputs:
%       r = random index (integer) between 1 and length(p)

if ( sum(p) <= 1-1e-8 || sum(p) >= 1+1e-8 )
   error('ERROR: Not valid probability distribution, sum(p) not equal to 1');
   sum(p)
   return;
end

x = rand; % random number uniformly distributed between 0 and 1

currentSum = 0;
for ii = 1:length(p)
   if ( x < currentSum+p(ii) )
       r = ii;
       return;
   else
       currentSum = currentSum+p(ii);
   end
end


end






