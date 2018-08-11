function [ RES ] = util_test_eq( X, Y, tol )
%util_test_eq( X, Y, tol )
% test if two matrices are equal, within tolerance
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

if sum(size(X) ~= size(Y))
    error('input dimensions must match');
end

if sum(sum(abs(X-Y)>tol))
    RES = 0;
else
    RES = 1;
end
