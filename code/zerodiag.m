function new = zerodiag(A, d)

% This function will set the diagnoal of an adjacency matrix to an 
% arbitrary value, d. A is the adjacency matrix.
% default is d = 0;

if nargin < 2
    d = 0;
end

n = size(A,1);
new = A;
new(1:n+1:end) = d;