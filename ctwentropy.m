function y=ctwentropy(x)
% Function CTWEntropy outputs the (random) conditional entropy of the
% sequential probability assignment given by CTW method. 
% x is a matrix, where for all i, x(:, i) is a probability vector

% eliminate elements in x that are either too large or small
x( find ( x<=10*eps | x>=1-10*eps) ) = 10*eps ;

% y is the row vector containing estimates of entropy of the corresponding
% sequential probability assignments.
y = sum(-x.*log2(x),1);