function prob = ctwprob(X,Y,Nx,D)
%'ctwprob' calculates the universal sequential probability assignment given
%by basic CTW method

% mapp the data pair (X,Y) into a single variable taking value with
% alphabet size |X||Y|
XY=X+Nx*Y;

%Calculate the universal sequential probability assignment
prob.pxy = ctwalgorithm(XY,Nx^2,D);
prob.px = ctwalgorithm(X,Nx,D);
prob.py = ctwalgorithm(Y,Nx,D);