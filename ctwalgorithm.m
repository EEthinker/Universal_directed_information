function Px_record = ctwalgorithm(x,Nx,D)
% Function CTWAlgorithm outputs the universal sequential probability
% assignments given by CTW method.
if size(x,1) ~= 1 
    error('The input vector must be a colum vector!');
end

n=length(x);
countTree = zeros(Nx, (Nx^(D+1) - 1) / (Nx-1)) ;
betaTree = ones(1,(Nx^(D+1) - 1 )/ (Nx-1))  ;
Px_record = zeros(Nx,n-D);
indexweight = Nx.^[0:D-1];
offset = (Nx^(D) - 1) / (Nx-1) + 1;

for i=D+1:n,
    context = x(i-D:i-1);
    leafindex = context*indexweight'+offset; 
    xt = x(i);
    eta = (countTree(1:Nx-1,leafindex)'+0.5)/(countTree(Nx,leafindex)+0.5);
    % update the leaf
    countTree(xt+1,leafindex) = countTree(xt+1,leafindex) + 1;    
    node =floor((leafindex+Nx-2)/Nx);  
    while ( node ~=0)
        [countTree, betaTree, eta] = ctwupdate(countTree,betaTree, eta, node, xt,1/2) ;
        node =floor((node+Nx-2)/Nx);
    end    
    eta_sum = sum(eta)+1;
    Px_record(:,i-D) = [eta 1]'/eta_sum ;
end
