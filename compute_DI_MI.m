function  [MI, DI, rev_DI]=compute_DI_MI(X,Y,Nx,D,alg,start_ratio,prob,flag)
% Function `compute_DI_MI' calculates the directed information I(X^n-->
% Y^n), mutual information I(X^n; Y^n) and reverse directed information I(Y^{n-1}-->X^n)
% for any positive integer n smaller than the length of X and Y. 

% X and Y: two input sequences;
% Nx:  the size of alphabet of X, assuming X and Y have the same size of
% alphabets;
% D:  the maximum depth of the context tree used in basic CTW algorithm,
% for references please see F. Willems, Y. Shtarkov and T. Tjalkens, 'The
% Context-Tree Weighting Method: Basic Properties', IEEE Transactions on
% Information Theory, 653-664, May 1995.
% alg:  indicates one of the four possible estimators proposed in J.
% Jiao. H. Permuter, L. Zhao, Y.-H. Kim and T. Weissman, 'Universal
% Estimation of Directed Information', http://arxiv.org/abs/1201.2334.
% Users can indicate strings 'E1','E2','E3' and 'E4' for corresponding
% estimators.
% start_ratio: indicates how large initial proportion of input data should be ignored when displaying
% the estimated results, for example, if start_ratio = 0.2, then the output DI
% only contains the estimate of I(X^n \to Y^n) for n larger than
% length(X)/5.

n_data=length(X);
% mapp the data pair (X,Y) into a single variable taking value with
% alphabet size |X||Y|
XY=X+Nx*Y;

if flag == 0
%Calculate the CTW probability assignment
    pxy = ctwalgorithm(XY,Nx^2,D);
    px = ctwalgorithm(X,Nx,D);
    py = ctwalgorithm(Y,Nx,D);
else 
    pxy = prob.pxy;
    px = prob.px;
    py = prob.py;
end

% px_xy is a Nx times n_data matrix, calculating p(x_i|x^{i-1},y^{i-1})
for i_x=1:Nx
    px_xy(i_x,:)=pxy(i_x,:);
    for j=2:Nx
        px_xy(i_x,:)=px_xy(i_x,:)+pxy(i_x+(j-1)*Nx,:); 
    end;
end;

%calculate P(y|x,X^{i-1},Y^{i-1})
temp= repmat(px_xy,Nx,1);
py_x_xy=pxy./temp;

if strcmp(alg,'E1')% use the trick that px(2,1)=px(2), px(2,2)=px(2+Nx)

temp_MI=-log2(px(X(D+1:end)+[1:Nx:end-Nx+1])) - log2(py(Y(D+1:end)+[1:Nx:end-Nx+1]))+  log2(pxy(XY(D+1:end)+[1:Nx^2:end-Nx^2+1]))  ;
temp_DI= - log2(py(Y(D+1:end)+[1:Nx:end-Nx+1]))+  log2(pxy(XY(D+1:end)+[1:Nx^2:end-Nx^2+1])) -log2(px_xy(X(D+1:end)+[1:Nx:end-Nx+1]))  ;
temp_rev_DI=-log2(px(X(D+1:end)+[1:Nx:end-Nx+1]))+log2(px_xy(X(D+1:end)+[1:Nx:end-Nx+1])) ;

elseif strcmp(alg,'E2')
temp_MI=ctwentropy(px)+ctwentropy(py)-ctwentropy(pxy);
temp_DI=ctwentropy(py)-ctwentropy(pxy)+ctwentropy(px_xy);
temp_rev_DI=ctwentropy(px)-ctwentropy(px_xy);

elseif strcmp(alg,'E3')
    temp_MI=zeros(1,size(px,2));
    temp_DI= temp_MI;
    temp_rev_DI=temp_MI;
    for iy=1:Nx
           temp_MI=temp_MI+py_x_xy(X(D+1:end)+(iy-1)*Nx+[1:Nx^2:end-Nx^2+1] ).*log2(pxy(X(D+1:end)+(iy-1)*Nx+[1:Nx^2:end-Nx^2+1] )./(py(iy,:).*px(X(D+1:end)+[1:Nx:end-Nx+1])));
            temp_DI=temp_DI+py_x_xy(X(D+1:end)+(iy-1)*Nx+[1:Nx^2:end-Nx^2+1] ).*log2(py_x_xy(X(D+1:end)+(iy-1)*Nx+[1:Nx^2:end-Nx^2+1] )./(py(iy,:)));
            temp_rev_DI=temp_rev_DI+px_xy(iy,: ).*log2(px_xy(iy,: )./px(iy,:));
    end;


elseif strcmp(alg,'E4')
    temp_MI=zeros(1,size(px,2));
    temp_DI= temp_MI;
    temp_rev_DI=temp_MI;
    for iy=1:Nx
       for ix=1:Nx
           temp_MI=temp_MI+pxy(ix+(iy-1)*Nx,:).*log2(pxy(ix+(iy-1)*Nx,:)./(py(iy,:).*px(ix,:)));
           temp_DI=temp_DI+pxy(ix+(iy-1)*Nx,:).*log2(pxy(ix+(iy-1)*Nx,:)./(py(iy,:).*px_xy(ix,:)));
           temp_rev_DI=temp_rev_DI+pxy(ix+(iy-1)*Nx,:).*log2(px_xy(ix,:)./px(ix,:));
       end;
    end;
end;

%MI= cumsum(temp_MI((floor(n_data*start_ratio)+1):end));
DI= cumsum(temp_DI((floor(n_data*start_ratio)+1):end));
rev_DI=cumsum(temp_rev_DI((floor(n_data*start_ratio)+1):end));
MI = DI+rev_DI;
