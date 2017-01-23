function y=ctwentropy1D(x)

if (x<=10*eps | x>=1-10*eps) 
    y= 0;
else
    y= - x.*log2(x)-(1-x).*log2(1-x) ;
end
