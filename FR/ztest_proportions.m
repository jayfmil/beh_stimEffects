function [z,p] = ztest_proportions(x1,n1,x2,n2)
% function [z,p] = ztest_proportions(x1,n2,x2,n2)
%
%

n1(n1==0) = NaN;
n2(n2==0) = NaN;
pooledP = (x1+x2) ./ (n1+n2);
se      = abs(sqrt(pooledP .* (1-pooledP) .* ((1./n1)+(1./n2))));
se(se==0) = NaN;
z       = ((x1./n1)-(x2./n2)) ./ se;

% z = (x1/n1-x2/n2)/(sqrt(phat*(1-phat)*(1/n1+1/n2)));