function [logOR,OR,SE] = computeLogOddsRatio(x1,n1,x2,n2)
% function logOdds = computeLogOddsRatio(x1,n1,x2,n2)
% 
% Compute log odds ratio number of successes in group 1, total n of group
% 1, number of successes in group 2, total n of group 2.

% odds = (((numData_crpStim{i}(:,:,r)+.5) ./ (totMinusActualStim{i}(:,:,r)+.5)) ./ ((numData_crpNonStim{i}(:,:,r)+.5) ./ (totMinusActualNonStim{i}(:,:,r)+.5)));


n1(n1==0) = NaN;
n2(n2==0) = NaN;
x12 = n1-x1;
x22 = n2-x2;
% x1  = x1+.5;
% x2  = x2+.5;
% x12 = x12+.5;
% x22 = x22+.5;

OR    = (x1./x12) ./ (x2./x22);
logOR = log(OR);
SE    = sqrt(1./x1 + 1./x12 + 1./x2 + 1./x22);
