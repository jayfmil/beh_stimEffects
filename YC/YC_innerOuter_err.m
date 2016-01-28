function [errInner,errOuter] = YC_innerOuter_err(events,trial_pres_mask)
% function [errInner,errOuter] = YC_innerOuter_err(events,trial_pres_mask)
%
%

if ~exist('trial_pres_mask','var') || isempty(trial_pres_mask)
    trial_pres_mask = strcmp({events.type},'NAV_TEST');
end

events = addExtraYCFields(events);

% filter to test events
inds            = strcmp({events.type},'NAV_TEST');
test_ev         = events(inds);
trial_pres_mask = trial_pres_mask(inds);

if length(test_ev) < 30
    fprintf('Not enough events for %s. Skipping\n',test_ev(1).subject);
    errInner = NaN;
    errOuter = NaN;
    return
end

% vectorize
sessVec  = [test_ev.session];
isInner  = [test_ev.inner];
distErrs = 1-[test_ev.testError];
sessions = unique(sessVec);

% just call grpstats
[m,sdev,n] = grpstats(distErrs(trial_pres_mask==1),isInner(trial_pres_mask==1));
errOuter = m(1);
errInner = m(2);