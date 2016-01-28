function [distErr_nBack,closestObj_nBack] = YC_nBack_err(events,trial_pres_mask)
% function [distErr_nBack,closestObj_nBack] = YC_nBack_err(events,trial_pres_mask);
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
    distErr_nBack    = NaN(1,4);
    closestObj_nBack = NaN(1,4);
    return
end

% vectorize
sessVec     = [test_ev.session];
objLocs     = vertcat(test_ev.objLocs);
respLocs    = vertcat(test_ev.respLocs);
sessions    = unique(sessVec);

% loop over each session
distErr_nBack_subj        = [];
closestErr_nBack_subj     = zeros(1,4);
closestErr_nBackNorm_subj = zeros(1,4);

for sess = 1:length(sessions)
    sessInds = sessVec == sessions(sess);
    
    % filter to session
    mask_sess     = trial_pres_mask(sessInds);
    respLocs_sess = respLocs(sessInds,:);
    objLocs_sess  = objLocs(sessInds,:);
    
    % loop over each list. maybe there is a better way but it's
    % unclear. calc error between response location on current list and
    % item locations on previous lists
    distErr_nBack_sess = [];
    for item = 4:sum(sessInds)
        
        % only compute for this trial if in mask
        if mask_sess(item)
            % calc error for current item and response locations, as wells
            % current response location and previous three item locations
            distErr_nBack_trial        = NaN(1,4);
            closestObj_nBack_trial     = NaN(1,4);
            for nBack = 0:3
                thisObj_loc  = objLocs_sess(item-nBack,:);
                thisResp_loc = respLocs_sess(item,:);
                distErr_nBack_trial(nBack+1) = 1-calc_YC_error(thisObj_loc,thisResp_loc);
                closestObj_nBack_trial(nBack+1) = sqrt(sum([thisObj_loc - thisResp_loc].^2));
            end
            distErr_nBack_sess = [distErr_nBack_sess;distErr_nBack_trial];
            [m,i] = min(closestObj_nBack_trial);
            closestErr_nBack_subj(i) = closestErr_nBack_subj(i) + 1;
            [m,i] = max(distErr_nBack_trial);
            closestErr_nBackNorm_subj(i) = closestErr_nBackNorm_subj(i) + 1;
        end
    end
    distErr_nBack_subj = [distErr_nBack_subj;distErr_nBack_sess];
end
distErr_nBack = mean(distErr_nBack_subj,1);
closestObj_nBack = closestErr_nBack_subj / sum(closestErr_nBack_subj);
closestObjNorm_nBack = closestErr_nBackNorm_subj / sum(closestErr_nBackNorm_subj);