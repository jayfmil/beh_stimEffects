function run_YC1_nBack(subjs)
%

if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end


% n subjs by 4 (same item, 1 back, 2 back, 3 back)
distErr_nBack        = NaN(length(subjs),4);
closestObj_nBack     = NaN(length(subjs),4);

% loop over each subject
for s = 1:length(subjs)    
    subj = subjs{s};
    fprintf('Processing %s.\n',subj);
    
    % load events
    events = get_sub_events('RAM_YC1',subj);
    
    % compute for each subject
    [distErr_nBack_subj,closestObj_nBack_subj] = YC_nBack_err(events);
    distErr_nBack(s,:) = distErr_nBack_subj;
    closestObj_nBack(s,:) = closestObj_nBack_subj;
    
end

keyboard