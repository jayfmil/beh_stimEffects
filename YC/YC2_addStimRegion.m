function events = YC2_addStimRegion(subj,events);

% this is stupid
if strcmp(subj,'R1047D')
    [events.stimAnodeTag] = deal('LOTD3');
    [events.stimCathodeTag] = deal('LOTD4');
end

% get bipolar electrode pairs for subject
tal     = getBipolarSubjElecs(subj,1);
if ~isfield(tal,'locTag'); [tal.locTag] = deal(''); end
monoTal = getBipolarSubjElecs(subj,0);
if ~isfield(monoTal,'locTag'); [monoTal.locTag] = deal(''); end

% figure out unique stimulations sites (prefer using tag string instead of
% electrode number, but not all events have the string)
if isfield(events,'stimAnodeTag')
    doTag = 1;
    sessStrs = cellfun(@num2str, num2cell([events.session]), 'UniformOutput', false)';
    an_cath_sess = horzcat({events.stimAnodeTag}', {events.stimCathodeTag}',sessStrs);
    an_cath_sess = an_cath_sess(sum(cellfun('isempty',an_cath_sess),2) == 0,:);
    
    % reduce to unqiue sites
    [~,i] = unique(cell2mat(an_cath_sess),'rows');
    stim_sess = an_cath_sess(i,:);
else
    doTag = 0;
    an_cath_sess = [[events.stimAnode]' [events.stimCathode]' [events.session]'];
    an_cath_sess = an_cath_sess(sum(isnan(an_cath_sess),2) == 0,:);
    
    % reduce to unqiue sites
    stim_sess = unique(an_cath_sess,'rows');
    stim_sess(sum(stim_sess(:,1:2),2) == 0,:) = [];    
end

% loop over each stim session in events
[events.stimRegion] = deal('');
[events.stimElec1] = deal('');
[events.stimElec2] = deal('');
for sess = 1:size(stim_sess,1)
   
    if ~doTag
        chans       = sort(stim_sess(sess,1:2));    
        stimElecInd = sum(ismember(vertcat(tal.channel),chans),2) == 2;
        stimElec1   = [monoTal.channel] == chans(1);
        stimElec2   = [monoTal.channel] == chans(2);
        sessInds    = [events.session]==stim_sess(sess,3);
    else
        thisTag1    = [stim_sess{sess,1},'-',stim_sess{sess,2}];
        thisTag2    = [stim_sess{sess,2},'-',stim_sess{sess,1}];
        stimElecInd = strcmpi(vertcat({tal.tagName}),thisTag1) | strcmpi(vertcat({tal.tagName}),thisTag2); 
        stimElec1   = strcmpi(vertcat({monoTal.tagName}),stim_sess{sess,1});
        stimElec2   = strcmpi(vertcat({monoTal.tagName}),stim_sess{sess,2});
        sessInds    = [events.session]==str2num(stim_sess{sess,3});        
    end
    
    % WHAT IS WITH THE EVENTS BEING SO POORLY LABELED
    stimLoc  = ''; if sum(stimElecInd) == 1; stimLoc = tal(stimElecInd).locTag; end
    stimLoc1 = ''; if sum(stimElec1) == 1; stimLoc1 = monoTal(stimElec1).locTag; end
    stimLoc2 = ''; if sum(stimElec2) == 1; stimLoc2 = monoTal(stimElec2).locTag; end
    if sum(stimElecInd) ~= 1 || ~isfield(tal,'locTag') || isempty(tal(stimElecInd).locTag) 
        fprintf('Manually adding stim location for %s\n',subj)       
        if strcmp(subj,'TJ082')            
            stimLoc = 'Left EC';                   
%         else
%             keyboard
        end
    end
    [events(sessInds).stimRegion] = deal(stimLoc);
    [events(sessInds).stimElec1] = deal(stimLoc1);
    [events(sessInds).stimElec2] = deal(stimLoc2);    
end
% if any(strcmp({events.stimRegion},''))
%     keyboard
% end

