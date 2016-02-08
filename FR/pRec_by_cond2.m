function pRec_by_cond(subjs,figDir)


if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_FR2');
end
subjs(strcmp(subjs,'UT009a')) = [];


excludeOPs = [];

% will hold prob of recall results
% num lists/items x num conditions x num subjs x num regions
pRec_list      = NaN(25,4,length(subjs),5);
pRec_serialPos = NaN(12,4,length(subjs),5);
pRec_binned    = NaN(2,4,length(subjs),5);
pRec_binnedNum = NaN(2,4,length(subjs),5);
pRec_binnedDen = NaN(2,4,length(subjs),5);

pRec_FR1       = NaN(2,length(subjs));
pRec_FR1_num   = NaN(2,length(subjs));
pRec_FR1_den   = NaN(2,length(subjs));

pRecDiff_serialPos_prob    = NaN(12,length(subjs),5);
pRecDiff_serialPos_logOdds = NaN(12,length(subjs),5);
pRecDiff_serialPos_zProp   = NaN(12,length(subjs),5);

pRecDiff_binned_prob    = NaN(2,length(subjs),5);
pRecDiff_binned_logOdds = NaN(2,length(subjs),5);
pRecDiff_binned_zProp   = NaN(2,length(subjs),5);

pRec_stimItem    = NaN(length(subjs),5);
pRec_nonStimItem = NaN(length(subjs),5);

% will hold crp results
% lags x subjs x regions
% list level
crp_stimLists           = NaN(23,length(subjs),5);
crp_nonStimLists        = NaN(23,length(subjs),5);
crp_stimNonStim_logOdds = NaN(23,length(subjs),5);
crp_stimNonStim_zProp   = NaN(23,length(subjs),5);

% also will bin -2,-1 and 1,2
crp_stimLists_bin           = NaN(3,length(subjs),5);
crp_nonStimLists_bin        = NaN(3,length(subjs),5);
crp_stimNonStim_logOdds_bin = NaN(3,length(subjs),5);
crp_stimNonStim_zProp_bin   = NaN(3,length(subjs),5);

% same to same item type
crp_stimItems            = NaN(23,length(subjs),5);
crp_nonStimItems         = NaN(23,length(subjs),5);
crp_stimNonStimItems_logOdds = NaN(23,length(subjs),5);

% item type to any
crp_stimItems_toAny        = NaN(23,length(subjs),5);
crp_nonStimItems_toAny     = NaN(23,length(subjs),5);
crp_stimNonStimItems_toAny_logOdds = NaN(23,length(subjs),5);

% train crps
crp_stimLists_train      = NaN(11,length(subjs),5);
crp_nonStimLists_train   = NaN(11,length(subjs),5);

% temporal factor scores
temp_fact_stim           = NaN(length(subjs),5);
temp_fact_nonStim        = NaN(length(subjs),5);
temp_fact_signed_stim    = NaN(length(subjs),5);
temp_fact_signed_nonStim = NaN(length(subjs),5);

% intrusion results
nPlis_stimLists = NaN(length(subjs),5);
nPlis_nonStimLists = NaN(length(subjs),5);
nXlis_stimLists = NaN(length(subjs),5);
nXlis_nonStimLists = NaN(length(subjs),5);

nPlis_stimLists_1Back = NaN(length(subjs),5);
nPlis_nonStimLists_1Back = NaN(length(subjs),5);

% plis
plis_counts_per_list = NaN(2,3,length(subjs),5);


% loop over each subject
subjs_in_region = zeros(1,5);
stimEff_in_region = NaN(length(subjs),5);
for s = 1:length(subjs)
    subj = subjs{s};
    fprintf('Processing %s.\n',subj);
    
    % load events and and stim region
    events = get_sub_events('RAM_FR2',subj);
    try
        events = addRegionToEvents(subj,events);
    catch e
        fprintf('Warning: cannot find stim location for %s. Skipping.\n',subj)
        continue
    end
    
    % if subjects was stimulated in our regions of interest filter to
    % those events and process
    hipp_events = ~cellfun('isempty',regexpi({events.stimRegion},['CA1|CA2|CA3|DG|SUB']));
    hipp_events = hipp_events | ~cellfun('isempty',regexpi({events.stimElec1},['CA1|CA2|CA3|DG|SUB']));
    hipp_events = hipp_events | ~cellfun('isempty',regexpi({events.stimElec2},['CA1|CA2|CA3|DG|SUB']));
    
    ec_events = ~cellfun('isempty',regexpi({events.stimRegion},['EC']));
    ec_events = ec_events | ~cellfun('isempty',regexpi({events.stimElec1},['EC']));
    ec_events = ec_events | ~cellfun('isempty',regexpi({events.stimElec2},['EC']));
    
    prc_events = ~cellfun('isempty',regexpi({events.stimRegion},['PRC|PHC']));
    prc_events = prc_events | ~cellfun('isempty',regexpi({events.stimElec1},['PRC|PHC|PHG']));
    prc_events = prc_events | ~cellfun('isempty',regexpi({events.stimElec2},['PRC|PHC|PHG']));
    eventInds = {hipp_events,ec_events,prc_events,hipp_events|ec_events};
    
    [pRec_FR1(:,s),pRec_FR1_num(:,s),pRec_FR1_den(:,s)] = FR1_binned_pRec(subj);
    
    % loop ove region
    for roi = 1:5
        
        if roi == 1 || (roi > 1 && any(eventInds{roi-1}));
            
            
            % filter to just this stimulation location
            if roi > 1
                eventsROI = events(eventInds{roi-1});
            else
                eventsROI = events;
            end
            
            %--------------------------------------------------------------
            % Prob Recall Analyses
            
            % filter to word events. will average the recalled field as a
            % function of various groupings
            encInds        = strcmp({eventsROI.type},'WORD');
            sessVec        = [eventsROI(encInds).session];
            recVec         = [eventsROI(encInds).recalled];
            serialVec      = [eventsROI(encInds).serialpos];
            listVec        = [eventsROI(encInds).list];
            nonStimListVec = [eventsROI(encInds).stimList]==0;
            stimItemVec    = [eventsROI(encInds).isStim]==1;
            nonStimItemStimList = nonStimListVec == 0 & stimItemVec==0;
            
            pRec_stimItem(s,roi) = sum(recVec & stimItemVec)/sum(stimItemVec);
            pRec_nonStimItem(s,roi) = sum(recVec & ~stimItemVec)/sum(~stimItemVec);
            if abs(pRec_stimItem(s,roi) - pRec_nonStimItem(s,roi)) > 0
                subjs_in_region(1,roi) = subjs_in_region(1,roi) + 1;
                stimEff_in_region(s,roi) = 1;
            else
                stimEff_in_region(s,roi) = 0;
                % continue
            end
                        
            % will mean the data as a function of list or serial
            % position
            binned_vec    = double(serialVec >= 6) + 1;
            group_vars    = {listVec serialVec binned_vec};
            
            % data to be averaged will be filtered by 1. stim items; 2.
            % non stim items on stim lists, and 3. non-stim items on
            % non-stim lists, 4. stim lists
            stimType_vars = {stimItemVec,nonStimItemStimList,nonStimListVec,~nonStimListVec};
            

            pRec_listTmp = NaN(1,25,4);
            pRec_listNum = NaN(1,25,4);
            pRec_listDen = NaN(1,25,4);

            pRec_serTmp  = NaN(1,12,4);
            pRec_serNum  = NaN(1,12,4);
            pRec_serDen  = NaN(1,12,4);
            
            pRec_binTmp  = NaN(1,2,4);
            pRec_binNum  = NaN(1,2,4);
            pRec_binDen  = NaN(1,2,4);            
            
            pRec_cell    = {pRec_listTmp,pRec_serTmp,pRec_binTmp};
            pRec_Nums    = {pRec_listNum,pRec_serNum,pRec_binNum};
            pRec_Dens    = {pRec_listDen,pRec_serDen,pRec_binDen};
            for stimType = 1:length(stimType_vars)
                for group = 1:length(group_vars)
                    [m,sdev,den,num]=grpstats(recVec(stimType_vars{stimType}),group_vars{group}(stimType_vars{stimType}),{'mean','std','numel','sum'});                    
                    uniq_group = unique(group_vars{group}(stimType_vars{stimType}));
                    pRec_cell{group}(1,uniq_group,stimType) = m;       
                    pRec_Nums{group}(1,uniq_group,stimType) = num;
                    pRec_Dens{group}(1,uniq_group,stimType) = den;                    
                end
            end
            
            pRec_list(:,:,s,roi)      = squeeze(nanmean(pRec_cell{1},1));
            pRec_serialPos(:,:,s,roi) = squeeze(nanmean(pRec_cell{2},1));
            pRec_binned(:,:,s,roi)    = squeeze(nanmean(pRec_cell{3},1));            
            pRec_binnedNum(:,:,s,roi) = squeeze(pRec_Nums{3});
            pRec_binnedDen(:,:,s,roi) = squeeze(pRec_Dens{3});
            
            % stim - non stim lists probability change
            x1 = squeeze(pRec_Nums{2}(1,:,4));
            n1 = squeeze(pRec_Dens{2}(1,:,4));
            x2 = squeeze(pRec_Nums{2}(1,:,3));
            n2 = squeeze(pRec_Dens{2}(1,:,3));            
            pRecDiff_serialPos_prob(:,s,roi)    = (x1./n1) - (x2./n2);
            pRecDiff_serialPos_logOdds(:,s,roi) = computeLogOddsRatio(x1,n1,x2,n2);
            pRecDiff_serialPos_zProp(:,s,roi)   = ztest_proportions(x1,n1,x2,n2);
            
            % stim - non stim lists probability change binned
            x1 = squeeze(pRec_Nums{3}(1,:,1));
            n1 = squeeze(pRec_Dens{3}(1,:,1));
            x2 = squeeze(pRec_Nums{3}(1,:,3));
            n2 = squeeze(pRec_Dens{3}(1,:,3));            
            pRecDiff_binned_prob(:,s,roi)    = (x1./n1) - (x2./n2);
            pRecDiff_binned_logOdds(:,s,roi) = computeLogOddsRatio(x1,n1,x2,n2);
            pRecDiff_binned_zProp(:,s,roi)   = ztest_proportions(x1,n1,x2,n2);            
            
            %--------------------------------------------------------------
            
            % CRP analyses
            %--------------------------------------------------------------
            % convert to data struct
            data          = FRdata(eventsROI,'list');
            cleanRecMat   = make_clean_recalls_mask2d(data.recalls);
            
            % vector of stim lists
            stimLists = all(data.pres.stimList,2);
            
            if any(stimLists) && any(~stimLists)
                
                %%%%%%%%%%%%%
                % TRY EXCLUDING FIRST 1 OR 2 OPs
                %%%%%%%%%%%%%%%%%%%%%%%
                % crp stim
                rec_mask = make_clean_recalls_mask2d(data.recalls(stimLists,:));                
                rec_mask(:,excludeOPs) = 0;
                [crp_stimLists_subj,numStim,denStim] = crp_jfm(data.recalls(stimLists,:),data.subject(stimLists),data.listLength,rec_mask);
                crp_stimLists(:,s,roi) = crp_stimLists_subj;
                
                numStim_bin = [sum(numStim(10:11)) NaN sum(numStim(13:14))];
                denStim_bin = [sum(denStim(10:11)) NaN sum(denStim(13:14))];                
                crp_stimLists_bin(:,s,roi) = numStim_bin./denStim_bin;
                               
                % crp nonstim
                rec_mask = make_clean_recalls_mask2d(data.recalls(~stimLists,:));                
                rec_mask(:,excludeOPs) = 0;                
                [crp_nonStimLists_subj,numNonStim,denNonStim] = crp_jfm(data.recalls(~stimLists,:),data.subject(~stimLists),data.listLength,rec_mask);
                crp_nonStimLists(:,s,roi) = crp_nonStimLists_subj;
                
                numNonStim_bin = [sum(numNonStim(10:11)) NaN sum(numNonStim(13:14))];
                denNonStim_bin = [sum(denNonStim(10:11)) NaN sum(denNonStim(13:14))];                
                crp_nonStimLists_bin(:,s,roi) = numNonStim_bin./denNonStim_bin;                
                
                % log odds stim / non stim
                crp_stimNonStim_logOdds(:,s,roi) = computeLogOddsRatio(numStim,denStim,numNonStim,denNonStim);
                crp_stimNonStim_logOdds_bin(:,s,roi) = computeLogOddsRatio(numStim_bin,denStim_bin,numNonStim_bin,denNonStim_bin);                                
                
                % z test for prop
                crp_stimNonStim_zProp(:,s,roi) = ztest_proportions(numStim,denStim,numNonStim,denNonStim);
                crp_stimNonStim_zProp_bin(:,s,roi) = ztest_proportions(numStim_bin,denStim_bin,numNonStim_bin,denNonStim_bin);                                
                
                % temp fact stim
                temp_fact_stim(s,roi) = temp_fact(data.recalls(stimLists,:),data.subject(stimLists),data.listLength);
                temp_fact_signed_stim(s,roi) = signed_temp_fact(data.recalls(stimLists,:),data.subject(stimLists),data.listLength);
                
                % temporal factor nonstim
                temp_fact_nonStim(s,roi) = temp_fact(data.recalls(~stimLists,:),data.subject(~stimLists),data.listLength);
                temp_fact_signed_nonStim(s,roi) = signed_temp_fact(data.recalls(~stimLists,:),data.subject(~stimLists),data.listLength);
                
                % train crp stim
                train_mask = repmat(1:6,2,1);
                train_mask = reshape(train_mask(:),12,1)';
                train_mask = repmat(train_mask,sum(stimLists),1);
                crp_stimLists_train_subj = train_crp(data.recalls(stimLists,:),train_mask,data.subject(stimLists));
                crp_stimLists_train(:,s,roi) = crp_stimLists_train_subj;
                
                % train crp nonsitm
                train_mask = repmat(1:6,2,1);
                train_mask = reshape(train_mask(:),12,1)';
                train_mask = repmat(train_mask,sum(~stimLists),1);
                crp_nonStimLists_train_subj = train_crp(data.recalls(~stimLists,:),train_mask,data.subject(~stimLists));
                crp_nonStimLists_train(:,s,roi) = crp_nonStimLists_train_subj;
                
                                
                % Now need to mask stim items specifically
                % recall stim mask
                stimItemMaskRec = study_mat2recall_mat(data.pres.isStim,data.recalls);
                
                % pres stim mask
                stimItemMaskPres = data.pres.isStim==1;
                
                % do crp for stim to stim transistion      
                rec_mask = stimItemMaskRec(stimLists,:)==1;
                rec_mask(:,excludeOPs) = 0;
                [crp_stimItems_subj,numStim,denStim] = crp_jfm(data.recalls(stimLists,:),data.subject(stimLists),data.listLength,rec_mask,rec_mask,stimItemMaskPres(stimLists,:),stimItemMaskPres(stimLists,:));
                crp_stimItems(:,s,roi) = crp_stimItems_subj;
                                
                % do crp for non-stim on stim lists to non-stim on stim list transistion        
                rec_mask = stimItemMaskRec(stimLists,:)==0;
                rec_mask(:,excludeOPs) = 0;
                [crp_nonStimItems_subj,numNonStim,denNonStim] = crp_jfm(data.recalls(stimLists,:),data.subject(stimLists),data.listLength,rec_mask,rec_mask,~stimItemMaskPres(stimLists,:),~stimItemMaskPres(stimLists,:));
                crp_nonStimItems(:,s,roi) = crp_nonStimItems_subj;
                crp_stimNonStimItems_logOdds(:,s,roi) = computeLogOddsRatio(numStim,denStim,numNonStim,denNonStim);
                                
                % do crp for stim to any transistion
                presAll = true(size(data.pres.itemno));                
                rec_mask = stimItemMaskRec(stimLists,:)==1;
                rec_mask(:,excludeOPs) = 0;
                to_rec_mask = cleanRecMat(stimLists,:);
                to_rec_mask(:,excludeOPs) = 0;
                [crp_stimItems_toAny_subj,numStim,denStim] = crp_jfm(data.recalls(stimLists,:),data.subject(stimLists),data.listLength,rec_mask,to_rec_mask,stimItemMaskPres(stimLists,:),presAll(stimLists,:));
                crp_stimItems_toAny(:,s,roi) = crp_stimItems_toAny_subj;
                               
                
                % do crp for non-stim on stim lists to any on stim list transistion         
                rec_mask = stimItemMaskRec(stimLists,:)==0;
                rec_mask(:,excludeOPs) = 0;
                to_rec_mask = cleanRecMat(stimLists,:);
                to_rec_mask(:,excludeOPs) = 0;                
                [crp_nonStimItems_toAny_subj,numNonStim,denNonStim] = crp_jfm(data.recalls(stimLists,:),data.subject(stimLists),data.listLength,rec_mask,to_rec_mask,~stimItemMaskPres(stimLists,:),presAll(stimLists,:));
                crp_nonStimItems_toAny(:,s,roi) = crp_nonStimItems_toAny_subj;
                crp_stimNonStimItems_toAny_logOdds(:,s,roi) = computeLogOddsRatio(numStim,denStim,numNonStim,denNonStim);
                
            end
            %--------------------------------------------------------------
            
            % Intrusion analysis
            %--------------------------------------------------------------
            plis_by_lists = sum(data.intrusions>0,2);
            nPlis_stimLists(s,roi)    = mean(plis_by_lists(stimLists));
            nPlis_nonStimLists(s,roi) = mean(plis_by_lists(~stimLists));
            
            plis_by_lists_1Back = sum(data.intrusions==1,2);
            nPlis_stimLists_1Back(s,roi) = mean(plis_by_lists_1Back(stimLists));
            nPlis_nonStimLists_1Back(s,roi) = mean(plis_by_lists_1Back(~stimLists));
            
            xlis_by_lists = sum(data.intrusions==-1,2);
            nXlis_stimLists(s,roi)    = mean(xlis_by_lists(stimLists));
            nXlis_nonStimLists(s,roi) = mean(xlis_by_lists(~stimLists));
            
            plis_count_subj = zeros(2,3);
            list_count_subj = zeros(2,3);
            
            maskAll = make_mask_exclude_repeats2d(data.rec_itemnos) & data.rec_itemnos > 0;
            
            sessions = unique(sessVec);
            for sess = 1:length(sessions)
                
                % just learned about trial subset and study_mat2recall_mat.
                % Neat.
                dataSess = trial_subset(data.session==sessions(sess),data);
                stimItemRecMat = study_mat2recall_mat(dataSess.pres.isStim,dataSess.recalls);
                stimListRecMat = study_mat2recall_mat(dataSess.pres.stimList,dataSess.recalls);
                stimListSess   = all(dataSess.pres.stimList,2);
                
                
                % plis
                for listNum = 2:size(dataSess.rec_itemnos,1)
                    currRecNos   = dataSess.rec_itemnos(listNum,:);
                    currRecalls  = dataSess.recalls(listNum,:);
                    currRecMask  = make_mask_exclude_repeats2d(currRecNos) & currRecNos > 0;
                    prevPresNos  = dataSess.pres_itemnos(listNum-1,:);
                    prevStimType = dataSess.pres.isStim(listNum-1,:);
                    
                    prevListPlis = ismember(prevPresNos,currRecNos(currRecMask));
                    
                    % categorized whether plis came from stim item (1),
                    % non-stim item on stim list (0), or non-stim list (-1)
                    prevListPliStimType = prevStimType(prevListPlis);
                    if stimListSess(listNum-1)==0
                        prevListPliStimType = prevListPliStimType - 1;
                    end
                    
                    if stimListSess(listNum)
                        plis_count_subj(1,:) = plis_count_subj(1,:) + histc(prevListPliStimType,-1:1);
                        if stimListSess(listNum-1)
                            list_count_subj(1,2:3) = list_count_subj(1,2:3) + 1;
                        else
                            list_count_subj(1,1) = list_count_subj(1,1) + 1;
                        end
                    else
                        plis_count_subj(2,:) = plis_count_subj(2,:) + histc(prevListPliStimType,-1:1);
                        if stimListSess(listNum-1)
                            list_count_subj(2,2:3) = list_count_subj(2,2:3) + 1;
                        else
                            list_count_subj(2,1) = list_count_subj(2,1) + 1;
                        end
                    end
                end
            end
            plis_counts_per_list(:,:,s,roi) = plis_count_subj./list_count_subj;
        end
    end
end


% make plots for each roi (plus all)
regions = {'All','Hipp','EC','PRC/PHC','Hipp/EC'};
figs = [];
if ~exist('figDir','dir')
    mkdir(figDir);
end


% prec data
plotData_pRec = {pRec_serialPos,pRec_list};
fields_pRec = {'spc','list'};
xlabels_pRec = {'Serial Position','List Number'};


% crp data
plotData_crpStim    = {crp_stimLists,crp_stimItems,crp_stimItems_toAny,crp_stimLists_train};
plotData_crpNonStim = {crp_nonStimLists,crp_nonStimItems,crp_nonStimItems_toAny,crp_nonStimLists_train};

x1=[pRec_FR1(1,:)',squeeze(pRec_binned(1,3,:,1)),squeeze(pRec_binned(1,4,:,1))];
x2=[pRec_FR1(2,:)',squeeze(pRec_binned(2,3,:,1)),squeeze(pRec_binned(2,4,:,1))];
x = [x1 x2];
bad=any(isnan(x),2);
dv = x(~bad,:);
iv1 = repmat([ones(size(dv,1),1) ones(size(dv,1),1)+1 ones(size(dv,1),1)+2],1,2);
iv2 = [ones(size(dv,1),3) ones(size(dv,1),3)+1];
s = repmat([1:size(dv,1)]',1,6);
x = [dv(:) iv1(:) iv2(:) s(:)];


x1=[pRec_FR1(1,:)',squeeze(pRec_binned(1,4,:,1))];
x2=[pRec_FR1(2,:)',squeeze(pRec_binned(2,4,:,1))];
x = [x1 x2];
bad=any(isnan(x),2);
dv = x(~bad,:);
iv1 = repmat([ones(size(dv,1),1) ones(size(dv,1),1)+1],1,2);
iv2 = [ones(size(dv,1),2) ones(size(dv,1),2)+1];
s = repmat([1:size(dv,1)]',1,4);
x = [dv(:) iv1(:) iv2(:) s(:)];

FR1_Fr2NS = (pRec_FR1_num+squeeze(pRec_binnedNum(:,3,:,1)))./(pRec_FR1_den+squeeze(pRec_binnedDen(:,3,:,1)));
x1=[FR1_Fr2NS(1,:)',squeeze(pRec_binned(1,4,:,1))];
x2=[FR1_Fr2NS(2,:)',squeeze(pRec_binned(2,4,:,1))];
x = [x1 x2];
bad=any(isnan(x),2);
dv = x(~bad,:);
iv1 = repmat([ones(size(dv,1),1) ones(size(dv,1),1)+1],1,2);
iv2 = [ones(size(dv,1),2) ones(size(dv,1),2)+1];
s = repmat([1:size(dv,1)]',1,4);
x = [dv(:) iv1(:) iv2(:) s(:)];


fields_crp = {'crp_list','crp_item','crp_item_any','crp_train'};
titles_crp = {'List Level','Item to (same) Item','Item to (any) Item','Train'};

for r = 1:length(regions)
    figs(r).region = regions{r};
    figSubDir = fullfile(figDir,regions{r});
    if ~exist(figSubDir,'dir');
        mkdir(figSubDir);
    end
    
    % plot p rec analyses
    for i = 1:2
        sdev  = nanstd(plotData_pRec{i},0,3);
        nSubj = sum(~isnan(plotData_pRec{i}),3);
        e     = squeeze(sdev ./ sqrt(nSubj-1));
        m     = squeeze(nanmean(plotData_pRec{i},3));
        
        figure(1)
        clf
        hold on
        errorbar(1:size(m,1),m(:,1,r),e(:,1,r),'linewidth',2)
        errorbar(1:size(m,1),m(:,2,r),e(:,2,r),'linewidth',2)
        errorbar(1:size(m,1),m(:,3,r),e(:,3,r),'linewidth',2)
        h=legend('Stim','Non-Stim Items','Non-Stim Lists');
        set(h,'fontsize',16);
        xlabel(xlabels_pRec{i},'fontsize',20);
        ylabel('Prob. Recall','fontsize',20)
        grid on
        set(gca,'fontsize',20);
        set(gca,'gridlinestyle',':')
        title(['Region ' regions{r}]);
        fname = fullfile(figSubDir,fields_pRec{i});
        print('-depsc2','-loose',fname);
        figs(r).(fields_pRec{i}) = fname;
    end          
            
    
    figure(7)
    clf
    FR1_Fr2NS = (pRec_FR1_num+squeeze(pRec_binnedNum(:,3,:,r)))./(pRec_FR1_den+squeeze(pRec_binnedDen(:,3,:,r)));
    x1=[FR1_Fr2NS(1,:)',squeeze(pRec_binned(1,4,:,r))];
    x2=[FR1_Fr2NS(2,:)',squeeze(pRec_binned(2,4,:,r))];
    x = [x1 x2];
    bad=any(isnan(x),2);
    x1 = x1(~bad,:);
    x2 = x2(~bad,:);
    dv = x(~bad,:);
    iv1 = repmat([ones(size(dv,1),1) ones(size(dv,1),1)+1],1,2);
    iv2 = [ones(size(dv,1),2) ones(size(dv,1),2)+1];
    s = repmat([1:size(dv,1)]',1,4);
    x = [dv(:) iv1(:) iv2(:) s(:)];
    
    pRec_bin1_region =[pRec_FR1(1,:)',squeeze(pRec_binned(1,3,:,r)),squeeze(pRec_binned(1,4,:,r))];
    pRec_bin2_region =[pRec_FR1(2,:)',squeeze(pRec_binned(2,3,:,r)),squeeze(pRec_binned(2,4,:,r))];
    bad = any(isnan([pRec_bin1_region pRec_bin2_region]),2);
    pRec_bin1_region = pRec_bin1_region(~bad,:);
    pRec_bin2_region = pRec_bin2_region(~bad,:);    
    
%     pRec_bin1_region = pRec_bin1_region - repmat(nanmean(pRec_bin1_region,2),[1 size(pRec_bin1_region,2)]);
%     pRec_bin2_region = pRec_bin2_region - repmat(nanmean(pRec_bin2_region,2),[1 size(pRec_bin2_region,2)]);
    c = {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0 0.4470 0.7410]};
    c = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
    plotData = [nanmean(pRec_bin1_region);nanmean(pRec_bin2_region)];
    plotData = [nanmean(x1);nanmean(x2)];
    e = [nanstd(pRec_bin1_region);nanstd(pRec_bin2_region)]./sqrt(size(pRec_bin1_region,1)-1);
    e = [loftus_masson(x1);loftus_masson(x2)];
    e = loftus_masson([x1 x2]);
    e=reshape(e,[],2)';
    h = bar(plotData);
    grid on
    set(gca,'gridlinestyle',':')
    
    xs = [(1-.8/3):.8/3:(1+.8/3);(2-.8/3):.8/3:(2+.8/3)];
    xs = [(1-1/4.5):1/4.5:(1+1/4.5);(2-1/4.5):1/4.5:(2+1/4.5)];
    xs = [(1-1/7):2/7:(1+1/7);(2-1/7):2/7:(2+1/7)];
    for i = 1:length(h)
       h(i).LineWidth = 2; 
       h(i).FaceColor = c{i};   
       hold on
       errorbar(xs(:,i),plotData(:,i),e(:,i),'k','linewidth',2,'linestyle','none')
    end
    set(gca,'xticklabel',{'Items 1-5','Items 6-12'})
    ylabel('Prob. Recall','fontsize',20)
    set(gca,'fontsize',20)
    h=legend('FR1+FR2 Non-Stim Lists','FR2 Stim Lists');
    h.FontSize = 16;
    p1=anova_rm(pRec_bin1_region,'off');
    p2=anova_rm(pRec_bin2_region,'off');
    p = RMAOV2_mod(x);
    titleStr = sprintf('Region %s: Bin p = %.3f, Exp p = %.3f, Int p = %.3f',regions{r},p(2),p(1),p(3));
    title(titleStr);
    set(gca,'TitleFontWeight','normal')
    set(gca,'titlefontsize',.9)
    fname = fullfile(figSubDir,'spcDiff_binned');
    print('-depsc2','-loose',fname);
    figs(r).spcDiff_binned = fname;     
     
    
    % plot spc diff analyses    
    spcDiff_region = {pRecDiff_serialPos_prob(:,:,r),pRecDiff_serialPos_logOdds(:,:,r),pRecDiff_serialPos_zProp(:,:,r)};
    labels = {'\Delta Prob.','Log Odds Ratio','Z Prop.'};
    
    % crp stim non stim change
    figure(2)
    clf    
    for i = 1:length(spcDiff_region)
        subplot(3,1,i)
        data = spcDiff_region{i};
        data(isinf(data)) = NaN;
        m    = squeeze(nanmean(data,2));
        n    = sum(~isnan(data),2);
        e    = squeeze(nanstd(data,[],2)) ./ sqrt(n-1);
        errorbar(1:12,m,e*1.96,'-k','linewidth',2);
        ylim = get(gca,'ylim');
        set(gca,'ylim',[-1 1]*max(abs(ylim)));
        hold on
        set(gca,'xlim',[0 13])
        x = get(gca,'xlim');
        plot(x,[0 0],'--k');
        grid on
        set(gca,'gridlinestyle',':')
        set(gca,'xtick',1:12);
        if i < 3
            set(gca,'xticklabel','');
        else
            xlabel('Serial Position','fontsize',20)
        end        
        [hSig,p,c,s] = ttest(data');        
        if any(hSig==1)
            inds = find(hSig==1);            
            for ind = 1:length(inds)
                errorbar(inds(ind),m(inds(ind)),e(inds(ind))*1.96,'r','linewidth',2);
            end
        end   
        ylabel(labels{i},'fontsize',20)
        set(gca,'fontsize',20)
        if i == 1
            title(['Region ' regions{r},': SPC Stim - Non Stim List']);
        end        
    end
    fname = fullfile(figSubDir,'spcDiff');
    print('-depsc2','-loose',fname);
    figs(r).spcDiff = fname;                   
        
    % plot crp diff analyses
    crpDiff_prop   = crp_stimLists(:,:,r) - crp_nonStimLists(:,:,r);
    crpDiff_region = {crpDiff_prop,crp_stimNonStim_logOdds(:,:,r),crp_stimNonStim_zProp(:,:,r)};
    labels = {'\Delta Prob.','Log Odds Ratio','Z Prop.'};
    
    % crp stim non stim change
    figure(3)
    clf    
    for i = 1:length(crpDiff_region)
        subplot(3,1,i)
        data = crpDiff_region{i};
        data(isinf(data)) = NaN;
        m    = squeeze(nanmean(data,2));
        n    = sum(~isnan(data),2);
        e    = squeeze(nanstd(data,[],2)) ./ sqrt(n-1);
        errorbar(1:7,m(9:15),e(9:15)*1.96,'-k','linewidth',2);
        ylim = get(gca,'ylim');
        set(gca,'ylim',[-1 1]*max(abs(ylim)));
        hold on
        x = get(gca,'xlim');
        plot(x,[0 0],'--k');
        grid on
        set(gca,'gridlinestyle',':')
        set(gca,'xtick',1:7);
        if i < 3
            set(gca,'xticklabel','');
        else
            set(gca,'xticklabel',-3:3);
            xlabel('Lag','fontsize',20)
        end
        
        [hSig,p,c,s] = ttest(data(9:15,:)');        
        if any(hSig==1)
            inds = find(hSig==1);            
            for ind = 1:length(inds)
                errorbar(inds(ind),m(8+inds(ind)),e(8+inds(ind))*1.96,'r','linewidth',2);
            end
        end   
        ylabel(labels{i},'fontsize',20)
        set(gca,'fontsize',20)
        if i == 1
            title(['Region ' regions{r},': List Level']);
        end            
    end
    fname = fullfile(figSubDir,'crpDiff');
    print('-depsc2','-loose',fname);
    figs(r).crpDiff = fname;
    
    % crp stim non-stim binned lags   
    crpDiff_prop_bin   = crp_stimLists_bin(:,:,r) - crp_nonStimLists_bin(:,:,r);
    crpDiff_region_bin = {crpDiff_prop_bin,crp_stimNonStim_logOdds_bin(:,:,r),crp_stimNonStim_zProp_bin(:,:,r)};
    labels = {'\Delta Prob.','Log Odds Ratio','Z Prop.'};
    
    % crp stim non stim change
    figure(4)
    clf    
    for i = 1:length(crpDiff_region_bin)
        subplot(3,1,i)
        data = crpDiff_region_bin{i};
        data(isinf(data)) = NaN;
        m    = squeeze(nanmean(data,2));
        n    = sum(~isnan(data),2);
        e    = squeeze(nanstd(data,[],2)) ./ sqrt(n-1);
        bar([1 3],m([1 3]),'w','linewidth',2)
        hold on
        errorbar(1:3,m,e*1.96,'-k','linewidth',2);
        ylim = get(gca,'ylim');
        set(gca,'ylim',[-1 1]*max(abs(ylim)));
        hold on
        x = get(gca,'xlim');
        plot(x,[0 0],'--k');
        grid on
        set(gca,'gridlinestyle',':')
        set(gca,'xtick',[1 3]);
        if i < 3
            set(gca,'xticklabel','');
        else
            set(gca,'xticklabel',{'-2,-1','1,2'});
            xlabel('Lag Bin','fontsize',20)
        end
        
        [hSig,p,c,s] = ttest(data');        
        if any(hSig==1)
            inds = find(hSig==1);            
            for ind = 1:length(inds)
                errorbar(inds(ind),m(inds(ind)),e(inds(ind))*1.96,'r','linewidth',2);
            end
        end   
        ylabel(labels{i},'fontsize',20)
        set(gca,'fontsize',20)
        if i == 1
            title(['Region ' regions{r},': List Level']);
        end        
    end
    fname = fullfile(figSubDir,'crpDiff_bin');
    print('-depsc2','-loose',fname);
    figs(r).crpDiff_bin = fname;
   
    
    for i = 1:4
                
        mStim     = squeeze(nanmean(plotData_crpStim{i}(:,:,r),2));
        nSubjStim = sum(~isnan(plotData_crpStim{i}(:,:,r)),2);
        eStim     = squeeze(nanstd(plotData_crpStim{i}(:,:,r),[],2)) ./ sqrt(nSubjStim-1);
                        
        mNonStim     = squeeze(nanmean(plotData_crpNonStim{i}(:,:,r),2));
        nSubjNonStim = sum(~isnan(plotData_crpNonStim{i}(:,:,r)),2);
        eSNonStim    = squeeze(nanstd(plotData_crpNonStim{i}(:,:,r),[],2)) ./ sqrt(nSubjNonStim-1);
        
        crpDiff      = plotData_crpStim{i}(:,:,r) - plotData_crpNonStim{i}(:,:,r);
        mDiff        = squeeze(nanmean(crpDiff,2));
        nSubjDiff    = sum(~isnan(crpDiff),2);
        eDiff        = squeeze(nanstd(crpDiff,[],2)) ./ sqrt(nSubjDiff-1);                        
       
        figure(5)
        clf
        h = [];
        pos = get(gca,'position');
        clf;
        if i < 4
            
            axes('position',[pos(1) pos(2) + .3 pos(3) pos(4) - .3])            
            h(1)=errorbar(1:7,mStim(9:15),eStim(9:15)*1.96,'-k','linewidth',2);
            set(h,'Color',[154, 51, 52]/255);
            hold on
            h(2)=errorbar(1:7,mNonStim(9:15),eSNonStim(9:15)*1.96,'-k','linewidth',2);
            set(gca,'xtick',1:7);
            set(gca,'xticklabel','');
            ylabel('Cond. Resp. Prob.','fontsize',20);
            legend(h,'S','NS')
            set(gca,'fontsize',20)
            grid on
            set(gca,'gridlinestyle',':')
            title(['Region ' regions{r},': ',titles_crp{i}]);
            
            axes('position',[pos(1) pos(2)+.03 pos(3) .22])
            errorbar(1:7,mDiff(9:15),eDiff(9:15)*1.96,'-k','linewidth',2);
            set(gca,'xtick',1:7);
            set(gca,'xticklabel',-3:3);
            xlabel('Lag','fontsize',20)
            set(gca,'fontsize',20)
            ylim = get(gca,'ylim');
            set(gca,'ylim',[-1 1]*max(abs(ylim)));
            hold on
            x = get(gca,'xlim');
            plot(xlim,[0 0],'--k');
            grid on
            set(gca,'gridlinestyle',':')
            [hSig,p,c,s] = ttest(crpDiff(9:15,:)');
            
            if any(hSig==1)
                inds = find(hSig==1);
                
                for ind = 1:length(inds)
                    errorbar(inds(ind),mDiff(8+inds(ind)),eDiff(8+inds(ind))*1.96,'r','linewidth',2);
                end
            end            
        else
            h(1)=errorbar(1:11,mStim,eStim*1.96,'-k','linewidth',2);
            set(h,'Color',[154, 51, 52]/255);
            hold on
            h(2)=errorbar(1:11,mNonStim,eSNonStim*1.96,'-k','linewidth',2);
            xlabel('Train Lag','fontsize',20)
            set(gca,'xtick',1:11);
            set(gca,'xlim',[0 12]);
            set(gca,'xticklabel',-5:5);
            title(['Region ' regions{r},': ','Train CRP']);
            grid on
            set(gca,'gridlinestyle',':')
            set(gca,'fontsize',20)
        end          
        fname = fullfile(figSubDir,fields_crp{i});        
        print('-depsc2','-loose',fname);
        figs(r).(fields_crp{i}) = fname;
        
    end
    
    % pli
    figure(6)
    clf
    m     = nanmean(plis_counts_per_list(:,:,:,r),3);
    nSubj = sum(~isnan(plis_counts_per_list(:,:,:,r)),3);
    e     = nanstd(plis_counts_per_list(:,:,:,r),[],3) ./ sqrt(nSubj-1);
    colors = {[0,98,139]/255,[230,230,220]/255,[129,165,148]/255};
    x1 = 1:3;
    x2 = 5:7;
    hAll = [];
    for i = 1:3
        hold on
        hAll(i)=bar(x1(i),m(1,i),'w','linewidth',2);
        set(hAll(i),'facecolor',colors{i});
        
        h=bar(x2(i),m(2,i),'w','linewidth',2);
        set(h,'facecolor',colors{i});
        
        set(gca,'xtick',[2 6])
        set(gca,'xticklabel',{'Stim','Non-Stim'})
        xlabel('Current List Type','fontsize',20)
        ylabel('PLIs Per List to List Type','fontsize',20)
        set(gca,'fontsize',20)
        
        
        
    end
    h=legend(hAll,'Non-stim lists','Non-stim items','Stim items','location','best');
    set(h,'fontsize',14)
    
    errorbar(x1,m(1,:),e(1,:)*1.96,'k','linewidth',2,'linestyle','none')
    errorbar(x2,m(2,:),e(2,:)*1.96,'k','linewidth',2,'linestyle','none')
    grid on
    set(gca,'gridlinestyle',':')
    title(['Region ' regions{r},': PLIs']);
    fname = fullfile(figSubDir,'plis');
    print('-depsc2','-loose',fname);
    figs(r).plis = fname;    
end

%%%% make report
texName = 'FR2_stimEff_diffPlots.tex';
write_texfile(figDir,texName,figs,subjs_in_region)

curr_dir = pwd;
cd(figDir);
fprintf('Compiling pdf...\n');
unix(['pdflatex -shell-escape ' fullfile(figDir, texName)]);
unix(['rm ' texName(1:end-3) 'aux']);
unix(['rm ' texName(1:end-3) 'log']);
fprintf('Done!\n');
cd(curr_dir);




% Start making the tex file
function write_texfile(saveDir,texName, figs,subjs_in_region)

% Write the document. If you do not have write permission, this will crash.
fid = fopen(fullfile(saveDir,texName),'w');

if fid==-1;
    error(sprintf('cannot open %s',texName))
end

% Write out the preamble to the tex doc. This is standard stuff and doesn't
% need to be changed
fprintf(fid,'\\documentclass[a4paper]{article} \n');
fprintf(fid,'\\usepackage[usenames,dvipsnames,svgnames,table]{xcolor}\n');
fprintf(fid,'\\usepackage{graphicx,multirow} \n');
fprintf(fid,'\\usepackage{epstopdf} \n');
fprintf(fid,'\\usepackage[small,bf,it]{caption}\n');
fprintf(fid,'\\usepackage{subfigure,amsmath} \n');
% fprintf(fid,'\\usepackage{wrapfig} \n');
% fprintf(fid,'\\usepackage{longtable} \n');
% fprintf(fid,'\\usepackage{pdfpages}\n');
% fprintf(fid,'\\usepackage{mathtools}\n');
% fprintf(fid,'\\usepackage{array}\n');
% fprintf(fid,'\\usepackage{enumitem}\n');
% fprintf(fid,'\\usepackage{sidecap} \\usepackage{soul}\n');

% fprintf(fid,'\\setlength\\belowcaptionskip{5pt}\n');
fprintf(fid,'\n');
fprintf(fid,'\\addtolength{\\oddsidemargin}{-.875in} \n');
fprintf(fid,'\\addtolength{\\evensidemargin}{-.875in} \n');
fprintf(fid,'\\addtolength{\\textwidth}{1.75in} \n');
fprintf(fid,'\\addtolength{\\topmargin}{-.75in} \n');
fprintf(fid,'\\addtolength{\\textheight}{1.75in} \n');
fprintf(fid,'\n');
fprintf(fid,'\\newcolumntype{C}[1]{>{\\centering\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}} \n');

fprintf(fid,'\\usepackage{fancyhdr}\n');
fprintf(fid,'\\pagestyle{fancy}\n');
fprintf(fid,'\\fancyhf{}\n');
fprintf(fid,'\\lhead{ }\n');
fprintf(fid,'\\rhead{Date created: %s}\n',date);

fprintf(fid,'\\usepackage{hyperref}\n');

% Start the document
fprintf(fid,'\\begin{document}\n\n\n');

for r = 1:length(figs)
    
    % error by trial
    fprintf(fid,'\\begin{figure}\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).spc);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).spcDiff);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).spcDiff_binned);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).crp_list);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).crpDiff);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).crpDiff_bin);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).crp_train);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).plis);
    fprintf(fid,['\\caption{\\textbf{%d subjects.  }\\textbf{A: } Probability of recalll by serial position. ',...
        '\\textbf{B: } SPC stim - non stim lists. ',...
        '\\textbf{C: } Probabilty of recall by bined serial position. ',...
        '\\textbf{D: } CRP for \\emph{stim lists} and \\emph{non-stim lists.} ',...
        '\\textbf{E: } CRP stim - non stim lists. ',...        
        '\\textbf{F: } CRP stim - non stim lists with binned first two positive and negative lags. ',...
        '\\textbf{G: } Train Lag CRP. ',...
        '\\textbf{H:} Number of prior list intrusions per list as a function of list type.}\n'],subjs_in_region(r));
    fprintf(fid,'\\end{figure}\n\n\n');
    fprintf(fid,'\\clearpage\n\n\n');
    
end
fprintf(fid,'\\end{document}\n\n\n');
fclose(fid);





function events = addRegionToEvents(subj,events);


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
        if strcmp(subj,'R1016M')
            stimLoc = 'Left DLPFC';
        elseif strcmp(subj,'R1028M')
            stimLoc = 'Right EC';
        elseif strcmp(subj,'R1035M')
            stimLoc = 'Left PRC';
        elseif strcmp(subj,'R1053M')
            stimLoc = 'Left PRC';
        elseif strcmp(subj,'R1074M')
            stimLoc = 'Left ACG';
        elseif strcmp(subj,'UT009a')
            stimLoc = 'Left PRC';
        elseif strcmp(subj,'R1115T')
            stimLoc = 'Right Anterior Insula';
        elseif strcmp(subj,'UT010')
            stimLoc = 'Left PRC';
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

