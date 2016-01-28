function pRec_by_cond(subjs,figDir)


if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_FR2');
end
subjs(strcmp(subjs,'UT009a')) = [];

% will hold prob of recall results
% num lists/items x num conditions x num subjs x num regions
pRec_list      = NaN(25,4,length(subjs),5);
pRec_serialPos = NaN(12,4,length(subjs),5);

pRec_stimItem    = NaN(length(subjs),5);
pRec_nonStimItem = NaN(length(subjs),5);

% will hold crp results
% lags x subjs x regions
% stim lists
crp_stimLists          = NaN(23,length(subjs),5);
crp_stimLists_num      = NaN(23,length(subjs),5);
crp_stimLists_den      = NaN(23,length(subjs),5);

% and non stim lists
crp_nonStimLists       = NaN(23,length(subjs),5);
crp_nonStimLists_num   = NaN(23,length(subjs),5);
crp_nonStimLists_den   = NaN(23,length(subjs),5);

% stim items
crp_stimItems          = NaN(23,length(subjs),5);
crp_stimItems_num      = NaN(23,length(subjs),5);
crp_stimItems_den      = NaN(23,length(subjs),5);

% non stim items
crp_nonStimItems       = NaN(23,length(subjs),5);
crp_nonStimItems_num   = NaN(23,length(subjs),5);
crp_nonStimItems_den   = NaN(23,length(subjs),5);

% stim items to any
crp_stimItems_toAny     = NaN(23,length(subjs),5);
crp_stimItems_toAny_num = NaN(23,length(subjs),5);
crp_stimItems_toAny_den = NaN(23,length(subjs),5);

% non stim items to any
crp_nonStimItems_toAny     = NaN(23,length(subjs),5);
crp_nonStimItems_toAny_num = NaN(23,length(subjs),5);
crp_nonStimItems_toAny_den = NaN(23,length(subjs),5);

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
            if abs(pRec_stimItem(s,roi) - pRec_nonStimItem(s,roi)) > .03
                subjs_in_region(1,roi) = subjs_in_region(1,roi) + 1;
                stimEff_in_region(s,roi) = 1;
            else
                stimEff_in_region(s,roi) = 0;
                % continue
            end
            
            
            % loop over each session. Maybe don't need to do this
            sessions      = unique(sessVec);
            pRec_listSess = NaN(length(sessions),25,3);
            pRec_serSess  = NaN(length(sessions),12,3);
            pRec_cell     = {pRec_listSess,pRec_serSess};
            for sess = 1:length(sessions);
                
                % pull out data for this session
                sessInds      = sessVec == sessions(sess);
                recSess       = recVec(sessInds);
                
                % will mean the data as a function of list or serial
                % position
                group_vars    = {listVec(sessInds) serialVec(sessInds)};
                
                % data to be averaged will be filtered by 1. stim items; 2.
                % non stim items on stim lists, and 3. non-stim items on
                % non-stim lists
                stimType_vars = {stimItemVec(sessInds),nonStimItemStimList(sessInds),nonStimListVec(sessInds),~nonStimListVec(sessInds)};
                
                for stimType = 1:length(stimType_vars)
                    for group = 1:length(group_vars)
                        [m,sdev,n]=grpstats(recSess(stimType_vars{stimType}),group_vars{group}(stimType_vars{stimType}),{'mean','std','numel'});
                        uniq_group = unique(group_vars{group}(stimType_vars{stimType}));
                        pRec_cell{group}(sess,uniq_group,stimType) = m;
                    end
                end
                
            end
            pRec_list(:,:,s,roi)      = squeeze(nanmean(pRec_cell{1},1));
            pRec_serialPos(:,:,s,roi) = squeeze(nanmean(pRec_cell{2},1));
            %--------------------------------------------------------------
            
            % CRP analyses
            %--------------------------------------------------------------
            % convert to data struct
            data          = FRdata(eventsROI,'list');
            cleanRecMat   = make_clean_recalls_mask2d(data.recalls);
            
            % vector of stim lists
            stimLists = all(data.pres.stimList,2);
            
            % Do crp, temporal factor, and train crp for sitm lists
            if any(stimLists)
                
                % crp
                [crp_stimLists_subj,num,den] = crp_jfm(data.recalls(stimLists,:),data.subject(stimLists),data.listLength);
                crp_stimLists(:,s,roi) = crp_stimLists_subj;
                crp_stimLists_num(:,s,roi) = num;
                crp_stimLists_den(:,s,roi) = den;
                
                % do temp fact
                temp_fact_stim(s,roi) = temp_fact(data.recalls(stimLists,:),data.subject(stimLists),data.listLength);
                temp_fact_signed_stim(s,roi) = signed_temp_fact(data.recalls(stimLists,:),data.subject(stimLists),data.listLength);
                
                % do train crp
                train_mask = repmat(1:6,2,1);
                train_mask = reshape(train_mask(:),12,1)';
                train_mask = repmat(train_mask,sum(stimLists),1);
                crp_stimLists_train_subj = train_crp(data.recalls(stimLists,:),train_mask,data.subject(stimLists));
                crp_stimLists_train(:,s,roi) = crp_stimLists_train_subj;
            end
            
            % Do crp, temporal factor, and train crp for NON stim lists
            if any(~stimLists)
                
                % train crp
                train_mask = repmat(1:6,2,1);
                train_mask = reshape(train_mask(:),12,1)';
                train_mask = repmat(train_mask,sum(~stimLists),1);
                crp_nonStimLists_train_subj = train_crp(data.recalls(~stimLists,:),train_mask,data.subject(~stimLists));
                crp_nonStimLists_train(:,s,roi) = crp_nonStimLists_train_subj;
                
                % temporal factor
                temp_fact_nonStim(s,roi) = temp_fact(data.recalls(~stimLists,:),data.subject(~stimLists),data.listLength);
                temp_fact_signed_nonStim(s,roi) = signed_temp_fact(data.recalls(~stimLists,:),data.subject(~stimLists),data.listLength);
                
                % crp
                [crp_nonStimLists_subj,num,den] = crp_jfm(data.recalls(~stimLists,:),data.subject(~stimLists),data.listLength);
                crp_nonStimLists(:,s,roi) = crp_nonStimLists_subj;
                crp_nonStimLists_num(:,s,roi) = num;
                crp_nonStimLists_den(:,s,roi) = den;
                
            end
            
            % Now need to mask stim items specifically
            % recall stim mask
            stimItemMaskRec = study_mat2recall_mat(data.pres.isStim,data.recalls);
            
            % pres stim mask
            stimItemMaskPres = data.pres.isStim==1;
            
            % do crp for stim to stim transistion
            if any(stimLists)
                [crp_stimItems_subj,num,den] = crp_jfm(data.recalls(stimLists,:),data.subject(stimLists),data.listLength,stimItemMaskRec(stimLists,:)==1,stimItemMaskRec(stimLists,:)==1,stimItemMaskPres(stimLists,:),stimItemMaskPres(stimLists,:));
                crp_stimItems(:,s,roi) = crp_stimItems_subj;
                crp_stimItems_num(:,s,roi) = num;
                crp_stimItems_den(:,s,roi) = den;
            end
            
            % do crp for non-stim on stim lists to non-stim on stim list transistion
            if any(~stimLists)
                [crp_nonStimItems_subj,num,den] = crp_jfm(data.recalls(stimLists,:),data.subject(stimLists),data.listLength,stimItemMaskRec(stimLists,:)==0,stimItemMaskRec(stimLists,:)==0,~stimItemMaskPres(stimLists,:),~stimItemMaskPres(stimLists,:));
                crp_nonStimItems(:,s,roi) = crp_nonStimItems_subj;
                crp_nonStimItems_num(:,s,roi) = num;
                crp_nonStimItems_den(:,s,roi) = den;
            end
            
            % do crp for stim to any transistion
            presAll = true(size(data.pres.itemno));
            if any(stimLists)
                [crp_stimItems_toAny_subj,num,den] = crp_jfm(data.recalls(stimLists,:),data.subject(stimLists),data.listLength,stimItemMaskRec(stimLists,:)==1,cleanRecMat(stimLists,:),stimItemMaskPres(stimLists,:),presAll(stimLists,:));
                crp_stimItems_toAny(:,s,roi) = crp_stimItems_toAny_subj;
                crp_stimItems_toAny_num(:,s,roi) = num;
                crp_stimItems_toAny_den(:,s,roi) = den;
                
            end
            
            % do crp for non-stim on stim lists to any on stim list transistion
            if any(~stimLists)
                [crp_nonStimItems_toAny_subj,num,den] = crp_jfm(data.recalls(stimLists,:),data.subject(stimLists),data.listLength,stimItemMaskRec(stimLists,:)==0,cleanRecMat(stimLists,:),~stimItemMaskPres(stimLists,:),presAll(stimLists,:));
                crp_nonStimItems_toAny(:,s,roi) = crp_nonStimItems_toAny_subj;
                crp_nonStimItems_toAny_num(:,s,roi) = num;
                crp_nonStimItems_toAny_den(:,s,roi) = den;
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

% compute standard error of log odds ratio. Need
% 2 x 2:     # transitions made | # poss. transitions - # transitions made
%       stim 
%    nonstim
% is that right?

numData_crpStim     = {crp_stimLists_num,crp_stimItems_num,crp_stimItems_toAny_num};
numData_crpNonStim     = {crp_nonStimLists_num,crp_nonStimItems_num,crp_nonStimItems_toAny_num};

denData_crpStim     = {crp_stimLists_den,crp_stimItems_den,crp_stimItems_toAny_den};

totMinusActualStim     = {crp_stimLists_den-crp_stimLists_num,...
                          crp_stimItems_den-crp_stimItems_num,...
                          crp_stimItems_toAny_den-crp_stimItems_toAny_num};

totMinusActualNonStim  = {crp_nonStimLists_den-crp_nonStimLists_num,...
                          crp_nonStimItems_den-crp_nonStimItems_num,...
                          crp_nonStimItems_toAny_den-crp_nonStimItems_toAny_num};
                      
                      

                      
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
        title(regions{r});
        fname = fullfile(figSubDir,fields_pRec{i});
        print('-depsc2','-loose',fname);
        figs(r).(fields_pRec{i}) = fname;
    end
    
    % crp
    for i = 1:1
        
        %         plotData_crpStim{i} = log10(plotData_crpStim{i}./(1-plotData_crpStim{i}));
        %         plotData_crpNonStim{i} = log10(plotData_crpNonStim{i}./(1-plotData_crpNonStim{i}));
        
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
        
%         crpDiffLogOddsRatio  = log((plotData_crpStim{i}(:,:,r))./(1-plotData_crpStim{i}(:,:,r)) ./ (plotData_crpNonStim{i}(:,:,r))./(1-plotData_crpNonStim{i}(:,:,r)));
        crpOddsRatio = (((numData_crpStim{i}(:,:,r)+.5) ./ (totMinusActualStim{i}(:,:,r)+.5)) ./ ((numData_crpNonStim{i}(:,:,r)+.5) ./ (totMinusActualNonStim{i}(:,:,r)+.5)));
        crpLogOddsRatio = log(crpOddsRatio);        
        se=sqrt(1./(numData_crpStim{i}(:,:,r)+.5) + 1./(numData_crpNonStim{i}(:,:,r)+.5) + 1./(totMinusActualStim{i}(:,:,r)+.5) + 1./(totMinusActualNonStim{i}(:,:,r)+.5));
        z = crpLogOddsRatio./se;
        z(isnan(crpDiff)) = NaN;
        meanZ = nanmean(z,2);
        nSubjZ = sum(~isnan(z),2);       
        eZ = nanstd(z,[],2)./sqrt(nSubjZ-1);
        
        figure(2)
        clf
        h = [];
        pos = get(gca,'position');
        clf;
        if i < 4
            
            axes('position',[pos(1) pos(2) + .3 pos(3) pos(4) - .3])
            title([regions{r},': ',titles_crp{i}]);
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
            
%             axes('position',[pos(1) pos(2)+.03 pos(3) .22])
%             errorbar(1:7,mDiff(9:15),eDiff(9:15)*1.96,'-k','linewidth',2);
%             set(gca,'xtick',1:7);
%             set(gca,'xticklabel',-3:3);
%             xlabel('Lag','fontsize',20)
%             set(gca,'fontsize',20)
%             ylim = get(gca,'ylim');
%             set(gca,'ylim',[-1 1]*max(abs(ylim)));
%             hold on
%             x = get(gca,'xlim');
%             plot(xlim,[0 0],'--k');
%             grid on
%             set(gca,'gridlinestyle',':')
%             [hSig,p,c,s] = ttest(crpDiff(9:15,:)');
%             
%             if any(hSig==1)
%                 inds = find(hSig==1);
%                 
%                 for ind = 1:length(inds)
%                     errorbar(inds(ind),mDiff(8+inds(ind)),eDiff(8+inds(ind))*1.96,'r','linewidth',2);
%                 end
%             end

            axes('position',[pos(1) pos(2)+.03 pos(3) .22])
            errorbar(1:7,meanZ(9:15),eZ(9:15)*1.96,'-k','linewidth',2);
            set(gca,'xtick',1:7);
            set(gca,'xticklabel',-3:3);
            xlabel('Lag','fontsize',20)
            set(gca,'fontsize',20)
            ylim = get(gca,'ylim');
%             set(gca,'ylim',[-1 1]*max(abs(ylim)));
            hold on
            x = get(gca,'xlim');
            plot(xlim,[0 0],'--k');
            grid on
            set(gca,'gridlinestyle',':')
            [hSig,p,c,s] = ttest(z(9:15,:)');
            
            if any(hSig==1)
                inds = find(hSig==1);
                
                for ind = 1:length(inds)
                    errorbar(inds(ind),meanZ(8+inds(ind)),eZ(8+inds(ind))*1.96,'r','linewidth',2);
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
        end
        keyboard
        
        
        
        
        
        fname = fullfile(figSubDir,fields_crp{i});
        
        print('-depsc2','-loose',fname);
        %         print('-depsc2',fname);
        figs(r).(fields_crp{i}) = fname;
        
    end
    % pli
    figure(3)
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
    title([regions{r},': PLIs']);
    fname = fullfile(figSubDir,'plis');
    print('-depsc2','-loose',fname);
    figs(r).plis = fname;
end

%%%% make report
texName = 'FR2_stimEff2.tex';
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
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).list);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).crp_list);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).crp_item);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).crp_item_any);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).crp_train);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.4\\textwidth]{%s}}}\n',figs(r).plis);
    fprintf(fid,['\\caption{\\textbf{%d subjects.  }\\textbf{A: } Probability of recalll by serial position. ',...
        '\\textbf{B: } Probability of recalll by list number. ',...
        '\\textbf{C: } CRP for \\emph{stim lists} and \\emph{non-stim lists.} ',...
        '\\textbf{D: } CRP for \\emph{stim items to stim items} and ',...
        '\\emph{non-stim items to non-stim items on stim lists}. ',...
        '\\textbf{E: } CRP for \\emph{stim items any} and ',...
        '\\emph{non-stim items to any item on stim lists}. ',...
        '\\textbf{F: } Train Lag CRP. ',...
        '\\textbf{G:} Number of prior list intrusions per list as a function of list type.}\n'],subjs_in_region(r));
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

