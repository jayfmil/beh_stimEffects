function spc_for_exp()

subjs_fr1 = get_subs('RAM_FR1');
subjs_fr2 = get_subs('RAM_FR2');
subjs = unique({subjs_fr1{:} subjs_fr2{:}});

for s = 1:length(subjs)
    
    if ismember(subjs{s},subjs_fr1) && ismember(subjs{s},subjs_fr2)        
        events_fr1 = get_sub_events('RAM_FR1',subjs{s});
        events_fr2 = get_sub_events('RAM_FR2',subjs{s});
        try
            data_fr1   = FRdata(events_fr1,'list');
            data_fr2   = FRdata(events_fr2,'list');
            
            if s == 1
                spcs_subj_fr1 = NaN(length(subjs),data_fr1.listLength);
                spcs_subj_fr2 = NaN(length(subjs),data_fr2.listLength);
                fr1_mean_spc = NaN(length(subjs),2);
                fr2_mean_spc = NaN(length(subjs),2);
                
                spcs_subj_fr2_stimLists = NaN(length(subjs),data_fr2.listLength);
                fr2_mean_spc_stimLists = NaN(length(subjs),2);
                spcs_subj_fr2_nonStimLists = NaN(length(subjs),data_fr2.listLength);
                fr2_mean_spc_nonStimLists = NaN(length(subjs),2);
                
            end
            spcs_subj_fr1(s,:) = spc(data_fr1.recalls,data_fr1.subject,data_fr1.listLength);                                    
            fr1_mean_spc(s,:) = [ mean(spcs_subj_fr1(s,1:5))  mean(spcs_subj_fr1(s,6:end))];
            
            spcs_subj_fr2(s,:) = spc(data_fr2.recalls,data_fr2.subject,data_fr2.listLength);
            fr2_mean_spc(s,:) = [ mean(spcs_subj_fr2(s,1:5))  mean(spcs_subj_fr2(s,6:end))];
            
            stimLists = all(data_fr2.pres.stimList,2);            
            spcs_subj_fr2_stimLists(s,:) = spc(data_fr2.recalls(stimLists,:),data_fr2.subject(stimLists,:),data_fr2.listLength);          
            fr2_mean_spc_stimLists(s,:) = [mean(spcs_subj_fr2_stimLists(s,1:5))  mean(spcs_subj_fr2_stimLists(s,6:end))];
            spcs_subj_fr2_nonStimLists(s,:) = spc(data_fr2.recalls(~stimLists,:),data_fr2.subject(~stimLists,:),data_fr2.listLength);
            fr2_mean_spc_nonStimLists(s,:) = [mean(spcs_subj_fr2_nonStimLists(s,1:5))  mean(spcs_subj_fr2_nonStimLists(s,6:end))];
               
            
        end
    end
end
keyboard









