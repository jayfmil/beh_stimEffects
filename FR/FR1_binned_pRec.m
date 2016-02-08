function [pRec,pRec_num,pRec_den] = FR1_binned_pRec(subj)

pRec     = NaN(2,1);
pRec_num = NaN(2,1);
pRec_den = NaN(2,1);
try
    events     = get_sub_events('RAM_FR1',subj);
    encInds    = strcmp({events.type},'WORD');
    recVec     = [events(encInds).recalled];
    serialVec  = [events(encInds).serialpos];
    binned_vec = double(serialVec >= 6) + 1;
    
    [pRec,sdev,pRec_den,pRec_num]=grpstats(recVec,binned_vec,{'mean','std','numel','sum'});                            

end