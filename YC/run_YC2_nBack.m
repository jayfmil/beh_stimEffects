function run_YC2_nBack(subjs,figDir)
%

if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC2');
end

% third dimension will be all regions, then hipp, ec, prc/phc
distErr_nBack_stim    = NaN(length(subjs),4,5);
closestObj_nBack_stim = NaN(length(subjs),4,5);

% and no stim
distErr_nBack_nonStim    = NaN(length(subjs),4,5);
closestObj_nBack_nonStim = NaN(length(subjs),4,5);

% loop over each subject
subjs_in_region = zeros(1,5);
for s = 1:length(subjs)
    subj = subjs{s};
    fprintf('Processing %s.\n',subj);
    
    % load events
    events = get_sub_events('RAM_YC2',subj);
    events = YC2_addStimRegion(subj,events);
    
    % first all stim
    stim_mask = [events.isStim];
    [distErr_nBack_subj,closestObj_nBack_subj] = YC_nBack_err(events,stim_mask);
    distErr_nBack_stim(s,:,1) = distErr_nBack_subj;
    closestObj_nBack_stim(s,:,1) = closestObj_nBack_subj;
    
    % and all nonstim
    [distErr_nBack_nonStim_subj,closestObj_nBack_nonStim_subj] = YC_nBack_err(events,stim_mask==0);
    distErr_nBack_nonStim(s,:,1) = distErr_nBack_nonStim_subj;
    closestObj_nBack_nonStim(s,:,1) = closestObj_nBack_nonStim_subj;
    subjs_in_region(1,1) = subjs_in_region(1,1) + 1;
    
    % if subject was stimulated in our regions of interest filter to
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
    
    
    % loop over region
    for roi = 1:4
        
        if any(eventInds{roi});
            subjs_in_region(1,roi+1) = subjs_in_region(1,roi+1) + 1;
            % filter to just this stimulation location
            eventsROI     = events(eventInds{roi});
            
            % first all stim in region
            stim_mask = [eventsROI.isStim];
            [distErr_nBack_subj,closestObj_nBack_subj] = YC_nBack_err(eventsROI,stim_mask);
            distErr_nBack_stim(s,:,roi+1) = distErr_nBack_subj;
            closestObj_nBack_stim(s,:,roi+1) = closestObj_nBack_subj;
            
            % and all nonstim in region
            [distErr_nBack_nonStim_subj,closestObj_nBack_nonStim_subj] = YC_nBack_err(eventsROI,stim_mask==0);
            distErr_nBack_nonStim(s,:,roi+1) = distErr_nBack_nonStim_subj;
            closestObj_nBack_nonStim(s,:,roi+1) = closestObj_nBack_nonStim_subj;
            
        end
    end
end


% make plots for each roi (plus all)
regions = {'All','Hipp','EC','PRC/PHC','Hipp/EC'};
figs = [];
if ~exist('figDir','dir')
    mkdir(figDir);
end

for r = 1:length(regions)
    figs(r).region = regions{r};
    figSubDir = fullfile(figDir,regions{r});
    if ~exist(figSubDir,'dir');
        mkdir(figSubDir);
    end
    
    % Figure 1 accuracy by nBack for stim and non-stim
    figure(1)
    clf
    stim_data    = distErr_nBack_stim(:,:,r);
    nonStim_data = distErr_nBack_nonStim(:,:,r);
    plotData = {stim_data,nonStim_data};
    titles = {'Stim','Non-stim'};
    for i = 1:2
        m = nanmean(plotData{i},1);
        s = nanstd(plotData{i});
        e = s./sqrt(sum(~isnan(plotData{i}))-1);
        
        subplot(2,1,i)
        h=bar(1:4,m,'w','linewidth',2);
        hold on
        errorbar(1:4,m,e*1.96,'k','linestyle','none','linewidth',2);
        %     set(gca,'ylim',[-.1 .1])
        %     set(gca,'ytick',[-.1:.05:.1]);
        set(gca,'xticklabel',{'0-Back','1-Back','2-Back','3-Back'});
        ylabel('Accuracy','fontsize',20);
        set(gca,'fontsize',20);
        grid on
        set(gca,'gridlinestyle',':')
        title([titles{i},' ',regions{r}]);
    end
    fname = fullfile(figSubDir,'nBack_accuracy');
    print('-depsc2','-loose',fname);
    figs(r).nBack_accuracy = fname;
    
    % Figure 2 accuracy by nBack for stim and non-stim
    figure(2)
    clf
    stim_data    = closestObj_nBack_stim(:,:,r)*100;
    nonStim_data = closestObj_nBack_nonStim(:,:,r)*100;
    plotData = {stim_data,nonStim_data};
    titles = {'Stim','Non-stim'};
    for i = 1:2
        m = nanmean(plotData{i},1);
        s = nanstd(plotData{i});
        e = s./sqrt(sum(~isnan(plotData{i}))-1);
        
        subplot(2,1,i)
        h=bar(1:4,m,'w','linewidth',2);
        hold on
        errorbar(1:4,m,e*1.96,'k','linestyle','none','linewidth',2);
        %     set(gca,'ylim',[-.1 .1])
        %     set(gca,'ytick',[-.1:.05:.1]);
        set(gca,'xticklabel',{'0-Back','1-Back','2-Back','3-Back'});
        ylabel('% Closest','fontsize',20);
        set(gca,'fontsize',20);
        grid on
        set(gca,'gridlinestyle',':')
        title([titles{i},' ',regions{r}]);
    end
    fname = fullfile(figSubDir,'nBack_closest');
    print('-depsc2','-loose',fname);
    figs(r).nBack_closest = fname;
    
    
    % Figure 3
    figure(3)
    clf
    stim_data    = distErr_nBack_stim(:,:,r);
    nonStim_data = distErr_nBack_nonStim(:,:,r);
    plotData = stim_data-nonStim_data;
    
    m = nanmean(plotData,1);
    s = nanstd(plotData);
    e = s./sqrt(sum(~isnan(plotData))-1);
    
    h=bar(1:4,m,'w','linewidth',2);
    hold on
    errorbar(1:4,m,e*1.96,'k','linestyle','none','linewidth',2);
    %     set(gca,'ylim',[-.1 .1])
    %     set(gca,'ytick',[-.1:.05:.1]);
    set(gca,'xticklabel',{'0-Back','1-Back','2-Back','3-Back'});
    ylabel('Acc. Stim - Acc. Non-Stim','fontsize',20);
    set(gca,'fontsize',20);
    grid on
    set(gca,'gridlinestyle',':')
    title(regions{r});
    
    fname = fullfile(figSubDir,'stimEff_nBack_accuracy');
    print('-depsc2','-loose',fname);
    figs(r).stimEff_nBack_accuracy = fname;
    
    
    
    %
    figure(4)
    clf
    stim_data    = closestObj_nBack_stim(:,:,r)*100;
    nonStim_data = closestObj_nBack_nonStim(:,:,r)*100;
    plotData = stim_data-nonStim_data;
    
    m = nanmean(plotData,1);
    s = nanstd(plotData);
    e = s./sqrt(sum(~isnan(plotData))-1);
    
    h=bar(1:4,m,'w','linewidth',2);
    hold on
    errorbar(1:4,m,e*1.96,'k','linestyle','none','linewidth',2);
    %     set(gca,'ylim',[-.1 .1])
    %     set(gca,'ytick',[-.1:.05:.1]);
    set(gca,'xticklabel',{'0-Back','1-Back','2-Back','3-Back'});
    ylabel('% Closest Stim - % Closest Non-Stim','fontsize',20);
    set(gca,'fontsize',20);
    grid on
    set(gca,'gridlinestyle',':')
    title(regions{r});
    
    fname = fullfile(figSubDir,'stimEff_nBack_closest');
    print('-depsc2','-loose',fname);
    figs(r).stimEff_nBack_closest = fname;
    
end

%%%% make report
texName = 'YC2_nBack.tex';
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
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.49\\textwidth]{%s}}}\n',figs(r).nBack_accuracy);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.49\\textwidth]{%s}}}\n',figs(r).nBack_closest);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.49\\textwidth]{%s}}}\n',figs(r).stimEff_nBack_accuracy);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.49\\textwidth]{%s}}}\n',figs(r).stimEff_nBack_closest);
    fprintf(fid,['\\caption{\\textbf{%d subjects.  }\\textbf{A: } Mean performance as a function of accuracy of ',...
                 'current response location to n-back object locations. ',...
                 '\\textbf{B: } Mean probability of a given n-back object ',...
                 'location being the closest to the current response location. ',...
                 '\\textbf{C: } Panel A Stim - Non-Stim. \\textbf{D: } Panel B Stim - Non-Stim.}\n'],subjs_in_region(r));
    fprintf(fid,'\\end{figure}\n\n\n');
    fprintf(fid,'\\clearpage\n\n\n');
    
end
fprintf(fid,'\\end{document}\n\n\n');
fclose(fid);








