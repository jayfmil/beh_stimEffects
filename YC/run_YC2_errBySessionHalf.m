function run_YC2_errBySessionHalf(subjs,figDir)
%

if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC2');
end

% third dimension will be all regions, then hipp, ec, prc/phc
distErr_sessionHalf_stim    = NaN(length(subjs),2,5);
distErr_sessionHalf_nonStim = NaN(length(subjs),2,5);

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
    [distErr_firstHalf_stim,distErr_secondHalf_stim] = YC_sessionHalf_err(events,stim_mask);
    distErr_sessionHalf_stim(s,1,1) = distErr_firstHalf_stim;
    distErr_sessionHalf_stim(s,2,1) = distErr_secondHalf_stim;
    
    % and all nonstim
    stim_mask = [events.isStim];
    [distErr_firstHalf_nonStim,distErr_secondHalf_nonStim] = YC_innerOuter_err(events,stim_mask==0);
    distErr_sessionHalf_nonStim(s,1,1) = distErr_firstHalf_nonStim;
    distErr_sessionHalf_nonStim(s,2,1) = distErr_secondHalf_nonStim;
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
            [distErr_firstHalf_stim,distErr_secondHalf_stim] = YC_innerOuter_err(eventsROI,stim_mask);
            distErr_sessionHalf_stim(s,1,roi+1) = distErr_firstHalf_stim;
            distErr_sessionHalf_stim(s,2,roi+1) = distErr_secondHalf_stim;
            
            % and all nonstim in region
            [distErr_firstHalf_nonStim,distErr_secondHalf_nonStim] = YC_innerOuter_err(eventsROI,stim_mask==0);
            distErr_sessionHalf_nonStim(s,1,roi+1) = distErr_firstHalf_nonStim;
            distErr_sessionHalf_nonStim(s,2,roi+1) = distErr_secondHalf_nonStim;
            
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
    
    % Figure 1. second half - first half diff for stim and non stim
    stim_data    = distErr_sessionHalf_stim(:,:,r);
    nonStim_data = distErr_sessionHalf_nonStim(:,:,r);
    
    halfDiff = [diff(stim_data,[],2) diff(nonStim_data,[],2)];
    m = nanmean(halfDiff,1);
    s = nanstd(halfDiff);
    e = s./sqrt(sum(~isnan(halfDiff))-1);
    figure(1)
    clf
    h=bar(1:2,m,'w','linewidth',2);
    hold on
    errorbar(1:2,m,e*1.96,'k','linestyle','none','linewidth',2);
    %     set(gca,'ylim',[-.1 .1])
    %     set(gca,'ytick',[-.1:.05:.1]);
    set(gca,'xticklabel',{'Stim', 'Non-stim'});
    ylabel('2nd Half - 1st Half (Raw)','fontsize',20);
    set(gca,'fontsize',20);
    grid on
    set(gca,'gridlinestyle',':')
    title(regions{r});
    fname = fullfile(figSubDir,'beh_byHalf');
    print('-depsc2','-loose',fname);
    figs(r).beh_byHalf = fname;
    
    [h,p,c,s] = ttest(halfDiff);
    
    % Figure 2. second half - first half diff for stim and non stim percent
    % change
    halfDiff = [diff(stim_data,[],2)./stim_data(:,1) diff(nonStim_data,[],2)./nonStim_data(:,1)]*100;
    m = nanmean(halfDiff,1);
    s = nanstd(halfDiff);
    e = s./sqrt(sum(~isnan(halfDiff))-1);
    figure(2)
    clf
    h=bar(1:2,m,'w','linewidth',2);
    hold on
    errorbar(1:2,m,e*1.96,'k','linestyle','none','linewidth',2);
    %     set(gca,'ylim',[-.15 .15]*100)
    %     set(gca,'ytick',[-.15:.05:.15]*100);
    set(gca,'xticklabel',{'Stim', 'Non-stim'});
    ylabel('2nd Half - 1st Half (%)','fontsize',20);
    set(gca,'fontsize',20);
    grid on
    set(gca,'gridlinestyle',':')
    title(regions{r});
    fname = fullfile(figSubDir,'beh_byHalf_percent');
    print('-depsc2','-loose',fname);
    figs(r).beh_byHalf_percent = fname;
    
    [h,p,c,s] = ttest(halfDiff);
    
    
    % Figure 3. change due to stim in each half
    stim_change = stim_data - nonStim_data;
    m = nanmean(stim_change,1);
    s = nanstd(stim_change);
    e = s./sqrt(sum(~isnan(stim_change))-1);
    figure(3)
    clf
    h=bar(1:2,m,'w','linewidth',2);
    hold on
    errorbar(1:2,m,e*1.96,'k','linestyle','none','linewidth',2);
    %     set(gca,'ylim',[-.15 .15])
    %     set(gca,'ytick',[-.15:.05:.15]);
    set(gca,'xticklabel',{'First Half', 'Second Half'});
    ylabel('Raw change due to stim','fontsize',20);
    set(gca,'fontsize',20);
    grid on
    set(gca,'gridlinestyle',':')
    title(regions{r});
    fname = fullfile(figSubDir,'stimEff_byHalf');
    print('-depsc2','-loose',fname);
    figs(r).stimEff_byHalf = fname;
    
    % Figure 4. change due to stim in each half percent
    stim_change = (stim_data - nonStim_data)./nonStim_data;
    m = nanmean(stim_change,1);
    s = nanstd(stim_change);
    e = s./sqrt(sum(~isnan(stim_change))-1);
    figure(4)
    clf
    h=bar(1:2,m,'w','linewidth',2);
    hold on
    errorbar(1:2,m,e*1.96,'k','linestyle','none','linewidth',2);
    %     set(gca,'ylim',[-.15 .15])
    %     set(gca,'ytick',[-.15:.05:.15]);
    set(gca,'xticklabel',{'First Half', 'Second Half'});
    ylabel('% Change due to stim','fontsize',20);
    set(gca,'fontsize',20);
    grid on
    set(gca,'gridlinestyle',':')
    title(regions{r});
    fname = fullfile(figSubDir,'stimEff_byHalf_percent');
    print('-depsc2','-loose',fname);
    figs(r).stimEff_byHalf_percent = fname;    
end
keyboard
%%%% make report
texName = 'YC2_errBySessionHalf.tex';
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
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.49\\textwidth]{%s}}}\n',figs(r).beh_byHalf);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.49\\textwidth]{%s}}}\n',figs(r).beh_byHalf_percent);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.49\\textwidth]{%s}}}\n',figs(r).stimEff_byHalf);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.49\\textwidth]{%s}}}\n',figs(r).stimEff_byHalf_percent);
    fprintf(fid,['\\caption{\\textbf{%d subjects.  }\\textbf{A: } Raw increase in accuracy for 2nd half of sessions. ',...
                 '\\textbf{B: } Percent increase in accuracy for 2nd half of sessions ',...
                 '\\textbf{C: } Raw increase in performance due to stim. \\textbf{D: } Percent increase in performance due to stim.}\n'],subjs_in_region(r));
    fprintf(fid,'\\end{figure}\n\n\n');
    fprintf(fid,'\\clearpage\n\n\n');
    
end
fprintf(fid,'\\end{document}\n\n\n');
fclose(fid);
