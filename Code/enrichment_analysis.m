%% enrichment analysis for metabolic subsystems

% read model file
tmp = load(fullfile('..','Data','ecYeastGEM_batch.mat'));
model = tmp.ecModel_batch; clear tmp

%% Clb2
% light
disp('----------------------------')
disp('Clb2')
disp('----------------------------')
disp('LIGHT')
fvaFile = '../Results/Clb2/fva_results_Clb2_light.csv';
[P_adj_Clb2_L,subSyst_Clb2_L] = getEnrichedSubSyst(fvaFile,model);
disp(table(P_adj_Clb2_L,subSyst_Clb2_L,'VariableNames',{'P_adj','KEGG map'}))
% dark
disp('----------------------------')
disp('DARK')
disp('----------------------------')
fvaFile = '../Results/Clb2/fva_results_Clb2_dark.csv';
[P_adj_Clb2_D,subSyst_Clb2_D] = getEnrichedSubSyst(fvaFile,model);
disp(table(P_adj_Clb2_D,subSyst_Clb2_D,'VariableNames',{'P_adj','KEGG map'}))

%% CDC48
disp('----------------------------')
disp('CDC48')
disp('----------------------------')
% light
disp('LIGHT')
fvaFile = '../Results/CDC48/fva_results_CDC48_light.csv';
[P_adj_CDC48_L,subSyst_CDC48_L] = getEnrichedSubSyst(fvaFile,model);
disp(table(P_adj_CDC48_L,subSyst_CDC48_L,'VariableNames',{'P_adj','KEGG map'}))
% dark
disp('DARK')
fvaFile = '../Results/CDC48/fva_results_CDC48_dark.csv';
[P_adj_CDC48_D,subSyst_CDC48_D] = getEnrichedSubSyst(fvaFile,model);
disp(table(P_adj_CDC48_D,subSyst_CDC48_D,'VariableNames',{'P_adj','KEGG map'}))

%% Combine results and plot heatmap
enr_subs_comb = {subSyst_Clb2_D;subSyst_Clb2_L;...
    subSyst_CDC48_D;subSyst_CDC48_L};
enr_subs_comb_uniq = unique(vertcat(enr_subs_comb{:}));

P_val_comb = {P_adj_Clb2_D;P_adj_Clb2_L;...
    P_adj_CDC48_D;P_adj_CDC48_L};

colors_rgb = [[100 100 100];[90 155 213]]/255;
cond_labels = strcat(repmat({['\color[rgb]{' regexprep(num2str(colors_rgb(1,:)),'\ +',',') '}'],...
    ['\color[rgb]{' regexprep(num2str(colors_rgb(2,:)),'\ +',',') '}']},1,2),...
    {'{\bf Clb2}','{\bf Clb2}','{\bf CDC48}','{\bf CDC48}'});

enr_matrix = repelem(0.05,numel(enr_subs_comb_uniq),numel(cond_labels));

for i=1:numel(cond_labels)
    row_idx = cellfun(@(x)find(ismember(enr_subs_comb_uniq,x)),enr_subs_comb{i});
    enr_matrix(row_idx,i) = P_val_comb{i};
end

h = imagesc(enr_matrix);

xticks(1:numel(cond_labels))
xticklabels(cond_labels)
xtickangle(90)
yticklabels(enr_subs_comb_uniq)

colorbar
ax = struct(gca);
ax.Colorbar.TickLabels(end) = {'> 0.05'};

set(gca,...
    'ticklength',[0,0],...
    'LineWidth',1.3,...
    'Colormap',colormap('summer'),...
    'FontSize',14,...
    'FontName','Arial')

hold on
g = zeros(2,1);
g(1) = plot(NaN,NaN,'markersize',150,'MarkerFaceColor',colors_rgb(1,:),'Visible','off',...
    'Marker','s','MarkerEdgeColor',colors_rgb(1,:),'LineStyle','none');
g(2) = plot(NaN,NaN,'markersize',150,'MarkerFaceColor',colors_rgb(2,:),'Visible','off',...
    'Marker','s','MarkerEdgeColor',colors_rgb(2,:),'LineStyle','none');

l = legend(g,{'\color{black}darkness','\color{black}blue light'},'location','northoutside','NumColumns',2,...
    'box','off','Fontsize',14);
l.FontName = 'Arial';

set(gcf,'OuterPosition',1000*[-1.0163    0.0097    0.6140    0.4973])

exportgraphics(gcf,'../Results/enrichment_non_overlapping_rxns.tiff','Resolution',300)

function [P_adj,subSystNames] = getEnrichedSubSyst(fvaFile,model,alpha)
if nargin<3;alpha=0.05;end

subsystems = cellfun(@(x)regexp(x,'sce\d{5}','match'),model.subSystems,'un',0);
subsyst_per_rxn = subsystems;
for i=1:numel(subsystems)
    if ~isempty(subsystems{i})
        subsyst_per_rxn{i} = strjoin([subsystems{i}{:}],',');
    else
        subsyst_per_rxn{i} = '';
    end
end
subsystems = [subsystems{:}]';
unique_subsyst = unique([subsystems{:}]');

fva_results = readtable(fvaFile);
diff_rxn_names = fva_results.reaction_name(fva_results.distinct_rxn_bool==1);
% exclude protein draw reactions, pseudoreactions, and growth sink reaction
exclude_idx = startsWith(diff_rxn_names,'draw_') | ...
    contains(diff_rxn_names,'pseudoreaction') | ismember(diff_rxn_names,'growth');

diff_rxn_ids = model.rxns(fva_results.distinct_rxn_bool==1);
diff_rxn_ids = diff_rxn_ids(~exclude_idx);

% get subsystems of diff rxns
diff_rxn_subsyst = model.subSystems(findRxnIDs(model,diff_rxn_ids));
diff_rxn_subsyst = cellfun(@(x)regexp(x,'sce\d{5}','match'),diff_rxn_subsyst,'un',0);
for i=1:numel(diff_rxn_subsyst)
    if ~isempty(diff_rxn_subsyst{i})
        diff_rxn_subsyst{i} = strjoin([diff_rxn_subsyst{i}{:}],',');
    else
        diff_rxn_subsyst{i} = '';
    end
end

% perform hypergeometric test
%                | in subsystem   | not in subsystem | Total
% ---------------------------------------------------|--------
%  diff rxn      |      n11       |     n12          | n1+
%  not diff rxn  |      n21       |     n22          | n2+
% ---------------------------------------------------|--------
%  Total         |      n+1       |     n+2          | n

n_1_plus = numel(diff_rxn_names);
P = nan(size(unique_subsyst));
n_test = 0;
for i=1:numel(unique_subsyst)
    x = sum(contains(diff_rxn_subsyst,unique_subsyst(i)));
    if x > 0
        n_test = n_test + 1;
        n_plus_1 = sum(contains(subsyst_per_rxn,unique_subsyst(i)));
        n_plus_2 = sum(~contains(subsyst_per_rxn,unique_subsyst(i)));
        n = n_plus_1 + n_plus_2;       
        % calculate one-sided mid-P-value
        P(i) = hygecdf(x-1,n,n_plus_1,n_1_plus,'upper') + 0.5*hygepdf(x,n,n_plus_1,n_1_plus);
        
    end
end

% Bonferroni correction
P_adj = P*n_test;
subSyst = unique_subsyst(P_adj<alpha);

subSystNames = cell(size(subSyst));
for i=1:numel(subSyst)
    subsyst_idx = find(contains(subsyst_per_rxn,subSyst{i}),1);
    tmp = regexp(model.subSystems{subsyst_idx},[subSyst{i} '.*'],'match','once');
    subSystNames(i) = strtrim(strrep(tmp(~cellfun(@isempty,tmp)),subSyst{i},''));
end

P_adj = P_adj(P_adj<alpha);
[P_adj,ia] = sort(P_adj,'ascend');
subSystNames = subSystNames(ia);
end
