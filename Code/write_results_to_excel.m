
% read model file
tmp = load(fullfile('..','Data','ecYeastGEM_batch.mat'));
model = tmp.ecModel_batch; clear tmp

% genotypes = {'bPAC','Clb2','CDC48'};
genotypes = {'Clb2','CDC48'};
conditions = {'light','dark'};

fva_base_name = 'fva_results_';

exclude_idx = startsWith(model.rxnNames,'draw_') | ...
    contains(model.rxnNames,'pseudoreaction') | ismember(model.rxnNames,'growth');


for i=1:numel(genotypes)
    for j=1:numel(conditions)
        fva_file = ['../Results/' genotypes{i} '/' fva_base_name genotypes{i} '_' conditions{j} '.csv'];
        fva_tab = readtable(fva_file);
        diff_idx = fva_tab.distinct_rxn_bool==1 & ~exclude_idx;
        
        diff_rxn_ids = model.rxns(diff_idx);
        
        % get subsystems of diff rxns
        diff_rxn_subsyst = model.subSystems(findRxnIDs(model,diff_rxn_ids));
        diff_rxn_subsyst = cellfun(@(x)regexp(x,'sce\d{5}.*','match'),diff_rxn_subsyst,'un',0);
        for k=1:numel(diff_rxn_subsyst)
            if ~isempty(diff_rxn_subsyst{k})
                diff_rxn_subsyst{k} = strjoin([diff_rxn_subsyst{k}{:}],',');
            else
                diff_rxn_subsyst{k} = '';
            end
        end
        
        conditions = strrep(conditions,'light','blue light');
        
        diff_rxn_subsyst = strtrim(regexprep(diff_rxn_subsyst,'sce\d{5}',''));
        writetable(table(model.rxnNames(diff_idx),diff_rxn_subsyst,...
            fva_tab{diff_idx,2},fva_tab{diff_idx,3},fva_tab{diff_idx,4},fva_tab{diff_idx,5},...
            'VariableNames', {'Reaction name','Subsystem','FVA min (mut) [mmol/gDW/h]',...
            'FVA max (mut) [mmol/gDW/h]', 'FVA min (wt) [mmol/gDW/h]', 'FVA max (wt) [mmol/gDW/h]'}),...
            '../Results/Suppl_Table_X_FVA_Results.xlsx','Sheet',[genotypes{i} ' ' conditions{j}],...
            'Range','A3')
    end
end