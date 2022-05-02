% Analysis for bPAC
clear

% add path to GECKO toolbox
GECKO_PATH = '../../GECKO';
% for server
addpath(genpath('/home/mpimp-golm.mpg.de/wendering1302/bin/gurobi911/linux64/matlab/'))

NCPU = 20;
% NCPU = 1;

% write data files to csv
data_file_base_name = 'bPAC_global';

%{
excel_file = ['../Data/' data_file_base_name '.ods'];
prot_abundance_sheet = 'bPAC_global';
diff_genes_sheet = 'D';

% abundance table
tab = readtable(excel_file,'Sheet',prot_abundance_sheet,'range','A:Y','readrownames',true);
tab = tab(~startsWith(tab.Properties.RowNames,'Row'),:);
writetable(tab,['../Data/' data_file_base_name '_abundances.csv'],'WriteRowNames',true);

% differential genes
tab = readtable(excel_file,'Sheet',diff_genes_sheet,'range','A:A');
tab = tab(1:end-1,:);
writetable(tab,['../Data/' data_file_base_name '_D.csv']);
%}

% if not existing, create output directory
out_dir = fullfile('..','Results',strtok(data_file_base_name,'_'));
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end


% read yeast batch model obtained from github.com/SysBioChalmers/ecModels (2022/03/10)
tmp = load(fullfile('..','Data','ecYeastGEM_batch.mat'));
model = tmp.ecModel_batch; clear tmp
% add rules field
model_cobra = ravenCobraWrapper(model);
model.rules = model_cobra.rules;

% define medium
medium = {'Adenine','L-Alanine','L-Arginine','L-Asparagine','L-Aspartate',...
    'L-Cysteine','L-Glutamine','L-Glutamate','L-Glycine','myo-Inositol','L-Isoleucine',...
    'L-Lysine','L-Methionine','4-aminobenzoate','L-Phenylalanine','L-Proline',...
    'L-Serine','L-Threonine','L-Tyrosine','L-Valine','L-Histidine','L-Leucine',...
    'Uracil','L-Tryptophan','Biotin','(R)-pantothenate','nicotinate',...
    'pyridoxine','thiamine(1+)'};

uptake_rxns = cellfun(@(x)model.rxns(startsWith(model.rxnNames,[x ' exchange (reversible)'],'IgnoreCase',true)),...
    medium,'un',0);
uptake_rxns = [uptake_rxns{:}]';

model.ub(findRxnIDs(model,uptake_rxns)) = Inf;

% set unwanted nutrients to zero (ammonia uptake)
model.ub(findRxnIDs(model,{'r_1654_REV'})) = 0;

% change default reaction flux limit to 1000
model.ub(isinf(model.ub)) = 1000;

% adapt biomass composition to create a wild type and a mutant model
biomass = find(contains(model.rxnNames,'biomass'));

DNA = find(strcmp(model.metNames,'DNA'));
RNA = find(strcmp(model.metNames,'RNA'));
carbohydrate = find(strcmp(model.metNames,'carbohydrate'));
protein = find(strcmp(model.metNames,'protein'));
cofactor = find(strcmp(model.metNames,'cofactor'));
ion = find(strcmp(model.metNames,'ion'));
lipid = find(strcmp(model.metNames,'lipid'));

% WT biomass composition
wild_protein = 0.236;
wild_carbohydrate = 0.437;
wild_DNA = 0.0013;
wild_RNA = 0.035;
wild_lipid = 0.064;

% (general?) mutant biomass composition
mut_protein = 0.317;
mut_carbohydrate = 0.341;
mut_DNA = 0.0011;
mut_RNA = 0.057;
mut_lipid = 0.076;

% get default biomass composition
orig_dir = pwd;
cd(fileparts(which('sumBioMass')))
[~,P_default,C_default,R_default,D_default,L_default] = sumBioMass(model);
cd(orig_dir)

sum_biomass_default = P_default + C_default + R_default + D_default + L_default;

% ceate WT-specific model

% re-scale measured biomass composition to match the same sum of biomass
sum_wt_biomass = sum([wild_protein,wild_carbohydrate,wild_DNA,wild_RNA,wild_lipid]);
wild_protein = sum_biomass_default * wild_protein/sum_wt_biomass;
wild_carbohydrate = sum_biomass_default * wild_carbohydrate/sum_wt_biomass;
wild_DNA = sum_biomass_default * wild_DNA/sum_wt_biomass;
wild_RNA = sum_biomass_default * wild_RNA/sum_wt_biomass;
wild_lipid = sum_biomass_default * wild_lipid/sum_wt_biomass;

model_wild = model;
model_wild.S(protein,biomass) = wild_protein/P_default * model_wild.S(protein,biomass);
model_wild.S(carbohydrate,biomass) = wild_carbohydrate/C_default * model_wild.S(carbohydrate,biomass);
model_wild.S(DNA,biomass) = wild_DNA/D_default * model_wild.S(DNA,biomass);
model_wild.S(RNA,biomass) = wild_RNA/R_default * model_wild.S(RNA,biomass);
model_wild.S(lipid,biomass) = wild_lipid/L_default * model_wild.S(lipid,biomass);

% create mutant-specific model
% re-scale measured biomass composition to match the same sum of biomass
sum_mut_biomass = sum([mut_protein,mut_carbohydrate,mut_DNA,mut_RNA,mut_lipid]);
mut_protein = sum_biomass_default * mut_protein/sum_mut_biomass;
mut_carbohydrate = sum_biomass_default * mut_carbohydrate/sum_mut_biomass;
mut_DNA = sum_biomass_default * mut_DNA/sum_mut_biomass;
mut_RNA = sum_biomass_default * mut_RNA/sum_mut_biomass;
mut_lipid = sum_biomass_default * mut_lipid/sum_mut_biomass;

model_mut = model;
model_mut.S(protein,biomass) = mut_protein/P_default * model_mut.S(protein,biomass);
model_mut.S(carbohydrate,biomass) = mut_carbohydrate/C_default * model_mut.S(carbohydrate,biomass);
model_mut.S(DNA,biomass) = mut_DNA/D_default * model_mut.S(DNA,biomass);
model_mut.S(RNA,biomass) = mut_RNA/R_default    * model_mut.S(RNA,biomass);
model_mut.S(lipid,biomass) = mut_lipid/L_default * model_mut.S(lipid,biomass);

% read proteomics data (only gene IDs)
load(fullfile(GECKO_PATH,'databases','ProtDatabase.mat'));
tab = readtable(['../Data/' data_file_base_name '_abundances.csv'],...
    'ReadRowNames',1);
msrd_proteins = tab.Properties.RowNames; clear tab

% match gene IDs from experiment to SwissProt IDs
msrd_proteins_swissprot = cell(length(msrd_proteins),1);

for i = 1:length(msrd_proteins)
    msrd_proteins{i} = strrep(msrd_proteins{i},'_SK1','');
end

for i = 1:length(msrd_proteins)
    if ~isempty(find(strcmp(kegg(:,3),msrd_proteins{i})))
        msrd_proteins_swissprot{i} = kegg{find(strcmp(kegg(:,3),msrd_proteins{i})),1};
    end
end
clear kegg swissprot

% map protein names to model
inx_msrd_proteins_swissprot = zeros(length(msrd_proteins_swissprot),1);
for i = 1:length(msrd_proteins_swissprot)
    for j = find(startsWith(model.rxns,'draw_prot_'))'
        prot_model = model.rxnNames{j}(11:end);
        if strcmp(prot_model,msrd_proteins_swissprot{i})
            inx_msrd_proteins_swissprot(i) = j;
        end
    end
end

% read proteomics data (abundances)
protein_abundance_raw = readmatrix(['../Data/' data_file_base_name '_abundances.csv'],...
    'Range',[2,2]);
% take mean of replicates
protein_abundance_wild_D = mean(protein_abundance_raw(:,1:6),2);
protein_abundance_wild_L = mean(protein_abundance_raw(:,7:12),2);
protein_abundance_mut_D = mean(protein_abundance_raw(:,13:18),2);
protein_abundance_mut_L = mean(protein_abundance_raw(:,19:24),2);

% for each protein calculate the ratio between abundances in mutant and WT
% (in light and dark conditions)
protein_abundance_ratio_L = protein_abundance_mut_L ./ protein_abundance_wild_L;
protein_abundance_ratio_D = protein_abundance_mut_D ./ protein_abundance_wild_D;

clear protein_abundance_raw protein_abundance_mut_L protein_abundance_wild_L ...
    protein_abundance_mut_D protein_abundance_wild_D

% ratio of biomass (mutant/WT)
biomass_ratio_L = 3.0/2.8;
biomass_ratio_D = 2.8/3.1;

% protein content (why no the same value as for biomass composition?)
protein_ratio_m_L = 21; % Light
protein_ratio_m_D = 18.3; % Dark

protein_ratio_w_L = 18.8; % Light
protein_ratio_w_D = 18.8; % Dark

% update protein pool constraint
% function_path = which('getModelParameters.m');
% status = system(...
%     ['sed -i ''s/parameters.Ptot = [0-9]*\.[0-9]*/parameters.Ptot = '...
%     num2str(protein_ratio_m/100) '/'' ' function_path]);

% fit GAM
% orig_dir = pwd;
% cd(fileparts(which('constrainEnzymes')))
% [~,~,~,GAM_mut_L] = constrainEnzymes(model_mut,0.4421,[],protein_ratio_m_L/100);
% [~,~,~,GAM_mut_D] = constrainEnzymes(model_mut,0.4421,[],protein_ratio_m_D/100);
% [~,~,~,GAM_wt_L] = constrainEnzymes(model_wild,0.4421,[],protein_ratio_w_L/100);
% [~,~,~,GAM_wt_D] = constrainEnzymes(model_wild,0.4421,[],protein_ratio_w_D/100);
% cd(orig_dir)

GAM_mut_L = 75.6;
GAM_mut_D = 75.2;
GAM_wt_L = 70.1;
GAM_wt_D = 70.1;

% light
fprintf('\n\nLIGHT\n\n')
GAM = struct('mut',GAM_mut_L,'wt',GAM_wt_L);
result = RENZO(model_mut,model_wild,msrd_proteins_swissprot,...
    protein_ratio_m_L, protein_ratio_w_L,...
    protein_abundance_ratio_L,biomass_ratio_L,GAM,NCPU);

fprintf('infeasible min: %d\n',result.infeasibility_min)
fprintf('infeasible max: %d\n',result.infeasibility_max)

results_light = {...
    result.FVA_min(1:length(model.rxns)),...
    result.FVA_max(1:length(model.rxns)),...
    result.FVA_min(length(model.rxns) + 1:2*length(model.rxns)),...
    result.FVA_max(length(model.rxns) + 1:2*length(model.rxns))};

writetable(...
    [cell2table(model.rxnNames,'VariableNames',{'reaction_name'}),...
    array2table([cell2mat(results_light) result.distinct_rxns],'VariableNames',...
        {'minFlux_mut','maxFlux_mut','minFlux_wt','maxFlux_wt','distinct_rxn_bool'})],...
    [out_dir filesep 'fva_results_' strtok(data_file_base_name,'_') '_light.csv'],...
    'QuoteStrings',false,'Delimiter','\t')

Genes_in_distinc_rxns_light = result.Genes_in_distinc_rxns;

clear results_light result

% dark
fprintf('\n\nDARK\n\n')
GAM = struct('mut',GAM_mut_D,'wt',GAM_wt_D);
result = RENZO(model_mut,model_wild,msrd_proteins_swissprot,...
    protein_ratio_m_D, protein_ratio_w_D,...
    protein_abundance_ratio_D,biomass_ratio_D,GAM,NCPU);

fprintf('infeasible min: %d\n',result.infeasibility_min)
fprintf('infeasible max: %d\n',result.infeasibility_max)

results_dark = {...
    result.FVA_min(1:length(model.rxns)),...
    result.FVA_max(1:length(model.rxns)),...
    result.FVA_min(length(model.rxns) + 1:2*length(model.rxns)),...
    result.FVA_max(length(model.rxns) + 1:2*length(model.rxns))};

writetable(...
    [cell2table(model.rxnNames,'VariableNames',{'reaction_name'}),...
    array2table([cell2mat(results_dark) result.distinct_rxns],'VariableNames',...
        {'minFlux_mut','maxFlux_mut','minFlux_wt','maxFlux_wt','distinct_rxn_bool'})],...
    [out_dir filesep 'fva_results_' strtok(data_file_base_name,'_') '_dark.csv'],...
    'QuoteStrings',false,'Delimiter','\t')

Genes_in_distinc_rxns_dark = result.Genes_in_distinc_rxns;

clear result
