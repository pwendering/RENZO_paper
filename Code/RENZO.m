function results = RENZO(model_mt,model_wt,msrd_protein_ids,Ptot_mt, Ptot_wt,...
    protein_abundance_ratio,biomass_ratio,GAM,NCPU,bio_name)
%% results = RENZO(model_mt,model_wt,msrd_protein_ids,Ptot_mt, Ptot_wt,protein_abundance_ratio,biomass_ratio,GAM,NCPU,bio_name)
% Impose protein abundance ratios on protein abundances in mutant and wild
% type ecModel (GECKO).
% Input:
%   struct model_mt:                        mutant-specific model
%   struct model_wt:                        wild type-specific model
%   cellstr msrd_protein_ids:               Swissprot IDs of measured proteins
%   double Ptot_mt:                         total protein content (mutant)
%   double Ptot_wt:                         total protein content (wild type)
%   double protein_abundance_ratio:         ratios of protein abundances (mutant/wild type)
%   double biomass_ratio:                   ratio between mutant and wild type biomass
%   struct GAM:                             growth-associated maintenance
%   double NCPU:                            (opt) number of workers to use
%   char bio_name:                          (opt) name of biomass reaction in model
% 
% Output:
%   struct results:                         results from flux variability
%                                           analysis
% 
% Dr. Zahra Razaghi Moghadam, Philipp Wendering (University of Potsdam)

% check input
if nargin < 9
    NCPU = 1;
elseif NCPU > 1
    p = gcp('nocreate');
    if isempty(p) || p.NumWorkers ~= NCPU
        delete(p);
        parpool(NCPU);
    end
end

% determine biomass index
if nargin < 10 || isempty(bio_name)
    bio_name = 'biomass';
end
bio_idx = find(contains(model_wt.rxnNames,bio_name));

if isempty(bio_idx)
    error('Biomass reaction index for %s could not be determined.',bio_name)
elseif numel(bio_idx) > 1
    error('Could not determine biomass reaction index: multiple matches for %s',bio_name)
end

% check if wild type and mutant models have the same set of reactions
if ~isequal(model_mt.rxns,model_wt.rxns)
    error('Reaction sets of wild type and mutant models are not identical.')
end

% number of reactions in the model
NRXNS = numel(model_wt.rxns);

% map protein names to model
inx_msrd_proteins_swissprot = zeros(length(msrd_protein_ids),1);
for i = 1:length(msrd_protein_ids)
    for j = find(startsWith(model_wt.rxns,'draw_prot_'))'
        prot_model = model_wt.rxnNames{j}(11:end);
        if strcmp(prot_model,msrd_protein_ids{i})
            inx_msrd_proteins_swissprot(i) = j;
        end
    end
end

% update protein pool for wild type and mutant model
orig_dir = pwd;
cd(fileparts(which('constrainEnzymes')))
model_mt = constrainEnzymes(model_mt,0.4421,GAM.mut,Ptot_mt/100);
model_wt = constrainEnzymes(model_wt,0.4421,GAM.wt,Ptot_wt/100);
cd(orig_dir)

% construct LP
num_added_variables = nnz(inx_msrd_proteins_swissprot);

Aeq = [...
    model_mt.S,zeros(size(model_wt.S)),zeros(size(model_mt.S,1),num_added_variables);...
    zeros(size(model_mt.S)),model_wt.S,zeros(size(model_mt.S,1),num_added_variables)];
beq = zeros(size(Aeq,1),1);

Aineq = [];
bineq = [];
for i = 1:num_added_variables
    idx = find(inx_msrd_proteins_swissprot ~= 0,i,'first');
    inx = idx(end);
    
    Aineq(end+1,inx_msrd_proteins_swissprot(inx)) = -1;
    Aineq(end,NRXNS + inx_msrd_proteins_swissprot(inx)) = protein_abundance_ratio(inx);
    Aineq(end,2*NRXNS + i) = -1;
    
    bineq(end+1,1) = 0;
    
    Aineq(end+1,inx_msrd_proteins_swissprot(inx)) = 1;
    Aineq(end,NRXNS + inx_msrd_proteins_swissprot(inx)) = -protein_abundance_ratio(inx);
    Aineq(end,2*NRXNS + i) = -1;
    
    bineq(end+1,1) = 0;
end

% biomass optimization
f_biomass = zeros(2*NRXNS + num_added_variables,1);
f_biomass(bio_idx + NRXNS) = -1;

m = struct();
m.obj = f_biomass';
m.A = [sparse(Aineq);sparse(Aeq)];
n = size(m.A, 2);
m.vtype = repmat('C', n, 1);
m.sense = [repmat('<',size(Aineq,1),1); repmat('=',size(Aeq,1),1)];
m.rhs = full([bineq(:); beq(:)]);
m.lb = [model_mt.lb;model_wt.lb;zeros(num_added_variables,1)];
m.ub = [model_mt.ub;model_wt.ub;1000*ones(num_added_variables,1)];

params = struct();
params.FeasibilityTol = 1e-9;
params.OutputFlag = 0;

x_biomass = gurobi(m,params);

% norm 1 optimization
f = zeros(2*NRXNS + num_added_variables,1);
f(2*NRXNS+1:2*NRXNS+num_added_variables) = 1;

Aeq(end+1,bio_idx) = 1;
Aeq(end,NRXNS + bio_idx) = -biomass_ratio;
beq(end+1) = 0;

m = struct();
m.obj = f';
m.A = [sparse(Aineq);sparse(Aeq)];
n = size(m.A, 2);
m.vtype = repmat('C', n, 1);
m.sense = [repmat('<',size(Aineq,1),1); repmat('=',size(Aeq,1),1)];
m.rhs = full([bineq(:); beq(:)]);
m.lb = [model_mt.lb;model_wt.lb;zeros(num_added_variables,1)];
m.ub = [model_mt.ub;model_wt.ub;1000*ones(num_added_variables,1)];
m.lb(bio_idx + NRXNS) = 0.85 * x_biomass.x(bio_idx + NRXNS);

x_D = gurobi(m,params);

% FVA

% add lower bounds for minimized difference of the ratio of predicted protein
% abundances for WT and mutant to the measured ratio
for i = 1:num_added_variables
    Aineq(end+1,2*NRXNS + i) = 1;
    bineq(end+1,1) = x_D.x(2*NRXNS + i) + 1e-5;
end

% update LP
m.A = [sparse(Aineq);sparse(Aeq)];
m.sense = [repmat('<',size(Aineq,1),1); repmat('=',size(Aeq,1),1)];
m.rhs = full([bineq(:); beq(:)]);

fprintf('\n\nStarting FVA for %d reactions...\n\n',2*NRXNS)
FVA_min = zeros(2*NRXNS,1);
FVA_max = zeros(2*NRXNS,1);
infeasibility_min = 0;
infeasibility_max = 0;

if NCPU > 1
    params.Threads = 1;
    parfor i = 1:2*NRXNS
        % minimization
        f_FVA = zeros(2*NRXNS + num_added_variables,1);
        f_FVA(i) = 1;
        tmpProblem = m;
        tmpProblem.obj = f_FVA';
        x_FVA = gurobi(tmpProblem,params);
        if any(strcmp(x_FVA.status,{'INFEASIBLE','INF_OR_UNBD'}))
            FVA_min(i) = 0;
            infeasibility_min = infeasibility_min + 1;
        else
            FVA_min(i) = x_FVA.x(i);
        end
        
        % maximization
        f_FVA = zeros(2*NRXNS + num_added_variables,1);
        f_FVA(i) = -1;
        tmpProblem.obj = f_FVA';
        x_FVA = gurobi(tmpProblem,params);
        if any(strcmp(x_FVA.status,{'INFEASIBLE','INF_OR_UNBD'}))
            FVA_max(i) = 1000;
            infeasibility_max = infeasibility_max + 1;
        else
            FVA_max(i) = x_FVA.x(i);
        end
    end
else
    for i = 1:2*length(model_wt.rxns)
        if mod(i,100)==0; fprintf('Processed %d reactions\n',i);end
        
        % minimization
        f_FVA = zeros(2*NRXNS + num_added_variables,1);
        f_FVA(i) = 1;
        tmpProblem = m;
        tmpProblem.obj = f_FVA';
        x_FVA = gurobi(tmpProblem,params);
        if any(strcmp(x_FVA.status,{'INFEASIBLE','INF_OR_UNBD'}))
            FVA_min(i) = 0;
            infeasibility_min = infeasibility_min + 1;
        else
            FVA_min(i) = x_FVA.x(i);
        end
        
        % maximization
        f_FVA = zeros(2*NRXNS + num_added_variables,1);
        f_FVA(i) = -1;
        tmpProblem.obj = f_FVA';
        x_FVA = gurobi(tmpProblem,params);
        if any(strcmp(x_FVA.status,{'INFEASIBLE','INF_OR_UNBD'}))
            FVA_max(i) = 1000;
            infeasibility_max = infeasibility_max + 1;
        else
            FVA_max(i) = x_FVA.x(i);
        end
    end
end

FVA_min(abs(FVA_min) < 1e-6) = 0;
FVA_max(abs(FVA_max) < 1e-6) = 0;

% find distinct flux ranges
distinct_rxns = zeros(NRXNS,1);
for i = 1:NRXNS
    if FVA_min(i) > FVA_max(i + NRXNS)
        distinct_rxns(i) = 1;
    end
    if FVA_min(i + NRXNS) > FVA_max(i)
        distinct_rxns(i) = 1;
    end
end

% find genes associated to reactions with non-overlapping flux ranges
Genes_in_distinc_rxns = [];
for i = 1:NRXNS
    if distinct_rxns(i) == 1
        if ~isempty(model_wt.rules(i))
            t = char(model_wt.rules{i});
            x_pos = strfind(t,'x');
            while ~isempty(x_pos)
                k = x_pos(1);
                while t(k) ~= ')'
                    k = k + 1;
                end
                Genes_in_distinc_rxns = union(str2num(t(x_pos(1)+2:k-1)),Genes_in_distinc_rxns);
                t = strrep(t,t(x_pos(1):k),model_wt.genes{str2num(t(x_pos(1)+2:k-1))});
                x_pos = [];
                x_pos = strfind(t,'x');
            end
        end
    end
end

results = struct('FVA_min',FVA_min,'FVA_max',FVA_max,...
    'infeasible_min',infeasibility_min,'infeasible_max',infeasibility_max,...
    'distinct_rxns',distinct_rxns,'genes_in_distinct_rxns',Genes_in_distinc_rxns);
end
