function [simtable guideLevelRates GeneLevelRates] = MEDUSA(ugr,dgr,ddr,onset,endpoint,L2FC_dataset,sgRNA_num)

% Full simulation %
% Rates for "nontargeting" sgRNAs
untreated_growth = log2(ugr); % untreated growth rate
treated_growth = log2(dgr/1)./onset; % growth rate with drug - first phase
treated_death = log2(dgr/ddr)./(endpoint-onset); %death rate with drug - second phase

% SGKO rate ranges
ko_growth = sort([linspace(0, 1.2, 498),.5,1]');  % changes in KO relative growth rate
ko_death = sort([linspace(0, 4, 495),.5,1,1.5,2.5,3.5]'); % changes in KO relative death rate

% Times
x = linspace(0,endpoint,500)'; % vector of timepoints
if isempty(find(x == onset))
    x = sort([linspace(0,endpoint,499),onset]'); % vector of timepoints
end

x0 = onset; % death 'onset'
start = 300; % arbitrary "cell" starting number

% Untreated - nontargetings
gdt = untreated_growth;
tauDG = x .* gdt;  
y = start .* 2.^(tauDG);
xuntreated = x;
yuntreated = y;
clear y tauDG gdt

% Treated - nontargetings
peak = find(x == x0);
% growth phase
gdt = treated_growth;
tauDG = x(1:peak) .* gdt;
y(1:peak) = start .* 2.^(tauDG);  
% death phase
ddt = treated_death;
tauDD = -(x(peak+1:end) - x0) * ddt;
y(peak+1:length(x)) = y(peak) * 2.^(tauDD);
xtreated = x;
ytreated = y';
clear y peak tauDG tauDD gdt

% Array for simulation
qmat = NaN(length(ko_death)*length(ko_growth),1);

w = 0;
for iii = 1:length(ko_death)

for ii = 1:length(ko_growth) 
    w = w+1;  

    % Untreated KOs
    gdt = untreated_growth .* ko_growth(ii);
    tauDG = x .* gdt;
    y = start .* 2.^(tauDG);   
    kountreated = y;
    clear y i tauDG gdt

    % Treated KOs (two phases, same scale)
    peak = find(x == onset);
    % growth phase
    gdt = treated_growth .* ko_growth(ii);
    tauDG = x(1:peak) .* gdt;
    y(1:peak) = start .* 2.^(tauDG);  
    % death phase  
    ddt = treated_death .* ko_death(iii);
    tauDD = -(x(peak+1:end)-onset) .* ddt;
    y(peak+1:length(x)) = y(peak) .* 2.^(tauDD);
    koetop = y';
    clear i peak tauDG tauDD
    
    log2utvt0 = log2((kountreated(end)./yuntreated(end)) ./ (kountreated(1)./yuntreated(1)));
    log2trvut = log2((koetop(end)./ytreated(end)) ./ (kountreated(end)./yuntreated(end)));

    qmat(w,1) = gdt./treated_growth; % relative growth rate
    qmat(w,2) = ddt./treated_death;  % relative death rate
    qmat(w,3) = log2utvt0;     % Untreated/T0
    qmat(w,4) = log2trvut;     % Etop/Untreated
    
    clear y gdt ddt kountreated koetop log2utvt0 log2trvut
end

clear i ii
end
clear iii w

clear untreated_growth treated_growth treated_death ko_growth ko_death start x x0 xuntreated yuntreated xtreated ytreated

simtable = [array2table(qmat)];
simtable.Properties.VariableNames = [{'GR_relative'},{'DR_relative'},{'sim_UTvT0'},{'sim_TRvUT'}];

clear qmat 

% Relative growth rate (untreated vs T0) %
q = L2FC_dataset.UTvT0;

GrowthTable = simtable(find(simtable.DR_relative == 1),[1,3]);
GrowthTable = sortrows(GrowthTable,'sim_UTvT0','ascend');

fc = GrowthTable.sim_UTvT0;
for i = 1:size(q,1)
    
    qfc = q(i);
    
    rGR(i,1) = GrowthTable.GR_relative(find(abs(fc-qfc) == min(abs(fc-qfc))));
  
    clear qfc
end
L2FC_dataset.Relative_GR = rGR;

clear q GrowthTable rGR fc i 

% Relative death rate (treated vs untreated) %
q = L2FC_dataset.TRvUT;

for i = 1:size(q,1)
    qg = L2FC_dataset.Relative_GR(i);
    qfc = q(i);
    
    qd = simtable((find(simtable.GR_relative == qg)),:);
    fc = qd.sim_TRvUT;
   
    rD(i,1) = qd.DR_relative(find(abs(fc-qfc) == min(abs(fc-qfc))));
  
    clear qg qfc qd fc
end
L2FC_dataset.Relative_DR = rD;

clear q rD i

% Collapse to gene-level scores %
% Mean gene-level score
unigene = unique(L2FC_dataset.Gene);

GeneRates = table();
GrowthRate = NaN(length(unigene),sgRNA_num);
DeathRate = NaN(length(unigene),sgRNA_num);
for i = 1:length(unigene)
    
    idx = find(strcmp(L2FC_dataset.Gene, unigene(i)));
    
    Gene(i,1) =  unigene(i);
    
    GrowthRate(i,1:length(idx)) = L2FC_dataset.Relative_GR(idx)';
    MeanGrowth(i,1) = mean(L2FC_dataset.Relative_GR(idx));
    
    DeathRate(i,1:length(idx)) = L2FC_dataset.Relative_DR(idx)';
    MeanDeath(i,1) = mean(L2FC_dataset.Relative_DR(idx));
    
    numgenes(i,1) = length(idx);

    clear idx  
end

GeneRates.Gene = Gene;
GeneRates.GrowthRate = GrowthRate;
GeneRates.MeanGrowth = MeanGrowth;

% zscore growth rate to nontargetings
nontrows = find(startsWith(GeneRates.Gene,'Nont'));
meannont1 = mean(GeneRates.MeanGrowth(nontrows));
stdnont1 = std(GeneRates.MeanGrowth(nontrows)); 
GeneRates.zMeanGrowth = (GeneRates.MeanGrowth - meannont1)./stdnont1;

GeneRates.DeathRate = DeathRate;
GeneRates.MeanDeath = MeanDeath;

% zscore death rate to nontargetings
nontrows = find(startsWith(GeneRates.Gene,'Nont'));
meannont1 = mean(GeneRates.MeanDeath(nontrows));
stdnont1 = std(GeneRates.MeanDeath(nontrows)); 
GeneRates.zMeanDeath = (GeneRates.MeanDeath - meannont1)./stdnont1;

clear nontrows meannont1 stdnont1 unigene i Gene GrowthRate MeanGrowth DeathRate MeanDeath


% Empiric pvalue and FDR correction (GROWTH) %
% Iterate for empiric pvalue
sgrnaFC = L2FC_dataset.Relative_GR;
RawGenes = GeneRates.Gene;

% do this 10000 times
iterations = 10000; %10000;
% pre-allocate
zpermMean=zeros(length(numgenes),iterations);
perm_mean=zeros(length(numgenes),1);

for j = 1:iterations
    x = randi(length(sgrnaFC),[(length(RawGenes)*sgRNA_num),1]); 
    perm_1 = reshape(sgrnaFC(x),[length(RawGenes),sgRNA_num]);

    secrettable2 = array2table([numgenes,perm_1]);
    secrettable2.Properties.VariableNames(1) = {'numgenes'};
    secrettable2.meanlfc = zeros(size(secrettable2,1),1);

    for i = 1:max(numgenes)
        posit = ismember(secrettable2.numgenes, i);
        secrettable2{posit,'meanlfc'} = mean(secrettable2{posit,2:(i+1)},2);
        
        clear posit
    end
    perm_mean(:) = secrettable2{:,'meanlfc'};

    % z score transform this (1 iteration of mean random mean)
    permmean1 = mean(perm_mean);    %mean of whole matrix
    permstd1 = std(perm_mean);      %std of whole matrix
    zpermMean(:,j) = (perm_mean - permmean1)./permstd1;

    clear x perm_1 secrettable2 i permmean1 permstd1 perm_mean 
end
clear i j jj q* iterations perm_mean perm_1 secret* 

% zscore
zMeanFC(:,1) = GeneRates.zMeanGrowth;

% calculate p-values    
pval_mean = [];
for i = 1:size(zMeanFC,1)
    qq = [zMeanFC(i,1),zpermMean(i,:)]';
    abz = abs(qq);
    [qsort,qidx] = sort(abz,'descend');
    qnum = find(qidx ==1);
    pval_mean(i,1) = qnum./length(qidx);
    
    clear qq abz qsort qidx qnum
end

% FDR correction
GeneRates.Growth_BHFDR = mafdr(pval_mean(:,1),'BHFDR',true);

clear sgrnaFC RawGenes iterations zpermMean perm_mean i pval_mean zMeanFC


% Empiric pvalue and FDR correction (DEATH) %
% Iterate for empiric pvalue
sgrnaFC = L2FC_dataset.Relative_DR;
RawGenes = GeneRates.Gene;

% do this 10000 times
iterations = 10000; %10000;
% pre-allocate
zpermMean=zeros(length(numgenes),iterations);
perm_mean=zeros(length(numgenes),1);

for j = 1:iterations
    x = randi(length(sgrnaFC),[(length(RawGenes)*sgRNA_num),1]); 
    perm_1 = reshape(sgrnaFC(x),[length(RawGenes),sgRNA_num]);

    secrettable2 = array2table([numgenes,perm_1]);
    secrettable2.Properties.VariableNames(1) = {'numgenes'};
    secrettable2.meanlfc = zeros(size(secrettable2,1),1);

    for i = 1:max(numgenes)
        posit = ismember(secrettable2.numgenes, i);
        secrettable2{posit,'meanlfc'} = mean(secrettable2{posit,2:(i+1)},2);
        
        clear posit
    end
    perm_mean(:) = secrettable2{:,'meanlfc'};

    % z score transform this (1 iteration of mean random mean)
    permmean1 = mean(perm_mean);    %mean of whole matrix
    permstd1 = std(perm_mean);      %std of whole matrix
    zpermMean(:,j) = (perm_mean - permmean1)./permstd1;

    clear x perm_1 secrettable2 i permmean1 permstd1 perm_mean 
end
clear i j jj q* iterations perm_mean perm_1 secret* 

% zscore
zMeanFC(:,1) = GeneRates.zMeanDeath;

% calculate p-values    
pval_mean = [];
for i = 1:size(zMeanFC,1)
    qq = [zMeanFC(i,1),zpermMean(i,:)]';
    abz = abs(qq);
    [qsort,qidx] = sort(abz,'descend');
    qnum = find(qidx ==1);
    pval_mean(i,1) = qnum./length(qidx);
    
    clear qq abz qsort qidx qnum
end

% FDR correction
GeneRates.Death_BHFDR = mafdr(pval_mean(:,1),'BHFDR',true);

clear sgrnaFC RawGenes iterations zpermMean perm_mean i pval_mean zMeanFC

guideLevelRates = L2FC_dataset;
GeneLevelRates = GeneRates(:,[1,3,6,8,9]);

clear GeneRates numgenes 
end