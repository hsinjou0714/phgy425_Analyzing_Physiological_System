load Assignment3.mat

%% What kind of function should you use
% using the numbers from the hygepdf example
M=100;
K=20;
N=10;
X=5;

hygepdf(X, M, K, N) % Chance of picking EXACTLY 5
hygecdf(X, M, K, N) % Chance of picking 0 to 5, also equivalent to
sum(hygepdf(0:X, M, K, N))

% We need probability of picking 5 to 10 (we can't pick more than the group
% size). There are several options
sum(hygepdf(X:N, M, K, N)) % sum of probability of picking 5, 6, 7, 8, 9, 10

% Probabilities sum up to 1, so the probability desired is 
% 1 - probability of getting less than X
% 1 - sum of probability of picking 0, 1, 2, 3, 4
% 1 - hygecdf(X-1, M, K, N)
1 - hygecdf(X-1, M, K, N)

%% Differential expression between breast and prostate
noiseThreshold = -2; % To keep most of the data. 0, 2, and even 5 might be acceptable
numOfVeryLowSamplesPerGene = sum(expressionLogged < noiseThreshold, 2);
filteredExpression = expressionLogged(numOfVeryLowSamplesPerGene <200, :);

expressionBreast = expressionLogged(:, strcmp(ccleSamples.tissueType, 'BREAST'));
expressionProstate =  expressionLogged(:, strcmp(ccleSamples.tissueType, 'PROSTATE'));
filteredBreastProstateExpression = [expressionBreast,  expressionProstate];
numberOfBreastSamples = sum(strcmp(ccleSamples.tissueType, 'BREAST'));
numberOfProstateSamples = sum(strcmp(ccleSamples.tissueType, 'PROSTATE'));
numVeryLowBrPr = sum(filteredBreastProstateExpression < noiseThreshold, 2);
hist(numVeryLowBrPr) % How many samples is each gene NOT expressed in
% Genes that are all 0 will have this number close to total number of
% samples (61)
numberOfSamplesThatAreVeryLow = (numberOfBreastSamples + numberOfProstateSamples)/2;
genesToKeep = numVeryLowBrPr < numberOfSamplesThatAreVeryLow; 
% Genes that are expressed in at least 30 samples will be kept
filteredBreastProstateExpression = filteredBreastProstateExpression(genesToKeep, :);
[~, p] = ttest2(filteredBreastProstateExpression(:, 1:numberOfBreastSamples), filteredBreastProstateExpression(:, (numberOfBreastSamples+1):end), 'Dim', 2, 'Vartype', 'unequal');
averageDifference = mean(filteredBreastProstateExpression(:, 1:numberOfBreastSamples), 2) - mean(filteredBreastProstateExpression(:, (numberOfBreastSamples+1):end), 2);

filteredBreastProstateGenes = genes(genesToKeep);
filteredBreastProstateGeneSymbols = geneSymbols(genesToKeep);

upInBreastInd = p < 0.05 & averageDifference > 0;
downInBreastInd = p < 0.05 & averageDifference < 0;

upInBreastGenes = filteredBreastProstateGenes(upInBreastInd);
downInBreastGenes = filteredBreastProstateGenes(downInBreastInd);
%% Set up group matrix
groupMatrixBreastProstate = [upInBreastInd, downInBreastInd];
for i=1:length(GOLoaded)
    [GOResults(i).p, GOResults(i).commonGenes] = hyperGeometricEnrichment(GOLoaded(i).Genes, GOLoaded(i).Terms, GOLoaded(i).Matrix, filteredBreastProstateGenes, groupMatrixBreastProstate);
end

%% Remove terms with less than 5 or more than 500 genes
GOFiltered = GOLoaded;
for i=1:length(GOFiltered)
    toKeep = sum(GOFiltered(i).Matrix) > 5 & sum(GOFiltered(i).Matrix) < 500;
    GOFiltered(i).Matrix = GOFiltered(i).Matrix(:, toKeep);
    GOFiltered(i).Terms = GOFiltered(i).Terms(toKeep);
end

for i=1:length(GOFiltered)
    [GOResultsFiltered(i).p, GOResultsFiltered(i).commonGenes] = hyperGeometricEnrichment(GOFiltered(i).Genes, GOFiltered(i).Terms, GOFiltered(i).Matrix, filteredBreastProstateGenes, groupMatrixBreastProstate);
end

% What are the Biological Processes enriched in genes that are up in breast
GOFiltered(1).Terms(GOResultsFiltered(1).p(:, 1) < 0.05)

%% How to read GOTerm2ID and use it