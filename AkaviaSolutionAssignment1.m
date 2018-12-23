%% Define files
ccleExpressionFile = '~/Documents/Data/CCLE/expression/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct';
ccleSampleTypeFile = '~/Documents/Data/CCLE/expression/sampleType.txt';
%% Read data
ccleExpression = readtable(ccleExpressionFile, 'FileType', 'text', 'Delimiter', '\t', 'HeaderLines', 2, 'ReadRowNames', 1, 'ReadVariableNames', 1);
genes = ccleExpression.Properties.RowNames;
geneSymbols = ccleExpression.Description;
expressionMatrix = ccleExpression{:, 2:end};
ccleSamples = readtable(ccleSampleTypeFile);
%% Identify noise by identical samples

% First, log the data
expressionLogged = log2(expressionMatrix); % This will give too many -Inf, which is problematic
% Add one tenth the minimal value that is above 0
smallValue = 0.1 * min(expressionMatrix(expressionMatrix > 0));
expressionLogged = log2(expressionMatrix + smallValue);

identicalSamples = {'253JBV', '253J'; 'HCC1588', 'LS513'; 'DLD1', 'HCT15'};
% 253JBV and 253J were modified when loaded. You can identify them by
% looking at ccleSamples, and seeing which ones match. Or you can use
% various code examples
index253JBV = find(contains(ccleSamples.cellLine, '253JBV'));
index253J = find(contains(ccleSamples.cellLine, '253J'), 1, 'last'); % because contains will match 253J and 253JBV
% 1, last forces the last one, which is the sample we want.
% To demonstrate this, run 
ccleSamples(contains(ccleSamples.cellLine, '253J'), :)

plot(expressionLogged(:, index253JBV), expressionLogged(:, index253J), '+') % Noise seems between 0 and 5, let's say 2
M = expressionLogged(:, index253JBV) - expressionLogged(:, index253J);
A = mean([expressionLogged(:, index253JBV), expressionLogged(:, index253J)], 2);
plot(A, M, '+') % Noise seems to be a bit unclear - I would estimate at least 0 in this graph

%% Explanation of M/A and noise
% When looking at your assignments, a lot of you selected noise threshold
% on the M axis. This seems to indicate that I did not explain it well.
% M is the difference between samples. A is the average expression.
% At some point of average expression A, M becomes a lot bigger, and does
% look like noise.
hold on
line([-20 15], [2 2])
hold off


%% The other samples
for i=2:size(identicalSamples, 1)
    sample1 = expressionLogged(:, contains(ccleSamples.cellLine, identicalSamples{i, 1}));
    sample2 = expressionLogged(:, contains(ccleSamples.cellLine, identicalSamples{i, 2}));
    figure(); plot(sample1, sample2, '+'); title(['Scatter Plot ', identicalSamples{i, 1}, ' ', identicalSamples{i, 2}])
    M = sample1 - sample2;
    A = mean([sample1, sample2], 2);
    figure(); plot(A, M, '+'); title(['M/A plot ', identicalSamples{i, 1}, ' ', identicalSamples{i, 2}]);
end

%% Calculating noise against all other samples
samplesFlattened = identicalSamples(:);
samplesFlattened([1, 4]) = strcat('x', samplesFlattened([1, 4]));
[samplesFlattened, indSamples, indCCLE] = intersect(samplesFlattened, ccleSamples.cellLine);
for i=1:length(samplesFlattened)
    currentSampleExp = expressionLogged(:, indCCLE(i));
    otherSamples = setdiff(1:length(ccleSamples.cellLine), indCCLE(i));
    otherSamplesExp = mean(expressionLogged(:, otherSamples), 2); 
    figure(); plot(otherSamplesExp, currentSampleExp, '+'); title(['Scatter Plot ', samplesFlattened{i}, ' vs average']);
    M = currentSampleExp - otherSamplesExp;
    A = mean([currentSampleExp, otherSamplesExp], 2);
    figure(); plot(A, M, '+'); title(['M/A plot ', samplesFlattened{i}, ' vs average']);
end
%% Pick noise threshold
% In Assignment 1, if you did not filter by noise, you did not lose points
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

image(filteredBreastProstateExpression([find(upInBreastInd); find(downInBreastInd)], :), 'CDataMapping', 'scaled')
upInBreastGenes = filteredBreastProstateGenes(upInBreastInd);
downInBreastGenes = filteredBreastProstateGenes(downInBreastInd);