function [outputArray] = mapIndicies(datasetOne,datasetTwo, ctrl)



forwardMap = calculateMinimizer(datasetOne, datasetTwo, ctrl);
backwardMap = calculateMinimizer(datasetTwo, datasetOne, ctrl);

[maxPossibleMatches, smallestDataset] = min([numel(datasetOne) numel(datasetTwo)]);
outputArray = [];

if smallestDataset == 2
    % If the smallest dataset is the second one, we need to permute the
    % variables to ensure that we do not try to compare against a VOID.
    tmp = backwardMap;
    backwardMap = forwardMap;
    forwardMap = tmp;
end
for aLoop = 1:maxPossibleMatches
    
   forwardPredict   = forwardMap(aLoop,2);
   idxInSetTwo      = backwardMap(:,1) == forwardPredict;  
   backwardPredict  = backwardMap(idxInSetTwo,2);   
   
   if backwardPredict == forwardMap(aLoop,1)  
      outputArray = [outputArray ; forwardMap(aLoop,:)];
   end
    
end

if smallestDataset == 2
    outputArray = [outputArray(:,2) outputArray(:,1)];
    % Flip order
end

fprintf('Number of forward && backward matches: %d \n',size(outputArray,1))

