function [outputArray] = mapIndicies(datasetOne,datasetTwo, ctrl)



forwardMap = calculateMinimizer(datasetOne, datasetTwo, ctrl);
backwardMap = calculateMinimizer(datasetTwo, datasetOne, ctrl);


% Now check if maps coincide.
% Coincide in this context means that the best fit from forwardMap
% coincides with the best fit from backwardMap

[maxPossibleMatches, smallestDataset] = min([numel(datasetOne) numel(datasetTwo)]);
outputArray = [];
for aLoop = 1:maxPossibleMatches
    
   forwardPredict = forwardMap(aLoop,2);
   
   idxInSetTwo = backwardMap(:,1) == forwardPredict;
   
   backwardPredict = backwardMap(idxInSetTwo,2);   
    
    
   findIdx(aLoop) =  backwardPredict == forwardMap(aLoop,1);
   
   if findIdx(aLoop)
       
      outputArray = [outputArray ; forwardMap(aLoop,:)];
       
       
   end
    
end


fprintf('Number of forward && backward matches: %d \n',size(outputArray,1))

