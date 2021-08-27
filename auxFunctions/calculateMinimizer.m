function idxMap = calculateMinimizer(setOne, setTwo, ctrl)

% Initialize
distanceToFirstFiber = nan(numel(setTwo),1);
volDiffToFirstFiber  = nan(numel(setTwo),1);

idxDist = nan(numel(setOne),1);
idxVol = nan(numel(setOne),1);
idxComb = nan(numel(setOne),1);

saveSmallestDist = nan(numel(setOne),1);
saveSmallestDiff = nan(numel(setOne),1);
saveSmallestComb = nan(numel(setOne),1);

idxComb2A = nan(numel(setOne),1);
idxComb2B = nan(numel(setOne),1);

for aLoop = 1:numel(setOne)
    posOfFiber = [setOne(aLoop).posX setOne(aLoop).posY setOne(aLoop).posZ];
    volOfFiber = [setOne(aLoop).numFlags];
    % Extract data for legibility
    
	for bLoop = 1:numel(setTwo)
        posOfOtherFiber = [setTwo(bLoop).posX setTwo(bLoop).posY setTwo(bLoop).posZ];
        volOfOtherFiber = [setTwo(bLoop).numFlags];
        % Extract data for legibility
        
       
       distanceToFirstFiber(bLoop) = sqrt(sum((posOfFiber - posOfOtherFiber).^2));
       volDiffToFirstFiber(bLoop)  = abs(volOfFiber - volOfOtherFiber)./volOfFiber;
       % Perform calculations
    end
    
    [saveSmallestDist(aLoop), idxDist(aLoop)] = min(distanceToFirstFiber);
    [saveSmallestDiff(aLoop), idxVol(aLoop)]  = min(volDiffToFirstFiber./(mean(volDiffToFirstFiber)));
    [saveSmallestComb(aLoop), idxComb(aLoop)] = min(volDiffToFirstFiber .*distanceToFirstFiber.^2);
    % Save best match
    
    if ctrl.plotMode
        subplot(1,3,1)
        semilogy(sort(distanceToFirstFiber),'ow','MarkerFaceColor','k','MarkerSize',4)
        title(num2str(idxDist(aLoop)))
        subplot(1,3,2)
        semilogy(sort(volDiffToFirstFiber),'ow','MarkerFaceColor','k','MarkerSize',4)
        title(num2str(idxVol(aLoop)))
        subplot(1,3,3)
        semilogy(sort(volDiffToFirstFiber.*distanceToFirstFiber.^2),'ow','MarkerFaceColor','k','MarkerSize',4)
        title(num2str(idxComb(aLoop)))
    end
    
    idxComb2A(aLoop) = setOne(aLoop).idx;
    idxComb2B(aLoop) = setTwo(idxComb(aLoop)).idx;
end

idxMap = [idxComb2A idxComb2B];