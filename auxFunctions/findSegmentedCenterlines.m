function fiberResult = findSegmentedCenterlines(segmentedField,ctrl,hyperParameters)
%function fiberResult = findSegmentedCenterlines(segmentedField,ctrl,hyperParameters)
%returns a STRUCT of individual fiber observations given a segmentedField
%and some constants.
%
% INPUTS
%           segmentedField              {A x B x C} uint16 array
%               Each value in A (as well as B, C) corresponds to the voxel
%               at coordinate [a , b , c] relative to one of the corners of
%               the tomographed volume where 
%                   a = length(1:A), 
%                   b = length(1:B), etc.
%
%           ctrl                        {} struct
%               .plotMode {Bool}
%                   Plot intermediate results, 0/1
%               .colorArray {N x 3} double, where N = number of colors.
%                   Colors for consistent plotting, if desired
%               .formatSpecMsgL1 & .formatSpecMsgL1 {char} 
%                   Formating rules for verbose output, if desired
%           
%           hyperParameters             {} struct
%               .discSteps
%                   Number of slices to divide the fiber into during
%                   calculation of the center line.
%               .voxelSize
%                   Micrometers per voxel-side, such that the volume of the
%                   voxel in um^3 equals voxelSize^3.
%               .slendernessAcceptRatio
%                   Lower limit for the proportion of variation in the 
%                   point cloud which should be explained by eigenvector 1.
%
% OUTPUTS
%           fiberResult                 {} struct
%               .idx {1} double
%                   Index (in the segmentedField array)
%               .numFlags {1} double
%                   Number of voxels associated with this index.
%               .centerline {hyperParameters.discSteps x 3} double
%                   Calculated centerline of the fiber.
%               .SvalOne {1} double
%                   Weight associated with the first eigenvector.
%               .posX {1} double
%                   Mean position of fiber along global X-axis.
%               .posY {1} double
%                   Mean position of fiber along global Y-axis.
%               .posZ {1} double
%                   Mean position of fiber along global Z-axis.
%               .localCoordinates {fiberResult.numFlags x 3} double
%                   Position [ X, Y, Z] of each voxel associated with this
%                   index, in the rotated base of V1, V2, V3.
%
%
% TO DO:
%   - Optimize selIdx call if possible.
%
%
% created by : August Brandberg
% date       : 2021-08-25


uniqueFibers    = unique(segmentedField);
numelFibers     = length(uniqueFibers);
dimOfInput      = size(segmentedField);

flagsToFind     = histcounts(segmentedField,0:(numelFibers+eps));
% Determine how many entries to find in "selIdx" below. Note that
% histcounts is highly optimized and the added time is easily saved.

for tLoop = 2:numelFibers 
    % Skip number 1 which is assumed to be the background
    
    fiberToSel = uniqueFibers(tLoop);
    selIdx = find(segmentedField == fiberToSel,flagsToFind(tLoop));
    % Called in this way, selIdx returns the indicies that are non-zero.
    %
    % This version is tested but 95% of the compute time is here so
    % optimizing it is highly worthwhile
    

    [I1,I2,I3] = ind2sub(dimOfInput,selIdx);
    % Convert the voxel indicator to a set of three-dimensional coordinates
    % where each point (entry in I1, I2, I3) is the center of a voxel.
    % Essentially, this is the sparse representation of segmentedField,
    % with the added complication of having to keep which segment is
    % "selected" somewhere else.

    
    [~,score,~,~,Svals,centerOfMass] = pca([I1 I2 I3],'Algorithm','svd','Economy',true,'NumComponents',3);
    % Extract the three eigenvectors that
    % describe the most of the observed variance in the data cloud given by
    % [I1 , I2 , I3] using the SVD algorithm.
    
    % COEFF         = Principal components of each eigenvector.
    % SCORE         = Coordinates in the base given by COEFF.
    % LATENT        = Principal component variances, corresponding to the eigenvalues of the
    %                 covariance matrix
    % EXPLAINED     = Fraction of variance explained by the variable.
    % MU            = Center of fiber
    
    
    if ctrl.plotMode
        plot3(score(:,1),score(:,2),score(:,3),'k.','DisplayName','Raw data')
        hold on
        title(['S vals = [' num2str(Svals(1)) ' / ' num2str(Svals(2)) ' / ' num2str(Svals(3)) ' ]'], ...
              'interpreter',ctrl.interpreter)
        xlabel('Global $x$', 'interpreter',ctrl.interpreter)
        ylabel('Global $y$', 'interpreter',ctrl.interpreter)
        zlabel('Global $z$', 'interpreter',ctrl.interpreter)
        axis equal
        set(gca,'TickLabelInterpreter',ctrl.interpreter)
    end
    
    rangeV1         = [min(score(:,1)) max(score(:,1))];
    % Range of values on axis one (the one explaining most of the variance,
    % hypothesied to roughly coincide with the axis of the fiber.
    
    discAxis        = linspace(rangeV1(1),rangeV1(2),hyperParameters.discSteps);
    discStepLength  = diff(discAxis);
    % Discretize along this axis

    %%% Bookkeeping
    fiberResult(tLoop).idx          = fiberToSel;
    fiberResult(tLoop).numFlags     = numel(selIdx);
    fiberResult(tLoop).centerline   = nan(hyperParameters.discSteps,3);
    fiberResult(tLoop).SvalOne      = Svals(1);
    fiberResult(tLoop).posX         = centerOfMass(1);
    fiberResult(tLoop).posY         = centerOfMass(2);
    fiberResult(tLoop).posZ         = centerOfMass(3);

    for cLoop = 1:hyperParameters.discSteps   
        
        if cLoop == 1
            lowerDiff = 0;
            upperDiff = 0.5*discStepLength(cLoop);
        elseif cLoop == hyperParameters.discSteps
            lowerDiff = 0.5*discStepLength(cLoop-1);
            upperDiff = 0;
        else
            lowerDiff = 0.5*discStepLength(cLoop-1);
            upperDiff = 0.5*discStepLength(cLoop);
        end
    
        scoreSelect = ( score(:,1) > (discAxis(cLoop)-lowerDiff) )  & ( score(:,1) < (discAxis(cLoop)+upperDiff) );

        meanTemp    = mean(score(scoreSelect,:),1);
        
        %%% Bookkeeping
        fiberResult(tLoop).centerline(cLoop,:) = meanTemp;
        fiberResult(tLoop).localCoordinates    = score;
        
%         if ctrl.plotMode % Diagnostic plot for each step inside loop
%             plot3(meanTemp(1),meanTemp(2)+100,meanTemp(3),'ks')
%         end
    end
    
    %%% Bookkeeping
    fiberResult(tLoop).centerline(isnan(fiberResult(tLoop).centerline(:,1)),:) = [];
    % Find and delete parts that did not contain enough points
    
    if ctrl.plotMode
        plot3(fiberResult(tLoop).centerline(:,1),fiberResult(tLoop).centerline(:,2)+100,fiberResult(tLoop).centerline(:,3), ...
              '-sk','DisplayName','Centerline')
%         legend('location','best','interpreter',ctrl.interpreter)
        pause(0.5)
        hold off
    end
end

