function fiberResult = findSegmentedCrossSections(fiberResult,ctrl,hyperParameters)
%function fiberResult = findSegmentedCrossSections(fiberResult,ctrl,hyperParameters)
%calculates cross-sectional properties based on previous determination of
%the fiber centerline.
%
% INPUTS:
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
% OUTPUTS:
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
%               .w { x 1} double
%                   Fiber width
%               .h { x 1} double
%                   Fiber height
%               .t { x 1} double
%                   Fiber wall thickness
%
%
% TO DO:
%
%
% created by : August Brandberg
% date : 2021-08-25
% 

if ctrl.plotMode
    D = figure('color','w','units','centimeters','OuterPosition',[10 10 2*16 16]);
end

numelFibers = numel(fiberResult);


for tLoop = 1:numelFibers
    
    score = fiberResult(tLoop).localCoordinates;
    % Load the fiber coordinates
    
    rangeV1 = [min(score(:,1)) max(score(:,1))];
    discAxis = linspace(rangeV1(1),rangeV1(2),hyperParameters.discSteps);
    discStepLength = diff(discAxis);
    % Range of values on axis one (the one explaining most of the variance,
    % hypothesied to roughly coincide with the axis of the fiber.
    
    cVec = [fiberResult(tLoop).centerline(1,:) ; diff(fiberResult(tLoop).centerline,1,1)];
    cVecNorm = cVec ./ vecnorm(cVec,2,2);
    % Generate a normalized directional vector from each point in the centerline array
    % to the next point.
    
    if sum(cVecNorm(:,3)) == 0
        oVecNorm1 = [cVecNorm(:,[1 3]) -(cVecNorm(:,1).^2+cVecNorm(:,3).^2)./cVecNorm(:,2)];
    else
        oVecNorm1 = [cVecNorm(:,[1 2]) -(cVecNorm(:,1).^2+cVecNorm(:,2).^2)./cVecNorm(:,3)];
    end
    oVecNorm1 = oVecNorm1./vecnorm(oVecNorm1,2,2);
    oVecNorm2 = cross(cVecNorm,oVecNorm1);
    % Generate two orthogonal vectors (orthogonal to cVecNorm) that form a plane which is
    % locally orthogonal to the fiber axis (fiber axis defined by the center line.)

    for kLoop = 2:size(oVecNorm1)-1
        % Skip the ends of the fiber, which are by nature somewhat less well-defined than the bulk.

        tempRotMatrix = [cVecNorm(kLoop,:)' oVecNorm1(kLoop,:)' oVecNorm2(kLoop,:)'];
        % Define the rotation matrix to the new system from the coordinate system obtained from the 
        % first PCA (== the coordinate system with eigenvector 1 aligned with the global fiber axis).

        % Select the points to be rotated
        lowerDiff = 0.5*discStepLength(kLoop-1);
        upperDiff = 0.5*discStepLength(kLoop);
        if kLoop == 1
            lowerDiff = 0;
        elseif kLoop == hyperParameters.discSteps
            upperDiff = 0;
        end
        
        scoreSelect = ( score(:,1) > (discAxis(kLoop)-lowerDiff) )  & ( score(:,1) < (discAxis(kLoop)+upperDiff) );
        
        if sum(scoreSelect) > 3
            rotatedScore = score(scoreSelect,:)*tempRotMatrix' - mean(score(scoreSelect,:)*tempRotMatrix',1);
            % Rotate and normalize the coordinates.

            [~,score2,~,~,~,~] = pca(rotatedScore(:,[2 3]),'Algorithm','svd','Economy',true,'NumComponents',2);
            % Perform a second PCA analysis, extracting the two vectors in the plane containing most variation.
            
            if range(score2(:,1)) > 2 && range(score2(:,2)) > 2
                % If there is sufficient variation along each axis.

                [inVecX,inVecY] = prepareCrossSection(score2);
                %

                [w,h,t] = CSprop(inVecX,inVecY);
                % Fit the cross-section to a solid and hollow rectangle. Return whichever seems
                % better.

                fiberResult(tLoop).A(kLoop-1) = sum(inVecX ~= 0.0);
                fiberResult(tLoop).w(kLoop-1) = w;
                fiberResult(tLoop).h(kLoop-1) = h;
                fiberResult(tLoop).t(kLoop-1) = t;

                if ctrl.plotMode
                    
                    subplot(1,3,1)
                    plot3(score(:,1),score(:,2),score(:,3),'sk')
                    hold on
                    plot3(score(scoreSelect,1),score(scoreSelect,2),score(scoreSelect,3),'sr')
                    axis equal
%                     xlabel('Fiber-centric $x$','interpreter',ctrl.interpreter)
%                     ylabel('Fiber-centric $y$','interpreter',ctrl.interpreter)
%                     zlabel('Fiber-centric $z$','interpreter',ctrl.interpreter)
                    hold off
                    
                    subplot(2,3,2)
                    plot(rotatedScore(:,2),rotatedScore(:,3),'sr')
                    axis equal
                    xlabel('Fiber-centric $S2$','interpreter',ctrl.interpreter)
                    ylabel('Fiber-centric $S3$','interpreter',ctrl.interpreter)
                    set(gca,'TickLabelInterpreter',ctrl.interpreter)
                    hold off
                    
                    subplot(2,3,5)
                    plot(score2(:,1),score2(:,2),'sr')
                    xlabel('Fiber-segment $S2$','interpreter',ctrl.interpreter)
                    ylabel('Fiber-segment $S3$','interpreter',ctrl.interpreter)
                    set(gca,'TickLabelInterpreter',ctrl.interpreter)
                    hold on

                    plot(0.5.*[-w w w -w -w],0.5.*[-h -h h h -h],'-b','linewidth',2)        
                    plot(0.5.*[-w+2*t w-2*t w-2*t -w+2*t -w+2*t],0.5.*[-h+2*t -h+2*t h-2*t h-2*t -h+2*t],'-b','linewidth',2)        
                    axis equal
                    hold off
                    % hold on
%                     subplot(2,3,6)
%                     plot(score2(:,1),score2(:,2),'sr')
%                     hold on

    %                 plot(0.5.*[-w0-t0 w0+t0 w0+t0 -w0-t0 -w0-t0],0.5.*[-h0-t0 -h0-t0 h0+t0 h0+t0 -h0-t0],'-b','linewidth',2)        
    %                 plot(0.5.*[-w0+t0 w0-t0 w0-t0 -w0+t0 -w0+t0],0.5.*[-h0+t0 -h0+t0 h0-t0 h0-t0 -h0+t0],'-b','linewidth',2)        
                    % plot([-w w w -w -w],[-h -h h h -h],'-','linewidth',2)    
                    axis equal
                    hold off

                    
                end
            else
                fiberResult(tLoop).w(kLoop-1) = nan;
                fiberResult(tLoop).h(kLoop-1) = nan;
                fiberResult(tLoop).t(kLoop-1) = nan;
            end
        else
            fiberResult(tLoop).w(kLoop-1) = nan;
            fiberResult(tLoop).h(kLoop-1) = nan;
            fiberResult(tLoop).t(kLoop-1) = nan;
        end

        if ctrl.plotMode
            subplot(2,3,3)
            plot([fiberResult(tLoop).w],'bo-','DisplayName','Width')
            hold on
            plot([fiberResult(tLoop).h],'sr-','DisplayName','Height')
            plot([fiberResult(tLoop).t],'kd-','DisplayName','Wall thickness')
            xlabel('Segment idx','interpreter',ctrl.interpreter)
            ylabel('Property value [voxel]','interpreter',ctrl.interpreter)
            set(gca,'TickLabelInterpreter',ctrl.interpreter)
            legend('location','best','interpreter',ctrl.interpreter)
            hold off
            
            
            subplot(2,3,6)
             plot([fiberResult(tLoop).A],'bo-','DisplayName','Area')
            hold on
            xlabel('Segment idx','interpreter',ctrl.interpreter)
            ylabel('Area [voxel]$^2$','interpreter',ctrl.interpreter)
            set(gca,'TickLabelInterpreter',ctrl.interpreter)
            legend('location','best','interpreter',ctrl.interpreter)
            hold off
            
            
                    
            rangeW2 = 10*[-1 1];
            rangeW3 = 10*[-1 1];

            leftBottom  = [discAxis(kLoop) 100+fiberResult(tLoop).centerline(kLoop,2) fiberResult(tLoop).centerline(kLoop,3)] +[0 rangeW2(1) rangeW3(1)]*tempRotMatrix';
            rightBottom = [discAxis(kLoop) 100+fiberResult(tLoop).centerline(kLoop,2) fiberResult(tLoop).centerline(kLoop,3)] +[0 rangeW2(2) rangeW3(1)]*tempRotMatrix';
            rightTop    = [discAxis(kLoop) 100+fiberResult(tLoop).centerline(kLoop,2) fiberResult(tLoop).centerline(kLoop,3)] +[0 rangeW2(2) rangeW3(2)]*tempRotMatrix';
            leftTop     = [discAxis(kLoop) 100+fiberResult(tLoop).centerline(kLoop,2) fiberResult(tLoop).centerline(kLoop,3)] +[0 rangeW2(1) rangeW3(2)]*tempRotMatrix';
            subplot(1,3,1)
            patch('Faces',[1 2 3 4],'Vertices',[leftBottom ; rightBottom ; rightTop ; leftTop],'FaceColor','red');
            hold off
            pause(1)
            
            if ctrl.exportPlots
               print([ctrl.saveDir filesep 'crossSectionResult_' sprintf('%d',tLoop)],'-dpng','-r800') 
            end
        end

    end
end

