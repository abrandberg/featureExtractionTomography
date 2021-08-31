%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mainFeatureExtractor.m
%
% Start here!
%
% created by : August Brandberg
% date : 2021-01-28
%

% Meta instructions
clear; close all; clc;
format;
format compact;
addpath('auxFunctions')

% Control struct
ctrl.plotMode               = 0;
% 0 - no plotting, 1 - plotting
ctrl.exportPlots            = 0;
% 0 - no saving of plots, 1 - saving of plots to saveDir
ctrl.colorArray             = lines(4);
ctrl.interpreter            = 'latex';
ctrl.histogramInstructions  = {'DisplayStyle','stairs','linewidth',1};
ctrl.formatSpecMsgL1        = '          -> %s\n';
ctrl.formatSpecMsgL2        = '                    - %s\n';
ctrl.saveDir                = ['plots' filesep];


% Hyper parameters
hyperParameters.discSteps               = 50;
% Number of midpoints to calculate along the length of the fiber.
hyperParameters.voxelSize               = 0.7; 
% Scaling factor from voxels to um.
hyperParameters.slendernessAcceptRatio  = 0.85;
% Accept criterion for fiber slenderness.
hyperParameters.volumeAccept            = 20000;
% Accept criterion for fiber size in number of voxels.

assert(hyperParameters.discSteps > 1.0-eps,'discSteps parameter must be an integer \geq 1.')
assert(hyperParameters.voxelSize > 0.0-eps,'voxelSize parameter must not be negative.')
assert(hyperParameters.slendernessAcceptRatio > 0.0 && hyperParameters.slendernessAcceptRatio < 1.0,'slendernessAcceptRatio should be \in (0.0 , 1.0)')
assert(hyperParameters.volumeAccept > 0.0-eps,'volumeAccept parameter must not be negative.')

fprintf(ctrl.formatSpecMsgL1,'Start of mainFeatureExtractor.m');
fprintf(ctrl.formatSpecMsgL2,['plotMode is ' num2str(ctrl.plotMode)]);
fprintf(ctrl.formatSpecMsgL2,['voxelSize is ' num2str(hyperParameters.voxelSize)]);


segmentedInputFieldFile = {'data\Sample_6_Third_Revision.nii'; % Sample_4.nii
                           'data\Sample_9.nii'}; %'data\Sample_6_Third_Revision.nii'
                       
% segmentedInputFieldFile = {'data\Sample_6_Third_Revision.nii'; % Sample_4.nii
%                            'data\Sample_4.nii'}; %'data\Sample_6_Third_Revision.nii'
%                        
% segmentedInputFieldFile = {'data\Sample_9.nii'; % Sample_4.nii
%                            'data\Sample_4.nii'}; %'data\Sample_6_Third_Revision.nii'

scanNames               =  {'Sample 4'; 'Sample 6'; 'Sample 9'};
% To be displayed in the legends of plots if ctrl.plotMode == true



% Open plots
if ctrl.plotMode
    if (exist(ctrl.saveDir,'dir') ~= 7)
        mkdir(ctrl.saveDir)
    end
    A = figure();
    B = figure('color','w','units','centimeters','OuterPosition',[10 10 2*16 16]);
    C = figure();
    D = figure();
end

for aLoop = 1:numel(segmentedInputFieldFile)
    
    fprintf(ctrl.formatSpecMsgL1,['Input file is ' segmentedInputFieldFile{aLoop}]);
    fprintf(ctrl.formatSpecMsgL2,'Importing file');
    segmentedField = importSegmentedData(segmentedInputFieldFile{aLoop});
    fprintf(ctrl.formatSpecMsgL2,['Field contains ' sprintf('%d',length(unique(segmentedField))) ' unique segments, inc. background']);
    

    unfilteredFiberPopulation = findSegmentedCenterlines(segmentedField,ctrl,hyperParameters);
    % Find center lines
    
    
    selIdx = [false [unfilteredFiberPopulation.numFlags] > hyperParameters.volumeAccept];
    fiberResult = unfilteredFiberPopulation(selIdx);
    % Filters on volume
    
    fiberResult = fiberResult([fiberResult.SvalOne] > hyperParameters.slendernessAcceptRatio);
    % Cleans out fibers which are not really "fiber-like" (i.e., not slender)
    
    if ctrl.plotMode
        figure(C)
        plot([unfilteredFiberPopulation.idx], [unfilteredFiberPopulation.numFlags], ...
             'o','color','w','MarkerFaceColor',ctrl.colorArray(aLoop,:),'displayname',scanNames{aLoop},'MarkerSize',5)
        hold on
        
        
        figure(A);
        plot3([fiberResult.posX],[fiberResult.posY],[fiberResult.posZ], ...
               'o','color','w','MarkerFaceColor',ctrl.colorArray(aLoop,:),'displayname',scanNames{aLoop},'MarkerSize',5)
        axis equal
        hold on
        pause(0.25)
    end
    
    
    fiberResult = findSegmentedCrossSections(fiberResult,ctrl,hyperParameters);
    % Extract cross-sectional properties
    
    for bLoop = 1:numel(fiberResult)
        fiberResult(bLoop).AMean = hyperParameters.voxelSize.^2 * median([fiberResult(bLoop).A(fiberResult(bLoop).A>0)]); % OBS OBS OBS
        fiberResult(bLoop).wMean = hyperParameters.voxelSize * mean([fiberResult(bLoop).w(fiberResult(bLoop).w>0)]);
        fiberResult(bLoop).hMean = hyperParameters.voxelSize * mean([fiberResult(bLoop).h(fiberResult(bLoop).h>0)]);
        fiberResult(bLoop).tMean = hyperParameters.voxelSize * mean([fiberResult(bLoop).t(fiberResult(bLoop).t>-eps)]);
    end
    
    if ctrl.plotMode
        
        figure(D)
        histogram([fiberResult.AMean],                             ...
                  'normalization','probability','edgecolor',ctrl.colorArray(aLoop,:), ...
                  ctrl.histogramInstructions{:})
        hold on
        xline(mean([fiberResult.AMean]),'--','linewidth',1.5,'color',ctrl.colorArray(aLoop,:))
        xline(median([fiberResult.AMean]),'-.','linewidth',1.5,'color',ctrl.colorArray(aLoop,:))
        
        
        figure(B)
        subplot(2,3,1)
        histogram([fiberResult.wMean],linspace(0,100,25),                             ...
                  'normalization','probability','edgecolor',ctrl.colorArray(aLoop,:), ...
                  ctrl.histogramInstructions{:})
        subplot(2,3,2) 
        histogram([fiberResult.hMean],linspace(0,50,15),                              ...
                  'normalization','probability','edgecolor',ctrl.colorArray(aLoop,:), ...
                  ctrl.histogramInstructions{:})
        subplot(2,3,3)
        histogram([fiberResult.tMean],linspace(0,20,15),                              ...
                  'normalization','probability','edgecolor',ctrl.colorArray(aLoop,:), ...
                  ctrl.histogramInstructions{:})
        subplot(2,3,4)
        histogram([fiberResult.wMean],linspace(0,100,25),'normalization','cdf','edgecolor',ctrl.colorArray(aLoop,:),ctrl.histogramInstructions{:})
        subplot(2,3,5)
        histogram([fiberResult.hMean],linspace(0,50,15),'normalization','cdf','edgecolor',ctrl.colorArray(aLoop,:),ctrl.histogramInstructions{:})
        subplot(2,3,6)
        histogram([fiberResult.tMean],linspace(0,20,15),'normalization','cdf','edgecolor',ctrl.colorArray(aLoop,:),ctrl.histogramInstructions{:})
        if aLoop == 1
            for vLoop = 1:6
                subplot(2,3,vLoop)
                hold on
            end
        end
    end

    datasetSave(aLoop).data = fiberResult;
    % Save results for this iteration of aLoop to a bigger struct for
    % comparison between files.
end

if ctrl.exportPlots
    figure(A)
    xlabel('$x$ [voxel]','interpreter',ctrl.interpreter)
    ylabel('$y$ [voxel]','interpreter',ctrl.interpreter)
    zlabel('$z$ [voxel]','interpreter',ctrl.interpreter)
    legend('location','northeastoutside','interpreter',ctrl.interpreter)
    set(gca,'TickLabelInterpreter',ctrl.interpreter)
    axis equal
    print([ctrl.saveDir filesep 'centerOfGravityPerFiberComparison'],'-dpng','-r800')

    figure(B)
    print([ctrl.saveDir filesep 'histogramsComparison'],'-dpng','-r800')

    figure(C)
    print([ctrl.saveDir filesep 'fiberVolumeComparison'],'-dpng','-r800')
end
% Second part: Mapping fibers in each file to each other.
% Note that this part currently only maps two different sets to eachother
% in a one-way comparison.
%
% TO DO:
%       - Arbitrary number of substeps
idxToCompareOne = 1;
idxToCompareTwo = 2;


[idxMapping] = mapIndicies(datasetSave(idxToCompareOne).data,datasetSave(idxToCompareTwo).data, ctrl);

segmentedFieldOne = importSegmentedData(segmentedInputFieldFile{idxToCompareOne});
segmentedFieldTwo = importSegmentedData(segmentedInputFieldFile{idxToCompareTwo});
dimOfInputOne = size(segmentedFieldOne);
dimOfInputTwo = size(segmentedFieldTwo);


if ctrl.plotMode
    figure;
end
for sLoop = 1:size(idxMapping,1)
    
    segOne = find(segmentedFieldOne == idxMapping(sLoop,1));
    segTwo = find(segmentedFieldTwo == idxMapping(sLoop,2));
    
    [I1_1,I2_1,I3_1] = ind2sub(dimOfInputOne,segOne);
    [I1_2,I2_2,I3_2] = ind2sub(dimOfInputTwo,segTwo);

    if ctrl.plotMode
        subplot(1,2,1)
        plot3(I1_1,I2_1,I3_1,'o','MarkerSize',2.5,'MarkerEdgeColor','none','MarkerFaceColor',ctrl.colorArray(idxToCompareOne,:))
        axis equal
        xlabel('$x$ [voxel]','interpreter',ctrl.interpreter)
        ylabel('$y$ [voxel]','interpreter',ctrl.interpreter)
        zlabel('$z$ [voxel]','interpreter',ctrl.interpreter)
        set(gca,'TickLabelInterpreter',ctrl.interpreter)
        
        subplot(1,2,2)
        plot3(I1_2,I2_2,I3_2,'o','MarkerSize',2.5,'MarkerEdgeColor','none','MarkerFaceColor',ctrl.colorArray(idxToCompareTwo,:))
        axis equal
    
        xlabel('$x$ [voxel]','interpreter',ctrl.interpreter)
        ylabel('$y$ [voxel]','interpreter',ctrl.interpreter)
        zlabel('$z$ [voxel]','interpreter',ctrl.interpreter)
        set(gca,'TickLabelInterpreter',ctrl.interpreter)
        pause(0.25)
        print([ctrl.saveDir filesep 'fiberMapQuality_fiber_'  num2str(sLoop)],'-dpng','-r800')
        hold off
    end
    
    saveVol(sLoop,1) = numel(segOne);
    saveVol(sLoop,2) = numel(segTwo);
        
end

if ctrl.plotMode
    figure;
    plot(saveVol(:,1),saveVol(:,2),'ow','MarkerFaceColor','k')
    hold on
    plot([0 max(saveVol(:,1))],[0 max(saveVol(:,1))],'k--')
    xlabel([scanNames{idxToCompareOne} ' Voxels per fiber'],'interpreter',ctrl.interpreter)
    ylabel([scanNames{idxToCompareTwo} ' Voxels per fiber'],'interpreter',ctrl.interpreter)
    set(gca,'TickLabelInterpreter',ctrl.interpreter)
    if ctrl.exportPlots
        print('systematicDiffVolumeQuestionMark','-dpng','-r1200')
    end
end

a = saveVol(:,1); b = saveVol(:,2);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AREA INVESTIGATION
ctrl.plotMode = 1;
if ctrl.plotMode
    figure;
end
for tLoop = 1:size(idxMapping,1)
    
    selAreaOne = [datasetSave(idxToCompareOne).data.idx] == idxMapping(tLoop,1);
    selAreaTwo = [datasetSave(idxToCompareTwo).data.idx] == idxMapping(tLoop,2);
    
    saveArea(tLoop,1) = datasetSave(idxToCompareOne).data(selAreaOne).AMean;
    saveArea(tLoop,2) = datasetSave(idxToCompareTwo).data(selAreaTwo).AMean;
        
    
    saveWidth(tLoop,1) = datasetSave(idxToCompareOne).data(selAreaOne).wMean;
    saveWidth(tLoop,2) = datasetSave(idxToCompareTwo).data(selAreaTwo).wMean;
end

if ctrl.plotMode
    figure;
    plot(saveArea(:,1),saveArea(:,2),'ow','MarkerFaceColor','k')
    hold on
    plot([0 max(saveArea(:,1))],[0 max(saveArea(:,1))],'k--')
    xlabel([scanNames{idxToCompareOne} ' Mean area per fiber'],'interpreter',ctrl.interpreter)
    ylabel([scanNames{idxToCompareTwo} ' Mean area per fiber'],'interpreter',ctrl.interpreter)
    set(gca,'TickLabelInterpreter',ctrl.interpreter)
    if ctrl.exportPlots
        print('systematicDiffAreaQuestionMark','-dpng','-r1200')
    end
end

cp = saveArea(:,1);
dp = saveArea(:,2);

if ctrl.plotMode
    figure;
    plot(saveWidth(:,1),saveWidth(:,2),'ow','MarkerFaceColor','k')
    hold on
    plot([0 max(saveWidth(:,1))],[0 max(saveWidth(:,1))],'k--')
    xlabel([scanNames{idxToCompareOne} ' Mean w per fiber'],'interpreter',ctrl.interpreter)
    ylabel([scanNames{idxToCompareTwo} ' Mean w per fiber'],'interpreter',ctrl.interpreter)
    set(gca,'TickLabelInterpreter',ctrl.interpreter)
    if ctrl.exportPlots
        print('systematicDiffWidthQuestionMark','-dpng','-r1200')
    end
end

ep = saveWidth(:,1);
fp = saveWidth(:,2);
