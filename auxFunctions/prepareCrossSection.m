function [inVecX,inVecY] = prepareCrossSection(score2)
%function [inVecX,inVecY] = prepareCrossSection(score2) works on
%the coordinates returned by the second PCA analysis, preparing the
%data for fitting of the cross-section.
%
% INPUTS:
%			score2 {N x 2} double where N is number of selected points.
%				Matrix of coordinates in the rotated coordinate system
%               which is local for the specific mid-point currently 
%               being examined.
%
% OUTPUTS:
%			inVecX {N} double
% 				Vector of values that are either a coordinate X
%               or 0. 0 means there was no mass at that position.
%
%			inVecY {N} double
% 				Vector of values that are either a coordinate Y
%               or 0. 0 means there was no mass at that position.
%
%
% TO DO:
%
%
% created by : August Brandberg
% date : 2021-08-26


edgesAlongOne = (-eps+min(score2(:,1))):1:(max(score2(:,1))+eps);
edgesAlongTwo = (-eps+min(score2(:,2))):1:(max(score2(:,2))+eps);
% Discretize along the two base vectors from smallest to biggest value seen
% on that axis.

[xa,Xedges,Yedges] = histcounts2(score2(:,1),score2(:,2),edgesAlongOne,edgesAlongTwo);
% Perform a 2D histogram count for each point spanned by the edges* vectors.
% Return number of hits per bin + bin edges.

XCenters = Xedges(2:end) - 0.5*diff(Xedges(1:2));
YCenters = Yedges(2:end) - 0.5*diff(Yedges(1:2));
% Obtain the bin centers by subtracting half a bin from each vector entry.

[X,Y] = meshgrid(XCenters,YCenters);
% Generate 2D mesh of the centers.

xa(xa > 0) = 1;
% Logical filtering: If more than one histcount hit (i.e., more than one voxel belonging to fiber in bin)
% then reset to 1 for "filled".
xa(xa < 1) = 0;
% Logical filtering: If less than one histcount hit (i.e, box empty) then set to 0.
% Of course, should anyway be zero.

inVecX = reshape(X.*xa',numel(X.*xa'),[]);
inVecY = reshape(Y.*xa',numel(Y.*xa'),[]);
% Multiply and reshape for output. Basically what this does is:
%
% 1. Multiply each coordinate in X with the binary indicator for whether that coordinate had fiber
%    in it.
% 2. Reshape the resulting matrix to a vector.
% 3. Do the same for Y coordinates, using the same binary indicator.