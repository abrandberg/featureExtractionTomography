function segmentedField = importSegmentedData(segmentedInputFieldFile)
%function segmentedField = importSegmentedData(segmentedInputFieldFile)
%imports a .NII file with segmented tomography data. Optionally, a file
%path ending in .MAT may be passed. The assumption is then that the .MAT
%file contains a previously imported .NII file which has been assigned to
%the variable A2. This use can sometimes speed up the loading a bit.
%
% INPUTS    
%           segmentedInputFieldFile     {String} 
%               Path to the file to import.
%
% OUTPUTS
%           segmentedField              {A x B x C} uint16 array
%               Each value in A (as well as B, C) corresponds to the voxel
%               at coordinate [a , b , c] relative to one of the corners of
%               the tomographed volume where 
%                   a = length(1:A), 
%                   b = length(1:B), etc.
%               
%               The value of A(i,j,k) ( == B(i,j,k), == C(i,j,k)) is the
%               unique identifier for the volume in that voxel.
%
%               Example:
%               A(432, 123, 987) == 31
%               ==>
%               "The voxel at coordinate [432 , 123 , 987] belongs to entity
%               31."
%               
%
% TO DO:
%
% created by : August Brandberg
% date       : 2021-08-25
%

switch lower(segmentedInputFieldFile(end-2:end))   
    case 'mat'
        load(segmentedInputFieldFile,'A2');
        segmentedField = A2;
    case 'nii'
        segmentedField = niftiread(segmentedInputFieldFile);
    otherwise
        disp('File type not implemented. Consider adding a case in importSegmentedData.m')
        disp(stop) % Poor man's error handling.
end