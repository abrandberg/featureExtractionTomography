function [w,h,t] = CSprop(x,y)
%function [w,h,t] = CSprop(x,y) takes coordinates in space as input and returns 
%the best fitting solid or hollow rectangle.
%
% The best fitting cross-section is assumed to be aligned such that the width of
% the cross-section is oriented along the direction of maximum variation, i.e.
% aligned with the first eigenvector.
%
% INPUTS:
%           x {N x 1} double
%              Vector of x-coordinates describing positions of mass.
%
%           y {N x 1} double
%              Vector of x-coordinates describing positions of mass.
%
% OUTPUTS:
%           w {1} double
%              Fiber cross-section width in voxels
%           
%           h {1} double
%              Fiber cross-section height in voxels
%
%           t {1} double
%              Fiber cross-section wall thickness, if cross-section is hollow, 
%              in voxels. Otherwise 0.
%
% TO DO:
%
% created by : August Brandberg
% date : 2021-08-26
%
A   = sum(x~=0);
Ixx = sum(x.^2);
Iyy = sum(y.^2);

opts    = optimset('Display', 'off','TolX',1e-12,'TolFun',1e-16,'Algorithm','levenberg-marquardt');
options = optimset('Display', 'off','TolX',1e-12,'TolFun',1e-16);

ini = [30; 10 ; 0+eps];
% Initial guess for width, height and wall thickness.


costFcnSolid  = @(x) [(x(1)*x(2)-A)./A            ; ...
                      (1/12*x(1)*x(2)^3-Iyy)./Iyy ; ...
                      (1/12*x(1)^3*x(2)-Ixx)./Ixx  ];

costFcnHollow = @(x) sqrt(1/length(ini) .* sum( ...
                     ( [(x(1)*x(2)-(x(2)-2*x(3))*(x(1)-2*x(3))-A)./A                         ;  ...
                        (1/12*x(1)*x(2)^3-(1/12)*(x(1)-2*x(3))*(x(2)-2*x(3))^3-Iyy)./Iyy     ;  ...
                        (1/12*x(1)^3*x(2)-(1/12)*(x(1)-2*x(3))^3*(x(2)-2*x(3))-Ixx)./Ixx      ] ...
                        ).^2) );

[Y,fvalSolid] = fsolve(costFcnSolid, ini(1:2), opts);
% Solve under the assumption that the cross-section is a solid rectangle

[xOut,fvalHollow] = fmincon(costFcnHollow,ini,[],[],[],[],[0 ; 0 ; 0],[],@mycon,options);
% Solve under the assumption that the cross-section is a hollow rectangle


errSolid = sqrt(1/length(fvalSolid).*sum(fvalSolid.^2));
% Create comparable error estimates

if errSolid < fvalHollow
   Y(3) = 0;
else
   Y = xOut; 
end

w = Y(1);                                % Fiber Width
h = Y(2);                                % Fiber Height
t = Y(3);                                % Fiber Thickness
% Assign outputs