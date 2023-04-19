% =========================================================================
% Dispersion Calculator
% Copyright (C) 2018-2023 DLR
% Created by Armin Huber
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% For more information contact Armin Huber at armin.huber@dlr.de
% =========================================================================
function [BulkWaves,X,Y] = Computer_Isotropic_SnellsLaw(Fluid,Solid,Theta,OutputWindowUI8)
X(:,1) = pi:pi/180:2*pi; % Theta for solid (rad)
X(:,2) = 1e3/Solid.LongitudinalVelocity; % slowness of longitudinal mode in the solid (ms/m)
X(:,3) = 1e3/Solid.TransverseVelocity; % slowness of transverse mode in the solid (ms/m)
Y(:,1) = 0:pi/180:pi; % Theta for fluid (rad)
Y(:,2) = 1e3/Fluid.Velocity; % slowness in the fluid (ms/m)
BulkWaves(4,1) = 0;
for i = 1:2
    if  real(asind(sind(Theta)*1e3/Fluid.Velocity/X(1,i+1))) ~= 90
        BulkWaves(1,i) = asind(sind(Theta)*1e3/Fluid.Velocity/X(1,i+1));
        BulkWaves(2,i) = X(1,i+1);
        BulkWaves(3,i) = sind(BulkWaves(1,i));
        BulkWaves(4,i) = cosd(BulkWaves(1,i));
    else
        BulkWaves(1,i) = 90;
        BulkWaves(2,i) = Y(1,2)*sind(Theta);
        BulkWaves(3,i) = 1;
        BulkWaves(4,i) = 0;
    end
end
% content of "BulkWaves":
% columns 1: L 2: S
% lines: 
% 1: polar angle
% 2: slowness
% 3: propagation direction component n1
% 4: propagation direction component n3
Theta2 = Theta;
if  Theta2 == 0
    Theta2 = 0.001;
elseif Theta2 == 90
    Theta2 = 89.999;
end
PhaseVelocity = Fluid.Velocity/sind(Theta2); % phase velocity component along x1 shared among all participating bulk waves, as required by Snell's law
Alpha(1) = sqrt(PhaseVelocity^2/Solid.LongitudinalVelocity^2-1);
Alpha(2) = sqrt(PhaseVelocity^2/Solid.TransverseVelocity^2-1);   
AlphaFluid = 1/Fluid.Velocity*PhaseVelocity*cosd(Theta2);
W = (PhaseVelocity^2-Solid.LongitudinalVelocity^2-Solid.TransverseVelocity^2*Alpha.^2)./((Solid.LongitudinalVelocity^2-Solid.TransverseVelocity^2)*Alpha);
D = 1i/PhaseVelocity*[0 0;0 0;Solid.Lambda+(Solid.Lambda+2*Solid.Mu)*Alpha.*W;0 0;Solid.Mu*(Alpha+W)]; % sigma33 sigma13
DFluid = 1i*Fluid.Density*PhaseVelocity; % sigmaFluid33
Z1 = [W AlphaFluid;D(3,:) -DFluid;D(5,:) 0];
Z2 = [AlphaFluid;DFluid;0];
U = Z1\Z2; % UT(L) UT(S) UR(L)
v = -1i*[1 AlphaFluid;U(1) W(1)*U(1);U(2) W(2)*U(2);U(3) AlphaFluid*U(3)]; % I(L) T(L) T(S) R(L)
sigma = [DFluid 0;D(3,1)*U(1) D(5,1)*U(1);D(3,2)*U(2) D(5,2)*U(2);DFluid*U(3) 0]; % I(L) T(L) T(S) R(L)
P = -.5*[real(sigma(1,1)*conj(v(1,2)));real(sigma(2,1)*conj(v(2,2))+sigma(2,2)*conj(v(2,1)));real(sigma(3,1)*conj(v(3,2))+sigma(3,2)*conj(v(3,1)));real(sigma(4,1)*conj(v(4,2)))]; % I(L) T(L) T(S) R(L)
E(1:3) = P(2:4)/P(1); % T(L) T(S) R(L) Sum
if  E(3) > 1.001
    E(:) = [0 0 1 0];
end
E(4) = sum(E,2);
if  abs(E(1)) < 1e-7
    E(1) = 0;
end
if  abs(E(2)) < 1e-7
    E(2) = 0;
end
if  abs(E(3)) < 1e-7
    E(3) = 0;
end
String = ['Fluid: ',Fluid.Name,newline,'Solid: ',Solid.Name,newline,newline,'Critical angles:'];
if  isreal(asind(Fluid.Velocity/Solid.LongitudinalVelocity))
    String = append(String,newline,'  Theta_L = ',num2str(asind(Fluid.Velocity/Solid.LongitudinalVelocity)),' deg');
else
    String = append(String,newline,'  Theta_L = -');
end
if  isreal(asind(Fluid.Velocity/Solid.TransverseVelocity))
    String = append(String,newline,'  Theta_SV = ',num2str(asind(Fluid.Velocity/Solid.TransverseVelocity)),' deg');
else
    String = append(String,newline,'  Theta_SV = -');
end
String = append(String,newline,newline,...
'Plane wave incidence:',newline,...
'  Theta = ',num2str(Theta),' deg',newline,...
newline,...
'LONGITUDINAL WAVE (L)',newline,...
'Refraction angle = ',num2str(BulkWaves(1,1)),' deg',newline,...
'Phase velocity = ',num2str(1/BulkWaves(2,1)),' m/ms',newline,...
'Phase velocity x1 = ',num2str(1/BulkWaves(2,1)*BulkWaves(3,1)),' m/ms',newline,...
'Phase velocity x3 = ',num2str(1/BulkWaves(2,1)*BulkWaves(4,1)),' m/ms');
if  BulkWaves(1,1) ~= 90
    String = append(String,newline,'Polarization angle = 0 deg');
else
    String = append(String,newline,'Evanescent surface wave');
end
String = append(String,newline,newline,...
'SHEAR VERTICAL WAVE (SV)',newline,...
'Refraction angle = ',num2str(BulkWaves(1,2)),' deg',newline,...
'Phase velocity = ',num2str(1/BulkWaves(2,2)),' m/ms',newline,...
'Phase velocity x1 = ',num2str(1/BulkWaves(2,2)*BulkWaves(3,2)),' m/ms',newline,...
'Phase velocity x3 = ',num2str(1/BulkWaves(2,2)*BulkWaves(4,2)),' m/ms');
if  BulkWaves(1,2) ~= 90
    String = append(String,newline,'Polarization angle = 90 deg');
else
    String = append(String,newline,'Evanescent surface wave');
end
String = append(String,newline,newline,...
'Energy scattering coefficients:',newline,...
'  R(L) = ',num2str(E(3)),newline,...
'  T(L) = ',num2str(E(1)),newline,...
'  T(SV) = ',num2str(E(2)),newline,...
'  Sum = ',num2str(E(4)));
OutputWindowUI8.String = String;
disp([String,newline,'---------------------------------'])