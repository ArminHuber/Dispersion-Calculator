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
function CLaminate = Computer_LaminateStiffness(Material,LayerOrientations,LayerThicknesses,PropagationAngle)
%#ok<*AGROW>
Phi = LayerOrientations-PropagationAngle;
for i = 1:length(Phi) % transformed layer stiffness matrices for orthotropic materials including transversely isotropy ones [1], p. 28
    s = sind(Phi(i));
    g = cosd(Phi(i));
    c{i}(1,1) = Material{i}.C(1,1)*g^4+Material{i}.C(2,2)*s^4+2*(Material{i}.C(1,2)+2*Material{i}.C(6,6))*s^2*g^2;
    c{i}(1,2) = (Material{i}.C(1,1)+Material{i}.C(2,2)-2*Material{i}.C(1,2)-4*Material{i}.C(6,6))*s^2*g^2+Material{i}.C(1,2);
    c{i}(1,3) = Material{i}.C(1,3)*g^2+Material{i}.C(2,3)*s^2;
    c{i}(1,6) = (Material{i}.C(1,2)+2*Material{i}.C(6,6)-Material{i}.C(1,1))*s*g^3+(Material{i}.C(2,2)-Material{i}.C(1,2)-2*Material{i}.C(6,6))*g*s^3;
    c{i}(2,2) = Material{i}.C(1,1)*s^4+Material{i}.C(2,2)*g^4+2*(Material{i}.C(1,2)+2*Material{i}.C(6,6))*s^2*g^2;
    c{i}(2,3) = Material{i}.C(2,3)*g^2+Material{i}.C(1,3)*s^2;  
    c{i}(2,6) = (Material{i}.C(1,2)+2*Material{i}.C(6,6)-Material{i}.C(1,1))*g*s^3+(Material{i}.C(2,2)-Material{i}.C(1,2)-2*Material{i}.C(6,6))*s*g^3;
    c{i}(3,3) = Material{i}.C(3,3);
    c{i}(3,6) = (Material{i}.C(2,3)-Material{i}.C(1,3))*s*g;
    c{i}(4,4) = Material{i}.C(4,4)*g^2+Material{i}.C(5,5)*s^2;
    c{i}(4,5) = (Material{i}.C(4,4)-Material{i}.C(5,5))*s*g;
    c{i}(5,5) = Material{i}.C(5,5)*g^2+Material{i}.C(4,4)*s^2;
    c{i}(6,6) = Material{i}.C(6,6)+(Material{i}.C(1,1)+Material{i}.C(2,2)-2*Material{i}.C(1,2)-4*Material{i}.C(6,6))*s^2*g^2;
end
UnitCellThickness = sum(LayerThicknesses);
CLaminate = c{1}*LayerThicknesses(1)/UnitCellThickness;
for i = 2:length(Phi) % multiplying the layer stiffness matrices with their repsective layer thickness-to-unit cell thickness-ratio and summing up
    CLaminate = CLaminate+c{i}*LayerThicknesses(i)/UnitCellThickness;
end
CLaminate(2,1) = CLaminate(1,2);
CLaminate(3,1) = CLaminate(1,3);
CLaminate(3,2) = CLaminate(2,3);
CLaminate(5,4) = CLaminate(4,5);
CLaminate(6,1) = CLaminate(1,6);
CLaminate(6,2) = CLaminate(2,6);
CLaminate(6,3) = CLaminate(3,6);