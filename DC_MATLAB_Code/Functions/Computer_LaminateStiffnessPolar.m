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
function CLP = Computer_LaminateStiffnessPolar(Material,LayerOrientations,LayerThicknesses,PropagationAngleResolution)
%#ok<*AGROW>
PropagationAngle = 0:PropagationAngleResolution:360;
UnitCellThickness = sum(LayerThicknesses);
for i = 1:length(PropagationAngle)
    Phi(i,:) = PropagationAngle(i)-LayerOrientations;
    for j = 1:length(LayerOrientations)
        s = sind(Phi(i,j));
        g = cosd(Phi(i,j));
        c{i,j}(1,1) = Material{j}.C(1,1)*g^4+Material{j}.C(2,2)*s^4+2*(Material{j}.C(1,2)+2*Material{j}.C(6,6))*s^2*g^2;
        c{i,j}(1,2) = (Material{j}.C(1,1)+Material{j}.C(2,2)-2*Material{j}.C(1,2)-4*Material{j}.C(6,6))*s^2*g^2+Material{j}.C(1,2);
        c{i,j}(1,3) = Material{j}.C(1,3)*g^2+Material{j}.C(2,3)*s^2;
        c{i,j}(1,6) = (Material{j}.C(1,2)+2*Material{j}.C(6,6)-Material{j}.C(1,1))*s*g^3+(Material{j}.C(2,2)-Material{j}.C(1,2)-2*Material{j}.C(6,6))*g*s^3;
        c{i,j}(2,2) = Material{j}.C(1,1)*s^4+Material{j}.C(2,2)*g^4+2*(Material{j}.C(1,2)+2*Material{j}.C(6,6))*s^2*g^2;
        c{i,j}(2,3) = Material{j}.C(2,3)*g^2+Material{j}.C(1,3)*s^2;  
        c{i,j}(2,6) = (Material{j}.C(1,2)+2*Material{j}.C(6,6)-Material{j}.C(1,1))*g*s^3+(Material{j}.C(2,2)-Material{j}.C(1,2)-2*Material{j}.C(6,6))*s*g^3;
        c{i,j}(3,3) = Material{j}.C(3,3);
        c{i,j}(3,6) = (Material{j}.C(2,3)-Material{j}.C(1,3))*s*g;
        c{i,j}(4,4) = Material{j}.C(4,4)*g^2+Material{j}.C(5,5)*s^2;
        c{i,j}(4,5) = (Material{j}.C(4,4)-Material{j}.C(5,5))*s*g;
        c{i,j}(5,5) = Material{j}.C(5,5)*g^2+Material{j}.C(4,4)*s^2;
        c{i,j}(6,6) = Material{j}.C(6,6)+(Material{j}.C(1,1)+Material{j}.C(2,2)-2*Material{j}.C(1,2)-4*Material{j}.C(6,6))*s^2*g^2;
    end
    CLaminatePolar{i,1} = c{i,1}*LayerThicknesses(1)/UnitCellThickness; 
    for j = 2:length(LayerOrientations)
        CLaminatePolar{i,1} = CLaminatePolar{i,1}+c{i,j}*LayerThicknesses(j)/UnitCellThickness;
    end
    CLP(i,1) = CLaminatePolar{i}(1,1);
    CLP(i,2) = CLaminatePolar{i}(1,2);
    CLP(i,3) = CLaminatePolar{i}(1,3);
    CLP(i,4) = CLaminatePolar{i}(1,6);
    CLP(i,5) = CLaminatePolar{i}(2,2);
    CLP(i,6) = CLaminatePolar{i}(2,3);
    CLP(i,7) = CLaminatePolar{i}(2,6);
    CLP(i,8) = CLaminatePolar{i}(3,3);
    CLP(i,9) = CLaminatePolar{i}(3,6);
    CLP(i,10) = CLaminatePolar{i}(4,4);
    CLP(i,11) = CLaminatePolar{i}(4,5);
    CLP(i,12) = CLaminatePolar{i}(5,5);
    CLP(i,13) = CLaminatePolar{i}(6,6);
end
CLP(:,14) = 0:2*pi/360:2*pi;