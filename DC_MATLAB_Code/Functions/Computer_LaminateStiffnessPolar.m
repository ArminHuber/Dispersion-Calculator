% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2024 DLR
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
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