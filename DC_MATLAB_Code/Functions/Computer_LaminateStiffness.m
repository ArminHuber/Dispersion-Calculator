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