% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2026 DLR
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
function CLP = Computer_LaminateStiffness(Mode,Material,LayerOrientations,LayerThicknesses,PropagationAngle)
if  Mode == 2
    PropagationAngle = 0:PropagationAngle:360;
end
c = complex(zeros(length(PropagationAngle),13,length(LayerOrientations)));
for i = 1:length(PropagationAngle)
    Phi = PropagationAngle(i)-LayerOrientations;
    for j = 1:length(Phi)
        s = sind(Phi(j));
        g = cosd(Phi(j));
        s2 = s^2;
        g2 = g^2;
        s3 = s^3;
        g3 = g^3;
        s4 = s2^2;
        g4 = g2^2;
        c(i,1,j) = Material{j}.C(1,1)*g4+Material{j}.C(2,2)*s4+2*(Material{j}.C(1,2)+2*Material{j}.C(6,6))*s2*g2; % C11
        c(i,2,j) = (Material{j}.C(1,1)+Material{j}.C(2,2)-2*Material{j}.C(1,2)-4*Material{j}.C(6,6))*s2*g2+Material{j}.C(1,2); % C12
        c(i,3,j) = Material{j}.C(1,3)*g2+Material{j}.C(2,3)*s2; % C13
        c(i,4,j) = (Material{j}.C(1,2)+2*Material{j}.C(6,6)-Material{j}.C(1,1))*s*g3+(Material{j}.C(2,2)-Material{j}.C(1,2)-2*Material{j}.C(6,6))*g*s3; % C16
        c(i,5,j) = Material{j}.C(1,1)*s4+Material{j}.C(2,2)*g4+2*(Material{j}.C(1,2)+2*Material{j}.C(6,6))*s2*g2; % C22
        c(i,6,j) = Material{j}.C(2,3)*g2+Material{j}.C(1,3)*s2; % C23
        c(i,7,j) = (Material{j}.C(1,2)+2*Material{j}.C(6,6)-Material{j}.C(1,1))*g*s3+(Material{j}.C(2,2)-Material{j}.C(1,2)-2*Material{j}.C(6,6))*s*g3; % C26
        c(i,8,j) = Material{j}.C(3,3); % C33
        c(i,9,j) = (Material{j}.C(2,3)-Material{j}.C(1,3))*s*g; % C36
        c(i,10,j) = Material{j}.C(4,4)*g2+Material{j}.C(5,5)*s2; % C44
        c(i,11,j) = (Material{j}.C(4,4)-Material{j}.C(5,5))*s*g; % C45
        c(i,12,j) = Material{j}.C(5,5)*g2+Material{j}.C(4,4)*s2; % C55
        c(i,13,j) = Material{j}.C(6,6)+(Material{j}.C(1,1)+Material{j}.C(2,2)-2*Material{j}.C(1,2)-4*Material{j}.C(6,6))*s2*g2; % C66
    end
end
R = reshape(LayerThicknesses/sum(LayerThicknesses),1,1,length(LayerOrientations)); % layer thickness-to-unit cell thickness-ratios
CLP = sum(c.*R,3); % multiply the layer stiffness matrices with their respective layer thickness-to-unit cell thickness-ratios and sum up
if  Mode == 1
    CLP = [CLP(1) CLP(2) CLP(3) 0 0 CLP(4);CLP(2) CLP(5) CLP(6) 0 0 CLP(7);CLP(3) CLP(6) CLP(8) 0 0 CLP(9);0 0 0 CLP(10) CLP(11) 0;0 0 0 CLP(11) CLP(12) 0;CLP(4) CLP(7) CLP(9) 0 0 CLP(13)];
elseif Mode == 2
    CLP(:,14) = 0:pi/180:2*pi;
end