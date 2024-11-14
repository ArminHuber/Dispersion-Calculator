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
function Materials = MaterialList
% fluids
fileID = fopen('MaterialList_Fluid.txt'); % open txt file containing the fluids' names, densities, and phase velocities
A = textscan(fileID,'%s %f %f'); % read the txt file
fclose(fileID); % close the txt file
Materials.Fluid = struct; % generate a structure where to store the materials' parameters
for i = 1:size(A{1},1) % enter all materials contained in A
    Materials.Fluid = setfield(Materials.Fluid,{1,1},cell2mat(A{1}(i)),{1,1},'Name',cell2mat(A{1}(i)));
    Materials.Fluid = setfield(Materials.Fluid,{1,1},cell2mat(A{1}(i)),{1,1},'Class','Fluid');
    Materials.Fluid = setfield(Materials.Fluid,{1,1},cell2mat(A{1}(i)),{1,1},'Density',A{2}(i));
    Materials.Fluid = setfield(Materials.Fluid,{1,1},cell2mat(A{1}(i)),{1,1},'Velocity',A{3}(i));
end
Materials.Fluid = orderfields(Materials.Fluid); % sort material entries into alphabetical order

% isotropic materials
fileID = fopen('MaterialList_Isotropic.txt'); % open txt file containing isotropic materials' names, densities, and engineering constants, i.e., Young's moduli E, and Poisson's ratios v, and longitudinal and transverse attenuation
A = textscan(fileID,'%s %f %f %f %f %f'); % read the txt file
fclose(fileID); % close the txt file
A{7} = A{3}.*A{4}./((1+A{4}).*(1-2*A{4})); % Lame's first parameter (Pa)
A{8} = A{3}./(2*(1+A{4})); % Lame's second parameter (Pa)
A{9} = sqrt((A{7}+2*A{8})./A{2}); % longitudinal bulk wave phase velocity (m/s)
A{10} = sqrt(A{8}./A{2}); % transverse bulk wave phase velocity (m/s)
A{11} = A{9}./(1+1i*A{5}/(2*pi)); % complex longitudinal bulk wave velocity; Disperse manual, p. 183, Eq. (A.17)
A{12} = A{10}./(1+1i*A{6}/(2*pi)); % complex transverse bulk wave velocity;
A{13} = sqrt(A{3}./(A{2}.*(1-A{4}.^2))); % plate wave velocity (m/s)
A{14} = sqrt(A{3}./A{2}); % cylinder wave velocity (m/s)
A{15} = A{10}.*((.87+1.12*A{4})./(1+A{4})); % Rayleigh wave velocity (m/s)
Mu_real = A{2}.*A{10}.^2*4*pi^2.*(4*pi^2-A{6}.^2)./(A{6}.^2+4*pi^2).^2; % complex Lame's parameters; Disperse manual, p. 183, Eq. (A.19)
Mu_imag = A{2}.*A{10}.^2*4*pi^2.*(4*pi*A{6})./(A{6}.^2+4*pi^2).^2;
Lambda_real = A{2}.*A{9}.^2*4*pi^2.*(4*pi^2-A{5}.^2)./(A{5}.^2+4*pi^2).^2-2*Mu_real;
Lambda_imag = A{2}.*A{9}.^2*4*pi^2.*(4*pi*A{5})./(A{5}.^2+4*pi^2).^2-2*Mu_imag;
A{16} = Lambda_real+1i*Lambda_imag; % complex Lame's first parameter, Lambda (Pa)
A{17} = Mu_real+1i*Mu_imag; % complex Lame's second parameter, Mu (Pa)
A{18} = A{17}.*((3*A{16}+2*A{17})./(A{16}+A{17})); % complex Young's modulus (GPa)
A{19} = A{16}./(2*(A{16}+A{17})); % complex Poisson's ratio
A{20}{1,1} = A{16}+2*A{17}; % complex stiffness component C11 (Pa) 
A{20}{1,2} = A{16}; % C12 (Pa)
A{20}{1,3} = A{16}; % C13 (Pa)
A{20}{2,2} = A{16}+2*A{17}; % C22 (Pa)
A{20}{2,3} = A{16}; % C23 (Pa)
A{20}{3,3} = A{16}+2*A{17}; % C33 (Pa)
A{20}{4,4} = A{17}; % C44 (Pa)
A{20}{5,5} = A{17}; % C55 (Pa)
A{20}{6,6} = A{17}; % C66 (Pa)
Materials.Isotropic = struct; % generate a structure where to store the materials' parameters
% disp('vL vP vT')
for i = 1:size(A{1},1) % enter all materials contained in A
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Name',cell2mat(A{1}(i)));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Class','Isotropic');
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Density',A{2}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'YoungsModulus',A{3}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'YoungsModulus_complex',A{18}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'PoissonsNumber',A{4}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'PoissonsNumber_complex',A{19}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'LongitudinalAttenuation',A{5}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'TransverseAttenuation',A{6}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Lambda',A{7}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Lambda_complex',A{16}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Mu',A{8}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Mu_complex',A{17}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'LongitudinalVelocity',A{9}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'LongitudinalVelocity_complex',A{11}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'TransverseVelocity',A{10}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'TransverseVelocity_complex',A{12}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'PlateVelocity',A{13}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'CylinderVelocity',A{14}(i));    
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'RayleighVelocity',A{15}(i));
    Materials.Isotropic = setfield(Materials.Isotropic,{1,1},cell2mat(A{1}(i)),{1,1},'C',[A{20}{1,1}(i) A{20}{1,2}(i) A{20}{1,3}(i) 0 0 0;0 A{20}{2,2}(i) A{20}{2,3}(i) 0 0 0;0 0 A{20}{3,3}(i) 0 0 0;0 0 0 A{20}{4,4}(i) 0 0;0 0 0 0 A{20}{5,5}(i) 0;0 0 0 0 0 A{20}{6,6}(i)]);    
%     if  A{10}(i) < 1500
%         disp([num2str(A{9}(i),'%.0f'),' ',num2str(A{13}(i),'%.0f'),' ',num2str(A{10}(i),'%.0f'),' ',A{1}{i}])
%     end
end
Materials.Isotropic = orderfields(Materials.Isotropic); % sort material entries into alphabetical order

% cubic materials
fileID = fopen('MaterialList_Cubic.txt'); % open txt file containing cubic materials' names, densities, and engineering constants
A = textscan(fileID,'%s %f %f %f %f'); % read the txt file
fclose(fileID); % close the txt file
D = 1-3*A{5}.^2-2*A{5}.^3;
A{6}{1,1} = (1-A{5}.^2)./D.*A{3}; % stiffness component C11 (Pa)
A{6}{1,2} = (A{5}.^2+A{5})./D.*A{3}; % C12 (Pa)
A{6}{1,3} = A{6}{1,2}; % C13 (Pa)
A{6}{2,2} = A{6}{1,1}; % C22 (Pa)
A{6}{2,3} = A{6}{1,2}; % C23 (Pa)
A{6}{3,3} = A{6}{1,1}; % C33 (Pa)
A{6}{4,4} = A{4}; % C44 (Pa)
A{6}{5,5} = A{4}; % C55 (Pa)
A{6}{6,6} = A{4}; % C66 (Pa)
A{7} = sqrt(real(A{6}{1,1})./A{2}); % longitudinal velocity (m/s)
A{8} = sqrt(real(A{6}{4,4})./A{2}); % shear velocity (m/s)
Materials.Cubic = struct; % generate a structure where to store the materials' parameters 
for i = 1:size(A{1},1) % enter all materials contained in A
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'Name',cell2mat(A{1}(i)));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'Class','Cubic');
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'Density',A{2}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'E1',A{3}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'E2',A{3}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'E3',A{3}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'G12',A{4}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'G13',A{4}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'G23',A{4}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'v12',A{5}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'v13',A{5}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'v23',A{5}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'LongitudinalVelocity_1',A{7}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'LongitudinalVelocity_2',A{7}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'LongitudinalVelocity_3',A{7}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'FastShearVelocity_1',A{8}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'FastShearVelocity_2',A{8}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'FastShearVelocity_3',A{8}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'SlowShearVelocity_1',A{8}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'SlowShearVelocity_2',A{8}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'SlowShearVelocity_3',A{8}(i));
    Materials.Cubic = setfield(Materials.Cubic,{1,1},cell2mat(A{1}(i)),{1,1},'C',[A{6}{1,1}(i) A{6}{1,2}(i) A{6}{1,3}(i) 0 0 0;0 A{6}{2,2}(i) A{6}{2,3}(i) 0 0 0;0 0 A{6}{3,3}(i) 0 0 0;0 0 0 A{6}{4,4}(i) 0 0;0 0 0 0 A{6}{5,5}(i) 0;0 0 0 0 0 A{6}{6,6}(i)]);
end
Materials.Cubic = orderfields(Materials.Cubic); % sort material entries into alphabetical order

% transversely isotropic materials
fileID = fopen('MaterialList_TransverselyIsotropic.txt'); % open txt file containing transversely isotropic materials' names, densities, and engineering constants, i.e., Young's moduli E1 and E2, and shear Moduli G12, as well as Poisson's ratios v12 and v23
A = textscan(fileID,'%s %f %f %f %f %f %f'); % read the txt file
fclose(fileID); % close the txt file
A{8} = .5*A{4}./(1+A{7}); % shear modulus G23 (Pa)
A{9} = A{4}.*A{6}./A{3}; % Poisson's ratio v21
D = 1-A{7}.^2-2*A{9}.*A{6}-2*A{6}.*A{7}.*A{9};
A{10}{1,1} = (1-A{7}.^2).*A{3}./D; % stiffness component C11 (Pa) 
A{10}{1,2} = (A{6}.*A{7}+A{6}).*A{4}./D; % C12 (Pa)
A{10}{1,3} = A{10}{1,2}; % C13 (Pa)
A{10}{2,2} = (1-A{6}.*A{9}).*A{4}./D; % C22 (Pa)
A{10}{2,3} = (A{9}.*A{6}+A{7}).*A{4}./D; % C23 (Pa)
A{10}{3,3} = A{10}{2,2}; % C33 (Pa)
A{10}{4,4} = .5*(A{10}{2,2}-A{10}{2,3}); % C44 (Pa)
A{10}{5,5} = A{5}; % C55 (Pa)
A{10}{6,6} = A{10}{5,5}; % C66 (Pa)
A{11} = sqrt(real(A{10}{1,1})./A{2}); % longitudinal velocity along the fibers (m/s)
A{12} = sqrt(real(A{10}{5,5})./A{2}); % fast shear velocity along the fibers (m/s)
A{13} = sqrt(real(A{10}{6,6})./A{2}); % slow shear velocity along the fibers (m/s)
A{14} = sqrt(real(A{10}{2,2})./A{2}); % longitudinal velocity normal to the fibers (m/s)
A{15} = sqrt(real(A{10}{6,6})./A{2}); % fast shear velocity normal to the fibers (m/s)
A{16} = sqrt(real(A{10}{4,4})./A{2}); % slow shear velocity normal to the fibers (m/s)
Materials.TransverselyIsotropic = struct; % generate a structure where to store the materials' parameters 
for i = 1:size(A{1},1) % enter all materials contained in A
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Name',cell2mat(A{1}(i)));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Class','Transversely isotropic');
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Density',A{2}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'E1',A{3}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'E2',A{4}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'E3',A{4}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'G12',A{5}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'G13',A{5}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'G23',A{8}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'v12',A{6}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'v13',A{6}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'v23',A{7}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'LongitudinalVelocity_1',A{11}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'LongitudinalVelocity_2',A{14}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'LongitudinalVelocity_3',A{14}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'FastShearVelocity_1',A{12}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'FastShearVelocity_2',A{15}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'FastShearVelocity_3',A{15}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'SlowShearVelocity_1',A{13}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'SlowShearVelocity_2',A{16}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'SlowShearVelocity_3',A{16}(i));
    Materials.TransverselyIsotropic = setfield(Materials.TransverselyIsotropic,{1,1},cell2mat(A{1}(i)),{1,1},'C',[A{10}{1,1}(i) A{10}{1,2}(i) A{10}{1,3}(i) 0 0 0;0 A{10}{2,2}(i) A{10}{2,3}(i) 0 0 0;0 0 A{10}{3,3}(i) 0 0 0;0 0 0 A{10}{4,4}(i) 0 0;0 0 0 0 A{10}{5,5}(i) 0;0 0 0 0 0 A{10}{6,6}(i)]);
end
Materials.TransverselyIsotropic = orderfields(Materials.TransverselyIsotropic); % sort material entries into alphabetical order

% orthotropic materials
fileID = fopen('MaterialList_Orthotropic.txt'); % open txt file containing orthotropic materials' names, densities, and engineering constants, i.e., Young's moduli E1, E2, E3, and shear Moduli G12, G13, G23, as well as Poisson's ratios v12, v13, v23
A = textscan(fileID,'%s %f %f %f %f %f %f %f %f %f %f'); % read the txt file
fclose(fileID); % close the txt file
A{12} = A{4}./A{3}.*A{9}; % Poisson's ratio v21
A{13} = A{5}./A{3}.*A{10}; % v31
A{14} = A{5}./A{4}.*A{11}; % v32
D = 1-A{9}.*A{12}-A{10}.*A{13}-A{11}.*A{14}-2.*A{9}.*A{11}.*A{13};
A{15}{1,1} = (1-A{11}.*A{14})./D.*A{3}; % stiffness component C11 (Pa)
A{15}{1,2} = (A{10}.*A{14}+A{9})./D.*A{4}; % C12 (Pa)
A{15}{1,3} = (A{9}.*A{11}+A{10})./D.*A{5}; % C13 (Pa)
A{15}{2,2} = (1-A{10}.*A{13})./D.*A{4}; % C22 (Pa)
A{15}{2,3} = (A{12}.*A{10}+A{11})./D.*A{5}; % C23 (Pa)
A{15}{3,3} = (1-A{9}.*A{12})./D.*A{5}; % C33 (Pa)
A{15}{4,4} = A{8}; % C44 (Pa)
A{15}{5,5} = A{7}; % C55 (Pa)
A{15}{6,6} = A{6}; % C66 (Pa)
A{16} = sqrt(real(A{15}{1,1})./A{2}); % longitudinal velocity along the fibers (1-direction) (m/s)
A{17} = sqrt(real(A{15}{5,5})./A{2}); % fast shear velocity along the fibers (1-direction) (m/s)
A{18} = sqrt(real(A{15}{6,6})./A{2}); % slow shear velocity along the fibers (1-direction) (m/s)
A{19} = sqrt(real(A{15}{2,2})./A{2}); % longitudinal velocity normal to the fibers (2-direction) (m/s)
A{20} = sqrt(real(A{15}{6,6})./A{2}); % fast shear velocity normal to the fibers (2-direction) (m/s)
A{21} = sqrt(real(A{15}{4,4})./A{2}); % slow shear velocity normal to the fibers (2-direction) (m/s)
A{22} = sqrt(real(A{15}{3,3})./A{2}); % longitudinal velocity normal to the fibers (3-direction) (m/s)
A{23} = sqrt(real(A{15}{5,5})./A{2}); % fast shear velocity normal to the fibers (3-direction) (m/s)
A{24} = sqrt(real(A{15}{4,4})./A{2}); % slow shear velocity normal to the fibers (3-direction) (m/s)
Materials.Orthotropic = struct; % generate a structure where to store the materials' parameters 
for i = 1:size(A{1},1) % enter all materials contained in A
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Name',cell2mat(A{1}(i)));
    if  A{15}{3,3}(i) == A{15}{1,1}(i) && A{15}{2,2}(i) == A{15}{1,1}(i) && A{15}{1,3}(i) == A{15}{1,2}(i) && A{15}{2,3}(i) == A{15}{1,2}(i) && A{15}{4,4}(i) == A{15}{6,6}(i) && A{15}{5,5}(i) == A{15}{6,6}(i)
        if  abs(A{15}{1,1}(i)-A{15}{1,2}(i)-2*A{15}{4,4}(i))/A{15}{1,1}(i) < 5e-3
            Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Class','Isotropic');
            Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'PlateVelocity',sqrt(A{1,3}(i)/(A{1,2}(i)*(1-A{1,9}(i)^2))));
        else
            Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Class','Cubic');
        end
    elseif A{15}{1,3}(i) == A{15}{1,2}(i) && A{15}{3,3}(i) == A{15}{2,2}(i) && A{15}{5,5}(i) == A{15}{6,6}(i)
        Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Class','Transversely isotropic');
    else
        Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Class','Orthotropic');
    end
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'Density',A{2}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'E1',A{3}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'E2',A{4}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'E3',A{5}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'G12',A{6}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'G13',A{7}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'G23',A{8}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'v12',A{9}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'v13',A{10}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'v23',A{11}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'LongitudinalVelocity_1',A{16}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'LongitudinalVelocity_2',A{19}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'LongitudinalVelocity_3',A{22}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'FastShearVelocity_1',A{17}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'FastShearVelocity_2',A{20}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'FastShearVelocity_3',A{23}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'SlowShearVelocity_1',A{18}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'SlowShearVelocity_2',A{21}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'SlowShearVelocity_3',A{24}(i));
    Materials.Orthotropic = setfield(Materials.Orthotropic,{1,1},cell2mat(A{1}(i)),{1,1},'C',[A{15}{1,1}(i) A{15}{1,2}(i) A{15}{1,3}(i) 0 0 0;0 A{15}{2,2}(i) A{15}{2,3}(i) 0 0 0;0 0 A{15}{3,3}(i) 0 0 0;0 0 0 A{15}{4,4}(i) 0 0;0 0 0 0 A{15}{5,5}(i) 0;0 0 0 0 0 A{15}{6,6}(i)]);
end
Materials.Orthotropic = orderfields(Materials.Orthotropic); % sort material entries into alphabetical order