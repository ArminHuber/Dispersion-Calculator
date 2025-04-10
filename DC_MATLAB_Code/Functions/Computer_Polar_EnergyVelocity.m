% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2025 DLR
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
function A = Computer_Polar_EnergyVelocity(Multithreading,A,c,a11,a12,a21,a22,a23,a31,a32,a33,a34,Material,PlateThickness,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,FrequencyRange,I,SamplesX3)
Layers = Repetitions*length(LayerThicknesses);
if  SymmetricSystem
    Layers = 2*Layers;
end
SamplesX3 = ceil(SamplesX3/Layers); % samples per layer

%#ok<*AGROW>
%#ok<*MINV>
%#ok<*PFBNS>
%#ok<*PFTUSW>
%#ok<*GVMIS>
global Stop 
Stop = 0;
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',size(A,1)),'Name','Calculating energy velocity...');
x3 = (0:PlateThickness/SamplesX3:PlateThickness)';
for n = 1:size(A,1)
    if  Stop
        close(h)
        return
    end
    if  SuperLayerSize == 1
        if  Multithreading
            parfor q = 1:size(A,2)
                EnergyVelocity = [0 0];
                PhaseVelocity = A{n,q}(:,1);
                AngularFrequency = 2*pi*FrequencyRange*1e3;
                for r = 1:height(A{n,q})
                    v = 0; 
                    sigma = 0;
                    epsilon = 0;
                    PowerFlowDensity = 0;
                    Wavenumber = AngularFrequency(r)/PhaseVelocity(r);
                    rc2 = Material{1}.Density*PhaseVelocity(r)^2; 
                    r2c4 = rc2^2;
                    A1 = a11(n,1)+a12(n,1)*rc2;
                    A2 = a21(n,1)+a22(n,1)*rc2+a23(n,1)*r2c4;
                    A3 = a31(n,1)+a32(n,1)*rc2+a33(n,1)*r2c4+a34(n,1)*rc2^3;
                    d1 = A2/3-A1^2/9;
                    d2 = A1^3/27-A1*A2/6+A3/2;
                    d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                    d4 = d1/(2*d3)-d3/2;
                    d5 = d1/d3;
                    d6 = (sqrt(3)*(d3+d5)*1i)/2;
                    Alpha = [sqrt(d4-d6-A1/3) sqrt(d4+d6-A1/3) -sqrt(d3-d5-A1/3)];
                    Alpha2 = Alpha.^2;
                    m11 = c{n,1}(1,1)-rc2+c{n,1}(5,5)*Alpha2;
                    m12 = c{n,1}(1,6)+c{n,1}(4,5)*Alpha2;
                    m13 = (c{n,1}(1,3)+c{n,1}(5,5))*Alpha;
                    m22 = c{n,1}(6,6)-rc2+c{n,1}(4,4)*Alpha2;
                    m23 = (c{n,1}(3,6)+c{n,1}(4,5))*Alpha;
                    V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                    W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                    D1 = 1i*Wavenumber*(c{n,1}(1,1)+c{n,1}(1,6)*V+c{n,1}(1,3)*Alpha.*W); % sigma11
                    D2 = 1i*Wavenumber*(c{n,1}(1,2)+c{n,1}(2,6)*V+c{n,1}(2,3)*Alpha.*W); % sigma22
                    D3 = 1i*Wavenumber*(c{n,1}(1,3)+c{n,1}(3,6)*V+c{n,1}(3,3)*Alpha.*W); % sigma33
                    D4 = 1i*Wavenumber*(c{n,1}(4,5)*(Alpha+W)+c{n,1}(4,4)*Alpha.*V); % sigma23
                    D5 = 1i*Wavenumber*(c{n,1}(5,5)*(Alpha+W)+c{n,1}(4,5)*Alpha.*V); % sigma13
                    D6 = 1i*Wavenumber*(c{n,1}(1,6)+c{n,1}(6,6)*V+c{n,1}(3,6)*Alpha.*W); % sigma12
                    E1 = 1i*Wavenumber; % epsilon11
                    E3 = 1i*Wavenumber*Alpha.*W; % epsilon33
                    E4 = 1i*Wavenumber*Alpha.*V; % epsilon23
                    E5 = 1i*Wavenumber*(Alpha+W); % epsilon13
                    E6 = 1i*Wavenumber*V; % epsilon12
                    E = exp(1i*Wavenumber*Alpha*PlateThickness);
                    L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                    U = L1(2:6,2:6)\-L1(2:6,1);
                    E = [exp(1i*Wavenumber*Alpha.*x3) exp(1i*Wavenumber*Alpha.*(PlateThickness-x3))];
                    v(1:height(E),1) = -1i*AngularFrequency(r)*(E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4)+U(4)*E(:,5)+U(5)*E(:,6));
                    v(1:height(E),2) = -1i*AngularFrequency(r)*(V(1)*E(:,1)+V(2)*U(1)*E(:,2)+V(3)*U(2)*E(:,3)+V(1)*U(3)*E(:,4)+V(2)*U(4)*E(:,5)+V(3)*U(5)*E(:,6));
                    v(1:height(E),3) = -1i*AngularFrequency(r)*(W(1)*E(:,1)+W(2)*U(1)*E(:,2)+W(3)*U(2)*E(:,3)-W(1)*U(3)*E(:,4)-W(2)*U(4)*E(:,5)-W(3)*U(5)*E(:,6));
                    sigma(1:height(E),1) = D1(1)*E(:,1)+D1(2)*U(1)*E(:,2)+D1(3)*U(2)*E(:,3)+D1(1)*U(3)*E(:,4)+D1(2)*U(4)*E(:,5)+D1(3)*U(5)*E(:,6); % sigma11
                    sigma(1:height(E),2) = D2(1)*E(:,1)+D2(2)*U(1)*E(:,2)+D2(3)*U(2)*E(:,3)+D2(1)*U(3)*E(:,4)+D2(2)*U(4)*E(:,5)+D2(3)*U(5)*E(:,6); % sigma22
                    sigma(1:height(E),3) = D3(1)*E(:,1)+D3(2)*U(1)*E(:,2)+D3(3)*U(2)*E(:,3)+D3(1)*U(3)*E(:,4)+D3(2)*U(4)*E(:,5)+D3(3)*U(5)*E(:,6); % sigma33
                    sigma(1:height(E),4) = D4(1)*E(:,1)+D4(2)*U(1)*E(:,2)+D4(3)*U(2)*E(:,3)-D4(1)*U(3)*E(:,4)-D4(2)*U(4)*E(:,5)-D4(3)*U(5)*E(:,6); % sigma23
                    sigma(1:height(E),5) = D5(1)*E(:,1)+D5(2)*U(1)*E(:,2)+D5(3)*U(2)*E(:,3)-D5(1)*U(3)*E(:,4)-D5(2)*U(4)*E(:,5)-D5(3)*U(5)*E(:,6); % sigma13
                    sigma(1:height(E),6) = D6(1)*E(:,1)+D6(2)*U(1)*E(:,2)+D6(3)*U(2)*E(:,3)+D6(1)*U(3)*E(:,4)+D6(2)*U(4)*E(:,5)+D6(3)*U(5)*E(:,6); % sigma12
                    epsilon(1:height(E),1) = E1(1)*(E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4)+U(4)*E(:,5)+U(5)*E(:,6)); % epsilon11                    
                    epsilon(1:height(E),3) = E3(1)*E(:,1)+E3(2)*U(1)*E(:,2)+E3(3)*U(2)*E(:,3)+E3(1)*U(3)*E(:,4)+E3(2)*U(4)*E(:,5)+E3(3)*U(5)*E(:,6); % epsilon33
                    epsilon(1:height(E),4) = E4(1)*E(:,1)+E4(2)*U(1)*E(:,2)+E4(3)*U(2)*E(:,3)-E4(1)*U(3)*E(:,4)-E4(2)*U(4)*E(:,5)-E4(3)*U(5)*E(:,6); % epsilon23
                    epsilon(1:height(E),5) = E5(1)*E(:,1)+E5(2)*U(1)*E(:,2)+E5(3)*U(2)*E(:,3)-E5(1)*U(3)*E(:,4)-E5(2)*U(4)*E(:,5)-E5(3)*U(5)*E(:,6); % epsilon13
                    epsilon(1:height(E),6) = E6(1)*E(:,1)+E6(2)*U(1)*E(:,2)+E6(3)*U(2)*E(:,3)+E6(1)*U(3)*E(:,4)+E6(2)*U(4)*E(:,5)+E6(3)*U(5)*E(:,6); % epsilon12
                    StrainEnergyDensity = .5*(imag(epsilon(:,1)).*imag(sigma(:,1))+imag(epsilon(:,3)).*imag(sigma(:,3))+real(epsilon(:,4)).*real(sigma(:,4))+real(epsilon(:,5)).*real(sigma(:,5))+imag(epsilon(:,6)).*imag(sigma(:,6)));
                    KineticEnergyDensity = .5*(Material{1}.Density*(imag(v(:,1)).^2+imag(v(:,2)).^2+real(v(:,3)).^2));
                    PowerFlowDensity(1:height(E),1) = -.5*(imag(sigma(:,1)).*imag(v(:,1))+imag(sigma(:,6)).*imag(v(:,2))+real(sigma(:,5)).*real(v(:,3)));
                    PowerFlowDensity(1:height(E),2) = -.5*(imag(sigma(:,6)).*imag(v(:,1))+imag(sigma(:,2)).*imag(v(:,2))+real(sigma(:,4)).*real(v(:,3)));
                    PowerFlow = sum(PowerFlowDensity(1:end-1,:)+PowerFlowDensity(2:end,:)); % PowerFlow = x3(2)*sum(PowerFlowDensity(1:end-1,:)+PowerFlowDensity(2:end,:))/2;
                    TotalEnergy = sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = x3(2)*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
                    EnergyVelocity(r,:) = PowerFlow/TotalEnergy; % ce1,ce2 (m/s)
                end
                A{n,q}(:,2:3) = EnergyVelocity;
                A{n,q}(:,2) = filloutliers(A{n,q}(:,2),'spline','movmedian',5,'ThresholdFactor',1);
                A{n,q}(:,3) = filloutliers(A{n,q}(:,3),'spline','movmedian',5,'ThresholdFactor',1);
            end
        else
            for q = 1:size(A,2)
                for r = 1:height(A{n,q})
                    PhaseVelocity = A{n,q}(r,1);
                    AngularFrequency = 2*pi*FrequencyRange(r)*1e3;
                    Wavenumber = AngularFrequency/PhaseVelocity;
                    rc2 = Material{1}.Density*PhaseVelocity^2;
                    r2c4 = rc2^2;
                    A1 = a11(n,1)+a12(n,1)*rc2;
                    A2 = a21(n,1)+a22(n,1)*rc2+a23(n,1)*r2c4;
                    A3 = a31(n,1)+a32(n,1)*rc2+a33(n,1)*r2c4+a34(n,1)*rc2^3;
                    d1 = A2/3-A1^2/9;
                    d2 = A1^3/27-A1*A2/6+A3/2;
                    d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                    d4 = d1/(2*d3)-d3/2;
                    d5 = d1/d3;
                    d6 = (sqrt(3)*(d3+d5)*1i)/2;
                    Alpha(1) = sqrt(d4-d6-A1/3);
                    Alpha(2) = sqrt(d4+d6-A1/3);
                    Alpha(3) = -sqrt(d3-d5-A1/3);
                    Alpha2 = Alpha.^2;
                    m11 = c{n,1}(1,1)-rc2+c{n,1}(5,5)*Alpha2;
                    m12 = c{n,1}(1,6)+c{n,1}(4,5)*Alpha2;
                    m13 = (c{n,1}(1,3)+c{n,1}(5,5))*Alpha;
                    m22 = c{n,1}(6,6)-rc2+c{n,1}(4,4)*Alpha2;
                    m23 = (c{n,1}(3,6)+c{n,1}(4,5))*Alpha;
                    V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                    W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                    D1 = 1i*Wavenumber*(c{n,1}(1,1)+c{n,1}(1,6)*V+c{n,1}(1,3)*Alpha.*W); % sigma11
                    D2 = 1i*Wavenumber*(c{n,1}(1,2)+c{n,1}(2,6)*V+c{n,1}(2,3)*Alpha.*W); % sigma22
                    D3 = 1i*Wavenumber*(c{n,1}(1,3)+c{n,1}(3,6)*V+c{n,1}(3,3)*Alpha.*W); % sigma33
                    D4 = 1i*Wavenumber*(c{n,1}(4,5)*(Alpha+W)+c{n,1}(4,4)*Alpha.*V); % sigma23
                    D5 = 1i*Wavenumber*(c{n,1}(5,5)*(Alpha+W)+c{n,1}(4,5)*Alpha.*V); % sigma13
                    D6 = 1i*Wavenumber*(c{n,1}(1,6)+c{n,1}(6,6)*V+c{n,1}(3,6)*Alpha.*W); % sigma12
                    E1 = 1i*Wavenumber; % epsilon11
                    E3 = 1i*Wavenumber*Alpha.*W; % epsilon33
                    E4 = 1i*Wavenumber*Alpha.*V; % epsilon23
                    E5 = 1i*Wavenumber*(Alpha+W); % epsilon13
                    E6 = 1i*Wavenumber*V; % epsilon12
                    E = exp(1i*Wavenumber*Alpha*PlateThickness);
                    L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                    U = L1(2:6,2:6)\-L1(2:6,1);
                    E = [exp(1i*Wavenumber*Alpha.*x3) exp(1i*Wavenumber*Alpha.*(PlateThickness-x3))];
                    v(:,1) = -1i*AngularFrequency*(E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4)+U(4)*E(:,5)+U(5)*E(:,6));
                    v(:,2) = -1i*AngularFrequency*(V(1)*E(:,1)+V(2)*U(1)*E(:,2)+V(3)*U(2)*E(:,3)+V(1)*U(3)*E(:,4)+V(2)*U(4)*E(:,5)+V(3)*U(5)*E(:,6));
                    v(:,3) = -1i*AngularFrequency*(W(1)*E(:,1)+W(2)*U(1)*E(:,2)+W(3)*U(2)*E(:,3)-W(1)*U(3)*E(:,4)-W(2)*U(4)*E(:,5)-W(3)*U(5)*E(:,6));
                    sigma(:,1) = D1(1)*E(:,1)+D1(2)*U(1)*E(:,2)+D1(3)*U(2)*E(:,3)+D1(1)*U(3)*E(:,4)+D1(2)*U(4)*E(:,5)+D1(3)*U(5)*E(:,6); % sigma11
                    sigma(:,2) = D2(1)*E(:,1)+D2(2)*U(1)*E(:,2)+D2(3)*U(2)*E(:,3)+D2(1)*U(3)*E(:,4)+D2(2)*U(4)*E(:,5)+D2(3)*U(5)*E(:,6); % sigma22
                    sigma(:,3) = D3(1)*E(:,1)+D3(2)*U(1)*E(:,2)+D3(3)*U(2)*E(:,3)+D3(1)*U(3)*E(:,4)+D3(2)*U(4)*E(:,5)+D3(3)*U(5)*E(:,6); % sigma33
                    sigma(:,4) = D4(1)*E(:,1)+D4(2)*U(1)*E(:,2)+D4(3)*U(2)*E(:,3)-D4(1)*U(3)*E(:,4)-D4(2)*U(4)*E(:,5)-D4(3)*U(5)*E(:,6); % sigma23
                    sigma(:,5) = D5(1)*E(:,1)+D5(2)*U(1)*E(:,2)+D5(3)*U(2)*E(:,3)-D5(1)*U(3)*E(:,4)-D5(2)*U(4)*E(:,5)-D5(3)*U(5)*E(:,6); % sigma13
                    sigma(:,6) = D6(1)*E(:,1)+D6(2)*U(1)*E(:,2)+D6(3)*U(2)*E(:,3)+D6(1)*U(3)*E(:,4)+D6(2)*U(4)*E(:,5)+D6(3)*U(5)*E(:,6); % sigma12
                    epsilon(:,1) = E1(1)*(E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4)+U(4)*E(:,5)+U(5)*E(:,6)); % epsilon11                    
                    epsilon(:,3) = E3(1)*E(:,1)+E3(2)*U(1)*E(:,2)+E3(3)*U(2)*E(:,3)+E3(1)*U(3)*E(:,4)+E3(2)*U(4)*E(:,5)+E3(3)*U(5)*E(:,6); % epsilon33
                    epsilon(:,4) = E4(1)*E(:,1)+E4(2)*U(1)*E(:,2)+E4(3)*U(2)*E(:,3)-E4(1)*U(3)*E(:,4)-E4(2)*U(4)*E(:,5)-E4(3)*U(5)*E(:,6); % epsilon23
                    epsilon(:,5) = E5(1)*E(:,1)+E5(2)*U(1)*E(:,2)+E5(3)*U(2)*E(:,3)-E5(1)*U(3)*E(:,4)-E5(2)*U(4)*E(:,5)-E5(3)*U(5)*E(:,6); % epsilon13
                    epsilon(:,6) = E6(1)*E(:,1)+E6(2)*U(1)*E(:,2)+E6(3)*U(2)*E(:,3)+E6(1)*U(3)*E(:,4)+E6(2)*U(4)*E(:,5)+E6(3)*U(5)*E(:,6); % epsilon12
                    StrainEnergyDensity = .5*(imag(epsilon(:,1)).*imag(sigma(:,1))+imag(epsilon(:,3)).*imag(sigma(:,3))+real(epsilon(:,4)).*real(sigma(:,4))+real(epsilon(:,5)).*real(sigma(:,5))+imag(epsilon(:,6)).*imag(sigma(:,6)));
                    KineticEnergyDensity = .5*(Material{1}.Density*(imag(v(:,1)).^2+imag(v(:,2)).^2+real(v(:,3)).^2));
                    PowerFlowDensity(:,1) = -.5*(imag(sigma(:,1)).*imag(v(:,1))+imag(sigma(:,6)).*imag(v(:,2))+real(sigma(:,5)).*real(v(:,3)));
                    PowerFlowDensity(:,2) = -.5*(imag(sigma(:,6)).*imag(v(:,1))+imag(sigma(:,2)).*imag(v(:,2))+real(sigma(:,4)).*real(v(:,3)));
                    PowerFlow = sum(PowerFlowDensity(1:end-1,:)+PowerFlowDensity(2:end,:)); % PowerFlow = x3(2)*sum(PowerFlowDensity(1:end-1,:)+PowerFlowDensity(2:end,:))/2;
                    TotalEnergy = sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = x3(2)*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
                    A{n,q}(r,2:3) = PowerFlow/TotalEnergy; % ce1,ce2 (m/s)
                end
                A{n,q}(:,2) = filloutliers(A{n,q}(:,2),'spline','movmedian',5,'ThresholdFactor',1);
                A{n,q}(:,3) = filloutliers(A{n,q}(:,3),'spline','movmedian',5,'ThresholdFactor',1);
            end
        end
    else
        if  Multithreading
            parfor q = 1:size(A,2)
                EnergyVelocity = [0 0];
                PhaseVelocity = A{n,q}(:,1);
                AngularFrequency = 2*pi*FrequencyRange*1e3;
                for r = 1:height(A{n,q})
                    v = 0; 
                    sigma = 0;
                    epsilon = 0;
                    KineticEnergyDensity = 0;
                    PowerFlowDensity = 0;
                    x3Total = 0;
                    Layup = cell(0);
                    Wavenumber = AngularFrequency(r)/PhaseVelocity(r);
                    for m = 1:SuperLayerSize
                        rc2 = Material{m}.Density*PhaseVelocity(r)^2;
                        r2c4 = rc2^2;
                        A1 = a11(n,m)+a12(n,m)*rc2;
                        A2 = a21(n,m)+a22(n,m)*rc2+a23(n,m)*r2c4;
                        A3 = a31(n,m)+a32(n,m)*rc2+a33(n,m)*r2c4+a34(n,m)*rc2^3;
                        d1 = A2/3-A1^2/9;
                        d2 = A1^3/27-A1*A2/6+A3/2;
                        d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                        d4 = d1/(2*d3)-d3/2;
                        d5 = d1/d3;
                        d6 = (sqrt(3)*(d3+d5)*1i)/2;
                        Alpha = [sqrt(d4-d6-A1/3) sqrt(d4+d6-A1/3) -sqrt(d3-d5-A1/3)];
                        Alpha2 = Alpha.^2;
                        m11 = c{n,m}(1,1)-rc2+c{n,m}(5,5)*Alpha2;
                        m12 = c{n,m}(1,6)+c{n,m}(4,5)*Alpha2;
                        m13 = (c{n,m}(1,3)+c{n,m}(5,5))*Alpha;
                        m22 = c{n,m}(6,6)-rc2+c{n,m}(4,4)*Alpha2;
                        m23 = (c{n,m}(3,6)+c{n,m}(4,5))*Alpha;
                        V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                        W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                        D3 = 1i*Wavenumber*(c{n,m}(1,3)+c{n,m}(3,6)*V+c{n,m}(3,3)*Alpha.*W);
                        D4 = 1i*Wavenumber*(c{n,m}(4,5)*(Alpha+W)+c{n,m}(4,4)*Alpha.*V);
                        D5 = 1i*Wavenumber*(c{n,m}(5,5)*(Alpha+W)+c{n,m}(4,5)*Alpha.*V);
                        E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                        L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                        Layup{m,7} = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                        Layup{m,1} = L1/Layup{m,7}; % L
                        Layup{m,3} = Alpha;
                        Layup{m,4} = V;
                        Layup{m,5} = W;
                        Layup{m,6} = LayerThicknesses(m);
                        Layup{m,8} = D3; % sigma33
                        Layup{m,9} = D5; % sigma13
                        Layup{m,10} = D4; % sigma23
                        Layup{m,11} = 1i*Wavenumber*(c{n,m}(1,1)+c{n,m}(1,6)*V+c{n,m}(1,3)*Alpha.*W); % sigma11
                        Layup{m,12} = 1i*Wavenumber*(c{n,m}(1,2)+c{n,m}(2,6)*V+c{n,m}(2,3)*Alpha.*W); % sigma22
                        Layup{m,13} = 1i*Wavenumber*(c{n,m}(1,6)+c{n,m}(6,6)*V+c{n,m}(3,6)*Alpha.*W); % sigma12
                        Layup{m,14} = 1i*Wavenumber*Alpha.*W; % epsilon33
                        Layup{m,15} = 1i*Wavenumber*(Alpha+W); % epsilon13
                        Layup{m,16} = 1i*Wavenumber*Alpha.*V; % epsilon23
                        Layup{m,17} = 1i*Wavenumber; % epsilon11
                        Layup{m,18} = 1i*Wavenumber*V; % epsilon12
                        Layup{m,19} = Material{m}.Density;
                    end
                    Layup = repmat(Layup,Repetitions,1);
                    if  SymmetricSystem
                        Layup = vertcat(Layup,flipud(Layup));
                    end
                    M = Layup{1};
                    for m = 2:SuperLayerSize
                        N = inv(Layup{m,1}(1:3,1:3)-M(4:6,4:6));
                        M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*Layup{m,1}(1:3,4:6);Layup{m,1}(4:6,1:3)*N*M(4:6,1:3) Layup{m,1}(4:6,4:6)-Layup{m,1}(4:6,1:3)*N*Layup{m,1}(1:3,4:6)];
                    end
                    MM = {M};
                    for m = 2:log2(Repetitions)+1
                        N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                        MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                    end
                    for m = m+1:length(Pattern)
                        N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                        MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                    end
                    if  SymmetricSystem
                        N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                        MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                    end
                    Z1 = [-MM{end}(1:3,1:3)*[1 1 1;Layup{1,4};-Layup{1,5}] -MM{end}(1:3,4:6)*[1 1 1;Layup{end,4};Layup{end,5}];MM{end}(4:6,1:3)*[1 1 1;Layup{1,4};-Layup{1,5}] MM{end}(4:6,4:6)*[1 1 1;Layup{end,4};Layup{end,5}]];
                    Z2 = [MM{end}(1:3,1:3)*[1 1 1;Layup{1,4};Layup{1,5}];-MM{end}(4:6,1:3)*[1 1 1;Layup{1,4};Layup{1,5}]];
                    RT = Z1\Z2(:,1);
                    uInterfaces = [1;Layup{1,4}(1);Layup{1,5}(1)]+[1 1 1;Layup{1,4};-Layup{1,5}]*RT(1:3);
                    if  SymmetricSystem
                        uInterfaces(:,2*Repetitions*SuperLayerSize+1) = [1 1 1;Layup{end,4};Layup{end,5}]*RT(4:6);
                    else
                        uInterfaces(:,Repetitions*SuperLayerSize+1) = [1 1 1;Layup{end,4};Layup{end,5}]*RT(4:6);
                    end
                    Layup{1,2} = Layup{1};
                    for m = 2:size(Layup,1)
                        N = inv(Layup{m,1}(1:3,1:3)-Layup{m-1,2}(4:6,4:6));
                        Layup{m,2} = [Layup{m-1,2}(1:3,1:3)+Layup{m-1,2}(1:3,4:6)*N*Layup{m-1,2}(4:6,1:3) -Layup{m-1,2}(1:3,4:6)*N*Layup{m,1}(1:3,4:6);Layup{m,1}(4:6,1:3)*N*Layup{m-1,2}(4:6,1:3) Layup{m,1}(4:6,4:6)-Layup{m,1}(4:6,1:3)*N*Layup{m,1}(1:3,4:6)];
                    end
                    for m = size(Layup,1):-1:2
                        N = inv(Layup{m,1}(1:3,1:3)-Layup{m-1,2}(4:6,4:6));
                        uInterfaces(:,m) = N*Layup{m-1,2}(4:6,1:3)*uInterfaces(:,1)-N*Layup{m,1}(1:3,4:6)*uInterfaces(:,m+1);
                    end
                    for m = 1:size(Layup,1)
                        U = Layup{m,7}\[uInterfaces(:,m);uInterfaces(:,m+1)];
                        x3 = (0:Layup{m,6}/SamplesX3:Layup{m,6})';
                        if  m == 1
                            x3Total = x3;
                        else
                            x3Total = vertcat(x3Total,x3+x3Total(end));
                        end
                        E = [exp(1i*Wavenumber*Layup{m,3}.*x3) exp(1i*Wavenumber*Layup{m,3}.*(Layup{m,6}-x3))];
                        v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = -1i*AngularFrequency(r)*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)+U(5)*E(:,5)+U(6)*E(:,6));
                        v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = -1i*AngularFrequency(r)*(Layup{m,4}(1)*U(1)*E(:,1)+Layup{m,4}(2)*U(2)*E(:,2)+Layup{m,4}(3)*U(3)*E(:,3)+Layup{m,4}(1)*U(4)*E(:,4)+Layup{m,4}(2)*U(5)*E(:,5)+Layup{m,4}(3)*U(6)*E(:,6));
                        v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = -1i*AngularFrequency(r)*(Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)+Layup{m,5}(3)*U(3)*E(:,3)-Layup{m,5}(1)*U(4)*E(:,4)-Layup{m,5}(2)*U(5)*E(:,5)-Layup{m,5}(3)*U(6)*E(:,6));
                        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,11}(1)*U(1)*E(:,1)+Layup{m,11}(2)*U(2)*E(:,2)+Layup{m,11}(3)*U(3)*E(:,3)+Layup{m,11}(1)*U(4)*E(:,4)+Layup{m,11}(2)*U(5)*E(:,5)+Layup{m,11}(3)*U(6)*E(:,6); % sigma11
                        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = Layup{m,12}(1)*U(1)*E(:,1)+Layup{m,12}(2)*U(2)*E(:,2)+Layup{m,12}(3)*U(3)*E(:,3)+Layup{m,12}(1)*U(4)*E(:,4)+Layup{m,12}(2)*U(5)*E(:,5)+Layup{m,12}(3)*U(6)*E(:,6); % sigma22
                        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,8}(1)*U(1)*E(:,1)+Layup{m,8}(2)*U(2)*E(:,2)+Layup{m,8}(3)*U(3)*E(:,3)+Layup{m,8}(1)*U(4)*E(:,4)+Layup{m,8}(2)*U(5)*E(:,5)+Layup{m,8}(3)*U(6)*E(:,6); % sigma33
                        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,10}(1)*U(1)*E(:,1)+Layup{m,10}(2)*U(2)*E(:,2)+Layup{m,10}(3)*U(3)*E(:,3)-Layup{m,10}(1)*U(4)*E(:,4)-Layup{m,10}(2)*U(5)*E(:,5)-Layup{m,10}(3)*U(6)*E(:,6); % sigma23
                        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(2)*U(2)*E(:,2)+Layup{m,9}(3)*U(3)*E(:,3)-Layup{m,9}(1)*U(4)*E(:,4)-Layup{m,9}(2)*U(5)*E(:,5)-Layup{m,9}(3)*U(6)*E(:,6); % sigma13
                        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,13}(1)*U(1)*E(:,1)+Layup{m,13}(2)*U(2)*E(:,2)+Layup{m,13}(3)*U(3)*E(:,3)+Layup{m,13}(1)*U(4)*E(:,4)+Layup{m,13}(2)*U(5)*E(:,5)+Layup{m,13}(3)*U(6)*E(:,6); % sigma12
                        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,17}*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)+U(5)*E(:,5)+U(6)*E(:,6)); % epsilon11
                        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,14}(1)*U(1)*E(:,1)+Layup{m,14}(2)*U(2)*E(:,2)+Layup{m,14}(3)*U(3)*E(:,3)+Layup{m,14}(1)*U(4)*E(:,4)+Layup{m,14}(2)*U(5)*E(:,5)+Layup{m,14}(3)*U(6)*E(:,6); % epsilon33
                        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,16}(1)*U(1)*E(:,1)+Layup{m,16}(2)*U(2)*E(:,2)+Layup{m,16}(3)*U(3)*E(:,3)-Layup{m,16}(1)*U(4)*E(:,4)-Layup{m,16}(2)*U(5)*E(:,5)-Layup{m,16}(3)*U(6)*E(:,6); % epsilon23
                        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,15}(1)*U(1)*E(:,1)+Layup{m,15}(2)*U(2)*E(:,2)+Layup{m,15}(3)*U(3)*E(:,3)-Layup{m,15}(1)*U(4)*E(:,4)-Layup{m,15}(2)*U(5)*E(:,5)-Layup{m,15}(3)*U(6)*E(:,6); % epsilon13
                        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,18}(1)*U(1)*E(:,1)+Layup{m,18}(2)*U(2)*E(:,2)+Layup{m,18}(3)*U(3)*E(:,3)+Layup{m,18}(1)*U(4)*E(:,4)+Layup{m,18}(2)*U(5)*E(:,5)+Layup{m,18}(3)*U(6)*E(:,6); % epsilon12  
                        KineticEnergyDensity((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = .5*(Layup{m,19}*(imag(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1)).^2+imag(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2)).^2+real(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3)).^2));
                    end
                    StrainEnergyDensity = .5*(imag(epsilon(:,1)).*imag(sigma(:,1))+imag(epsilon(:,3)).*imag(sigma(:,3))+real(epsilon(:,4)).*real(sigma(:,4))+real(epsilon(:,5)).*real(sigma(:,5))+imag(epsilon(:,6)).*imag(sigma(:,6)));
                    PowerFlowDensity(1:height(v),1) = -.5*(imag(sigma(:,1)).*imag(v(:,1))+imag(sigma(:,6)).*imag(v(:,2))+real(sigma(:,5)).*real(v(:,3)));
                    PowerFlowDensity(1:height(v),2) = -.5*(imag(sigma(:,6)).*imag(v(:,1))+imag(sigma(:,2)).*imag(v(:,2))+real(sigma(:,4)).*real(v(:,3)));
                    PowerFlow = trapz(x3Total,PowerFlowDensity);
                    TotalEnergy = trapz(x3Total,StrainEnergyDensity+KineticEnergyDensity)/2;
                    EnergyVelocity(r,:) = PowerFlow/TotalEnergy; % ce1,ce2 (m/s)
                end
                A{n,q}(:,2:3) = EnergyVelocity;
                A{n,q}(:,2) = filloutliers(A{n,q}(:,2),'spline','movmedian',5,'ThresholdFactor',1);
                A{n,q}(:,3) = filloutliers(A{n,q}(:,3),'spline','movmedian',5,'ThresholdFactor',1);
            end
        else
            for q = 1:size(A,2)
                for r = 1:height(A{n,q})
                    Layup = cell(0);
                    PhaseVelocity = A{n,q}(r,1);
                    AngularFrequency = 2*pi*FrequencyRange(r)*1e3;
                    Wavenumber = AngularFrequency/PhaseVelocity;
                    for m = 1:SuperLayerSize
                        rc2 = Material{m}.Density*PhaseVelocity^2;
                        r2c4 = rc2^2;
                        A1 = a11(n,m)+a12(n,m)*rc2;
                        A2 = a21(n,m)+a22(n,m)*rc2+a23(n,m)*r2c4;
                        A3 = a31(n,m)+a32(n,m)*rc2+a33(n,m)*r2c4+a34(n,m)*rc2^3;
                        d1 = A2/3-A1^2/9;
                        d2 = A1^3/27-A1*A2/6+A3/2;
                        d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                        d4 = d1/(2*d3)-d3/2;
                        d5 = d1/d3;
                        d6 = (sqrt(3)*(d3+d5)*1i)/2;
                        Alpha(1) = sqrt(d4-d6-A1/3);
                        Alpha(2) = sqrt(d4+d6-A1/3);
                        Alpha(3) = -sqrt(d3-d5-A1/3);
                        Alpha2 = Alpha.^2;
                        m11 = c{n,m}(1,1)-rc2+c{n,m}(5,5)*Alpha2;
                        m12 = c{n,m}(1,6)+c{n,m}(4,5)*Alpha2;
                        m13 = (c{n,m}(1,3)+c{n,m}(5,5))*Alpha;
                        m22 = c{n,m}(6,6)-rc2+c{n,m}(4,4)*Alpha2;
                        m23 = (c{n,m}(3,6)+c{n,m}(4,5))*Alpha;
                        V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                        W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                        D3 = 1i*Wavenumber*(c{n,m}(1,3)+c{n,m}(3,6)*V+c{n,m}(3,3)*Alpha.*W);
                        D4 = 1i*Wavenumber*(c{n,m}(4,5)*(Alpha+W)+c{n,m}(4,4)*Alpha.*V);
                        D5 = 1i*Wavenumber*(c{n,m}(5,5)*(Alpha+W)+c{n,m}(4,5)*Alpha.*V);
                        E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                        L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                        Layup{m,7} = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                        Layup{m,1} = L1/Layup{m,7}; % L
                        Layup{m,3} = Alpha;
                        Layup{m,4} = V;
                        Layup{m,5} = W;
                        Layup{m,6} = LayerThicknesses(m);
                        Layup{m,8} = D3; % sigma33
                        Layup{m,9} = D5; % sigma13
                        Layup{m,10} = D4; % sigma23
                        Layup{m,11} = 1i*Wavenumber*(c{n,m}(1,1)+c{n,m}(1,6)*V+c{n,m}(1,3)*Alpha.*W); % sigma11
                        Layup{m,12} = 1i*Wavenumber*(c{n,m}(1,2)+c{n,m}(2,6)*V+c{n,m}(2,3)*Alpha.*W); % sigma22
                        Layup{m,13} = 1i*Wavenumber*(c{n,m}(1,6)+c{n,m}(6,6)*V+c{n,m}(3,6)*Alpha.*W); % sigma12
                        Layup{m,14} = 1i*Wavenumber*Alpha.*W; % epsilon33
                        Layup{m,15} = 1i*Wavenumber*(Alpha+W); % epsilon13
                        Layup{m,16} = 1i*Wavenumber*Alpha.*V; % epsilon23
                        Layup{m,17} = 1i*Wavenumber; % epsilon11
                        Layup{m,18} = 1i*Wavenumber*V; % epsilon12
                        Layup{m,19} = Material{m}.Density;
                    end
                    Layup = repmat(Layup,Repetitions,1);
                    if  SymmetricSystem
                        Layup = vertcat(Layup,flipud(Layup));
                    end
                    M = Layup{1};
                    for m = 2:SuperLayerSize
                        N = inv(Layup{m,1}(1:3,1:3)-M(4:6,4:6));
                        M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*Layup{m,1}(1:3,4:6);Layup{m,1}(4:6,1:3)*N*M(4:6,1:3) Layup{m,1}(4:6,4:6)-Layup{m,1}(4:6,1:3)*N*Layup{m,1}(1:3,4:6)];
                    end
                    MM{1} = M;
                    for m = 2:log2(Repetitions)+1
                        N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                        MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                    end
                    for m = m+1:length(Pattern)
                        N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                        MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                    end
                    if  SymmetricSystem
                        N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                        MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                    end
                    Z1 = [-MM{end}(1:3,1:3)*[1 1 1;Layup{1,4};-Layup{1,5}] -MM{end}(1:3,4:6)*[1 1 1;Layup{end,4};Layup{end,5}];MM{end}(4:6,1:3)*[1 1 1;Layup{1,4};-Layup{1,5}] MM{end}(4:6,4:6)*[1 1 1;Layup{end,4};Layup{end,5}]];
                    Z2 = [MM{end}(1:3,1:3)*[1 1 1;Layup{1,4};Layup{1,5}];-MM{end}(4:6,1:3)*[1 1 1;Layup{1,4};Layup{1,5}]];
                    RT = Z1\Z2(:,1);
                    uInterfaces = [1;Layup{1,4}(1);Layup{1,5}(1)]+[1 1 1;Layup{1,4};-Layup{1,5}]*RT(1:3);
                    if  SymmetricSystem
                        uInterfaces(:,2*Repetitions*SuperLayerSize+1) = [1 1 1;Layup{end,4};Layup{end,5}]*RT(4:6);
                    else
                        uInterfaces(:,Repetitions*SuperLayerSize+1) = [1 1 1;Layup{end,4};Layup{end,5}]*RT(4:6);
                    end
                    Layup{1,2} = Layup{1};
                    for m = 2:size(Layup,1)
                        N = inv(Layup{m,1}(1:3,1:3)-Layup{m-1,2}(4:6,4:6));
                        Layup{m,2} = [Layup{m-1,2}(1:3,1:3)+Layup{m-1,2}(1:3,4:6)*N*Layup{m-1,2}(4:6,1:3) -Layup{m-1,2}(1:3,4:6)*N*Layup{m,1}(1:3,4:6);Layup{m,1}(4:6,1:3)*N*Layup{m-1,2}(4:6,1:3) Layup{m,1}(4:6,4:6)-Layup{m,1}(4:6,1:3)*N*Layup{m,1}(1:3,4:6)];
                    end
                    for m = size(Layup,1):-1:2
                        N = inv(Layup{m,1}(1:3,1:3)-Layup{m-1,2}(4:6,4:6));
                        uInterfaces(:,m) = N*Layup{m-1,2}(4:6,1:3)*uInterfaces(:,1)-N*Layup{m,1}(1:3,4:6)*uInterfaces(:,m+1);
                    end
                    for m = 1:size(Layup,1)
                        U = Layup{m,7}\[uInterfaces(:,m);uInterfaces(:,m+1)];
                        x3 = (0:Layup{m,6}/SamplesX3:Layup{m,6})';
                        if  m == 1
                            x3Total = x3;
                        else
                            x3Total = vertcat(x3Total,x3+x3Total(end));
                        end
                        E = [exp(1i*Wavenumber*Layup{m,3}.*x3) exp(1i*Wavenumber*Layup{m,3}.*(Layup{m,6}-x3))];
                        v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = -1i*AngularFrequency*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)+U(5)*E(:,5)+U(6)*E(:,6));
                        v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = -1i*AngularFrequency*(Layup{m,4}(1)*U(1)*E(:,1)+Layup{m,4}(2)*U(2)*E(:,2)+Layup{m,4}(3)*U(3)*E(:,3)+Layup{m,4}(1)*U(4)*E(:,4)+Layup{m,4}(2)*U(5)*E(:,5)+Layup{m,4}(3)*U(6)*E(:,6));
                        v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = -1i*AngularFrequency*(Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)+Layup{m,5}(3)*U(3)*E(:,3)-Layup{m,5}(1)*U(4)*E(:,4)-Layup{m,5}(2)*U(5)*E(:,5)-Layup{m,5}(3)*U(6)*E(:,6));
                        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,11}(1)*U(1)*E(:,1)+Layup{m,11}(2)*U(2)*E(:,2)+Layup{m,11}(3)*U(3)*E(:,3)+Layup{m,11}(1)*U(4)*E(:,4)+Layup{m,11}(2)*U(5)*E(:,5)+Layup{m,11}(3)*U(6)*E(:,6); % sigma11
                        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = Layup{m,12}(1)*U(1)*E(:,1)+Layup{m,12}(2)*U(2)*E(:,2)+Layup{m,12}(3)*U(3)*E(:,3)+Layup{m,12}(1)*U(4)*E(:,4)+Layup{m,12}(2)*U(5)*E(:,5)+Layup{m,12}(3)*U(6)*E(:,6); % sigma22
                        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,8}(1)*U(1)*E(:,1)+Layup{m,8}(2)*U(2)*E(:,2)+Layup{m,8}(3)*U(3)*E(:,3)+Layup{m,8}(1)*U(4)*E(:,4)+Layup{m,8}(2)*U(5)*E(:,5)+Layup{m,8}(3)*U(6)*E(:,6); % sigma33
                        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,10}(1)*U(1)*E(:,1)+Layup{m,10}(2)*U(2)*E(:,2)+Layup{m,10}(3)*U(3)*E(:,3)-Layup{m,10}(1)*U(4)*E(:,4)-Layup{m,10}(2)*U(5)*E(:,5)-Layup{m,10}(3)*U(6)*E(:,6); % sigma23
                        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(2)*U(2)*E(:,2)+Layup{m,9}(3)*U(3)*E(:,3)-Layup{m,9}(1)*U(4)*E(:,4)-Layup{m,9}(2)*U(5)*E(:,5)-Layup{m,9}(3)*U(6)*E(:,6); % sigma13
                        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,13}(1)*U(1)*E(:,1)+Layup{m,13}(2)*U(2)*E(:,2)+Layup{m,13}(3)*U(3)*E(:,3)+Layup{m,13}(1)*U(4)*E(:,4)+Layup{m,13}(2)*U(5)*E(:,5)+Layup{m,13}(3)*U(6)*E(:,6); % sigma12
                        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,17}*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)+U(5)*E(:,5)+U(6)*E(:,6)); % epsilon11
                        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,14}(1)*U(1)*E(:,1)+Layup{m,14}(2)*U(2)*E(:,2)+Layup{m,14}(3)*U(3)*E(:,3)+Layup{m,14}(1)*U(4)*E(:,4)+Layup{m,14}(2)*U(5)*E(:,5)+Layup{m,14}(3)*U(6)*E(:,6); % epsilon33
                        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,16}(1)*U(1)*E(:,1)+Layup{m,16}(2)*U(2)*E(:,2)+Layup{m,16}(3)*U(3)*E(:,3)-Layup{m,16}(1)*U(4)*E(:,4)-Layup{m,16}(2)*U(5)*E(:,5)-Layup{m,16}(3)*U(6)*E(:,6); % epsilon23
                        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,15}(1)*U(1)*E(:,1)+Layup{m,15}(2)*U(2)*E(:,2)+Layup{m,15}(3)*U(3)*E(:,3)-Layup{m,15}(1)*U(4)*E(:,4)-Layup{m,15}(2)*U(5)*E(:,5)-Layup{m,15}(3)*U(6)*E(:,6); % epsilon13
                        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,18}(1)*U(1)*E(:,1)+Layup{m,18}(2)*U(2)*E(:,2)+Layup{m,18}(3)*U(3)*E(:,3)+Layup{m,18}(1)*U(4)*E(:,4)+Layup{m,18}(2)*U(5)*E(:,5)+Layup{m,18}(3)*U(6)*E(:,6); % epsilon12  
                        KineticEnergyDensity((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = .5*(Layup{m,19}*(imag(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1)).^2+imag(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2)).^2+real(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3)).^2));
                    end
                    StrainEnergyDensity = .5*(imag(epsilon(:,1)).*imag(sigma(:,1))+imag(epsilon(:,3)).*imag(sigma(:,3))+real(epsilon(:,4)).*real(sigma(:,4))+real(epsilon(:,5)).*real(sigma(:,5))+imag(epsilon(:,6)).*imag(sigma(:,6)));
                    PowerFlowDensity(:,1) = -.5*(imag(sigma(:,1)).*imag(v(:,1))+imag(sigma(:,6)).*imag(v(:,2))+real(sigma(:,5)).*real(v(:,3)));
                    PowerFlowDensity(:,2) = -.5*(imag(sigma(:,6)).*imag(v(:,1))+imag(sigma(:,2)).*imag(v(:,2))+real(sigma(:,4)).*real(v(:,3)));
                    PowerFlow = trapz(x3Total,PowerFlowDensity);
                    TotalEnergy = trapz(x3Total,StrainEnergyDensity+KineticEnergyDensity)/2;
                    A{n,q}(r,2:3) = PowerFlow/TotalEnergy; % ce1,ce2 (m/s)
                end
                A{n,q}(:,2) = filloutliers(A{n,q}(:,2),'spline','movmedian',5,'ThresholdFactor',1);
                A{n,q}(:,3) = filloutliers(A{n,q}(:,3),'spline','movmedian',5,'ThresholdFactor',1);
            end
        end
    end
    waitbar(n/size(A,1),h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',n,size(A,1),100*n/size(A,1),toc))
end
close(h)