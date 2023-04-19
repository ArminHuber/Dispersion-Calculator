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
function ModeShapeLines_Anisotropic(FunctionMode,MakeFigure,ExportData,XSLX,TXT,MAT,Hybrid,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Crop,LayupString,Plot,PNGresolution,Color1,Color2,Color3,Color4,Color5,Color6,A,ALamb,AShear,AScholte,B,BLamb,BShear,BScholte,S,SLamb,SShear,SScholte,BoxLineWidth,c,Material,Directory,FileName,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FontSizeLegend,Frequency,HeadLine,LegendLocation,LineWidth,Mode,PDF,PNG,PlateThickness,PropagationAngle,SamplesX3,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,ShowHalfSpaces,HalfSpaces,Phase)
%#ok<*AGROW>
%#ok<*MINV>
%#ok<*FNDSB>
I = [1 1 -1;-1 -1 1;-1 -1 1];
I1 = [1 -1;-1 1];
if  ~ToggleUpperFluid
    UpperFluid.Velocity = 1e-10;
    UpperFluid.Density = 1e-10;
end
if  ~ToggleLowerFluid
    LowerFluid.Velocity = 1e-10;
    LowerFluid.Density = 1e-10;
end
p = str2double(regexp(Mode,'\d*','match'))+1;
if  ~Decoupled
    if  ~contains(Mode,'Scholte') && Mode(1) == 'S'
        q = find(S{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(S{p}(:,1))) || Frequency < min(S{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(S{p}(:,1))),' and ',num2str(ceil(max(S{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(S{p}(:,1)-Frequency) == min(abs(S{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = S{p}(q,1);
            end
        end
        PhaseVelocity = S{p}(q,4)*1e3;
        Attenuation = S{p}(q,7)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif ~contains(Mode,'Scholte') && Mode(1) == 'A'
        q = find(A{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(A{p}(:,1))) || Frequency < min(A{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(A{p}(:,1))),' and ',num2str(ceil(max(A{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(A{p}(:,1)-Frequency) == min(abs(A{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = A{p}(q,1);
            end
        end
        PhaseVelocity = A{p}(q,4)*1e3;
        Attenuation = A{p}(q,7)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif ~contains(Mode,'Scholte') && Mode(1) == 'B'
        q = find(B{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(B{p}(:,1))) || Frequency < min(B{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(B{p}(:,1))),' and ',num2str(ceil(max(B{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(B{p}(:,1)-Frequency) == min(abs(B{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = B{p}(q,1);
            end
        end
        PhaseVelocity = B{p}(q,4)*1e3;
        Attenuation = B{p}(q,7)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif contains(Mode,'Scholte') && Mode(1) == 'S' 
        q = find(SScholte{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SScholte{p}(:,1))) || Frequency < min(SScholte{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SScholte{p}(:,1))),' and ',num2str(ceil(max(SScholte{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(SScholte{p}(:,1)-Frequency) == min(abs(SScholte{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = SScholte{p}(q,1);
            end
        end
        PhaseVelocity = SScholte{p}(q,4)*1e3;
        Attenuation = SScholte{p}(q,7)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif contains(Mode,'Scholte') && Mode(1) == 'A'
        q = find(AScholte{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(AScholte{p}(:,1))) || Frequency < min(AScholte{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(AScholte{p}(:,1))),' and ',num2str(ceil(max(AScholte{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(AScholte{p}(:,1)-Frequency) == min(abs(AScholte{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = AScholte{p}(q,1);
            end
        end
        PhaseVelocity = AScholte{p}(q,4)*1e3;
        Attenuation = AScholte{p}(q,7)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif contains(Mode,'Scholte') && Mode(1) == 'B'
        q = find(BScholte{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(BScholte{p}(:,1))) || Frequency < min(BScholte{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(BScholte{p}(:,1))),' and ',num2str(ceil(max(BScholte{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(BScholte{p}(:,1)-Frequency) == min(abs(BScholte{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = BScholte{p}(q,1);
            end
        end
        PhaseVelocity = BScholte{p}(q,4)*1e3;
        Attenuation = BScholte{p}(q,7)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    end
    AngularFrequency = 2*pi*Frequency*1e3;
    WaveNumber = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
    WaveNumber2 = WaveNumber^2;
    WaveNumber4 = WaveNumber^4;
    WaveNumber6 = WaveNumber^6;
    for m = 1:SuperLayerSize
        Delta = c{m}(3,3)*c{m}(4,4)*c{m}(5,5)-c{m}(3,3)*c{m}(4,5)^2;
        rw2 = Material{m}.Density*AngularFrequency^2;
        r2w4 = Material{m}.Density^2*AngularFrequency^4;
        A1 = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta*WaveNumber2+(c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta*rw2;
        A2 = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta*WaveNumber4+(c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta*rw2*WaveNumber2+(c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta*r2w4;
        A3 = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta*WaveNumber6+(c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta*rw2*WaveNumber4+(c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta*r2w4*WaveNumber2-Material{m}.Density^3*AngularFrequency^6/Delta;
        k3a = A2/3-A1^2/9;
        k3b = A1^3/27-A1*A2/6+A3/2;
        k3c = (sqrt(k3b^2+k3a^3)-k3b)^(1/3);
        k3d = k3a/(2*k3c)-k3c/2;
        k3e = k3a/k3c;
        k3f = (sqrt(3)*(k3c+k3e)*1i)/2;
        k3(1) = sqrt(k3d-k3f-A1/3);
        k3(2) = sqrt(k3d+k3f-A1/3);
        k3(3) = -sqrt(k3c-k3e-A1/3);
        k32 = k3.^2;
        m11 = c{m}(1,1)*WaveNumber2-rw2+c{m}(5,5)*k32;
        m12 = c{m}(1,6)*WaveNumber2+c{m}(4,5)*k32;
        m13 = (c{m}(1,3)+c{m}(5,5))*WaveNumber*k3;
        m22 = c{m}(6,6)*WaveNumber2-rw2+c{m}(4,4)*k32;
        m23 = (c{m}(3,6)+c{m}(4,5))*WaveNumber*k3;
        V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
        W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
        D3 = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber*V+c{m}(3,3)*k3.*W); % sigma33
        D4 = 1i*(c{m}(4,5)*(k3+WaveNumber*W)+c{m}(4,4)*k3.*V); % sigma23
        D5 = 1i*(c{m}(5,5)*(k3+WaveNumber*W)+c{m}(4,5)*k3.*V); % sigma13
        E = exp(1i*k3*LayerThicknesses(m));
        L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
        Layup{m,7} = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W]; % L2
        Layup{m,1} = L1/Layup{m,7}; % L
        Layup{m,3} = k3;
        Layup{m,4} = V;
        Layup{m,5} = W;
        Layup{m,6} = LayerThicknesses(m);
        Layup{m,8} = Material{m}.Density;
        Layup{m,9} = 1i*(c{m}(1,1)*WaveNumber+c{m}(1,6)*WaveNumber*V+c{m}(1,3)*k3.*W); % sigma11
        Layup{m,10} = 1i*(c{m}(1,2)*WaveNumber+c{m}(2,6)*WaveNumber*V+c{m}(2,3)*k3.*W); % sigma22
        Layup{m,11} = D3; % sigma33
        Layup{m,12} = D4; % sigma23
        Layup{m,13} = D5; % sigma13
        Layup{m,14} = 1i*(c{m}(1,6)*WaveNumber+c{m}(6,6)*WaveNumber*V+c{m}(3,6)*k3.*W); % sigma12
        Layup{m,15} = 1i*WaveNumber; % epsilon11
        Layup{m,16} = 1i*k3.*W; % epsilon33
        Layup{m,17} = 1i*k3.*V; % epsilon23
        Layup{m,18} = 1i*(k3+WaveNumber*W); % epsilon13
        Layup{m,19} = 1i*WaveNumber*V; % epsilon12
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
    for m = 2:Repetitions
        N = inv(MM{1}(1:3,1:3)-MM{m-1}(4:6,4:6));
        MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{1}(1:3,4:6);MM{1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{1}(4:6,4:6)-MM{1}(4:6,1:3)*N*MM{1}(1:3,4:6)];
    end
    if  SymmetricSystem
        N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
        MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
    end
    if  FluidLoading
        G = inv(MM{end});
        k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2);
        k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2);
        if  Viscoelastic && contains(Mode,'Scholte')
            k3UpperFluid = -k3UpperFluid;
            k3LowerFluid = -k3LowerFluid;
        end
        WUpperFluid = k3UpperFluid/WaveNumber;
        WLowerFluid = k3LowerFluid/WaveNumber;
        DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber; % sigma11, sigma22, sigma33 in the upper fluid
        DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber; % in the lower fluid
        ULowerFluid = -(WUpperFluid+G(3,1)*DUpperFluid)/(G(3,4)*DLowerFluid);
        uInterfaces = G(1:3,1)*DUpperFluid+G(1:3,4)*DLowerFluid*ULowerFluid;
        if  SymmetricSystem
            uInterfaces(:,2*Repetitions*SuperLayerSize+1) = G(4:6,1)*DUpperFluid+G(4:6,4)*DLowerFluid*ULowerFluid;
        else
            uInterfaces(:,Repetitions*SuperLayerSize+1) = G(4:6,1)*DUpperFluid+G(4:6,4)*DLowerFluid*ULowerFluid;
        end
    else
        Z1 = [-MM{end}(1:3,1:3)*[1 1 1;Layup{1,4};-Layup{1,5}] -MM{end}(1:3,4:6)*[1 1 1;Layup{end,4};Layup{end,5}];MM{end}(4:6,1:3)*[1 1 1;Layup{1,4};-Layup{1,5}] MM{end}(4:6,4:6)*[1 1 1;Layup{end,4};Layup{end,5}]];
        Z2 = [MM{end}(1:3,1:3)*[1 1 1;Layup{1,4};Layup{1,5}];-MM{end}(4:6,1:3)*[1 1 1;Layup{1,4};Layup{1,5}]];
        RT = Z1\Z2(:,1);
        uInterfaces = [1;Layup{1,4}(1);Layup{1,5}(1)]+[1 1 1;Layup{1,4};-Layup{1,5}]*RT(1:3);
        if  SymmetricSystem
            uInterfaces(:,2*Repetitions*SuperLayerSize+1) = [1 1 1;Layup{end,4};Layup{end,5}]*RT(4:6);
        else
            uInterfaces(:,Repetitions*SuperLayerSize+1) = [1 1 1;Layup{end,4};Layup{end,5}]*RT(4:6);
        end
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
        x3 = 0:Layup{m,6}/SamplesX3:Layup{m,6};
        if  m == 1
            x3Total = x3;
        else
            x3Total = horzcat(x3Total,x3+x3Total(end));
        end
        E = [exp(1i*Layup{m,3}.*x3') exp(1i*Layup{m,3}.*(Layup{m,6}-x3)')];
        u((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)+U(5)*E(:,5)+U(6)*E(:,6);
        u((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = Layup{m,4}(1)*U(1)*E(:,1)+Layup{m,4}(2)*U(2)*E(:,2)+Layup{m,4}(3)*U(3)*E(:,3)+Layup{m,4}(1)*U(4)*E(:,4)+Layup{m,4}(2)*U(5)*E(:,5)+Layup{m,4}(3)*U(6)*E(:,6);
        u((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)+Layup{m,5}(3)*U(3)*E(:,3)-Layup{m,5}(1)*U(4)*E(:,4)-Layup{m,5}(2)*U(5)*E(:,5)-Layup{m,5}(3)*U(6)*E(:,6);
        v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),:) = -1i*AngularFrequency*u((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),:);
        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(2)*U(2)*E(:,2)+Layup{m,9}(3)*U(3)*E(:,3)+Layup{m,9}(1)*U(4)*E(:,4)+Layup{m,9}(2)*U(5)*E(:,5)+Layup{m,9}(3)*U(6)*E(:,6); % sigma11
        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = Layup{m,10}(1)*U(1)*E(:,1)+Layup{m,10}(2)*U(2)*E(:,2)+Layup{m,10}(3)*U(3)*E(:,3)+Layup{m,10}(1)*U(4)*E(:,4)+Layup{m,10}(2)*U(5)*E(:,5)+Layup{m,10}(3)*U(6)*E(:,6); % sigma22
        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,11}(1)*U(1)*E(:,1)+Layup{m,11}(2)*U(2)*E(:,2)+Layup{m,11}(3)*U(3)*E(:,3)+Layup{m,11}(1)*U(4)*E(:,4)+Layup{m,11}(2)*U(5)*E(:,5)+Layup{m,11}(3)*U(6)*E(:,6); % sigma33
        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,12}(1)*U(1)*E(:,1)+Layup{m,12}(2)*U(2)*E(:,2)+Layup{m,12}(3)*U(3)*E(:,3)-Layup{m,12}(1)*U(4)*E(:,4)-Layup{m,12}(2)*U(5)*E(:,5)-Layup{m,12}(3)*U(6)*E(:,6); % sigma23
        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,13}(1)*U(1)*E(:,1)+Layup{m,13}(2)*U(2)*E(:,2)+Layup{m,13}(3)*U(3)*E(:,3)-Layup{m,13}(1)*U(4)*E(:,4)-Layup{m,13}(2)*U(5)*E(:,5)-Layup{m,13}(3)*U(6)*E(:,6); % sigma13
        sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,14}(1)*U(1)*E(:,1)+Layup{m,14}(2)*U(2)*E(:,2)+Layup{m,14}(3)*U(3)*E(:,3)+Layup{m,14}(1)*U(4)*E(:,4)+Layup{m,14}(2)*U(5)*E(:,5)+Layup{m,14}(3)*U(6)*E(:,6); % sigma12
        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,15}*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)+U(5)*E(:,5)+U(6)*E(:,6)); % epsilon11
        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,16}(1)*U(1)*E(:,1)+Layup{m,16}(2)*U(2)*E(:,2)+Layup{m,16}(3)*U(3)*E(:,3)+Layup{m,16}(1)*U(4)*E(:,4)+Layup{m,16}(2)*U(5)*E(:,5)+Layup{m,16}(3)*U(6)*E(:,6); % epsilon33
        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,17}(1)*U(1)*E(:,1)+Layup{m,17}(2)*U(2)*E(:,2)+Layup{m,17}(3)*U(3)*E(:,3)-Layup{m,17}(1)*U(4)*E(:,4)-Layup{m,17}(2)*U(5)*E(:,5)-Layup{m,17}(3)*U(6)*E(:,6); % epsilon23
        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,18}(1)*U(1)*E(:,1)+Layup{m,18}(2)*U(2)*E(:,2)+Layup{m,18}(3)*U(3)*E(:,3)-Layup{m,18}(1)*U(4)*E(:,4)-Layup{m,18}(2)*U(5)*E(:,5)-Layup{m,18}(3)*U(6)*E(:,6); % epsilon13
        epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,19}(1)*U(1)*E(:,1)+Layup{m,19}(2)*U(2)*E(:,2)+Layup{m,19}(3)*U(3)*E(:,3)+Layup{m,19}(1)*U(4)*E(:,4)+Layup{m,19}(2)*U(5)*E(:,5)+Layup{m,19}(3)*U(6)*E(:,6); % epsilon12  
        KineticEnergyDensity((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = .5*Layup{m,8}*(abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1)).^2+abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2)).^2+abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3)).^2);      
    end
    StrainEnergyDensity = .5*(imag(epsilon(:,1)).*imag(sigma(:,1))+real(epsilon(:,1)).*real(sigma(:,1))+imag(epsilon(:,3)).*imag(sigma(:,3))+real(epsilon(:,3)).*real(sigma(:,3))+imag(epsilon(:,4)).*imag(sigma(:,4))+real(epsilon(:,4)).*real(sigma(:,4))+imag(epsilon(:,5)).*imag(sigma(:,5))+real(epsilon(:,5)).*real(sigma(:,5))+imag(epsilon(:,6)).*imag(sigma(:,6))+real(epsilon(:,6)).*real(sigma(:,6)));
    PowerFlowDensity(:,1) = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,6).*conj(v(:,2))+sigma(:,5).*conj(v(:,3)));
    PowerFlowDensity(:,2) = -.5*real(sigma(:,6).*conj(v(:,1))+sigma(:,2).*conj(v(:,2))+sigma(:,4).*conj(v(:,3)));
    PowerFlowDensity(:,3) = -.5*real(sigma(:,5).*conj(v(:,1))+sigma(:,4).*conj(v(:,2))+sigma(:,3).*conj(v(:,3)));
    if  Phase
        uPhase2 = rad2deg(angle(u(:,2)*exp(-1i*angle(u(1,1)))));
        sigmaPhase456 = rad2deg(angle(sigma(:,4:6)*exp(-1i*angle(u(1,1)))));
        epsilonPhase456 = rad2deg(angle(epsilon(:,4:6)*exp(-1i*angle(u(1,1)))));
    end
    u(:,2) = u(:,2)*exp(-1i*angle(u(2,2))); % zero in the fluid, so we make phase shift now
    sigma(:,4) = sigma(:,4)*exp(-1i*angle(sigma(2,4)));
    sigma(:,5) = sigma(:,5)*exp(-1i*angle(sigma(2,5)));
    sigma(:,6) = sigma(:,6)*exp(-1i*angle(sigma(2,6)));
    epsilon(:,4) = epsilon(:,4)*exp(-1i*angle(epsilon(2,4)));
    epsilon(:,5) = epsilon(:,5)*exp(-1i*angle(epsilon(2,5)));
    epsilon(:,6) = epsilon(:,6)*exp(-1i*angle(epsilon(2,6)));
else
    if  ~contains(Mode,'SH')
        if  ~contains(Mode,'Scholte') && Mode(1) == 'S'
            q = find(SLamb{p}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(SLamb{p}(:,1))) || Frequency < min(SLamb{p}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SLamb{p}(:,1))),' and ',num2str(ceil(max(SLamb{p}(:,1)))),' kHz.'],'Error');
                    return
                else
                    q = find(abs(SLamb{p}(:,1)-Frequency) == min(abs(SLamb{p}(:,1)-Frequency)));
                    q = q(1);
                    Frequency = SLamb{p}(q,1);
                end
            end
            PhaseVelocity = SLamb{p}(q,4)*1e3;
            Attenuation = SLamb{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        elseif ~contains(Mode,'Scholte') && Mode(1) == 'A'
            q = find(ALamb{p}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(ALamb{p}(:,1))) || Frequency < min(ALamb{p}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(ALamb{p}(:,1))),' and ',num2str(ceil(max(ALamb{p}(:,1)))),' kHz.'],'Error');
                    return
                else
                    q = find(abs(ALamb{p}(:,1)-Frequency) == min(abs(ALamb{p}(:,1)-Frequency)));
                    q = q(1);
                    Frequency = ALamb{p}(q,1);
                end
            end
            PhaseVelocity = ALamb{p}(q,4)*1e3;
            Attenuation = ALamb{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        elseif ~contains(Mode,'Scholte') && Mode(1) == 'B'
            q = find(BLamb{p}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(BLamb{p}(:,1))) || Frequency < min(BLamb{p}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(BLamb{p}(:,1))),' and ',num2str(ceil(max(BLamb{p}(:,1)))),' kHz.'],'Error');
                    return
                else
                    q = find(abs(BLamb{p}(:,1)-Frequency) == min(abs(BLamb{p}(:,1)-Frequency)));
                    q = q(1);
                    Frequency = BLamb{p}(q,1);
                end
            end
            PhaseVelocity = BLamb{p}(q,4)*1e3;
            Attenuation = BLamb{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        elseif contains(Mode,'Scholte') && Mode(1) == 'S' 
            q = find(SScholte{p}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(SScholte{p}(:,1))) || Frequency < min(SScholte{p}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SScholte{p}(:,1))),' and ',num2str(ceil(max(SScholte{p}(:,1)))),' kHz.'],'Error');
                    return
                else
                    q = find(abs(SScholte{p}(:,1)-Frequency) == min(abs(SScholte{p}(:,1)-Frequency)));
                    q = q(1);
                    Frequency = SScholte{p}(q,1);
                end
            end
            PhaseVelocity = SScholte{p}(q,4)*1e3;
            Attenuation = SScholte{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        elseif contains(Mode,'Scholte') && Mode(1) == 'A'
            q = find(AScholte{p}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(AScholte{p}(:,1))) || Frequency < min(AScholte{p}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(AScholte{p}(:,1))),' and ',num2str(ceil(max(AScholte{p}(:,1)))),' kHz.'],'Error');
                    return
                else
                    q = find(abs(AScholte{p}(:,1)-Frequency) == min(abs(AScholte{p}(:,1)-Frequency)));
                    q = q(1);
                    Frequency = AScholte{p}(q,1);
                end
            end
            PhaseVelocity = AScholte{p}(q,4)*1e3;
            Attenuation = AScholte{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        elseif contains(Mode,'Scholte') && Mode(1) == 'B'
            q = find(BScholte{p}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(BScholte{p}(:,1))) || Frequency < min(BScholte{p}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(BScholte{p}(:,1))),' and ',num2str(ceil(max(BScholte{p}(:,1)))),' kHz.'],'Error');
                    return
                else
                    q = find(abs(BScholte{p}(:,1)-Frequency) == min(abs(BScholte{p}(:,1)-Frequency)));
                    q = q(1);
                    Frequency = BScholte{p}(q,1);
                end
            end
            PhaseVelocity = BScholte{p}(q,4)*1e3;
            Attenuation = BScholte{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        end
        AngularFrequency = 2*pi*Frequency*1e3;
        WaveNumber = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        WaveNumber2 = WaveNumber^2;
        WaveNumber4 = WaveNumber^4;
        for m = 1:SuperLayerSize
            rw2 = Material{m}.Density*AngularFrequency^2;
            A1(m) = 2*c{m}(3,3)*c{m}(5,5);
            A2 = (c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2)*WaveNumber2-(c{m}(3,3)+c{m}(5,5))*rw2;
            A3 = c{m}(1,1)*c{m}(5,5)*WaveNumber4-(c{m}(1,1)+c{m}(5,5))*rw2*WaveNumber2+Material{m}.Density^2*AngularFrequency^4;
            k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
            k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
            W = (rw2-c{m}(1,1)*WaveNumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*WaveNumber*k3);
            D3 = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,3)*k3.*W); % sigma33
            D5 = 1i*(c{m}(5,5)*(k3+WaveNumber*W)); % sigma13
            E = exp(1i*k3*LayerThicknesses(m));
            L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
            Layup{m,7} = [1 1 E;W -W.*E;E 1 1;W.*E -W]; % L2
            Layup{m,1} = L1/Layup{m,7}; % L
            Layup{m,3} = k3;
            Layup{m,5} = W;
            Layup{m,6} = LayerThicknesses(m);
            Layup{m,8} = Material{m}.Density;
            Layup{m,9} = 1i*(c{m}(1,1)*WaveNumber+c{m}(1,3)*k3.*W); % sigma11
            Layup{m,10} = 1i*(c{m}(1,2)*WaveNumber+c{m}(2,3)*k3.*W); % sigma22
            Layup{m,11} = D3; % sigma33
            Layup{m,12} = D5; % sigma13
            Layup{m,13} = 1i*WaveNumber; % epsilon11
            Layup{m,14} = 1i*k3.*W; % epsilon33
            Layup{m,15} = 1i*(k3+WaveNumber*W); % epsilon13
        end
        Layup = repmat(Layup,Repetitions,1);
        if  SymmetricSystem
            Layup = vertcat(Layup,flipud(Layup));
        end
        M = Layup{1};
        for m = 2:SuperLayerSize
            N = inv(Layup{m,1}(1:2,1:2)-M(3:4,3:4));
            M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*Layup{m,1}(1:2,3:4);Layup{m,1}(3:4,1:2)*N*M(3:4,1:2) Layup{m,1}(3:4,3:4)-Layup{m,1}(3:4,1:2)*N*Layup{m,1}(1:2,3:4)];
        end
        MM{1} = M;
        for m = 2:Repetitions
            N = inv(MM{1}(1:2,1:2)-MM{m-1}(3:4,3:4));
            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{1}(1:2,3:4);MM{1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{1}(3:4,3:4)-MM{1}(3:4,1:2)*N*MM{1}(1:2,3:4)];
        end
        if  SymmetricSystem
            N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
            MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
        end
        if  FluidLoading
            G = inv(MM{end});
            k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2);
            k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2);
            if  Viscoelastic && contains(Mode,'Scholte')
                k3UpperFluid = -k3UpperFluid;
                k3LowerFluid = -k3LowerFluid;
            end
            WUpperFluid = k3UpperFluid/WaveNumber;
            WLowerFluid = k3LowerFluid/WaveNumber;
            DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber; % sigma11, sigma22, sigma33 in the upper fluid
            DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber; % in the lower fluid
            ULowerFluid = -(WUpperFluid+G(2,1)*DUpperFluid)/(G(2,3)*DLowerFluid);
            uInterfaces = G(1:2,1)*DUpperFluid+G(1:2,3)*DLowerFluid*ULowerFluid;
            if  SymmetricSystem
                uInterfaces(:,2*Repetitions*SuperLayerSize+1) = G(3:4,1)*DUpperFluid+G(3:4,3)*DLowerFluid*ULowerFluid;
            else
                uInterfaces(:,Repetitions*SuperLayerSize+1) = G(3:4,1)*DUpperFluid+G(3:4,3)*DLowerFluid*ULowerFluid;
            end
        else
            Z1 = [-MM{end}(1:2,1:2)*[1 1;-Layup{1,5}] -MM{end}(1:2,3:4)*[1 1;Layup{end,5}];MM{end}(3:4,1:2)*[1 1;-Layup{1,5}] MM{end}(3:4,3:4)*[1 1;Layup{end,5}]];
            Z2 = [MM{end}(1:2,1:2)*[1 1;Layup{1,5}];-MM{end}(3:4,1:2)*[1 1;Layup{1,5}]];
            RT = Z1\Z2(:,1);
            uInterfaces = [1;Layup{1,5}(1)]+[1 1;-Layup{1,5}]*RT(1:2);
            if  SymmetricSystem
                uInterfaces(:,2*Repetitions*SuperLayerSize+1) = [1 1;Layup{end,5}]*RT(3:4);
            else
                uInterfaces(:,Repetitions*SuperLayerSize+1) = [1 1;Layup{end,5}]*RT(3:4);
            end
        end
        Layup{1,2} = Layup{1};
        for m = 2:size(Layup,1)
            N = inv(Layup{m,1}(1:2,1:2)-Layup{m-1,2}(3:4,3:4));
            Layup{m,2} = [Layup{m-1,2}(1:2,1:2)+Layup{m-1,2}(1:2,3:4)*N*Layup{m-1,2}(3:4,1:2) -Layup{m-1,2}(1:2,3:4)*N*Layup{m,1}(1:2,3:4);Layup{m,1}(3:4,1:2)*N*Layup{m-1,2}(3:4,1:2) Layup{m,1}(3:4,3:4)-Layup{m,1}(3:4,1:2)*N*Layup{m,1}(1:2,3:4)];
        end
        for m = size(Layup,1):-1:2
            N = inv(Layup{m,1}(1:2,1:2)-Layup{m-1,2}(3:4,3:4));
            uInterfaces(:,m) = N*Layup{m-1,2}(3:4,1:2)*uInterfaces(:,1)-N*Layup{m,1}(1:2,3:4)*uInterfaces(:,m+1);
        end
        for m = 1:size(Layup,1)
            U = Layup{m,7}\[uInterfaces(:,m);uInterfaces(:,m+1)];
            x3 = 0:Layup{m,6}/SamplesX3:Layup{m,6};
            if  m == 1
                x3Total = x3;
            else
                x3Total = horzcat(x3Total,x3+x3Total(end));
            end
            E = [exp(1i*Layup{m,3}.*x3') exp(1i*Layup{m,3}.*(Layup{m,6}-x3)')];
            u((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4);
            u((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)-Layup{m,5}(1)*U(3)*E(:,3)-Layup{m,5}(2)*U(4)*E(:,4);
            v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),[1 3]) = -1i*AngularFrequency*u((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),[1 3]);
            sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(2)*U(2)*E(:,2)+Layup{m,9}(1)*U(3)*E(:,3)+Layup{m,9}(2)*U(4)*E(:,4); % sigma11
            sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = Layup{m,10}(1)*U(1)*E(:,1)+Layup{m,10}(2)*U(2)*E(:,2)+Layup{m,10}(1)*U(3)*E(:,3)+Layup{m,10}(2)*U(4)*E(:,4); % sigma22
            sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,11}(1)*U(1)*E(:,1)+Layup{m,11}(2)*U(2)*E(:,2)+Layup{m,11}(1)*U(3)*E(:,3)+Layup{m,11}(2)*U(4)*E(:,4); % sigma33
            sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,12}(1)*U(1)*E(:,1)+Layup{m,12}(2)*U(2)*E(:,2)-Layup{m,12}(1)*U(3)*E(:,3)-Layup{m,12}(2)*U(4)*E(:,4); % sigma13
            epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,13}*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)); % epsilon11
            epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,14}(1)*U(1)*E(:,1)+Layup{m,14}(2)*U(2)*E(:,2)+Layup{m,14}(1)*U(3)*E(:,3)+Layup{m,14}(2)*U(4)*E(:,4); % epsilon33
            epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,15}(1)*U(1)*E(:,1)+Layup{m,15}(2)*U(2)*E(:,2)-Layup{m,15}(1)*U(3)*E(:,3)-Layup{m,15}(2)*U(4)*E(:,4); % epsilon13
            KineticEnergyDensity((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = .5*Layup{m,8}*(abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1)).^2+abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3)).^2);
        end
        sigma(:,6) = 0;
        epsilon(:,6) = 0;
        StrainEnergyDensity = .5*(imag(epsilon(:,1)).*imag(sigma(:,1))+real(epsilon(:,1)).*real(sigma(:,1))+imag(epsilon(:,3)).*imag(sigma(:,3))+real(epsilon(:,3)).*real(sigma(:,3))+imag(epsilon(:,5)).*imag(sigma(:,5))+real(epsilon(:,5)).*real(sigma(:,5)));
        PowerFlowDensity(:,1) = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)));
        PowerFlowDensity(:,3) = -.5*real(sigma(:,5).*conj(v(:,1))+sigma(:,3).*conj(v(:,3)));
        if  Phase
            sigmaPhase5 = rad2deg(angle(sigma(:,5)*exp(-1i*angle(u(1,1)))));
            epsilonPhase5 = rad2deg(angle(epsilon(:,5)*exp(-1i*angle(u(1,1)))));
        end
        sigma(:,5) = sigma(:,5)*exp(-1i*angle(sigma(2,5))); % zero in the fluid, so we make phase shift now        
        epsilon(:,5) = epsilon(:,5)*exp(-1i*angle(epsilon(2,5)));
    else
        if  Mode(1) == 'S'
            q = find(SShear{p}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(SShear{p}(:,1))) || Frequency < min(SShear{p}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SShear{p}(:,1))),' and ',num2str(ceil(max(SShear{p}(:,1)))),' kHz.'],'Error');
                    return
                else
                    q = find(abs(SShear{p}(:,1)-Frequency) == min(abs(SShear{p}(:,1)-Frequency)));
                    q = q(1);
                    Frequency = SShear{p}(q,1);
                end
            end
            PhaseVelocity = SShear{p}(q,4)*1e3;
            Attenuation = SShear{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        elseif Mode(1) == 'A'
            q = find(AShear{p-1}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(AShear{p-1}(:,1))) || Frequency < min(AShear{p-1}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(AShear{p-1}(:,1))),' and ',num2str(ceil(max(AShear{p-1}(:,1)))),' kHz.'],'Error');
                    return
                else
                    q = find(abs(AShear{p-1}(:,1)-Frequency) == min(abs(AShear{p-1}(:,1)-Frequency)));
                    q = q(1);
                    Frequency = AShear{p-1}(q,1);
                end
            end
            PhaseVelocity = AShear{p-1}(q,4)*1e3;
            Attenuation = AShear{p-1}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        elseif Mode(1) == 'B'
            q = find(BShear{p}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(BShear{p}(:,1))) || Frequency < min(BShear{p}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(BShear{p}(:,1))),' and ',num2str(ceil(max(BShear{p}(:,1)))),' kHz.'],'Error');
                    return
                else
                    q = find(abs(BShear{p}(:,1)-Frequency) == min(abs(BShear{p}(:,1)-Frequency)));
                    q = q(1);
                    Frequency = BShear{p}(q,1);
                end
            end
            PhaseVelocity = BShear{p}(q,4)*1e3;
            Attenuation = BShear{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        end
        AngularFrequency = 2*pi*Frequency*1e3;
        WaveNumber = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        for m = 1:SuperLayerSize
            k3 = sqrt((Material{m}.Density*AngularFrequency^2-WaveNumber^2*c{m}(6,6))/c{m}(4,4));
            if  k3 == 0
                k3 = 1e-10;
            end
            E = exp(1i*k3*LayerThicknesses(m));
            Layup{m,1} = 1i*k3*c{m}(4,4)/(E^2-1)*[-1-E^2 2*E;-2*E 1+E^2];
            Layup{m,3} = k3;
            Layup{m,6} = LayerThicknesses(m);
            Layup{m,7} = [1 E;E 1];
            Layup{m,8} = 1i*k3*c{m}(4,4); % sigma23
            Layup{m,9} = 1i*WaveNumber*c{m}(6,6); % sigma12
            Layup{m,10} = 1i*k3; % epsilon23
            Layup{m,11} = 1i*WaveNumber; % epsilon12
            Layup{m,12} = Material{m}.Density;
        end
        Layup = repmat(Layup,Repetitions,1);
        if  SymmetricSystem
            Layup = vertcat(Layup,flipud(Layup));
        end
        M = Layup{1};
        for m = 2:SuperLayerSize
            N = inv(Layup{m,1}(1,1)-M(2,2));
            M = [M(1,1)+M(1,2)*N*M(2,1) -M(1,2)*N*Layup{m,1}(1,2);Layup{m,1}(2,1)*N*M(2,1) Layup{m,1}(2,2)-Layup{m,1}(2,1)*N*Layup{m,1}(1,2)];
        end
        MM{1} = M;
        for m = 2:Repetitions
            N = inv(MM{1}(1,1)-MM{m-1}(2,2));
            MM{m} = [MM{m-1}(1,1)+MM{m-1}(1,2)*N*MM{m-1}(2,1) -MM{m-1}(1,2)*N*MM{1}(1,2);MM{1}(2,1)*N*MM{m-1}(2,1) MM{1}(2,2)-MM{1}(2,1)*N*MM{1}(1,2)];
        end
        if  SymmetricSystem
            N = inv(-2*MM{end}(2,2));
            MM{end} = [MM{end}(1,1)+MM{end}(1,2)*N*MM{end}(2,1) MM{end}(1,2)*N*MM{end}(2,1);-MM{end}(1,2)*N*MM{end}(2,1) -MM{end}(1,1)-MM{end}(1,2)*N*MM{end}(2,1)];
        end
        Z1 = [MM{end}(1,1) -MM{end}(1,2);MM{end}(2,1) MM{end}(2,2)]; % changed sign in Z1(1,1)!
        Z2 = [MM{end}(1,1);-MM{end}(2,1)];
        RT = Z1\Z2;
        uInterfaces = 1+RT(1);
        if  SymmetricSystem
            uInterfaces(:,2*Repetitions*SuperLayerSize+1) = RT(2);
        else
            uInterfaces(:,Repetitions*SuperLayerSize+1) = RT(2);
        end
        Layup{1,2} = Layup{1};
        for m = 2:size(Layup,1)
            N = inv(Layup{m,1}(1,1)-Layup{m-1,2}(2,2));
            Layup{m,2} = [Layup{m-1,2}(1,1)+Layup{m-1,2}(1,2)*N*Layup{m-1,2}(2,1) -Layup{m-1,2}(1,2)*N*Layup{m,1}(1,2);Layup{m,1}(2,1)*N*Layup{m-1,2}(2,1) Layup{m,1}(2,2)-Layup{m,1}(2,1)*N*Layup{m,1}(1,2)];
        end
        for m = size(Layup,1):-1:2
            N = inv(Layup{m,1}(1,1)-Layup{m-1,2}(2,2));
            uInterfaces(:,m) = N*Layup{m-1,2}(2,1)*uInterfaces(:,1)-N*Layup{m,1}(1,2)*uInterfaces(:,m+1);
        end
        for m = 1:size(Layup,1)
            U = Layup{m,7}\[uInterfaces(:,m);uInterfaces(:,m+1)];
            x3 = 0:Layup{m,6}/SamplesX3:Layup{m,6};
            if  m == 1
                x3Total = x3;
            else
                x3Total = horzcat(x3Total,x3+x3Total(end));
            end
            E = [exp(1i*Layup{m,3}*x3') exp(1i*Layup{m,3}*(Layup{m,6}-x3)')];
            u((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = U(1)*E(:,1)+U(2)*E(:,2);
            v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = -1i*AngularFrequency*u((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2);
            sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,8}(1)*U(1)*E(:,1)-Layup{m,8}(1)*U(2)*E(:,2); % sigma23
            sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(1)*U(2)*E(:,2); % sigma12
            epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,10}(1)*U(1)*E(:,1)-Layup{m,10}(1)*U(2)*E(:,2); % epsilon23
            epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,11}(1)*U(1)*E(:,1)+Layup{m,11}(1)*U(2)*E(:,2); % epsilon12
            KineticEnergyDensity((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = .5*Layup{m,12}*abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2)).^2;
        end
        u(:,3) = 0;
        StrainEnergyDensity = .5*(imag(epsilon(:,4)).*imag(sigma(:,4))+real(epsilon(:,4)).*real(sigma(:,4))+imag(epsilon(:,6)).*imag(sigma(:,6))+real(epsilon(:,6)).*real(sigma(:,6)));
        PowerFlowDensity(:,1) = -.5*real(sigma(:,6).*conj(v(:,2)));
        PowerFlowDensity(:,3) = -.5*real(sigma(:,4).*conj(v(:,2)));
        if  Phase
            uPhase = rad2deg(angle(u));
            sigmaPhase = rad2deg(angle(sigma));
            epsilonPhase = rad2deg(angle(epsilon));
        end
        u(:,2) = u(:,2)*exp(-1i*angle(u(2,2)));
        sigma(:,4) = sigma(:,4)*exp(-1i*angle(sigma(2,4)));
        sigma(:,6) = sigma(:,6)*exp(-1i*angle(sigma(2,6)));
        epsilon(:,4) = epsilon(:,4)*exp(-1i*angle(epsilon(2,4)));
        epsilon(:,6) = epsilon(:,6)*exp(-1i*angle(epsilon(2,6)));
    end
end
PowerFlow = trapz(x3Total,PowerFlowDensity(:,1));
if  FluidLoading && ShowHalfSpaces
    x3Fluid = 0:PlateThickness/SamplesX3:HalfSpaces*PlateThickness;
    if  ~contains(Mode,'SH')
        if  ToggleUpperFluid
            EUpperFluid = exp(1i*k3UpperFluid*x3Fluid);
        end
        if  ToggleLowerFluid
            ELowerFluid = exp(1i*k3LowerFluid*x3Fluid);
        end
    end
    E1Fluid = 1i*WaveNumber; % epsilon11
    if  ToggleUpperFluid
        if  ~contains(Mode,'SH')
            E3UpperFluid = 1i*k3UpperFluid*WUpperFluid; % epsilon33
            uUpperFluid(:,1) = EUpperFluid;
            uUpperFluid(:,3) = -WUpperFluid*EUpperFluid;
            vUpperFluid = -1i*AngularFrequency*uUpperFluid;
            sigmaUpperFluid(:,1) = DUpperFluid*EUpperFluid; % sigma11
            sigmaUpperFluid(:,2) = sigmaUpperFluid(:,1); % sigma22
            sigmaUpperFluid(:,3) = sigmaUpperFluid(:,1); % sigma33
            sigmaUpperFluid(:,6) = 0;
            epsilonUpperFluid(:,1) = E1Fluid*EUpperFluid; % epsilon11
            epsilonUpperFluid(:,3) = E3UpperFluid*EUpperFluid; % epsilon33
            epsilonUpperFluid(:,6) = 0;
            StrainEnergyDensityUpperFluid = .5*(imag(epsilonUpperFluid(:,1)).*imag(sigmaUpperFluid(:,1))+real(epsilonUpperFluid(:,1)).*real(sigmaUpperFluid(:,1))+imag(epsilonUpperFluid(:,3)).*imag(sigmaUpperFluid(:,3))+real(epsilonUpperFluid(:,3)).*real(sigmaUpperFluid(:,3)));
            KineticEnergyDensityUpperFluid = .5*(UpperFluid.Density*(abs(vUpperFluid(:,1)).^2+abs(vUpperFluid(:,3)).^2));
            PowerFlowDensityUpperFluid(:,1) = -.5*real(sigmaUpperFluid(:,1).*conj(vUpperFluid(:,1)));
            PowerFlowDensityUpperFluid(:,3) = -.5*real(sigmaUpperFluid(:,3).*conj(vUpperFluid(:,3)));
        elseif contains(Mode,'SH')
            uUpperFluid(length(x3Fluid),3) = 0;
            sigmaUpperFluid(length(x3Fluid),6) = 0;
            epsilonUpperFluid(length(x3Fluid),6) = 0;
            StrainEnergyDensityUpperFluid(length(x3Fluid),1) = 0;
            KineticEnergyDensityUpperFluid(length(x3Fluid),1) = 0;
            PowerFlowDensityUpperFluid(length(x3Fluid),3) = 0;
            if  Phase
                uPhase = vertcat(zeros(length(x3Fluid),3),uPhase);
                sigmaPhase = vertcat(zeros(length(x3Fluid),6),sigmaPhase);
                epsilonPhase = vertcat(zeros(length(x3Fluid),6),epsilonPhase);
            end
        end
        x3Total = horzcat(-fliplr(x3Fluid),x3Total);
        u = vertcat(flipud(uUpperFluid),u);
        sigma = vertcat(flipud(sigmaUpperFluid),sigma);
        epsilon = vertcat(flipud(epsilonUpperFluid),epsilon);
        StrainEnergyDensity = vertcat(flipud(StrainEnergyDensityUpperFluid),StrainEnergyDensity);
        KineticEnergyDensity = vertcat(flipud(KineticEnergyDensityUpperFluid),KineticEnergyDensity);
        PowerFlowDensity = vertcat(flipud(PowerFlowDensityUpperFluid),PowerFlowDensity);
    end
    if  ToggleLowerFluid
        if  ~contains(Mode,'SH')
            E3LowerFluid = 1i*k3LowerFluid*WLowerFluid; % epsilon33
            uLowerFluid(:,1) = ULowerFluid*ELowerFluid;
            uLowerFluid(:,3) = WLowerFluid*ULowerFluid*ELowerFluid;
            vLowerFluid = -1i*AngularFrequency*uLowerFluid;
            sigmaLowerFluid(:,1) = DLowerFluid*ULowerFluid*ELowerFluid; % sigma11
            sigmaLowerFluid(:,2) = sigmaLowerFluid(:,1); % sigma22
            sigmaLowerFluid(:,3) = sigmaLowerFluid(:,1); % sigma33
            sigmaLowerFluid(:,6) = 0;
            epsilonLowerFluid(:,1) = E1Fluid*ULowerFluid*ELowerFluid; % epsilon11
            epsilonLowerFluid(:,3) = E3LowerFluid*ULowerFluid*ELowerFluid; % epsilon33
            epsilonLowerFluid(:,6) = 0;
            StrainEnergyDensityLowerFluid = .5*(imag(epsilonLowerFluid(:,1)).*imag(sigmaLowerFluid(:,1))+real(epsilonLowerFluid(:,1)).*real(sigmaLowerFluid(:,1))+imag(epsilonLowerFluid(:,3)).*imag(sigmaLowerFluid(:,3))+real(epsilonLowerFluid(:,3)).*real(sigmaLowerFluid(:,3)));
            KineticEnergyDensityLowerFluid = .5*(LowerFluid.Density*(abs(vLowerFluid(:,1)).^2+abs(vLowerFluid(:,3)).^2));
            PowerFlowDensityLowerFluid(:,1) = -.5*real(sigmaLowerFluid(:,1).*conj(vLowerFluid(:,1)));
            PowerFlowDensityLowerFluid(:,3) = -.5*real(sigmaLowerFluid(:,3).*conj(vLowerFluid(:,3)));
        elseif contains(Mode,'SH')
            uLowerFluid(length(x3Fluid),3) = 0;
            sigmaLowerFluid(length(x3Fluid),6) = 0;
            epsilonLowerFluid(length(x3Fluid),6) = 0;
            StrainEnergyDensityLowerFluid(length(x3Fluid),1) = 0;
            KineticEnergyDensityLowerFluid(length(x3Fluid),1) = 0;
            PowerFlowDensityLowerFluid(length(x3Fluid),3) = 0;
            if  Phase
                uPhase = vertcat(uPhase,zeros(length(x3Fluid),3));
                sigmaPhase = vertcat(sigmaPhase,zeros(length(x3Fluid),6));
                epsilonPhase = vertcat(epsilonPhase,zeros(length(x3Fluid),6));
            end
        end
        x3Total = horzcat(x3Total,x3Total(end)+x3Fluid);
        u = vertcat(u,uLowerFluid);
        sigma = vertcat(sigma,sigmaLowerFluid);
        epsilon = vertcat(epsilon,epsilonLowerFluid);
        StrainEnergyDensity = vertcat(StrainEnergyDensity,StrainEnergyDensityLowerFluid);
        KineticEnergyDensity = vertcat(KineticEnergyDensity,KineticEnergyDensityLowerFluid);
        PowerFlowDensity = vertcat(PowerFlowDensity,PowerFlowDensityLowerFluid);        
    end
end
if  ~contains(Mode,'SH')
    if  Phase
        if  FluidLoading && ShowHalfSpaces
            if  ToggleUpperFluid
                uPhase = rad2deg(angle(u*exp(-1i*angle(u(length(x3Fluid)+1,1)))));
                sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(length(x3Fluid)+1,1)))));
                epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(length(x3Fluid)+1,1)))));
                if  ~Decoupled
                    uPhase(length(x3Fluid)+1:length(x3Fluid)+length(uPhase2),2) = uPhase2;
                    sigmaPhase(length(x3Fluid)+1:length(x3Fluid)+length(sigmaPhase456),4:6) = sigmaPhase456;
                    epsilonPhase(length(x3Fluid)+1:length(x3Fluid)+length(epsilonPhase456),4:6) = epsilonPhase456;
                else
                    sigmaPhase(length(x3Fluid)+1:length(x3Fluid)+length(sigmaPhase5),5) = sigmaPhase5;
                    epsilonPhase(length(x3Fluid)+1:length(x3Fluid)+length(epsilonPhase5),5) = epsilonPhase5;
                end
            else
                uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
                sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
                epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
                if  ~Decoupled
                    uPhase(1:length(uPhase2),2) = uPhase2;
                    sigmaPhase(1:length(sigmaPhase456),4:6) = sigmaPhase456;
                    epsilonPhase(1:length(epsilonPhase456),4:6) = epsilonPhase456;
                else
                    sigmaPhase(1:length(sigmaPhase5),5) = sigmaPhase5;
                    epsilonPhase(1:length(epsilonPhase5),5) = epsilonPhase5;
                end
            end
        else
            uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
            sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
            epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
            if  ~Decoupled
                uPhase(:,2) = uPhase2;
                sigmaPhase(:,4:6) = sigmaPhase456;
                epsilonPhase(:,4:6) = epsilonPhase456;
            else
                sigmaPhase(:,5) = sigmaPhase5;
                epsilonPhase(:,5) = epsilonPhase5;
            end
        end
    end
    u(:,1) = u(:,1)*exp(-1i*angle(u(2,1)));
    u(:,3) = u(:,3)*exp(-1i*angle(u(2,3)));
    sigma(:,1) = sigma(:,1)*exp(-1i*angle(sigma(2,1)));
    sigma(:,2) = sigma(:,2)*exp(-1i*angle(sigma(2,2)));
    sigma(:,3) = sigma(:,3)*exp(-1i*angle(sigma(2,3)));
    epsilon(:,1) = epsilon(:,1)*exp(-1i*angle(epsilon(2,1)));
    epsilon(:,3) = epsilon(:,3)*exp(-1i*angle(epsilon(2,3)));
end
u = u/sqrt(PowerFlow);
sigma = sigma/sqrt(PowerFlow);
epsilon = epsilon/sqrt(PowerFlow);
if  real(u(2,1)) < 0
    u = -u;
end
if  real(sigma(2,1)) < 0
    sigma = -sigma;
end
if  real(epsilon(2,1)) < 0
    epsilon = -epsilon;
end
if  Phase
    uPhase(find(round(uPhase) == -180)) = 180;
    sigmaPhase(find(round(sigmaPhase) == -180)) = 180;
    epsilonPhase(find(round(epsilonPhase) == -180)) = 180;
end
StrainEnergyDensity = .5*StrainEnergyDensity/PowerFlow;
KineticEnergyDensity = .5*KineticEnergyDensity/PowerFlow;    
TotalEnergyDensity = StrainEnergyDensity+KineticEnergyDensity;
PowerFlowDensity = PowerFlowDensity/PowerFlow;
if  ExportData
    if  FunctionMode == 1
        if  ~Phase
            Table = table('Size',[length(x3) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'x3 (mm)','u1 (nm)','u2 (nm)','u3 (nm)'});
            Table(:,1) = num2cell(-1e3*x3');
            Table(:,2) = num2cell(real(u(:,1))*1e9);
            Table(:,3) = num2cell(real(u(:,2))*1e9);
            Table(:,4) = num2cell(real(u(:,3))*1e9);
        else
            Table = table('Size',[length(x3) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','u1 (nm)','u2 (nm)','u3 (nm)','Phase u1 (deg)','Phase u2 (deg)','Phase u3 (deg)'});
            Table(:,1) = num2cell(-1e3*x3');
            Table(:,2) = num2cell(real(u(:,1))*1e9);
            Table(:,3) = num2cell(real(u(:,2))*1e9);
            Table(:,4) = num2cell(real(u(:,3))*1e9);
            Table(:,5) = num2cell(uPhase(:,1));
            Table(:,6) = num2cell(uPhase(:,2));
            Table(:,7) = num2cell(uPhase(:,3));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Displacement_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Displacement_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.txt']));
            end
            if  MAT
                Displacement = Table;
                save(fullfile(Directory,[FileName,'_Displacement_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.mat']),'Displacement')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end    
    elseif FunctionMode == 2
        if  ~Phase
            Table = table('Size',[length(x3) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','sigma11 (kPa)','sigma22 (kPa)','sigma33 (kPa)','sigma23 (kPa)','sigma13 (kPa)','sigma12 (kPa)'});
            Table(:,1) = num2cell(-1e3*x3');
            Table(:,2) = num2cell(real(sigma(:,1))/1e3);
            Table(:,3) = num2cell(real(sigma(:,2))/1e3);
            Table(:,4) = num2cell(real(sigma(:,3))/1e3);
            Table(:,5) = num2cell(real(sigma(:,4))/1e3);
            Table(:,6) = num2cell(real(sigma(:,5))/1e3);
            Table(:,7) = num2cell(real(sigma(:,6))/1e3);
        else
            Table = table('Size',[length(x3) 13],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','sigma11 (kPa)','sigma22 (kPa)','sigma33 (kPa)','sigma23 (kPa)','sigma13 (kPa)','sigma12 (kPa)','Phase sigma11 (deg)','Phase sigma22 (deg)','Phase sigma33 (deg)','Phase sigma23 (deg)','Phase sigma13 (deg)','Phase sigma12 (deg)'});
            Table(:,1) = num2cell(-1e3*x3');
            Table(:,2) = num2cell(real(sigma(:,1))/1e3);
            Table(:,3) = num2cell(real(sigma(:,2))/1e3);
            Table(:,4) = num2cell(real(sigma(:,3))/1e3);
            Table(:,5) = num2cell(real(sigma(:,4))/1e3);
            Table(:,6) = num2cell(real(sigma(:,5))/1e3);
            Table(:,7) = num2cell(real(sigma(:,6))/1e3);
            Table(:,8) = num2cell(sigmaPhase(:,1));
            Table(:,9) = num2cell(sigmaPhase(:,2));
            Table(:,10) = num2cell(sigmaPhase(:,3));
            Table(:,11) = num2cell(sigmaPhase(:,4));
            Table(:,12) = num2cell(sigmaPhase(:,5));
            Table(:,13) = num2cell(sigmaPhase(:,6));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Stress_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Stress_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.txt']));
            end
            if  MAT
                Stress = Table;
                save(fullfile(Directory,[FileName,'_Stress_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.mat']),'Stress')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end        
    elseif FunctionMode == 3 
        if  ~Phase
            Table = table('Size',[length(x3) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','epsilon11','epsilon22','epsilon33','epsilon23','epsilon13','epsilon12'});
            Table(:,1) = num2cell(-1e3*x3');
            Table(:,2) = num2cell(real(epsilon(:,1)));
            Table(:,3) = num2cell(real(epsilon(:,2)));
            Table(:,4) = num2cell(real(epsilon(:,3)));
            Table(:,5) = num2cell(real(epsilon(:,4)));
            Table(:,6) = num2cell(real(epsilon(:,5)));
            Table(:,7) = num2cell(real(epsilon(:,6)));
        else
            Table = table('Size',[length(x3) 13],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','epsilon11','epsilon22','epsilon33','epsilon23','epsilon13','epsilon12','Phase epsilon11 (deg)','Phase epsilon22 (deg)','Phase epsilon33 (deg)','Phase epsilon23 (deg)','Phase epsilon13 (deg)','Phase epsilon12 (deg)'});
            Table(:,1) = num2cell(-1e3*x3');
            Table(:,2) = num2cell(real(epsilon(:,1)));
            Table(:,3) = num2cell(real(epsilon(:,2)));
            Table(:,4) = num2cell(real(epsilon(:,3)));
            Table(:,5) = num2cell(real(epsilon(:,4)));
            Table(:,6) = num2cell(real(epsilon(:,5)));
            Table(:,7) = num2cell(real(epsilon(:,6)));
            Table(:,8) = num2cell(epsilonPhase(:,1));
            Table(:,9) = num2cell(epsilonPhase(:,2));
            Table(:,10) = num2cell(epsilonPhase(:,3));
            Table(:,11) = num2cell(epsilonPhase(:,4));
            Table(:,12) = num2cell(epsilonPhase(:,5));
            Table(:,13) = num2cell(epsilonPhase(:,6));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Strain_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Strain_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.txt']));
            end
            if  MAT
                Strain = Table;
                save(fullfile(Directory,[FileName,'_Strain_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.mat']),'Strain')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end 
    elseif FunctionMode == 4   
        Table = table('Size',[length(x3Total) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'x3 (mm)','Estrain (J/m2)','Ekin (J/m2)','Etotal (J/m2)'});
        Table(:,1) = num2cell(-1e3*x3Total');
        Table(:,2) = num2cell(StrainEnergyDensity);
        Table(:,3) = num2cell(KineticEnergyDensity);
        Table(:,4) = num2cell(TotalEnergyDensity);
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.txt']));
            end
            if  MAT
                EnergyDensity = Table;
                save(fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.mat']),'EnergyDensity')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end
    elseif FunctionMode == 5
        Table = table('Size',[length(x3Total) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'x3 (mm)','P1 (W/m)','P2 (W/m)','P3 (W/m)'});
        Table(:,1) = num2cell(-1e3*x3Total');
        Table(:,2) = num2cell(PowerFlowDensity(:,1));
        Table(:,3) = num2cell(PowerFlowDensity(:,2));
        Table(:,4) = num2cell(PowerFlowDensity(:,3));
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.txt']));
            end
            if  MAT
                PowerFlowDensity = Table;
                save(fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.mat']),'PowerFlowDensity')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end
    end
end
if  MakeFigure
    if  Phase
        if  FunctionMode == 1 
            if  strcmp(LegendLocation,'in')
                f = figure('Name','Displacement phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 .72 1],'color','w');
            else
                f = figure('Name','Displacement phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 .905 1],'color','w');
            end
        elseif FunctionMode == 2
            if  strcmp(LegendLocation,'in')
                f = figure('Name','Stress phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 .72 1],'color','w');
            else
                f = figure('Name','Stress phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 .88 1],'color','w');
            end
        elseif FunctionMode == 3
            if  strcmp(LegendLocation,'in')
                f = figure('Name','Strain phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 .72 1],'color','w');
            else
                f = figure('Name','Strain phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 .88 1],'color','w');
            end
        end    
        datacursormode on
        jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
        jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));    
        x = xline(0,'Color',[.6 .6 .6]);
        hasbehavior(x,'legend',false);   
        hold on
        if  FunctionMode == 1
            if  Plot(3)
                plot(uPhase(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;
            end
            if  Plot(5)
                plot(uPhase(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;
            end
            if  Plot(4)
                plot(uPhase(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;
            end
        elseif FunctionMode == 2
            if  Plot(3)
                plot(sigmaPhase(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(sigmaPhase(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(sigmaPhase(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(sigmaPhase(:,5),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(sigmaPhase(:,4),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(sigmaPhase(:,6),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end 
        elseif FunctionMode == 3
            if  Plot(3)
                plot(epsilonPhase(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(epsilonPhase(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(epsilonPhase(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(epsilonPhase(:,5),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(epsilonPhase(:,4),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(epsilonPhase(:,6),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end
        end
        ax = gca;
        ax.Box = 'on';
        ax.LineWidth = BoxLineWidth;
        ax.FontSize = FontSizeAxes;
        ax.Title.Interpreter = 'latex';
        ax.Title.FontSize = FontSizeHeadLine;
        if  Hybrid
            Material{1}.Name = 'hybrid';
        end
        if  HeadLine > 0
            if  ~contains(Mode,'Scholte') && ~contains(Mode,'SH')
                ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
            elseif contains(Mode,'SH')
                ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
            elseif contains(Mode,'Scholte')
                ModeName = [Mode(1),'$^{\mathrm{Scholte}}_{',num2str(p-1),'}$'];
            end
            if  HeadLine == 1
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_'))];
            elseif ~SymmetricSystem && HeadLine == 2
                if  Repetitions == 1
                    String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']'];
                else
                    String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
                end
            elseif SymmetricSystem && HeadLine == 2
                if  Repetitions == 1
                    String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{\mathrm s}$'];
                else
                    String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
                end
            end
            if  FluidLoading
                if  ToggleUpperFluid && ToggleLowerFluid
                    String = append(String,' in ',char(join(split(UpperFluid.Name,'_'),'\_')),'/',char(join(split(LowerFluid.Name,'_'),'\_')));
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    String = append(String,' in ',char(join(split(UpperFluid.Name,'_'),'\_')),'/vacuum');
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    String = append(String,' in vacuum/',char(join(split(LowerFluid.Name,'_'),'\_')));
                end
            end
            ax.Title.String = String;
        end
        ax.XLabel.Interpreter = 'latex';
        ax.XLabel.FontSize = FontSizeAxesLabels;
        ax.XLim = [-max(abs(ax.XLim)) max(abs(ax.XLim))];
        ax.YLabel.Interpreter = 'latex';
        ax.YLabel.FontSize = FontSizeAxesLabels;
        ax.YLabel.String = '$x_3$ (mm)';    
        ax.TickLabelInterpreter = 'latex';    
        ax.YLim = -1e3*[x3Total(end) x3Total(1)];
        if  FluidLoading && ShowHalfSpaces
            if  ToggleUpperFluid && ToggleLowerFluid
                if  strcmp(UpperFluid.Name,LowerFluid.Name)
                    UpperFaceAlpha = .2;
                    LowerFaceAlpha = .2;
                else
                    if  UpperFluid.Density*UpperFluid.Velocity > LowerFluid.Density*LowerFluid.Velocity
                        UpperFaceAlpha = .2;
                        LowerFaceAlpha = .1;
                    else
                        UpperFaceAlpha = .1;
                        LowerFaceAlpha = .2;
                    end
                end
                patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',UpperFaceAlpha,'EdgeColor','none')
                patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[ax.YLim(1) ax.YLim(1) -PlateThickness*1e3 -PlateThickness*1e3],'b','FaceAlpha',LowerFaceAlpha,'EdgeColor','none')
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',.2,'EdgeColor','none')
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[ax.YLim(1) ax.YLim(1) -PlateThickness*1e3 -PlateThickness*1e3],'b','FaceAlpha',.2,'EdgeColor','none')
            end
        end
        if  FunctionMode == 1
            if  Plot(3)
                plot(uPhase(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;
            end
            if  Plot(5)
                plot(uPhase(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;
            end
            if  Plot(4)
                plot(uPhase(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;
            end
        elseif FunctionMode == 2
            if  Plot(3)
                plot(sigmaPhase(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(sigmaPhase(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(sigmaPhase(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(sigmaPhase(:,5),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(sigmaPhase(:,4),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(sigmaPhase(:,6),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end 
        elseif FunctionMode == 3
            if  Plot(3)
                plot(epsilonPhase(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(epsilonPhase(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(epsilonPhase(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(epsilonPhase(:,5),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(epsilonPhase(:,4),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(epsilonPhase(:,6),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end
        end
        if  FunctionMode == 1
            ax.XLabel.String = 'Displacement phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$u_3$','$u_1$','$u_2$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')
            else
                LegendNames = {'Out-of-plane ($u_3$)','In-plane ($u_1$)','Shear horizontal ($u_2$)'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')
            end
            if  Export
                try
                    if  PDF
                        if  ~Crop
                            if  strcmp(LegendLocation,'in')
                                set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[38 30])
                            else
                                set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[45 30])
                            end
                            print(f,fullfile(Directory,[FileName,'_Phase']),'-dpdf');
                        else
                            exportgraphics(f,fullfile(Directory,[FileName,'_Phase.pdf']))
                        end
                    end
                    if  PNG
                        if  ~Crop
                            print(f,fullfile(Directory,[FileName,'_Phase']),'-dpng',['-r',num2str(PNGresolution)])
                        else
                            exportgraphics(f,fullfile(Directory,[FileName,'_Phase.png']),'Resolution',PNGresolution)
                        end
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                    return
                end
            end
        elseif FunctionMode == 2
            ax.XLabel.String = 'Stress phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$\sigma_{33}$','$\sigma_{11}$','$\sigma_{22}$','$\sigma_{13}$','$\sigma_{23}$','$\sigma_{12}$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
            else
                LegendNames = {'Out-of-plane ($\sigma_{33}$)','In-plane ($\sigma_{11}$)','In-plane ($\sigma_{22}$)','Shear ($\sigma_{13}$)','Shear ($\sigma_{23}$)','Shear ($\sigma_{12}$)'}; 
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
            end
            if  Export
                try
                    if  PDF
                        if  ~Crop
                            if  strcmp(LegendLocation,'in')
                                set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[38 30])
                            else
                                set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[45 30])
                            end
                            print(f,fullfile(Directory,[FileName,'_Phase']),'-dpdf');
                        else
                            exportgraphics(f,fullfile(Directory,[FileName,'_Phase.pdf']))
                        end
                    end
                    if  PNG
                        if  ~Crop
                            print(f,fullfile(Directory,[FileName,'_Phase']),'-dpng',['-r',num2str(PNGresolution)])
                        else
                            exportgraphics(f,fullfile(Directory,[FileName,'_Phase.png']),'Resolution',PNGresolution)
                        end
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                    return
                end
            end
        elseif FunctionMode == 3
            ax.XLabel.String = 'Strain phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$\varepsilon_{33}$','$\varepsilon_{11}$','$\varepsilon_{22}$','$\varepsilon_{13}$','$\varepsilon_{23}$','$\varepsilon_{12}$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
            else
                LegendNames = {'Out-of-plane ($\varepsilon_{33}$)','In-plane ($\varepsilon_{11}$)','In-plane ($\varepsilon_{22}$)','Shear ($\varepsilon_{13}$)','Shear ($\varepsilon_{23}$)','Shear ($\varepsilon_{12}$)'}; 
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
            end 
            if  Export
                try
                    if  PDF
                        if  ~Crop
                            if  strcmp(LegendLocation,'in')
                                set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[38 30])
                            else
                                set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[45 30])
                            end
                            print(f,fullfile(Directory,[FileName,'_Phase']),'-dpdf');
                        else
                            exportgraphics(f,fullfile(Directory,[FileName,'_Phase.pdf']))
                        end
                    end
                    if  PNG
                        if  ~Crop
                            print(f,fullfile(Directory,[FileName,'_Phase']),'-dpng',['-r',num2str(PNGresolution)])
                        else
                            exportgraphics(f,fullfile(Directory,[FileName,'_Phase.png']),'Resolution',PNGresolution)
                        end
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                    return
                end
            end
        end
        tb = axtoolbar('default');
        tb.Visible = 'on';
        d = datacursormode(f);
        d.Interpreter = 'latex';
        d.UpdateFcn = @CursorPhase;
    end    
    if  FunctionMode == 1 
        if  strcmp(LegendLocation,'in')
            f = figure('Name','Displacement','Toolbar','none','Units','normalized','OuterPosition',[0 0 .72 1],'color','w');
        else
            f = figure('Name','Displacement','Toolbar','none','Units','normalized','OuterPosition',[0 0 .905 1],'color','w');
        end
    elseif FunctionMode == 2
        if  strcmp(LegendLocation,'in')
            f = figure('Name','Stress','Toolbar','none','Units','normalized','OuterPosition',[0 0 .72 1],'color','w');
        else
            f = figure('Name','Stress','Toolbar','none','Units','normalized','OuterPosition',[0 0 .88 1],'color','w');
        end
    elseif FunctionMode == 3
        if  strcmp(LegendLocation,'in')
            f = figure('Name','Strain','Toolbar','none','Units','normalized','OuterPosition',[0 0 .72 1],'color','w');
        else
            f = figure('Name','Strain','Toolbar','none','Units','normalized','OuterPosition',[0 0 .88 1],'color','w');
        end
    elseif FunctionMode == 4
        if  strcmp(LegendLocation,'in')
            f = figure('Name','Energy density','Toolbar','none','Units','normalized','OuterPosition',[0 0 .72 1],'colorPlot(3) == 1','w');
        else
            f = figure('Name','Energy density','Toolbar','none','Units','normalized','OuterPosition',[0 0 .88 1],'color','w');
        end
    elseif FunctionMode == 5 
        if  strcmp(LegendLocation,'in')
            f = figure('Name','Power flow density','Toolbar','none','Units','normalized','OuterPosition',[0 0 .72 1],'color','w');
        else
            f = figure('Name','Power flow density','Toolbar','none','Units','normalized','OuterPosition',[0 0 .88 1],'color','w');
        end
    end    
    datacursormode on
    jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));    
    x = xline(0,'Color',[.6 .6 .6]);
    hasbehavior(x,'legend',false);   
    hold on
    if  FunctionMode == 1
        if  Plot(3)
            plot(real(u(:,3))*1e9,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(real(u(:,1))*1e9,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;
        end
        if  Plot(4)
            plot(real(u(:,2))*1e9,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;
        end
    elseif FunctionMode == 2
        if  Plot(3)
            plot(real(sigma(:,3))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(sigma(:,1))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(sigma(:,2))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(sigma(:,5))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(sigma(:,4))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(sigma(:,6))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end 
    elseif FunctionMode == 3
        if  Plot(3)
            plot(real(epsilon(:,3)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(epsilon(:,1)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(epsilon(:,2)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(epsilon(:,5)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(epsilon(:,4)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(epsilon(:,6)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end
    elseif FunctionMode == 4
        if  Plot(5)
            plot(StrainEnergyDensity,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z = 1;
        else
            z = 0;
        end
        if  Plot(4)
            plot(KineticEnergyDensity,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(3)
            plot(TotalEnergyDensity,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z(3) = 1;
        else
            z(3) = 0;        
        end        
    elseif FunctionMode == 5
        if  Plot(3)
            plot(PowerFlowDensity(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(PowerFlowDensity(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(4)
            plot(PowerFlowDensity(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
    end
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = BoxLineWidth;
    ax.FontSize = FontSizeAxes;
    ax.Title.Interpreter = 'latex';
    ax.Title.FontSize = FontSizeHeadLine;
    if  Hybrid
        Material{1}.Name = 'hybrid';
    end
    if  HeadLine > 0
        if  ~contains(Mode,'Scholte') && ~contains(Mode,'SH')
            ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
        elseif contains(Mode,'SH')
            ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
        elseif contains(Mode,'Scholte')
            ModeName = [Mode(1),'$^{\mathrm{Scholte}}_{',num2str(p-1),'}$'];
        end
        if  HeadLine == 1
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_'))];
        elseif ~SymmetricSystem && HeadLine == 2
            if  Repetitions == 1
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']'];
            else
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
            end
        elseif SymmetricSystem && HeadLine == 2
            if  Repetitions == 1
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{\mathrm s}$'];
            else
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
            end
        end
        if  FluidLoading
            if  ToggleUpperFluid && ToggleLowerFluid
                String = append(String,' in ',char(join(split(UpperFluid.Name,'_'),'\_')),'/',char(join(split(LowerFluid.Name,'_'),'\_')));
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                String = append(String,' in ',char(join(split(UpperFluid.Name,'_'),'\_')),'/vacuum');
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                String = append(String,' in vacuum/',char(join(split(LowerFluid.Name,'_'),'\_')));
            end
        end
        ax.Title.String = String;
    end
    ax.XLabel.Interpreter = 'latex';
    ax.XLabel.FontSize = FontSizeAxesLabels;
    ax.XLim = [-max(abs(ax.XLim)) max(abs(ax.XLim))];
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.FontSize = FontSizeAxesLabels;
    ax.YLabel.String = '$x_3$ (mm)';    
    ax.TickLabelInterpreter = 'latex';    
    ax.YLim = -1e3*[x3Total(end) x3Total(1)];
    if  FluidLoading && ShowHalfSpaces
        if  ToggleUpperFluid && ToggleLowerFluid
            if  strcmp(UpperFluid.Name,LowerFluid.Name)
                UpperFaceAlpha = .2;
                LowerFaceAlpha = .2;
            else
                if  UpperFluid.Density*UpperFluid.Velocity > LowerFluid.Density*LowerFluid.Velocity
                    UpperFaceAlpha = .2;
                    LowerFaceAlpha = .1;
                else
                    UpperFaceAlpha = .1;
                    LowerFaceAlpha = .2;
                end
            end
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',UpperFaceAlpha,'EdgeColor','none')
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[ax.YLim(1) ax.YLim(1) -PlateThickness*1e3 -PlateThickness*1e3],'b','FaceAlpha',LowerFaceAlpha,'EdgeColor','none')
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',.2,'EdgeColor','none')
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[ax.YLim(1) ax.YLim(1) -PlateThickness*1e3 -PlateThickness*1e3],'b','FaceAlpha',.2,'EdgeColor','none')
        end
    end
    if  FunctionMode == 1
        if  Plot(3)
            plot(real(u(:,3))*1e9,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(real(u(:,1))*1e9,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;
        end
        if  Plot(4)
            plot(real(u(:,2))*1e9,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;
        end 
    elseif FunctionMode == 2
        if  Plot(3)
            plot(real(sigma(:,3))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(sigma(:,1))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(sigma(:,2))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(sigma(:,5))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(sigma(:,4))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(sigma(:,6))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end 
    elseif FunctionMode == 3
        if  Plot(3)
            plot(real(epsilon(:,3)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(epsilon(:,1)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(epsilon(:,2)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(epsilon(:,5)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(epsilon(:,4)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(epsilon(:,6)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end
    elseif FunctionMode == 4
        if  Plot(5)
            plot(StrainEnergyDensity,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z = 1;
        else
            z = 0;
        end
        if  Plot(4)
            plot(KineticEnergyDensity,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(3)
            plot(TotalEnergyDensity,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z(3) = 1;
        else
            z(3) = 0;        
        end        
    elseif FunctionMode == 5
        if  Plot(3)
            plot(PowerFlowDensity(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(PowerFlowDensity(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(4)
            plot(PowerFlowDensity(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
    end
    if  FunctionMode == 1
        ax.XLabel.String = 'Displacement (nm)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$u_3$','$u_1$','$u_2$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')
        else
            LegendNames = {'Out-of-plane ($u_3$)','In-plane ($u_1$)','Shear horizontal ($u_2$)'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')
        end
        if  Export
            try
                if  PDF
                    if  ~Crop
                        if  strcmp(LegendLocation,'in')
                            set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[38 30])
                        else
                            set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[45 30])
                        end
                        print(f,fullfile(Directory,FileName),'-dpdf');
                    else
                        exportgraphics(f,fullfile(Directory,[FileName,'.pdf']))
                    end
                end
                if  PNG
                    if  ~Crop
                        print(f,fullfile(Directory,FileName),'-dpng',['-r',num2str(PNGresolution)])
                    else
                        exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                    end
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end
    elseif FunctionMode == 2
        ax.XLabel.String = 'Stress (kPa)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$\sigma_{33}$','$\sigma_{11}$','$\sigma_{22}$','$\sigma_{13}$','$\sigma_{23}$','$\sigma_{12}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Out-of-plane ($\sigma_{33}$)','In-plane ($\sigma_{11}$)','In-plane ($\sigma_{22}$)','Shear ($\sigma_{13}$)','Shear ($\sigma_{23}$)','Shear ($\sigma_{12}$)'}; 
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
        end
        if  Export
            try
                if  PDF
                    if  ~Crop
                        if  strcmp(LegendLocation,'in')
                            set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[38 30])
                        else
                            set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[45 30])
                        end
                        print(f,fullfile(Directory,FileName),'-dpdf');
                    else
                        exportgraphics(f,fullfile(Directory,[FileName,'.pdf']))
                    end
                end
                if  PNG
                    if  ~Crop
                        print(f,fullfile(Directory,FileName),'-dpng',['-r',num2str(PNGresolution)])
                    else
                        exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                    end
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end
    elseif FunctionMode == 3
        ax.XLabel.String = 'Strain';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$\varepsilon_{33}$','$\varepsilon_{11}$','$\varepsilon_{22}$','$\varepsilon_{13}$','$\varepsilon_{23}$','$\varepsilon_{12}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Out-of-plane ($\varepsilon_{33}$)','In-plane ($\varepsilon_{11}$)','In-plane ($\varepsilon_{22}$)','Shear ($\varepsilon_{13}$)','Shear ($\varepsilon_{23}$)','Shear ($\varepsilon_{12}$)'}; 
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
        end 
        if  Export
            try
                if  PDF
                    if  ~Crop
                        if  strcmp(LegendLocation,'in')
                            set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[38 30])
                        else
                            set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[45 30])
                        end
                        print(f,fullfile(Directory,FileName),'-dpdf');
                    else
                        exportgraphics(f,fullfile(Directory,[FileName,'.pdf']))
                    end
                end
                if  PNG
                    if  ~Crop
                        print(f,fullfile(Directory,FileName),'-dpng',['-r',num2str(PNGresolution)])
                    else
                        exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                    end
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end
    elseif FunctionMode == 4
        ax.XLabel.String = 'Energy density (J/m$^2$)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$E_\mathrm{strain}$','$E_\mathrm{kin}$','$E_\mathrm{total}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Strain energy density','Kinetic energy density','Total energy density'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')        
        end
        if  Export
            try
                if  PDF
                    if  ~Crop
                        if  strcmp(LegendLocation,'in')
                            set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[38 30])
                        else
                            set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[45 30])
                        end
                        print(f,fullfile(Directory,FileName),'-dpdf');
                    else
                        exportgraphics(f,fullfile(Directory,[FileName,'.pdf']))
                    end
                end
                if  PNG
                    if  ~Crop
                        print(f,fullfile(Directory,FileName),'-dpng',['-r',num2str(PNGresolution)])
                    else
                        exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                    end
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end         
    elseif FunctionMode == 5
        ax.XLabel.String = 'Power flow density (W/m)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$p_3$','$p_1$','$p_2$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Out-of-plane ($p_3$)','In-plane ($p_1$)','Shear horizontal ($p_2$)'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')        
        end
        if  Export
            try
                if  PDF
                    if  ~Crop
                        if  strcmp(LegendLocation,'in')
                            set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[38 30])
                        else
                            set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[45 30])
                        end
                        print(f,fullfile(Directory,FileName),'-dpdf');
                    else
                        exportgraphics(f,fullfile(Directory,[FileName,'.pdf']))
                    end
                end
                if  PNG
                    if  ~Crop
                        print(f,fullfile(Directory,FileName),'-dpng',['-r',num2str(PNGresolution)])
                    else
                        exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                    end
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end
    end
    tb = axtoolbar('default');
    tb.Visible = 'on';
    d = datacursormode(f);
    d.Interpreter = 'latex';
    d.UpdateFcn = @Cursor;
end
function output_txt = Cursor(~,event_obj)
    if  length(event_obj.Target.XData) == 4
        if  event_obj.Target.YData(1) == 0
            output_txt = char(join(split(UpperFluid.Name,'_'),'\_'));
        else
            output_txt = char(join(split(LowerFluid.Name,'_'),'\_'));
        end
    else
        if  FunctionMode == 1 
            if  event_obj.Target.Color == Color1
                output_txt = {['$u_1$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$u_2$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$u_3$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 2
            if  event_obj.Target.Color == Color3
                output_txt = {['$\sigma_{33}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\sigma_{13}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\sigma_{23}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\sigma_{11}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\sigma_{22}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\sigma_{12}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 3
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varepsilon_{33}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varepsilon_{13}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varepsilon_{23}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varepsilon_{11}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varepsilon_{22}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varepsilon_{12}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 4
            if  event_obj.Target.Color == Color1
                output_txt = {['$E_\mathrm{strain}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$E_\mathrm{kin}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$E_\mathrm{total}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 5 
            if  event_obj.Target.Color == Color1
                output_txt = {['$p_1$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$p_2$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$p_3$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        end
    end
end    
function output_txt = CursorPhase(~,event_obj)
    if  length(event_obj.Target.XData) == 4
        if  event_obj.Target.YData(1) == 0
            output_txt = char(join(split(UpperFluid.Name,'_'),'\_'));
        else
            output_txt = char(join(split(LowerFluid.Name,'_'),'\_'));
        end
    else
        if  FunctionMode == 1 
            if  event_obj.Target.Color == Color1
                output_txt = {['$\varphi(u_1)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(u_2)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$\varphi(u_3)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 2
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varphi(\sigma_{33})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varphi(\sigma_{13})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varphi(\sigma_{23})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varphi(\sigma_{11})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(\sigma_{22})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varphi(\sigma_{12})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 3
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varphi(\varepsilon_{33})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varphi(\varepsilon_{13})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varphi(\varepsilon_{23})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varphi(\varepsilon_{11})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(\varepsilon_{22})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varphi(\varepsilon_{12})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        end
    end
end
end