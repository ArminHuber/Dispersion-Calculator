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
function [Time,u,x1,x3,x3Total,p] = ModeShapeAnimationComputer_Anisotropic(Material,Layers,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,c,A,ALamb,AShear,AScholte,B,BLamb,BShear,BScholte,S,SLamb,SShear,SScholte,CycleDuration,Cycles,FrameRate,Frequency,Length,Mode,PlateThickness,SamplesX1,SamplesX3,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,ShowHalfSpaces,HalfSpaces)
%#ok<*AGROW>
%#ok<*MINV> 
Time = 0:Cycles/(Frequency*1e3*Cycles*CycleDuration*FrameRate):Cycles/(Frequency*1e3);
disp(['Movie duration: ',num2str(Cycles*CycleDuration),' s',newline,'Frames @ ',num2str(FrameRate),'/s: ',num2str(Cycles*CycleDuration*FrameRate)]);
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
    x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
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
        for i = 1:length(x1) % displacements u1,u3 at locations x1 and x3
            for j = 1:length(x3)
                E = [exp(1i*(WaveNumber*x1(i)+Layup{m,3}.*x3(j)-AngularFrequency*Time')) exp(1i*(WaveNumber*x1(i)+Layup{m,3}.*(Layup{m,6}-x3(j))-AngularFrequency*Time'))];
                u{(m-1)*length(x3)+j,i}(:,1) = U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)+U(5)*E(:,5)+U(6)*E(:,6);
                u{(m-1)*length(x3)+j,i}(:,2) = Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)+Layup{m,5}(3)*U(3)*E(:,3)-Layup{m,5}(1)*U(4)*E(:,4)-Layup{m,5}(2)*U(5)*E(:,5)-Layup{m,5}(3)*U(6)*E(:,6);
            end
        end
    end
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
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
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
            for i = 1:length(x1) % displacements u1,u3 at locations x1 and x3
                for j = 1:length(x3)
                    E = [exp(1i*(WaveNumber*x1(i)+Layup{m,3}.*x3(j)-AngularFrequency*Time')) exp(1i*(WaveNumber*x1(i)+Layup{m,3}.*(Layup{m,6}-x3(j))-AngularFrequency*Time'))];
                    u{(m-1)*length(x3)+j,i}(:,1) = U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4);
                    u{(m-1)*length(x3)+j,i}(:,2) = Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)-Layup{m,5}(1)*U(3)*E(:,3)-Layup{m,5}(2)*U(4)*E(:,4);
                end
            end
        end
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
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
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
            for i = 1:length(x1) % displacements u2,u3 at locations x2 and x3
                for j = 1:length(x3)
                    E = [exp(1i*(WaveNumber*x1(i)+Layup{m,3}*x3(j)-AngularFrequency*Time')) exp(1i*(WaveNumber*x1(i)+Layup{m,3}*(Layup{m,6}-x3(j))-AngularFrequency*Time'))];
                    u{(m-1)*length(x3)+j,i}(:,1) = U(1)*E(:,1)+U(2)*E(:,2);
                    u{(m-1)*length(x3)+j,i}(:,2) = 0;
                end
            end
        end
    end
end 
if  SuperLayerSize == 1
    x3Total = x3;
else 
    for i = 2:length(x3Total)
        if  x3Total(i) == x3Total(i-1)
            Ind(i) = i;
        end
    end
    Ind(Ind == 0) = [];
    x3Total(Ind) = [];
    u(Ind,:) = [];
end
if  FluidLoading && ShowHalfSpaces
    x3Fluid = 0:PlateThickness/SamplesX3/Layers:HalfSpaces*PlateThickness;
    if  ToggleUpperFluid
        if  ~contains(Mode,'SH')
            for i = 1:length(x1)
                for j = 1:length(x3Fluid)
                    EUpperFluid = exp(1i*(WaveNumber*x1(i)+k3UpperFluid*x3Fluid(j)-AngularFrequency*Time));
                    uUpperFluid{j,i}(:,1) = EUpperFluid;
                    uUpperFluid{j,i}(:,2) = -WUpperFluid*EUpperFluid;
                end
            end
        elseif contains(Mode,'SH')
            for i = 1:length(x1)
                for j = 1:length(x3Fluid)
                    uUpperFluid{j,i}(length(Time),2) = 0;
                end
            end
        end
        x3Total = horzcat(-fliplr(x3Fluid),x3Total);
        u = vertcat(flipud(uUpperFluid),u);
    end
    if  ToggleLowerFluid
        if  ~contains(Mode,'SH')
            for i = 1:length(x1)
                for j = 1:length(x3Fluid)
                    ELowerFluid = exp(1i*(WaveNumber*x1(i)+k3LowerFluid*x3Fluid(j)-AngularFrequency*Time));
                    uLowerFluid{j,i}(:,1) = ULowerFluid*ELowerFluid;
                    uLowerFluid{j,i}(:,2) = WLowerFluid*ULowerFluid*ELowerFluid;
                end
            end
        elseif contains(Mode,'SH')
            for i = 1:length(x1)
                for j = 1:length(x3Fluid)
                    uLowerFluid{j,i}(length(Time),2) = 0;
                end
            end              
        end
        x3Total = horzcat(x3Total,x3Total(end)+x3Fluid);        
        u = vertcat(u,uLowerFluid);
    end  
end