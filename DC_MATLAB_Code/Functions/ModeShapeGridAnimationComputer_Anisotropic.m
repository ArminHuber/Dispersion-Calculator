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
function [Time,u,x1,x3Total,p] = ModeShapeGridAnimationComputer_Anisotropic(Material,Layers,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,c,A,ALamb,AShear,AScholte,B,BLamb,BShear,BScholte,S,SLamb,SShear,SScholte,CycleDuration,FrameRate,Frequency,Length,Mode,PlateThickness,SamplesX1,SamplesX3,I,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,ShowHalfSpaces,HalfSpaces)
%#ok<*AGROW>
%#ok<*MINV> 
Time = (0:1/(Frequency*1e3*CycleDuration*FrameRate):1/(Frequency*1e3))';
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
                [~,q] = min(abs(S{p}(:,1)-Frequency));
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
                [~,q] = min(abs(A{p}(:,1)-Frequency));
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
                [~,q] = min(abs(B{p}(:,1)-Frequency));
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
                [~,q] = min(abs(SScholte{p}(:,1)-Frequency));
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
                [~,q] = min(abs(AScholte{p}(:,1)-Frequency));
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
                [~,q] = min(abs(BScholte{p}(:,1)-Frequency));
                Frequency = BScholte{p}(q,1);
            end
        end
        PhaseVelocity = BScholte{p}(q,4)*1e3;
        Attenuation = BScholte{p}(q,7)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    end
    x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
    AngularFrequency = 2*pi*Frequency*1e3;
    AngularFrequency2 = AngularFrequency^2;
    Wavenumber = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
    Wavenumber2 = Wavenumber^2;
    Wavenumber4 = Wavenumber2^2;
    Wavenumber6 = Wavenumber2^3;
    for m = 1:SuperLayerSize
        Delta = c{m}(3,3)*c{m}(4,4)*c{m}(5,5)-c{m}(3,3)*c{m}(4,5)^2;
        rw2 = Material{m}.Density*AngularFrequency2;
        r2w4 = rw2^2;
        A1 = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta*Wavenumber2+(c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta*rw2;
        A2 = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta*Wavenumber4+(c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta*rw2*Wavenumber2+(c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta*r2w4;
        A3 = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta*Wavenumber6+(c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta*rw2*Wavenumber4+(c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta*r2w4*Wavenumber2-rw2^3/Delta;
        d1 = A2/3-A1^2/9;
        d2 = A1^3/27-A1*A2/6+A3/2;
        d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
        d4 = d1/(2*d3)-d3/2;
        d5 = d1/d3;
        d6 = (sqrt(3)*(d3+d5)*1i)/2;
        k3(1) = sqrt(d4-d6-A1/3);
        k3(2) = sqrt(d4+d6-A1/3);
        k3(3) = -sqrt(d3-d5-A1/3);
        k32 = k3.^2;
        m11 = c{m}(1,1)*Wavenumber2-rw2+c{m}(5,5)*k32;
        m12 = c{m}(1,6)*Wavenumber2+c{m}(4,5)*k32;
        m13 = (c{m}(1,3)+c{m}(5,5))*Wavenumber*k3;
        m22 = c{m}(6,6)*Wavenumber2-rw2+c{m}(4,4)*k32;
        m23 = (c{m}(3,6)+c{m}(4,5))*Wavenumber*k3;
        V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
        W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
        D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber*V+c{m}(3,3)*k3.*W); % sigma33
        D4 = 1i*(c{m}(4,5)*(k3+Wavenumber*W)+c{m}(4,4)*k3.*V); % sigma23
        D5 = 1i*(c{m}(5,5)*(k3+Wavenumber*W)+c{m}(4,5)*k3.*V); % sigma13
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
    if  FluidLoading
        G = inv(MM{end});
        k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-Wavenumber2);
        k3LowerFluid = sqrt(AngularFrequency2/LowerFluid.Velocity^2-Wavenumber2);
        if  contains(Mode,'Scholte') && Attenuation ~= 0
            if  PhaseVelocity < UpperFluid.Velocity
                k3UpperFluid = -k3UpperFluid;
            end
            if  PhaseVelocity < LowerFluid.Velocity
                k3LowerFluid = -k3LowerFluid;
            end
        end
        WUpperFluid = k3UpperFluid/Wavenumber;
        WLowerFluid = k3LowerFluid/Wavenumber;
        DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2/Wavenumber; % sigma11, sigma22, sigma33 in the upper fluid
        DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2/Wavenumber; % in the lower fluid
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
                E = [exp(1i*(Wavenumber*x1(i)+Layup{m,3}.*x3(j)-AngularFrequency*Time)) exp(1i*(Wavenumber*x1(i)+Layup{m,3}.*(Layup{m,6}-x3(j))-AngularFrequency*Time))];
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
                    [~,q] = min(abs(SLamb{p}(:,1)-Frequency));
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
                    [~,q] = min(abs(ALamb{p}(:,1)-Frequency));
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
                    [~,q] = min(abs(BLamb{p}(:,1)-Frequency));
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
                    [~,q] = min(abs(SScholte{p}(:,1)-Frequency));
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
                    [~,q] = min(abs(AScholte{p}(:,1)-Frequency));
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
                    [~,q] = min(abs(BScholte{p}(:,1)-Frequency));
                    Frequency = BScholte{p}(q,1);
                end
            end
            PhaseVelocity = BScholte{p}(q,4)*1e3;
            Attenuation = BScholte{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        end
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
        AngularFrequency = 2*pi*Frequency*1e3;
        AngularFrequency2 = AngularFrequency^2;
        Wavenumber = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        Wavenumber2 = Wavenumber^2;
        Wavenumber4 = Wavenumber2^2;
        for m = 1:SuperLayerSize
            rw2 = Material{m}.Density*AngularFrequency2;
            A1(m) = 2*c{m}(3,3)*c{m}(5,5);
            A2 = (c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2)*Wavenumber2-(c{m}(3,3)+c{m}(5,5))*rw2;
            A3 = c{m}(1,1)*c{m}(5,5)*Wavenumber4-(c{m}(1,1)+c{m}(5,5))*rw2*Wavenumber2+rw2^2;
            k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
            k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
            W = (rw2-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber*k3);
            D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*k3.*W); % sigma33
            D5 = 1i*(c{m}(5,5)*(k3+Wavenumber*W)); % sigma13
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
        for m = 2:log2(Repetitions)+1
            N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
        end
        for m = m+1:length(Pattern)
            N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
        end
        if  SymmetricSystem
            N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
            MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
        end
        if  FluidLoading
            G = inv(MM{end});
            k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-Wavenumber2);
            k3LowerFluid = sqrt(AngularFrequency2/LowerFluid.Velocity^2-Wavenumber2);
            if  contains(Mode,'Scholte') && Attenuation ~= 0
                if  PhaseVelocity < UpperFluid.Velocity
                    k3UpperFluid = -k3UpperFluid;
                end
                if  PhaseVelocity < LowerFluid.Velocity
                    k3LowerFluid = -k3LowerFluid;
                end
            end
            WUpperFluid = k3UpperFluid/Wavenumber;
            WLowerFluid = k3LowerFluid/Wavenumber;
            DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2/Wavenumber; % sigma11, sigma22, sigma33 in the upper fluid
            DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2/Wavenumber; % in the lower fluid
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
                    E = [exp(1i*(Wavenumber*x1(i)+Layup{m,3}.*x3(j)-AngularFrequency*Time)) exp(1i*(Wavenumber*x1(i)+Layup{m,3}.*(Layup{m,6}-x3(j))-AngularFrequency*Time))];
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
                    [~,q] = min(abs(SShear{p}(:,1)-Frequency));
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
                    [~,q] = min(abs(AShear{p-1}(:,1)-Frequency));
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
                    [~,q] = min(abs(BShear{p}(:,1)-Frequency));
                    Frequency = BShear{p}(q,1);
                end
            end
            PhaseVelocity = BShear{p}(q,4)*1e3;
            Attenuation = BShear{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        end
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
        AngularFrequency = 2*pi*Frequency*1e3;
        AngularFrequency2 = AngularFrequency^2;
        Wavenumber = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        for m = 1:SuperLayerSize
            k3 = sqrt((Material{m}.Density*AngularFrequency2-Wavenumber^2*c{m}(6,6))/c{m}(4,4));
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
        for m = 2:log2(Repetitions)+1
            N = inv(MM{m-1}(1,1)-MM{m-1}(2,2));
            MM{m} = [MM{m-1}(1,1)+MM{m-1}(1,2)*N*MM{m-1}(2,1) -MM{m-1}(1,2)*N*MM{m-1}(1,2);MM{m-1}(2,1)*N*MM{m-1}(2,1) MM{m-1}(2,2)-MM{m-1}(2,1)*N*MM{m-1}(1,2)];
        end
        for m = m+1:length(Pattern)
            N = inv(MM{Pattern(m)}(1,1)-MM{m-1}(2,2));
            MM{m} = [MM{m-1}(1,1)+MM{m-1}(1,2)*N*MM{m-1}(2,1) -MM{m-1}(1,2)*N*MM{Pattern(m)}(1,2);MM{Pattern(m)}(2,1)*N*MM{m-1}(2,1) MM{Pattern(m)}(2,2)-MM{Pattern(m)}(2,1)*N*MM{Pattern(m)}(1,2)];
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
                    E = [exp(1i*(Wavenumber*x1(i)+Layup{m,3}*x3(j)-AngularFrequency*Time)) exp(1i*(Wavenumber*x1(i)+Layup{m,3}*(Layup{m,6}-x3(j))-AngularFrequency*Time))];
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
                    EUpperFluid = exp(1i*(Wavenumber*x1(i)+k3UpperFluid*x3Fluid(j)-AngularFrequency*Time));
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
                    ELowerFluid = exp(1i*(Wavenumber*x1(i)+k3LowerFluid*x3Fluid(j)-AngularFrequency*Time));
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