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
function ModeShapeGrid_Anisotropic(Hybrid,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,LayupString,Layers,Material,c,PNGresolution,A,ALamb,AShear,AScholte,B,BLamb,BShear,BScholte,S,SLamb,SShear,SScholte,Directory,Export,FontSizeHeadLine,FontSizeAxesLabels,Frequency,GridLine,HeadLine,Length,LineWidth,Mode,FileName,PDF,PNG,PlateThickness,PropagationAngle,SamplesX1,SamplesX3,Gain,I,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Undistorted,Decoupled,ShowHalfSpaces,HalfSpaces)
%#ok<*AGROW>
%#ok<*MINV>
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
        x3 = (0:Layup{m,6}/SamplesX3:Layup{m,6})';
        if  m == 1
            x3Total = x3;
        else
            x3Total = vertcat(x3Total,x3+x3Total(end));
        end
        for i = 1:length(x1) % displacements u1,u3 at locations x1 and x3
            E = [exp(1i*(Wavenumber*x1(i)+Layup{m,3}.*x3)) exp(1i*(Wavenumber*x1(i)+Layup{m,3}.*(Layup{m,6}-x3)))];
            u{i,1}((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)+U(5)*E(:,5)+U(6)*E(:,6);
            u{i,1}((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)+Layup{m,5}(3)*U(3)*E(:,3)-Layup{m,5}(1)*U(4)*E(:,4)-Layup{m,5}(2)*U(5)*E(:,5)-Layup{m,5}(3)*U(6)*E(:,6);
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
            x3 = (0:Layup{m,6}/SamplesX3:Layup{m,6})';
            if  m == 1
                x3Total = x3;
            else
                x3Total = vertcat(x3Total,x3+x3Total(end));
            end
            for i = 1:length(x1) % displacements u1,u3 at locations x1 and x3
                E = [exp(1i*(Wavenumber*x1(i)+Layup{m,3}.*x3)) exp(1i*(Wavenumber*x1(i)+Layup{m,3}.*(Layup{m,6}-x3)))];
                u{i,1}((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4);
                u{i,1}((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)-Layup{m,5}(1)*U(3)*E(:,3)-Layup{m,5}(2)*U(4)*E(:,4);
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
            x3 = (0:Layup{m,6}/SamplesX3:Layup{m,6})';
            if  m == 1
                x3Total = x3;
            else
                x3Total = vertcat(x3Total,x3+x3Total(end));
            end
            for i = 1:length(x1) % displacements u2,u3 at locations x2 and x3
                E = [exp(1i*(Wavenumber*x1(i)+Layup{m,3}*x3)) exp(1i*(Wavenumber*x1(i)+Layup{m,3}*(Layup{m,6}-x3)))];
                u{i,1}((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = U(1)*E(:,1)+U(2)*E(:,2);
                u{i,1}((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = 0;
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
    for i = 1:size(u,1)
        u{i}(Ind,:) = [];
    end
end
if  FluidLoading && ShowHalfSpaces
    x3Fluid = (0:PlateThickness/SamplesX3/Layers:HalfSpaces*PlateThickness)';
    if  ToggleUpperFluid
        if  ~contains(Mode,'SH')
            for i = 1:length(x1)
                EUpperFluid = exp(1i*(Wavenumber*x1(i)+k3UpperFluid*x3Fluid));
                uUpperFluid{i,1}(:,1) = EUpperFluid;
                uUpperFluid{i,1}(:,2) = -WUpperFluid*EUpperFluid;
            end
        elseif contains(Mode,'SH')
            for i = 1:length(x1)
                uUpperFluid{i,1}(length(x3Fluid),2) = 0;
            end
        end
        x3Total = vertcat(-flipud(x3Fluid),x3Total);
        for i = 1:length(x1)
            u{i,:} = vertcat(flipud(uUpperFluid{i,:}),u{i,:});
        end
    end
    if  ToggleLowerFluid
        if  ~contains(Mode,'SH')
            for i = 1:length(x1)
                ELowerFluid = exp(1i*(Wavenumber*x1(i)+k3LowerFluid*x3Fluid));
                uLowerFluid{i,1}(:,1) = ULowerFluid*ELowerFluid;
                uLowerFluid{i,1}(:,2) = WLowerFluid*ULowerFluid*ELowerFluid;
            end
        elseif contains(Mode,'SH')
            for i = 1:length(x1)
                uLowerFluid{i,1}(length(x3Fluid),2) = 0;
            end               
        end
        x3Total = vertcat(x3Total,x3Total(end)+x3Fluid);
        for i = 1:length(x1)
            u{i,:} = vertcat(u{i,:},uLowerFluid{i,:});
        end
    end  
end
[X1,X3] = meshgrid(x1,x3Total);
Ratio = (abs(x3Total(1))+abs(x3Total(end)))/x1(end);
u1Max = max(abs(real(u{1}(:,1))));
u3Max = max(abs(real(u{1}(:,2))));
if  FluidLoading && ShowHalfSpaces
    k = 1+2*HalfSpaces;
else
    k = 1;
end
if  u1Max > u3Max
    Compensation = Gain*x1(end)/80/k/u1Max;
else
    Compensation = Gain*(abs(x3Total(1))+abs(x3Total(end)))/40/k/u3Max;
end
for i = 1:size(X1,2)
    X1Distorted(:,i) = X1(:,i)+real(u{i}(:,1))*Compensation;
    X3Distorted(:,i) = X3(:,i)+real(u{i}(:,2))*Compensation*Ratio;
end
f = figure('Name','Mode shape','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
if  ~FluidLoading || ~ShowHalfSpaces
    if  Undistorted
        line(X1(:,1:GridLine:end),X3(:,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(X1(1:GridLine:end,:)',X3(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
    end
    line(X1Distorted(:,1:GridLine:end),X3Distorted(:,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
    line(X1Distorted(1:GridLine:end,:)',X3Distorted(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k')
elseif FluidLoading && ShowHalfSpaces
    n = HalfSpaces*SamplesX3*Layers;
    if  Undistorted
        line(X1(1:n+1,1:GridLine:end),X3(1:n+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(X1(1:GridLine:n+1,:)',X3(1:GridLine:n+1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        line(X1(end-n:end,1:GridLine:end),X3(end-n:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(X1(end-n:GridLine:end,:)',X3(end-n:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        line(X1(n+2:end-n-1,1:GridLine:end),X3(n+2:end-n-1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(X1(n+2:GridLine:end-n-1,:)',X3(n+2:GridLine:end-n-1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
    end
    if  ToggleUpperFluid && ToggleLowerFluid
        if  strcmp(UpperFluid.Name,LowerFluid.Name)
            UpperFluidColor = [0 .5 1];
            LowerFluidColor = [0 .5 1];
        else
            if  UpperFluid.Density*UpperFluid.Velocity > LowerFluid.Density*LowerFluid.Velocity
                UpperFluidColor = [0 .5 1];
                LowerFluidColor = [0 .7 1];
            else
                UpperFluidColor = [0 .7 1];
                LowerFluidColor = [0 .5 1];
            end
        end
        line(X1Distorted(1:n+1,1:GridLine:end),X3Distorted(1:n+1,1:GridLine:end),'LineWidth',LineWidth,'Color',UpperFluidColor)
        line(X1Distorted(1:GridLine:n+1,:)',X3Distorted(1:GridLine:n+1,:)','LineWidth',LineWidth,'Color',UpperFluidColor)
        line(X1Distorted(end-n:end,1:GridLine:end),X3Distorted(end-n:end,1:GridLine:end),'LineWidth',LineWidth,'Color',LowerFluidColor)
        line(X1Distorted(end-n:GridLine:end,:)',X3Distorted(end-n:GridLine:end,:)','LineWidth',LineWidth,'Color',LowerFluidColor)
        line(X1Distorted(n+2:end-n-1,1:GridLine:end),X3Distorted(n+2:end-n-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X1Distorted(n+2:GridLine:end-n-1,:)',X3Distorted(n+2:GridLine:end-n-1,:)','LineWidth',LineWidth,'Color','k')
    elseif ToggleUpperFluid && ~ToggleLowerFluid
        line(X1Distorted(1:n+1,1:GridLine:end),X3Distorted(1:n+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1])
        line(X1Distorted(1:GridLine:n+1,:)',X3Distorted(1:GridLine:n+1,:)','LineWidth',LineWidth,'Color',[0 .5 1])
        line(X1Distorted(n+2:end,1:GridLine:end),X3Distorted(n+2:end,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X1Distorted(n+2:GridLine:end,:)',X3Distorted(n+2:GridLine:end,:)','LineWidth',LineWidth,'Color','k')
    elseif ~ToggleUpperFluid && ToggleLowerFluid
        line(X1Distorted(end-n:end,1:GridLine:end),X3Distorted(end-n:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1])
        line(X1Distorted(end-n:GridLine:end,:)',X3Distorted(end-n:GridLine:end,:)','LineWidth',LineWidth,'Color',[0 .5 1])
        line(X1Distorted(1:end-n-1,1:GridLine:end),X3Distorted(1:end-n-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X1Distorted(1:GridLine:end-n-1,:)',X3Distorted(1:GridLine:end-n-1,:)','LineWidth',LineWidth,'Color','k')
    end
end
ax = gca;
axis off
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
        String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_')];
    elseif ~SymmetricSystem && HeadLine == 2
        if  Repetitions == 1
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
        else
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
        end
    elseif SymmetricSystem && HeadLine == 2
        if  Repetitions == 1
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
        else
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
        end
    end
    if  FluidLoading
        if  ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/',replace(LowerFluid.Name,'_','\_'));
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/vacuum');
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in vacuum/',replace(LowerFluid.Name,'_','\_'));
        end
    end
    ax.Title.String = String;
end
if  ~contains(Mode,'SH')
    text(.5-.155*FontSizeAxesLabels/30,.05,'Propagation direction ($x_1)$','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels)
    text(-.02,.5-.14*FontSizeAxesLabels/30,'Thickness ($x_3$)','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels,'rotation',90);
else
    ax.Position = [0.13,0.21,0.775,0.615]; % default [.13 .11 .775 .815]
    text(.5-.125*FontSizeAxesLabels/30,-.073,'Shear horizontal ($x_2$)','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels)
    text(-.02,.5-.18*FontSizeAxesLabels/30,'Thickness ($x_3$)','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels,'rotation',90);
end
ax.XLim = [-.075*x1(end) x1(end)+.075*x1(end)];
ax.YDir = 'reverse';
if  Export
    try
        if  PDF
            exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
        end
        if  PNG
            exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
        return
    end
end
tb = axtoolbar('default');
tb.Visible = 'on';
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @Cursor;
function output_txt = Cursor(~,event_obj) %#ok<*INUSD>
    output_txt = {[]};
end     
end