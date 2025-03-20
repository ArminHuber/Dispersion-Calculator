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
function [u,epsilon,sigma,StrainEnergyDensity,KineticEnergyDensity,TotalEnergyDensity,PowerFlowDensity,uPhase,epsilonPhase,sigmaPhase,x3Total,p] = ModeShapeLinesComputer_Anisotropic(FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,ALamb,AShear,AScholte,BLamb,BShear,BScholte,SLamb,SShear,SScholte,c,Delta,Material,Layers,Frequency,Mode,PlateThickness,SamplesX3,I,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,ShowHalfSpaces,HalfSpaces,Phase)
%#ok<*AGROW>
uPhase=0;epsilonPhase=0;sigmaPhase=0;
p = str2double(regexp(Mode,'\d*','match'))+1;
if  ~contains(Mode,'Scholte') && ~contains(Mode,'SH') && Mode(1) == 'S'
    if  Frequency > ceil(max(SLamb{p}(:,1))) || Frequency < min(SLamb{p}(:,1))
        errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SLamb{p}(:,1))),' and ',num2str(ceil(max(SLamb{p}(:,1)))),' kHz.'],'Error');
        return
    else
        [~,q] = min(abs(SLamb{p}(:,1)-Frequency));
        Frequency = SLamb{p}(q,1);
    end
    PhaseVelocity = SLamb{p}(q,4)*1e3;
    Attenuation = SLamb{p}(q,7)*PhaseVelocity/Frequency/1e3; % Np/wavelength
elseif ~contains(Mode,'Scholte') && ~contains(Mode,'SH') && Mode(1) == 'A'
    if  Frequency > ceil(max(ALamb{p}(:,1))) || Frequency < min(ALamb{p}(:,1))
        errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(ALamb{p}(:,1))),' and ',num2str(ceil(max(ALamb{p}(:,1)))),' kHz.'],'Error');
        return
    else
        [~,q] = min(abs(ALamb{p}(:,1)-Frequency));
        Frequency = ALamb{p}(q,1);
    end
    PhaseVelocity = ALamb{p}(q,4)*1e3;
    Attenuation = ALamb{p}(q,7)*PhaseVelocity/Frequency/1e3; % Np/wavelength
elseif ~contains(Mode,'Scholte') && ~contains(Mode,'SH') && Mode(1) == 'B'
    if  Frequency > ceil(max(BLamb{p}(:,1))) || Frequency < min(BLamb{p}(:,1))
        errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(BLamb{p}(:,1))),' and ',num2str(ceil(max(BLamb{p}(:,1)))),' kHz.'],'Error');
        return
    else
        [~,q] = min(abs(BLamb{p}(:,1)-Frequency));
        Frequency = BLamb{p}(q,1);
    end
    PhaseVelocity = BLamb{p}(q,4)*1e3;
    Attenuation = BLamb{p}(q,7)*PhaseVelocity/Frequency/1e3; % Np/wavelength
elseif contains(Mode,'Scholte') && Mode(1) == 'S' 
    if  Frequency > ceil(max(SScholte{p}(:,1))) || Frequency < min(SScholte{p}(:,1))
        errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SScholte{p}(:,1))),' and ',num2str(ceil(max(SScholte{p}(:,1)))),' kHz.'],'Error');
        return
    else
        [~,q] = min(abs(SScholte{p}(:,1)-Frequency));
        Frequency = SScholte{p}(q,1);
    end
    PhaseVelocity = SScholte{p}(q,4)*1e3;
    Attenuation = SScholte{p}(q,7)*PhaseVelocity/Frequency/1e3; % Np/wavelength
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
    Attenuation = AScholte{p}(q,7)*PhaseVelocity/Frequency/1e3; % Np/wavelength
elseif contains(Mode,'Scholte') && Mode(1) == 'B'
    if  Frequency > ceil(max(BScholte{p}(:,1))) || Frequency < min(BScholte{p}(:,1))
        errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(BScholte{p}(:,1))),' and ',num2str(ceil(max(BScholte{p}(:,1)))),' kHz.'],'Error');
        return
    else
        [~,q] = min(abs(BScholte{p}(:,1)-Frequency));
        Frequency = BScholte{p}(q,1);
    end
    PhaseVelocity = BScholte{p}(q,4)*1e3;
    Attenuation = BScholte{p}(q,7)*PhaseVelocity/Frequency/1e3; % Np/wavelength
elseif  contains(Mode,'SH') && Mode(1) == 'S'
    if  Frequency > ceil(max(SShear{p}(:,1))) || Frequency < min(SShear{p}(:,1))
        errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SShear{p}(:,1))),' and ',num2str(ceil(max(SShear{p}(:,1)))),' kHz.'],'Error');
        return
    else
        [~,q] = min(abs(SShear{p}(:,1)-Frequency));
        Frequency = SShear{p}(q,1);
    end
    PhaseVelocity = SShear{p}(q,4)*1e3;
    Attenuation = SShear{p}(q,7)*PhaseVelocity/Frequency/1e3; % Np/wavelength
elseif contains(Mode,'SH') && Mode(1) == 'A'
    if  Frequency > ceil(max(AShear{p-1}(:,1))) || Frequency < min(AShear{p-1}(:,1))
        errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(AShear{p-1}(:,1))),' and ',num2str(ceil(max(AShear{p-1}(:,1)))),' kHz.'],'Error');
        return
    else
        [~,q] = min(abs(AShear{p-1}(:,1)-Frequency));
        Frequency = AShear{p-1}(q,1);
    end
    PhaseVelocity = AShear{p-1}(q,4)*1e3;
    Attenuation = AShear{p-1}(q,7)*PhaseVelocity/Frequency/1e3; % Np/wavelength
elseif contains(Mode,'SH') && Mode(1) == 'B'
    if  Frequency > ceil(max(BShear{p}(:,1))) || Frequency < min(BShear{p}(:,1))
        errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(BShear{p}(:,1))),' and ',num2str(ceil(max(BShear{p}(:,1)))),' kHz.'],'Error');
        return
    else
        [~,q] = min(abs(BShear{p}(:,1)-Frequency));
        Frequency = BShear{p}(q,1);
    end
    PhaseVelocity = BShear{p}(q,4)*1e3;
    Attenuation = BShear{p}(q,7)*PhaseVelocity/Frequency/1e3; % Np/wavelength
end
if  ~ToggleUpperFluid
    UpperFluid.Velocity = 1e-10;
    UpperFluid.Density = 1e-10;
end
if  ~ToggleLowerFluid
    LowerFluid.Velocity = 1e-10;
    LowerFluid.Density = 1e-10;
end
AngularFrequency = Frequency*pi*2e3;
AngularFrequency2 = AngularFrequency^2;
Wavenumber = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/2/pi);
Wavenumber2 = Wavenumber^2;
if  ~Decoupled
    Wavenumber4 = Wavenumber2^2;
    Wavenumber6 = Wavenumber2^3;
    for m = 1:SuperLayerSize
        rw2 = Material{m}.Density*AngularFrequency2;
        r2w4 = rw2^2;
        A1 = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m)*Wavenumber2+(c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m)*rw2;
        A2 = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m)*Wavenumber4+(c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m)*rw2*Wavenumber2+(c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m)*r2w4;
        A3 = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m)*Wavenumber6+(c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m)*rw2*Wavenumber4+(c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m)*r2w4*Wavenumber2-rw2^3/Delta(m);
        d1 = A1/3;
        d2 = A2/3-d1^2;
        d3 = d1^3-d1*A2/2+A3/2;
        d4 = (sqrt(d2^3+d3^2)-d3)^(1/3);
        d5 = d2/d4;
        d6 = (d5-d4)/2-d1;
        d7 = (d5+d4)/2i*sqrt(3);
        k3(1) = sqrt(d6+d7);
        k3(2) = sqrt(d6-d7);
        k3(3) = sqrt(d4-d5-d1);
        k32 = k3.^2;
        k3k = k3*Wavenumber;
        m11 = c{m}(1,1)*Wavenumber2+c{m}(5,5)*k32-rw2;
        m22 = c{m}(6,6)*Wavenumber2+c{m}(4,4)*k32-rw2;
        m12 = c{m}(1,6)*Wavenumber2+c{m}(4,5)*k32;
        m13 = (c{m}(1,3)+c{m}(5,5))*k3k;
        m23 = (c{m}(3,6)+c{m}(4,5))*k3k;
        m1 = m13.*m22-m12.*m23;
        V = (m11.*m23-m13.*m12)./m1;
        W = (m11.*m22-m12.^2)./-m1;
        e1 = Wavenumber*W+k3;
        e2 = k3.*V;
        e3 = k3.*W;
        D3 = 1i*((c{m}(1,3)+c{m}(3,6)*V)*Wavenumber+c{m}(3,3)*e3);
        D4 = 1i*(c{m}(4,5)*e1+c{m}(4,4)*e2);
        D5 = 1i*(c{m}(5,5)*e1+c{m}(4,5)*e2);
        E = exp(1i*k3*LayerThicknesses(m));
        L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
        A{m,3} = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W]; % L2
        A{m,1} = L1/A{m,3}; % L
        A{m,4} = LayerThicknesses(m);
        A{m,5} = Material{m}.Density;
        A{m,6} = k3;
        A{m,7} = V;
        A{m,8} = W;
        A{m,9} = 1i*e3; % epsilon33
        A{m,10} = 1i*e2; % epsilon23
        A{m,11} = 1i*e1; % epsilon13
        A{m,12} = 1i*Wavenumber*V; % epsilon12
        A{m,13} = 1i*((c{m}(1,1)+c{m}(1,6)*V)*Wavenumber+c{m}(1,3)*e3); % sigma11
        A{m,14} = 1i*((c{m}(1,2)+c{m}(2,6)*V)*Wavenumber+c{m}(2,3)*e3); % sigma22
        A{m,15} = D3; % sigma33
        A{m,16} = D4; % sigma23
        A{m,17} = D5; % sigma13
        A{m,18} = 1i*((c{m}(1,6)+c{m}(6,6)*V)*Wavenumber+c{m}(3,6)*e3); % sigma12
    end
    A = repmat(A,Repetitions,1);
    M{1} = A{1};
    for m = 2:SuperLayerSize
        M0 = A{m,1}(1:3,1:3)-M{1}(4:6,4:6);
        M1 = M{1}(1:3,4:6)/M0;
        M2 = A{m,1}(4:6,1:3)/M0;
        M{1} = [M{1}(1:3,1:3)+M1*M{1}(4:6,1:3) -M1*A{m,1}(1:3,4:6);M2*M{1}(4:6,1:3) A{m,1}(4:6,4:6)-M2*A{m,1}(1:3,4:6)];
    end
    for m = 1:length(Pattern)
        M0 = M{Pattern(m)}(1:3,1:3)-M{m}(4:6,4:6);
        M1 = M{m}(1:3,4:6)/M0;
        M2 = M{Pattern(m)}(4:6,1:3)/M0;
        M{m+1} = [M{m}(1:3,1:3)+M1*M{m}(4:6,1:3) -M1*M{Pattern(m)}(1:3,4:6);M2*M{m}(4:6,1:3) M{Pattern(m)}(4:6,4:6)-M2*M{Pattern(m)}(1:3,4:6)];
    end
    if  SymmetricSystem
        A = [A;flipud(A)];
        M0 = M{end}(4:6,4:6).*I-M{end}(4:6,4:6);
        M1 = M{end}(1:3,4:6)/M0;
        M2 = M{end}(1:3,4:6).*I/M0;
        M{end} = [M{end}(1:3,1:3,:)+M1*M{end}(4:6,1:3) -M1*(M{end}(4:6,1:3).*I);M2*M{end}(4:6,1:3,:) M{end}(1:3,1:3,:).*I-M2*(M{end}(4:6,1:3,:).*I)];
    end
    if  FluidLoading
        M{end} = inv(M{end});
        k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-Wavenumber2);
        k3LowerFluid = sqrt(AngularFrequency2/LowerFluid.Velocity^2-Wavenumber2);
        if  contains(Mode,'Scholte')
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
        ULowerFluid = -(WUpperFluid+M{end}(3,1)*DUpperFluid)/(M{end}(3,4)*DLowerFluid);
        uInterfaces = M{end}(1:3,1)*DUpperFluid+M{end}(1:3,4)*DLowerFluid*ULowerFluid;
        uInterfaces(:,Layers+1) = M{end}(4:6,1)*DUpperFluid+M{end}(4:6,4)*DLowerFluid*ULowerFluid;
    else
        Z1 = [-M{end}(1:3,1:3)*[1 1 1;A{1,7};-A{1,8}] -M{end}(1:3,4:6)*[1 1 1;A{end,7};A{end,8}];M{end}(4:6,1:3)*[1 1 1;A{1,7};-A{1,8}] M{end}(4:6,4:6)*[1 1 1;A{end,7};A{end,8}]];
        Z2 = [M{end}(1:3,1:3)*[1 1 1;A{1,7};A{1,8}];-M{end}(4:6,1:3)*[1 1 1;A{1,7};A{1,8}]];
        RT = Z1\Z2(:,1);
        uInterfaces = [1;A{1,7}(1);A{1,8}(1)]+[1 1 1;A{1,7};-A{1,8}]*RT(1:3);
        uInterfaces(:,Layers+1) = [1 1 1;A{end,7};A{end,8}]*RT(4:6);
    end
    A{1,2} = A{1};
    for m = 2:Layers
        M0 = A{m,1}(1:3,1:3)-A{m-1,2}(4:6,4:6);
        M1 = A{m-1,2}(1:3,4:6)/M0;
        M2 = A{m,1}(4:6,1:3)/M0;
        A{m,2} = [A{m-1,2}(1:3,1:3)+M1*A{m-1,2}(4:6,1:3) -M1*A{m,1}(1:3,4:6);M2*A{m-1,2}(4:6,1:3) A{m,1}(4:6,4:6)-M2*A{m,1}(1:3,4:6)];
    end
    for m = Layers:-1:2
        M0 = A{m,1}(1:3,1:3)-A{m-1,2}(4:6,4:6);
        uInterfaces(:,m) = M0\A{m-1,2}(4:6,1:3)*uInterfaces(:,1)-M0\A{m,1}(1:3,4:6)*uInterfaces(:,m+1);
    end
    for m = 1:Layers
        x3 = (0:A{m,4}/SamplesX3:A{m,4})';
        if  m == 1
            x3Total = x3;
        else
            x3Total = [x3Total;x3Total(end)+x3];
        end
        r = (m-1)*length(x3)+1:m*length(x3);
        E = exp(1i*[A{m,6}.*x3 A{m,6}.*(A{m,4}-x3)]);
        U = A{m,3}\[uInterfaces(:,m);uInterfaces(:,m+1)];
        u(r,1) = E*U; % u1
        u(r,2) = [A{m,7} A{m,7}].*E*U; % u2
        u(r,3) = [A{m,8} -A{m,8}].*E*U; % u3
        v(r,:) = -1i*AngularFrequency*u(r,:); % v1,v2,v3
        epsilon(r,1) = 1i*Wavenumber*u(r,1); % epsilon11
        epsilon(r,3) = [A{m,9} A{m,9}].*E*U; % epsilon33
        epsilon(r,4) = [A{m,10} -A{m,10}].*E*U; % epsilon23
        epsilon(r,5) = [A{m,11} -A{m,11}].*E*U; % epsilon13
        epsilon(r,6) = [A{m,12} A{m,12}].*E*U; % epsilon12
        sigma(r,1) = [A{m,13} A{m,13}].*E*U; % sigma11
        sigma(r,2) = [A{m,14} A{m,14}].*E*U; % sigma22
        sigma(r,3) = [A{m,15} A{m,15}].*E*U; % sigma33
        sigma(r,4) = [A{m,16} -A{m,16}].*E*U; % sigma23
        sigma(r,5) = [A{m,17} -A{m,17}].*E*U; % sigma13
        sigma(r,6) = [A{m,18} A{m,18}].*E*U; % sigma12
        KineticEnergyDensity(r,1) = A{m,5}*sum(abs(v(r,:)).^2,2)/2;
    end
    StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
    PowerFlowDensity(:,1) = -real(sigma(:,1).*conj(v(:,1))+sigma(:,6).*conj(v(:,2))+sigma(:,5).*conj(v(:,3)))/2;
    PowerFlowDensity(:,2) = -real(sigma(:,6).*conj(v(:,1))+sigma(:,2).*conj(v(:,2))+sigma(:,4).*conj(v(:,3)))/2;
    PowerFlowDensity(:,3) = -real(sigma(:,5).*conj(v(:,1))+sigma(:,4).*conj(v(:,2))+sigma(:,3).*conj(v(:,3)))/2;
    if  Phase
        uPhase2 = rad2deg(angle(u(:,2)*exp(-1i*angle(u(1,1)))));
        epsilonPhase456 = rad2deg(angle(epsilon(:,4:6)*exp(-1i*angle(u(1,1)))));
        sigmaPhase456 = rad2deg(angle(sigma(:,4:6)*exp(-1i*angle(u(1,1)))));
    end
    u(:,2) = u(:,2)*exp(-1i*angle(u(2,2))); % zero in the fluid, so we make phase shift now
    epsilon(:,4) = epsilon(:,4)*exp(-1i*angle(epsilon(2,4)));
    epsilon(:,5) = epsilon(:,5)*exp(-1i*angle(epsilon(2,5)));
    epsilon(:,6) = epsilon(:,6)*exp(-1i*angle(epsilon(2,6)));
    sigma(:,4) = sigma(:,4)*exp(-1i*angle(sigma(2,4)));
    sigma(:,5) = sigma(:,5)*exp(-1i*angle(sigma(2,5)));
    sigma(:,6) = sigma(:,6)*exp(-1i*angle(sigma(2,6)));
elseif Decoupled && ~contains(Mode,'SH')
    Wavenumber4 = Wavenumber2^2;
    for m = 1:SuperLayerSize
        rw2 = Material{m}.Density*AngularFrequency2;
        A1 = 2*c{m}(3,3)*c{m}(5,5);
        A2 = (c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2)*Wavenumber2-(c{m}(3,3)+c{m}(5,5))*rw2;
        A3 = c{m}(1,1)*c{m}(5,5)*Wavenumber4-(c{m}(1,1)+c{m}(5,5))*rw2*Wavenumber2+rw2^2;
        d1 = sqrt(A2^2-2*A1*A3);
        k3(1) = sqrt((-A2+d1)/A1);
        k3(2) = sqrt((-A2-d1)/A1);
        W = (rw2-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber*k3);
        e1 = Wavenumber*W+k3;
        e3 = k3.*W;
        D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*e3);
        D5 = 1i*c{m}(5,5)*e1;
        E = exp(1i*k3*LayerThicknesses(m));
        L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
        A{m,3} = [1 1 E;W -W.*E;E 1 1;W.*E -W]; % L2
        A{m,1} = L1/A{m,3}; % L
        A{m,4} = LayerThicknesses(m);
        A{m,5} = Material{m}.Density;
        A{m,6} = k3;
        A{m,7} = W;
        A{m,8} = 1i*e3; % epsilon33
        A{m,9} = 1i*e1; % epsilon13
        A{m,10} = 1i*(c{m}(1,1)*Wavenumber+c{m}(1,3)*e3); % sigma11
        A{m,11} = 1i*(c{m}(1,2)*Wavenumber+c{m}(2,3)*e3); % sigma22
        A{m,12} = D3; % sigma33
        A{m,13} = D5; % sigma13
    end
    A = repmat(A,Repetitions,1);
    M{1} = A{1};
    for m = 2:SuperLayerSize
        M0 = A{m,1}(1:2,1:2)-M{1}(3:4,3:4);
        M1 = M{1}(1:2,3:4)/M0;
        M2 = A{m,1}(3:4,1:2)/M0;
        M{1} = [M{1}(1:2,1:2)+M1*M{1}(3:4,1:2) -M1*A{m,1}(1:2,3:4);M2*M{1}(3:4,1:2) A{m,1}(3:4,3:4)-M2*A{m,1}(1:2,3:4)];
    end
    for m = 1:length(Pattern)
        M0 = M{Pattern(m)}(1:2,1:2)-M{m}(3:4,3:4);
        M1 = M{m}(1:2,3:4)/M0;
        M2 = M{Pattern(m)}(3:4,1:2)/M0;
        M{m+1} = [M{m}(1:2,1:2)+M1*M{m}(3:4,1:2) -M1*M{Pattern(m)}(1:2,3:4);M2*M{m}(3:4,1:2) M{Pattern(m)}(3:4,3:4)-M2*M{Pattern(m)}(1:2,3:4)];
    end
    if  SymmetricSystem
        A = [A;flipud(A)];
        M0 = M{end}(3:4,3:4).*I1-M{end}(3:4,3:4);
        M1 = M{end}(1:2,3:4)/M0;
        M2 = M{end}(1:2,3:4).*I1/M0;
        M{end} = [M{end}(1:2,1:2,:)+M1*M{end}(3:4,1:2) -M1*(M{end}(3:4,1:2).*I1);M2*M{end}(3:4,1:2,:) M{end}(1:2,1:2,:).*I1-M2*(M{end}(3:4,1:2,:).*I1)];
    end
    if  FluidLoading
        M{end} = inv(M{end});
        k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-Wavenumber2);
        k3LowerFluid = sqrt(AngularFrequency2/LowerFluid.Velocity^2-Wavenumber2);
        if  contains(Mode,'Scholte')
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
        ULowerFluid = -(WUpperFluid+M{end}(2,1)*DUpperFluid)/(M{end}(2,3)*DLowerFluid);
        uInterfaces = M{end}(1:2,1)*DUpperFluid+M{end}(1:2,3)*DLowerFluid*ULowerFluid;
        uInterfaces(:,Layers+1) = M{end}(3:4,1)*DUpperFluid+M{end}(3:4,3)*DLowerFluid*ULowerFluid;
    else
        Z1 = [-M{end}(1:2,1:2)*[1 1;-A{1,7}] -M{end}(1:2,3:4)*[1 1;A{end,7}];M{end}(3:4,1:2)*[1 1;-A{1,7}] M{end}(3:4,3:4)*[1 1;A{end,7}]];
        Z2 = [M{end}(1:2,1:2)*[1 1;A{1,7}];-M{end}(3:4,1:2)*[1 1;A{1,7}]];
        RT = Z1\Z2(:,1);
        uInterfaces = [1;A{1,7}(1)]+[1 1;-A{1,7}]*RT(1:2);
        uInterfaces(:,Layers+1) = [1 1;A{end,7}]*RT(3:4);
    end
    A{1,2} = A{1};
    for m = 2:Layers
        M0 = A{m,1}(1:2,1:2)-A{m-1,2}(3:4,3:4);
        M1 = A{m-1,2}(1:2,3:4)/M0;
        M2 = A{m,1}(3:4,1:2)/M0;
        A{m,2} = [A{m-1,2}(1:2,1:2)+M1*A{m-1,2}(3:4,1:2) -M1*A{m,1}(1:2,3:4);M2*A{m-1,2}(3:4,1:2) A{m,1}(3:4,3:4)-M2*A{m,1}(1:2,3:4)];
    end
    for m = Layers:-1:2
        M0 = A{m,1}(1:2,1:2)-A{m-1,2}(3:4,3:4);
        uInterfaces(:,m) = M0\A{m-1,2}(3:4,1:2)*uInterfaces(:,1)-M0\A{m,1}(1:2,3:4)*uInterfaces(:,m+1);
    end
    for m = 1:Layers
        x3 = (0:A{m,4}/SamplesX3:A{m,4})';
        if  m == 1
            x3Total = x3;
        else
            x3Total = [x3Total;x3Total(end)+x3];
        end
        r = (m-1)*length(x3)+1:m*length(x3);
        E = exp(1i*[A{m,6}.*x3 A{m,6}.*(A{m,4}-x3)]);
        U = A{m,3}\[uInterfaces(:,m);uInterfaces(:,m+1)];
        u(r,1) = E*U; % u1
        u(r,3) = [A{m,7} -A{m,7}].*E*U; % u3
        v(r,[1 3]) = -1i*AngularFrequency*u(r,[1 3]); % v1,v3
        epsilon(r,1) = 1i*Wavenumber*u(r,1); % epsilon11
        epsilon(r,3) = [A{m,8} A{m,8}].*E*U; % epsilon33
        epsilon(r,5) = [A{m,9} -A{m,9}].*E*U; % epsilon13
        sigma(r,1) = [A{m,10} A{m,10}].*E*U; % sigma11
        sigma(r,2) = [A{m,11} A{m,11}].*E*U; % sigma22
        sigma(r,3) = [A{m,12} A{m,12}].*E*U; % sigma33
        sigma(r,5) = [A{m,13} -A{m,13}].*E*U; % sigma13
        KineticEnergyDensity(r,1) = A{m,5}*sum(abs(v(r,:)).^2,2)/2;
    end
    epsilon(:,6) = 0;
    sigma(:,6) = 0;
    StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
    PowerFlowDensity(:,1) = -real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)))/2;
    PowerFlowDensity(:,3) = -real(sigma(:,5).*conj(v(:,1))+sigma(:,3).*conj(v(:,3)))/2;
    if  Phase
        epsilonPhase5 = rad2deg(angle(epsilon(:,5)*exp(-1i*angle(u(1,1)))));
        sigmaPhase5 = rad2deg(angle(sigma(:,5)*exp(-1i*angle(u(1,1)))));
    end
    epsilon(:,5) = epsilon(:,5)*exp(-1i*angle(epsilon(2,5))); % zero in the fluid, so we make phase shift now
    sigma(:,5) = sigma(:,5)*exp(-1i*angle(sigma(2,5)));
elseif Decoupled && contains(Mode,'SH')
    for m = 1:SuperLayerSize
        k3 = sqrt((Material{m}.Density*AngularFrequency2-Wavenumber2*c{m}(6,6))/c{m}(4,4));
        if  k3 == 0
            k3 = 1e-10;
        end
        D4 = 1i*k3*c{m}(4,4);
        E = exp(1i*k3*LayerThicknesses(m));
        E2 = E^2;
        A{m,1} = D4/(E2-1)*[-1-E2 2*E;-2*E 1+E2]; % L
        A{m,3} = [1 E;E 1]; % L2
        A{m,4} = LayerThicknesses(m);
        A{m,5} = Material{m}.Density;
        A{m,6} = 1i*k3; % epsilon23
        A{m,7} = D4; % sigma23
        A{m,8} = 1i*Wavenumber*c{m}(6,6); % sigma12
    end
    A = repmat(A,Repetitions,1);
    M{1} = A{1};
    for m = 2:SuperLayerSize
        M0 = A{m,1}(1,1)-M{1}(2,2);
        M1 = M{1}(1,2)/M0;
        M2 = A{m,1}(2,1)/M0;
        M{1} = [M{1}(1,1)+M1*M{1}(2,1) -M1*A{m,1}(1,2);M2*M{1}(2,1) A{m,1}(2,2)-M2*A{m,1}(1,2)];
    end
    for m = 1:length(Pattern)
        M0 = M{Pattern(m)}(1,1)-M{m}(2,2);
        M1 = M{m}(1,2)/M0;
        M2 = M{Pattern(m)}(2,1)/M0;
        M{m+1} = [M{m}(1,1)+M1*M{m}(2,1) -M1*M{Pattern(m)}(1,2);M2*M{m}(2,1) M{Pattern(m)}(2,2)-M2*M{Pattern(m)}(1,2)];
    end
    if  SymmetricSystem
        A = [A;flipud(A)];
        M1 = -M{end}(1,2)/(2*M{end}(2,2));
        M{end} = [M{end}(1,1)+M1*M{end}(2,1) M1*M{end}(2,1);-M1*M{end}(2,1) -M{end}(1,1)-M1*M{end}(2,1)];
    end
    Z1 = [M{end}(1,1) -M{end}(1,2);M{end}(2,:)]; % changed sign in Z1(1,1)!
    Z2 = [M{end}(1,1);-M{end}(2,1)];
    RT = Z1\Z2;
    uInterfaces = 1+RT(1);
    uInterfaces(:,Layers+1) = RT(2);
    A{1,2} = A{1};
    for m = 2:Layers
        M0 = A{m,1}(1,1)-A{m-1,2}(2,2);
        M1 = A{m-1,2}(1,2)/M0;
        M2 = A{m,1}(2,1)/M0;
        A{m,2} = [A{m-1,2}(1,1)+M1*A{m-1,2}(2,1) -M1*A{m,1}(1,2);M2*A{m-1,2}(2,1) A{m,1}(2,2)-M2*A{m,1}(1,2)];
    end
    for m = Layers:-1:2
        M0 = A{m,1}(1,1)-A{m-1,2}(2,2);
        uInterfaces(:,m) = M0\A{m-1,2}(2,1)*uInterfaces(:,1)-M0\A{m,1}(1,2)*uInterfaces(:,m+1);
    end
    for m = 1:Layers
        x3 = (0:A{m,4}/SamplesX3:A{m,4})';
        if  m == 1
            x3Total = x3;
        else
            x3Total = [x3Total;x3Total(end)+x3];
        end
        r = (m-1)*length(x3)+1:m*length(x3);
        E = exp(A{m,6}*[x3 A{m,4}-x3]);
        U = A{m,3}\[uInterfaces(:,m);uInterfaces(:,m+1)];
        u(r,2) = E*U; % u2
        v(r,2) = -1i*AngularFrequency*u(r,2); % v2
        epsilon(r,4) = [A{m,6} -A{m,6}].*E*U; % epsilon23
        epsilon(r,6) = 1i*Wavenumber*u(r,2); % epsilon12
        sigma(r,4) = [A{m,7} -A{m,7}].*E*U; % sigma23
        sigma(r,6) = [A{m,8} A{m,8}].*E*U; % sigma12
        KineticEnergyDensity(r,1) = A{m,5}*abs(v(r,2)).^2/2;
    end
    u(:,3) = 0;
    StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
    PowerFlowDensity(:,1) = -real(sigma(:,6).*conj(v(:,2)))/2;
    PowerFlowDensity(:,3) = -real(sigma(:,4).*conj(v(:,2)))/2;
    if  Phase
        uPhase = rad2deg(angle(u));
        epsilonPhase = rad2deg(angle(epsilon));
        sigmaPhase = rad2deg(angle(sigma));
    end
    u(:,2) = u(:,2)*exp(-1i*angle(u(2,2)));
    epsilon(:,4) = epsilon(:,4)*exp(-1i*angle(epsilon(2,4)));
    epsilon(:,6) = epsilon(:,6)*exp(-1i*angle(epsilon(2,6)));
    sigma(:,4) = sigma(:,4)*exp(-1i*angle(sigma(2,4)));
    sigma(:,6) = sigma(:,6)*exp(-1i*angle(sigma(2,6)));
end
PowerFlow = trapz(x3Total,PowerFlowDensity(:,1));
% disp(['ce: ',num2str(2*PowerFlow/trapz(x3Total,StrainEnergyDensity+KineticEnergyDensity)),' m/s'])
if  FluidLoading && ShowHalfSpaces
    x3Fluid = (0:PlateThickness/SamplesX3:HalfSpaces*PlateThickness)';
    E1Fluid = 1i*Wavenumber; % epsilon11
    if  ToggleUpperFluid
        if  ~contains(Mode,'SH')
            EUpperFluid = exp(1i*k3UpperFluid*x3Fluid);
            E3UpperFluid = 1i*k3UpperFluid*WUpperFluid; % epsilon33
            uUpperFluid(:,1) = EUpperFluid;
            uUpperFluid(:,3) = -WUpperFluid*EUpperFluid;
            vUpperFluid = -1i*AngularFrequency*uUpperFluid;
            epsilonUpperFluid(:,1) = E1Fluid*EUpperFluid; % epsilon11
            epsilonUpperFluid(:,3) = E3UpperFluid*EUpperFluid; % epsilon33
            epsilonUpperFluid(:,6) = 0;
            sigmaUpperFluid(:,1) = DUpperFluid*EUpperFluid; % sigma11
            sigmaUpperFluid(:,2) = sigmaUpperFluid(:,1); % sigma22
            sigmaUpperFluid(:,3) = sigmaUpperFluid(:,1); % sigma33
            sigmaUpperFluid(:,6) = 0;
            StrainEnergyDensityUpperFluid = real(sum(epsilonUpperFluid.*conj(sigmaUpperFluid),2))/2;
            KineticEnergyDensityUpperFluid = UpperFluid.Density*sum(abs(vUpperFluid).^2,2)/2;
            PowerFlowDensityUpperFluid(:,1) = -real(sigmaUpperFluid(:,1).*conj(vUpperFluid(:,1)))/2;
            PowerFlowDensityUpperFluid(:,3) = -real(sigmaUpperFluid(:,3).*conj(vUpperFluid(:,3)))/2;
        elseif contains(Mode,'SH')
            uUpperFluid(length(x3Fluid),3) = 0;
            epsilonUpperFluid(length(x3Fluid),6) = 0;
            sigmaUpperFluid(length(x3Fluid),6) = 0;
            StrainEnergyDensityUpperFluid(length(x3Fluid),1) = 0;
            KineticEnergyDensityUpperFluid(length(x3Fluid),1) = 0;
            PowerFlowDensityUpperFluid(length(x3Fluid),3) = 0;
            if  Phase
                uPhase = [zeros(length(x3Fluid),3);uPhase];
                sigmaPhase = [zeros(length(x3Fluid),6);sigmaPhase];
                epsilonPhase = [zeros(length(x3Fluid),6);epsilonPhase];
            end
        end
        x3Total = [-flipud(x3Fluid);x3Total];
        u = [flipud(uUpperFluid);u];
        epsilon = [flipud(epsilonUpperFluid);epsilon];
        sigma = [flipud(sigmaUpperFluid);sigma];
        StrainEnergyDensity = [flipud(StrainEnergyDensityUpperFluid);StrainEnergyDensity];
        KineticEnergyDensity = [flipud(KineticEnergyDensityUpperFluid);KineticEnergyDensity];
        PowerFlowDensity = [flipud(PowerFlowDensityUpperFluid);PowerFlowDensity];
    end
    if  ToggleLowerFluid
        if  ~contains(Mode,'SH')
            ELowerFluid = exp(1i*k3LowerFluid*x3Fluid);
            E3LowerFluid = 1i*k3LowerFluid*WLowerFluid; % epsilon33
            uLowerFluid(:,1) = ULowerFluid*ELowerFluid;
            uLowerFluid(:,3) = WLowerFluid*ULowerFluid*ELowerFluid;
            vLowerFluid = -1i*AngularFrequency*uLowerFluid;
            epsilonLowerFluid(:,1) = E1Fluid*ULowerFluid*ELowerFluid; % epsilon11
            epsilonLowerFluid(:,3) = E3LowerFluid*ULowerFluid*ELowerFluid; % epsilon33
            epsilonLowerFluid(:,6) = 0;
            sigmaLowerFluid(:,1) = DLowerFluid*ULowerFluid*ELowerFluid; % sigma11
            sigmaLowerFluid(:,2) = sigmaLowerFluid(:,1); % sigma22
            sigmaLowerFluid(:,3) = sigmaLowerFluid(:,1); % sigma33
            sigmaLowerFluid(:,6) = 0;
            StrainEnergyDensityLowerFluid = real(sum(epsilonLowerFluid.*conj(sigmaLowerFluid),2))/2;
            KineticEnergyDensityLowerFluid = LowerFluid.Density*sum(abs(vLowerFluid).^2,2)/2;
            PowerFlowDensityLowerFluid(:,1) = -real(sigmaLowerFluid(:,1).*conj(vLowerFluid(:,1)))/2;
            PowerFlowDensityLowerFluid(:,3) = -real(sigmaLowerFluid(:,3).*conj(vLowerFluid(:,3)))/2;
        elseif contains(Mode,'SH')
            uLowerFluid(length(x3Fluid),3) = 0;
            epsilonLowerFluid(length(x3Fluid),6) = 0;
            sigmaLowerFluid(length(x3Fluid),6) = 0;
            StrainEnergyDensityLowerFluid(length(x3Fluid),1) = 0;
            KineticEnergyDensityLowerFluid(length(x3Fluid),1) = 0;
            PowerFlowDensityLowerFluid(length(x3Fluid),3) = 0;
            if  Phase
                uPhase = [uPhase;zeros(length(x3Fluid),3)];
                epsilonPhase = [epsilonPhase;zeros(length(x3Fluid),6)];
                sigmaPhase = [sigmaPhase;zeros(length(x3Fluid),6)];
            end
        end
        x3Total = [x3Total;x3Total(end)+x3Fluid];
        u = [u;uLowerFluid];
        epsilon = [epsilon;epsilonLowerFluid];
        sigma = [sigma;sigmaLowerFluid];
        StrainEnergyDensity = [StrainEnergyDensity;StrainEnergyDensityLowerFluid];
        KineticEnergyDensity = [KineticEnergyDensity;KineticEnergyDensityLowerFluid];
        PowerFlowDensity = [PowerFlowDensity;PowerFlowDensityLowerFluid];
    end
end
if  ~contains(Mode,'SH')
    if  Phase
        if  FluidLoading && ShowHalfSpaces
            if  ToggleUpperFluid
                uPhase = rad2deg(angle(u*exp(-1i*angle(u(length(x3Fluid)+1,1)))));
                epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(length(x3Fluid)+1,1)))));
                sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(length(x3Fluid)+1,1)))));
                if  ~Decoupled
                    uPhase(length(x3Fluid)+1:length(x3Fluid)+length(uPhase2),2) = uPhase2;
                    epsilonPhase(length(x3Fluid)+1:length(x3Fluid)+length(epsilonPhase456),4:6) = epsilonPhase456;
                    sigmaPhase(length(x3Fluid)+1:length(x3Fluid)+length(sigmaPhase456),4:6) = sigmaPhase456;
                else
                    epsilonPhase(length(x3Fluid)+1:length(x3Fluid)+length(epsilonPhase5),5) = epsilonPhase5;
                    sigmaPhase(length(x3Fluid)+1:length(x3Fluid)+length(sigmaPhase5),5) = sigmaPhase5;
                end
            else
                uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
                epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
                sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
                if  ~Decoupled
                    uPhase(1:length(uPhase2),2) = uPhase2;
                    epsilonPhase(1:length(epsilonPhase456),4:6) = epsilonPhase456;
                    sigmaPhase(1:length(sigmaPhase456),4:6) = sigmaPhase456;
                else
                    epsilonPhase(1:length(epsilonPhase5),5) = epsilonPhase5;
                    sigmaPhase(1:length(sigmaPhase5),5) = sigmaPhase5;
                end
            end
        else
            uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
            epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
            sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
            if  ~Decoupled
                uPhase(:,2) = uPhase2;
                epsilonPhase(:,4:6) = epsilonPhase456;
                sigmaPhase(:,4:6) = sigmaPhase456;
            else
                epsilonPhase(:,5) = epsilonPhase5;
                sigmaPhase(:,5) = sigmaPhase5;
            end
        end
    end
    u(:,1) = u(:,1)*exp(-1i*angle(u(2,1)));
    u(:,3) = u(:,3)*exp(-1i*angle(u(2,3)));
    epsilon(:,1) = epsilon(:,1)*exp(-1i*angle(epsilon(2,1)));
    epsilon(:,3) = epsilon(:,3)*exp(-1i*angle(epsilon(2,3)));
    sigma(:,1) = sigma(:,1)*exp(-1i*angle(sigma(2,1)));
    sigma(:,2) = sigma(:,2)*exp(-1i*angle(sigma(2,2)));
    sigma(:,3) = sigma(:,3)*exp(-1i*angle(sigma(2,3)));
end
u = u/sqrt(PowerFlow);
epsilon = epsilon/sqrt(PowerFlow);
sigma = sigma/sqrt(PowerFlow);
if  real(u(2,1)) < 0
    u = -u;
end
if  real(epsilon(2,1)) < 0
    epsilon = -epsilon;
end
if  real(sigma(2,1)) < 0
    sigma = -sigma;
end
if  Phase
    uPhase(round(uPhase) == -180) = 180;
    epsilonPhase(round(epsilonPhase) == -180) = 180;
    sigmaPhase(round(sigmaPhase) == -180) = 180;
end
StrainEnergyDensity = StrainEnergyDensity/PowerFlow/2;
KineticEnergyDensity = KineticEnergyDensity/PowerFlow/2;
TotalEnergyDensity = StrainEnergyDensity+KineticEnergyDensity;
PowerFlowDensity = PowerFlowDensity/PowerFlow;