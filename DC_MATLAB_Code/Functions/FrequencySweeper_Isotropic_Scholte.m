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
function [HSScholte,HAScholte,HBScholte] = FrequencySweeper_Isotropic_Scholte(Material,Half,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,SweepRange,Symmetric)
Resolution = 1e-5; % (kHz)
PhaseVelocityOffset = 1-1e-4; % (m/s)

%#ok<*AGROW>
HSScholte = [];
HAScholte = [];
HBScholte = [];
if  Resolution < abs(SweepRange(1)-SweepRange(2))
    Bisections = ceil(log2(Resolution/(abs(SweepRange(1)-SweepRange(2))))/log2(2*.25));
else
    Bisections = 1;
end
AngularFrequency = 2*pi*SweepRange*1e3;
if  Symmetric
    Fluid = UpperFluid;
    PhaseVelocity = PhaseVelocityOffset*Fluid.Velocity;
    k2 = (AngularFrequency/PhaseVelocity).^2;
    kT2 = (AngularFrequency/Material.TransverseVelocity).^2;
    x = sqrt((AngularFrequency/Material.LongitudinalVelocity).^2-k2);
    y = sqrt(kT2-k2);
    a1 = (y.^2-k2).^2./y;
    a2 = 4*k2.*x;
    a3 = tan(x*Half);
    a4 = tan(y*Half);
    F = Fluid.Density*kT2.^2.*x./(y*Material.Density.*sqrt((AngularFrequency/Fluid.Velocity).^2-k2));
    YSScholte = abs(a1./a3+a2./a4-1i*F);
    YAScholte = abs(a1.*a3+a2.*a4+1i*F);
    XRoughSScholte = [];
    XRoughAScholte = [];
    for i = 2:length(SweepRange)-1
        if  YSScholte(i) < YSScholte(i-1) && YSScholte(i) < YSScholte(i+1)
            XRoughSScholte(end+1) = SweepRange(i);
        end
        if  YAScholte(i) < YAScholte(i-1) && YAScholte(i) < YAScholte(i+1)
            XRoughAScholte(end+1)  = SweepRange(i);
        end
    end
    if  isempty(XRoughAScholte)
        XRoughSScholte = [];
    else
        XRoughSScholte(XRoughSScholte < XRoughAScholte(1)) = [];
    end
    if  isempty(XRoughSScholte) && isempty(XRoughAScholte)
        disp('No higher order Scholte modes found!')
        return
    end
    for j = 1:length(XRoughSScholte)
        Frequency = [XRoughSScholte(j)-(SweepRange(2)-SweepRange(1)) XRoughSScholte(j)+(SweepRange(2)-SweepRange(1))];
        for o = 1:Bisections
            Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
            AngularFrequency = 2*pi*Frequency*1e3;
            for i = 1:length(Frequency)
                k2 = (AngularFrequency(i)/PhaseVelocity)^2;
                kT2 = (AngularFrequency(i)/Material.TransverseVelocity)^2;
                x = sqrt((AngularFrequency(i)/Material.LongitudinalVelocity)^2-k2);
                y = sqrt(kT2-k2);
                Y(i) = abs((y^2-k2)^2/(y*tan(x*Half))+4*k2*x/tan(y*Half)-1i*Fluid.Density*kT2^2*x/(y*Material.Density*sqrt((AngularFrequency(i)/Fluid.Velocity)^2-k2)));
                if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                    if  o == Bisections
                        HSScholte(end+1) = Frequency(i-1);
                    end
                    Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                    break
                end
            end
        end
    end
    for j = 1:length(XRoughAScholte)
        Frequency = [XRoughAScholte(j)-(SweepRange(2)-SweepRange(1)) XRoughAScholte(j)+(SweepRange(2)-SweepRange(1))];
        for o = 1:Bisections
            Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
            AngularFrequency = 2*pi*Frequency*1e3;
            for i = 1:length(Frequency)
                k2 = (AngularFrequency(i)/PhaseVelocity)^2;
                kT2 = (AngularFrequency(i)/Material.TransverseVelocity)^2;
                x = sqrt((AngularFrequency(i)/Material.LongitudinalVelocity)^2-k2);
                y = sqrt(kT2-k2);
                Y(i) = abs((y^2-k2)^2/y*tan(x*Half)+4*k2*x*tan(y*Half)+1i*Fluid.Density*kT2^2*x/(y*Material.Density*sqrt((AngularFrequency(i)/Fluid.Velocity)^2-k2)));
                if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                    if  o == Bisections
                        HAScholte(end+1) = Frequency(i-1);
                    end
                    Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                    break
                end
            end
        end
    end
    String = ['Frq. @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode         Frq.(kHz)'];
    if  any(HAScholte)
        for i = 1:length(HAScholte)
            if  i < 10
                if  HAScholte(i) < 1e2
                    String = append(String,newline,'AScholte',num2str(i),'      ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e2 && HAScholte(i) < 1e3
                    String = append(String,newline,'AScholte',num2str(i),'     ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e3 && HAScholte(i) < 1e4
                    String = append(String,newline,'AScholte',num2str(i),'    ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e4
                    String = append(String,newline,'AScholte',num2str(i),'   ',num2str(HAScholte(i),'%.3f'));
                end
            else
                if  HAScholte(i) < 1e2
                    String = append(String,newline,'AScholte',num2str(i),'     ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e2 && HAScholte(i) < 1e3
                    String = append(String,newline,'AScholte',num2str(i),'    ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e3 && HAScholte(i) < 1e4
                    String = append(String,newline,'AScholte',num2str(i),'   ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e4
                    String = append(String,newline,'AScholte',num2str(i),'  ',num2str(HAScholte(i),'%.3f'));
                end
            end
        end
    end
    String = append(String,newline);
    if  any(HSScholte)
        for i = 1:length(HSScholte)
            if  i < 10
                if  HSScholte(i) < 1e2
                    String = append(String,newline,'SScholte',num2str(i),'      ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e2 && HSScholte(i) < 1e3
                    String = append(String,newline,'SScholte',num2str(i),'     ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e3 && HSScholte(i) < 1e4
                    String = append(String,newline,'SScholte',num2str(i),'    ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e4
                    String = append(String,newline,'SScholte',num2str(i),'   ',num2str(HSScholte(i),'%.3f'));
                end
            else
                if  HSScholte(i) < 1e2
                    String = append(String,newline,'SScholte',num2str(i),'     ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e2 && HSScholte(i) < 1e3
                    String = append(String,newline,'SScholte',num2str(i),'    ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e3 && HSScholte(i) < 1e4
                    String = append(String,newline,'SScholte',num2str(i),'   ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e4
                    String = append(String,newline,'SScholte',num2str(i),'  ',num2str(HSScholte(i),'%.3f'));
                end
            end
        end
    end
    String = append(String,newline,newline);
    if  any(HAScholte)
        String = append(String,'AScholte: ',num2str(length(HAScholte)));
    end
    String = append(String,newline);
    if  any(HSScholte)
        String = append(String,'SScholte: ',num2str(length(HSScholte)));
    end
else
    if  ToggleUpperFluid && ToggleLowerFluid
        if  UpperFluid.Velocity > LowerFluid.Velocity
            PhaseVelocity = PhaseVelocityOffset*LowerFluid.Velocity;
        else
            PhaseVelocity = PhaseVelocityOffset*UpperFluid.Velocity;
        end
    elseif ToggleUpperFluid && ~ToggleLowerFluid
        PhaseVelocity = PhaseVelocityOffset*UpperFluid.Velocity;
        Fluid = UpperFluid;
    elseif ~ToggleUpperFluid && ToggleLowerFluid
        PhaseVelocity = PhaseVelocityOffset*LowerFluid.Velocity;
        Fluid = LowerFluid;
    end
    XRoughScholte = [];
    k = AngularFrequency/PhaseVelocity;
    k2 = k.^2;
    kL2 = (AngularFrequency/Material.LongitudinalVelocity).^2;
    x = sqrt(kL2-k2);
    y = sqrt((AngularFrequency/Material.TransverseVelocity).^2-k2);
    xH = x*Half;
    yH = y*Half;
    SinxH = sin(xH);
    SinyH = sin(yH);
    CosxH = cos(xH);
    CosyH = cos(yH);
    Sin_xH = sin(-xH);
    Sin_yH = sin(-yH);
    Cos_xH = cos(-xH);
    Cos_yH = cos(-yH);
    a1 = -(Material.Lambda*kL2+2*Material.Mu*x.^2);
    a2 = Material.Mu*(k2-y.^2);
    a3 = 2i*Material.Mu*k.*x;
    a4 = 2i*Material.Mu*k.*y;
    a5 = -1i*k;
    if  ToggleUpperFluid && ToggleLowerFluid
        kUpperFluid2 = (AngularFrequency/UpperFluid.Velocity).^2;
        kLowerFluid2 = (AngularFrequency/LowerFluid.Velocity).^2;
        k3UpperFluid = sqrt(kUpperFluid2-k2);
        k3LowerFluid = sqrt(kLowerFluid2-k2);
        e = exp(1i*k3UpperFluid*Half);
        e_ = exp(-1i*k3LowerFluid*Half); 
        a6 = -UpperFluid.Velocity^2*UpperFluid.Density*kUpperFluid2;
        a7 = LowerFluid.Velocity^2*LowerFluid.Density*kLowerFluid2;
        a8 = 1i*k3UpperFluid;
        a9 = 1i*k3LowerFluid;
        for i = 1:length(SweepRange)
            M(1,1) = a1(i)*SinxH(i);
            M(1,2) = a1(i)*CosxH(i);
            M(1,3) = -a4(i)*CosyH(i);
            M(1,4) = a4(i)*SinyH(i);
            M(1,5) = a6(i)*e(i);
            M(2,1) = a3(i)*CosxH(i);
            M(2,2) = -a3(i)*SinxH(i);
            M(2,3) = a2(i)*SinyH(i);
            M(2,4) = a2(i)*CosyH(i);
            M(3,1) = x(i)*CosxH(i);
            M(3,2) = -x(i)*SinxH(i);
            M(3,3) = a5(i)*SinyH(i);
            M(3,4) = a5(i)*CosyH(i);
            M(3,5) = a8(i)*e(i);
            M(4,1) = a1(i)*Sin_xH(i);
            M(4,2) = a1(i)*Cos_xH(i);
            M(4,3) = -a4(i)*Cos_yH(i);
            M(4,4) = a4(i)*Sin_yH(i);
            M(4,6) = a7(i)*e_(i);
            M(5,1) = a3(i)*Cos_xH(i);
            M(5,2) = -a3(i)*Sin_xH(i);
            M(5,3) = a2(i)*Sin_yH(i);
            M(5,4) = a2(i)*Cos_yH(i);
            M(6,1) = x(i)*Cos_xH(i);
            M(6,2) = -x(i)*Sin_xH(i);
            M(6,3) = a5(i)*Sin_yH(i);
            M(6,4) = a5(i)*Cos_yH(i);
            M(6,6) = a9(i)*e_(i);
            Y(i) = abs(det(M));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                XRoughScholte(end+1) = SweepRange(i-1);
            end
        end
    else
        kFluid2 = (AngularFrequency/Fluid.Velocity).^2;
        k3Fluid = sqrt(kFluid2-k2);
        e = exp(1i*k3Fluid*Half);
        a6 = -Fluid.Velocity^2*Fluid.Density*kFluid2;
        a8 = 1i*k3Fluid;
        for i = 1:length(SweepRange)
            M(1,1) = a1(i)*SinxH(i);
            M(1,2) = a1(i)*CosxH(i);
            M(1,3) = -a4(i)*CosyH(i);
            M(1,4) = a4(i)*SinyH(i);
            M(1,5) = a6(i)*e(i);
            M(2,1) = a3(i)*CosxH(i);
            M(2,2) = -a3(i)*SinxH(i);
            M(2,3) = a2(i)*SinyH(i);
            M(2,4) = a2(i)*CosyH(i);
            M(3,1) = x(i)*CosxH(i);
            M(3,2) = -x(i)*SinxH(i);
            M(3,3) = a5(i)*SinyH(i);
            M(3,4) = a5(i)*CosyH(i);
            M(3,5) = a8(i)*e(i);
            M(4,1) = a1(i)*Sin_xH(i);
            M(4,2) = a1(i)*Cos_xH(i);
            M(4,3) = -a4(i)*Cos_yH(i);
            M(4,4) = a4(i)*Sin_yH(i);
            M(5,1) = a3(i)*Cos_xH(i);
            M(5,2) = -a3(i)*Sin_xH(i);
            M(5,3) = a2(i)*Sin_yH(i);
            M(5,4) = a2(i)*Cos_yH(i);
            Y(i) = abs(det(M));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                XRoughScholte(end+1) = SweepRange(i-1);
            end
        end
    end
    if  isempty(XRoughScholte)
        disp('No higher order Scholte modes found!')
        return
    end
    for j = 1:length(XRoughScholte)
        Frequency = [XRoughScholte(j)-(SweepRange(2)-SweepRange(1)) XRoughScholte(j)+(SweepRange(2)-SweepRange(1))];
        for o = 1:Bisections
            Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
            AngularFrequency = 2*pi*Frequency*1e3;
            for i = 1:length(Frequency)
                k = AngularFrequency(i)/PhaseVelocity;
                k2 = k^2;
                kL2 = (AngularFrequency(i)/Material.LongitudinalVelocity)^2;
                x = sqrt(kL2-k2);
                y = sqrt((AngularFrequency(i)/Material.TransverseVelocity)^2-k2);
                xH = x*Half;
                yH = y*Half;
                SinxH = sin(xH);
                SinyH = sin(yH);
                CosxH = cos(xH);
                CosyH = cos(yH);
                Sin_xH = sin(-xH);
                Sin_yH = sin(-yH);
                Cos_xH = cos(-xH);
                Cos_yH = cos(-yH);
                a1 = -(Material.Lambda*kL2+2*Material.Mu*x^2);
                a2 = Material.Mu*(k2-y^2);
                a3 = 2i*Material.Mu*k*x;
                a4 = 2i*Material.Mu*k*y;
                a5 = -1i*k;
                if  ToggleUpperFluid && ToggleLowerFluid
                    kUpperFluid2 = (AngularFrequency(i)/UpperFluid.Velocity)^2;
                    kLowerFluid2 = (AngularFrequency(i)/LowerFluid.Velocity)^2;
                    k3UpperFluid = sqrt(kUpperFluid2-k2);
                    k3LowerFluid = sqrt(kLowerFluid2-k2);
                    e = exp(1i*k3UpperFluid*Half);
                    e_ = exp(-1i*k3LowerFluid*Half); 
                    a6 = -UpperFluid.Velocity^2*UpperFluid.Density*kUpperFluid2;
                    a7 = LowerFluid.Velocity^2*LowerFluid.Density*kLowerFluid2;
                    a8 = 1i*k3UpperFluid;
                    a9 = 1i*k3LowerFluid;
                    M(1,1) = a1*SinxH;
                    M(1,2) = a1*CosxH;
                    M(1,3) = -a4*CosyH;
                    M(1,4) = a4*SinyH;
                    M(1,5) = a6*e;
                    M(2,1) = a3*CosxH;
                    M(2,2) = -a3*SinxH;
                    M(2,3) = a2*SinyH;
                    M(2,4) = a2*CosyH;
                    M(3,1) = x*CosxH;
                    M(3,2) = -x*SinxH;
                    M(3,3) = a5*SinyH;
                    M(3,4) = a5*CosyH;
                    M(3,5) = a8*e;
                    M(4,1) = a1*Sin_xH;
                    M(4,2) = a1*Cos_xH;
                    M(4,3) = -a4*Cos_yH;
                    M(4,4) = a4*Sin_yH;
                    M(4,6) = a7*e_;
                    M(5,1) = a3*Cos_xH;
                    M(5,2) = -a3*Sin_xH;
                    M(5,3) = a2*Sin_yH;
                    M(5,4) = a2*Cos_yH;
                    M(6,1) = x*Cos_xH;
                    M(6,2) = -x*Sin_xH;
                    M(6,3) = a5*Sin_yH;
                    M(6,4) = a5*Cos_yH;
                    M(6,6) = a9*e_;
                    Y(i) = abs(det(M));
                else
                    kFluid2 = (AngularFrequency(i)/Fluid.Velocity)^2;
                    k3Fluid = sqrt(kFluid2-k2);
                    e = exp(1i*k3Fluid*Half);
                    a6 = -Fluid.Velocity^2*Fluid.Density*kFluid2;
                    a8 = 1i*k3Fluid;
                    M(1,1) = a1*SinxH;
                    M(1,2) = a1*CosxH;
                    M(1,3) = -a4*CosyH;
                    M(1,4) = a4*SinyH;
                    M(1,5) = a6*e;
                    M(2,1) = a3*CosxH;
                    M(2,2) = -a3*SinxH;
                    M(2,3) = a2*SinyH;
                    M(2,4) = a2*CosyH;
                    M(3,1) = x*CosxH;
                    M(3,2) = -x*SinxH;
                    M(3,3) = a5*SinyH;
                    M(3,4) = a5*CosyH;
                    M(3,5) = a8*e;
                    M(4,1) = a1*Sin_xH;
                    M(4,2) = a1*Cos_xH;
                    M(4,3) = -a4*Cos_yH;
                    M(4,4) = a4*Sin_yH;
                    M(5,1) = a3*Cos_xH;
                    M(5,2) = -a3*Sin_xH;
                    M(5,3) = a2*Sin_yH;
                    M(5,4) = a2*Cos_yH;
                    Y(i) = abs(det(M));
                end
                if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                    if  o == Bisections
                        HBScholte(end+1) = Frequency(i-1);
                    end
                    Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                    break
                end
            end
        end
    end
    String = ['Frq. @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode         Frq.(kHz)'];
    if  any(HBScholte)
        for i = 1:length(HBScholte)
            if  i < 10
                if  HBScholte(i) < 1e2
                    String = append(String,newline,'BScholte',num2str(i),'      ',num2str(HBScholte(i),'%.3f'));
                elseif HBScholte(i) >= 1e2 && HBScholte(i) < 1e3
                    String = append(String,newline,'BScholte',num2str(i),'     ',num2str(HBScholte(i),'%.3f'));
                elseif HBScholte(i) >= 1e3 && HBScholte(i) < 1e4
                    String = append(String,newline,'BScholte',num2str(i),'    ',num2str(HBScholte(i),'%.3f'));
                elseif HBScholte(i) >= 1e4
                    String = append(String,newline,'BScholte',num2str(i),'   ',num2str(HBScholte(i),'%.3f'));
                end
            else
                if  HBScholte(i) < 1e2
                    String = append(String,newline,'BScholte',num2str(i),'     ',num2str(HBScholte(i),'%.3f'));
                elseif HBScholte(i) >= 1e2 && HBScholte(i) < 1e3
                    String = append(String,newline,'BScholte',num2str(i),'    ',num2str(HBScholte(i),'%.3f'));
                elseif HBScholte(i) >= 1e3 && HBScholte(i) < 1e4
                    String = append(String,newline,'BScholte',num2str(i),'   ',num2str(HBScholte(i),'%.3f'));
                elseif HBScholte(i) >= 1e4
                    String = append(String,newline,'BScholte',num2str(i),'  ',num2str(HBScholte(i),'%.3f'));
                end
            end
        end
    end
    String = append(String,newline,newline);
    if  any(HBScholte)
        String = append(String,'BScholte: ',num2str(length(HBScholte)));
    end
end
disp([String,newline,'----------------'])