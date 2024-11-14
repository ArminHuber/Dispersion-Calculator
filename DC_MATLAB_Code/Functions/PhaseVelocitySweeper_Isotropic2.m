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
function HBScholte = PhaseVelocitySweeper_Isotropic2(Material,Viscoelastic,Half,UpperFluid,LowerFluid,Frequency)
Resolution = 1e-6; % (m/s)
VelocityStep = 1; % (m/s); for coarse sweep
VelocityStepFine = .25; % (m/s); for first iteration in the fine sweeps
EnlargementReal = 3; % SweepRangeReal = XRough(p,1)+EnlargementReal*VelocityStep*[-1 1];
EnlargementImag = 2; % SweepRangeImag = XRough(p,2)+EnlargementImag*VelocityStep*[-1 1];
ImagRealThreshold1 = .1; % exclude wavenumbers with ratios of imaginary part/real part greater than this threshold
ImagRealThreshold2 = .5; % for above slow fluid velocity
DensityThreshold = .01; % rho_fluid/rho_solid
% disp(1+2*VelocityStep/VelocityStepFine*[EnlargementReal EnlargementImag]) % array size of first iteration in the fine sweeps

%#ok<*AGROW>
if  UpperFluid.Velocity > LowerFluid.Velocity
    FastFluidVelocity = UpperFluid.Velocity;
    SlowFluidVelocity = LowerFluid.Velocity;
    Density_ = LowerFluid.Density/Material.Density;
else
    FastFluidVelocity = LowerFluid.Velocity;
    SlowFluidVelocity = UpperFluid.Velocity;
    Density_ = UpperFluid.Density/Material.Density;
end
if  Density_ > DensityThreshold
    RangeImag = .2;
else
    RangeImag = .04;
end
H = cell(1,2);
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
AngularFrequency = 2*pi*Frequency*1e3;
kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
kUpperFluid2 = (AngularFrequency/UpperFluid.Velocity)^2;
kLowerFluid2 = (AngularFrequency/LowerFluid.Velocity)^2;
for r = 1:2
    XRough = [];
    if  r == 1
        SweepRangeReal = 50:VelocityStep:SlowFluidVelocity+VelocityStep;
        if  Viscoelastic
            SweepRangeImag = -1e-10:-VelocityStep:-RangeImag*SlowFluidVelocity;
        else
            SweepRangeImag = -1e-10;
        end
        SweepRangeImag(abs(SweepRangeImag/SweepRangeReal(end)) > ImagRealThreshold1) = [];
    elseif r == 2
        SweepRangeReal = SlowFluidVelocity-VelocityStep:VelocityStep:FastFluidVelocity+VelocityStep;
        SweepRangeImag = -1e-10:-VelocityStep:-RangeImag*FastFluidVelocity;
        SweepRangeImag(abs(SweepRangeImag/SweepRangeReal(end)) > ImagRealThreshold2) = [];
    end
    Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
    if  r == 1
        Wavenumber(imag(Wavenumber)./real(Wavenumber) > ImagRealThreshold1) = NaN;
    elseif r == 2
        Wavenumber(imag(Wavenumber)./real(Wavenumber) > ImagRealThreshold2) = NaN;
    end
    Y = Computer(r,Wavenumber,kL2,kT2,Half,Lambda,Mu,LowerFluid,UpperFluid,kLowerFluid2,kUpperFluid2);
    Y(isinf(Y)) = NaN;
    if  ~Viscoelastic && r == 1
        for j = 2:size(Y,1)-1
            if  Y(j) < Y(j-1) && Y(j) < Y(j+1)
                XRough(end+1,:) = [SweepRangeReal(j) SweepRangeImag]; % real velocity (m/s), imaginary velocity (m/s)
            end
        end
% figure,plot(SweepRangeReal,20*log10(Y))
    else
        Y = [max(Y,[],'all')*ones(height(Y),1) Y]; % CAUTION: XRough is at SweepRangeImag(l-1) instead of SweepRangeImag(l) if wall is added at the START of SweepRangeImag instead of at the END!
        for l = 2:size(Y,2)-1
            for j = 2:size(Y,1)-1
                if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                    XRough(end+1,:) = [SweepRangeReal(j) SweepRangeImag(l-1)]; % real velocity (m/s), imaginary velocity (m/s)
                end
            end
        end
% figure;surf([-SweepRangeImag(2) SweepRangeImag],SweepRangeReal,20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
    end
    if  ~isempty(XRough)
        if  r == 1
            XRough(end+1,:) = [SlowFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        elseif r == 2
            XRough(end+1,:) = [FastFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        end
    else
        if  r == 1
            XRough = [SlowFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        elseif r == 2
            XRough = [FastFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        end
    end
% disp(XRough)
    p = 1;
    while p <= height(XRough)
        SweepRangeReal = XRough(p,1)+EnlargementReal*VelocityStep*[-1 1];
        SweepRangeImag = XRough(p,2)+EnlargementImag*VelocityStep*[-1 1];
        for k = 1:1e2
            if  k == 1
                SweepRangeReal = SweepRangeReal(1):VelocityStepFine:SweepRangeReal(end);
                SweepRangeImag = SweepRangeImag(1):VelocityStepFine:SweepRangeImag(end);
                if  r == 1
                    SweepRangeReal(SweepRangeReal > SlowFluidVelocity) = [];
                elseif r == 2
                    SweepRangeReal(SweepRangeReal > FastFluidVelocity) = [];
                    SweepRangeReal(SweepRangeReal < SlowFluidVelocity) = [];
                end
            else
                SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
            end
            SweepRangeImag(SweepRangeImag > -1e-10) = [];
            Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
            cp = AngularFrequency./real(Wavenumber);
            for l = 1:width(Wavenumber)
                for j = 2:height(Wavenumber)-1
                    if  cp(j-1,l) < Material.TransverseVelocity && cp(j+1,l) > Material.TransverseVelocity
                        Wavenumber(j,l) = NaN;
                    end
                end
            end
            if  ~isempty(H{r})
                for j = 2:height(Wavenumber)-1
                    for m = 1:height(H{r})
                        if  SweepRangeReal(j-1) < H{r}(m,3) && SweepRangeReal(j+1) > H{r}(m,3)
                            Wavenumber(j,:) = NaN;
                        end
                    end
                end
            end
            Y = Computer(r,Wavenumber,kL2,kT2,Half,Lambda,Mu,LowerFluid,UpperFluid,kLowerFluid2,kUpperFluid2);
            Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
% if  k <= 1
% f = figure;surf([SweepRangeImag SweepRangeImag(end)-(SweepRangeImag(1)-SweepRangeImag(2))],SweepRangeReal,20*log10(Y))
% close(f)
% end
            Min = zeros(size(Y));
            for l = 2:size(Y,2)-1
                for j = 2:size(Y,1)-1
                    if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                        Min(j,l) = 1;
                    end
                end
            end
            [b1,b2] = find(Min);
            if  ~isempty(b1)
                for l = 1:length(b1)
                    if  l == 1
                        MIN = [b1(1) b2(1)];
                    else
                        XRough(end+1,:) = [SweepRangeReal(b1(l)) SweepRangeImag(b2(l))]; % real velocity (m/s), imaginary velocity (m/s)
                    end
                end
            else
                MIN = 0;
                break
            end
            if  k == 100 || Resolution > abs(SweepRangeReal(1)-SweepRangeReal(2))
                break
            end
            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
        end
        if  any(MIN)
            H{r}(end+1,:) = [AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) imag(Wavenumber(MIN(1),MIN(2))) SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))]; % phase velocity (m/s), attenuation (Np/m), real velocity (m/s), imaginary velocity (m/s)
        end
        p = p+1;
    end
    if  ~isempty(H{r})
        H{r} = sortrows(H{r},1);
    end
end
HBScholte = [H{1};H{2}];
if  isempty(HBScholte)
    disp('No higher order Scholte modes found!')
    return
end
String = ['cp @ ',num2str(Frequency),' kHz:',newline,'Mode        cp(m/ms)  alpha(Np/m)'];
for p = 1:height(HBScholte)
    if  p < 11
        String = append(String,newline,'BScholte',num2str(p-1),'    ',num2str(HBScholte(p,1)/1e3,'%.5f'),'   ',num2str(HBScholte(p,2)));
    else
        String = append(String,newline,'BScholte',num2str(p-1),'   ',num2str(HBScholte(p,1)/1e3,'%.5f'),'   ',num2str(HBScholte(p,2)));
    end
end
disp(append(String,newline,newline,'BScholte: ',num2str(height(HBScholte)),newline,'----------------'));
end
function Y = Computer(r,Wavenumber,kL2,kT2,Half,Lambda,Mu,LowerFluid,UpperFluid,kLowerFluid2,kUpperFluid2)
    k2 = Wavenumber.^2;
    x = sqrt(kL2-k2);
    y = sqrt(kT2-k2);
    if  r == 1
        k3UpperFluid = -sqrt(kUpperFluid2-k2);
        k3LowerFluid = -sqrt(kLowerFluid2-k2);
    else
        if  UpperFluid.Velocity > LowerFluid.Velocity
            k3UpperFluid = -sqrt(kUpperFluid2-k2);
            k3LowerFluid = sqrt(kLowerFluid2-k2);
        else
            k3UpperFluid = sqrt(kUpperFluid2-k2);
            k3LowerFluid = -sqrt(kLowerFluid2-k2);
        end
    end
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
    a1 = -(Lambda*kL2+2*Mu*x.^2);
    a2 = Mu*(k2-y.^2);
    a3 = 2i*Mu*Wavenumber.*x;
    a4 = 2i*Mu*Wavenumber.*y;
    a5 = -1i*Wavenumber;
    a6 = -UpperFluid.Velocity^2*UpperFluid.Density*kUpperFluid2;
    a7 = LowerFluid.Velocity^2*LowerFluid.Density*kLowerFluid2;
    a8 = 1i*k3UpperFluid;
    a9 = 1i*k3LowerFluid;
    e = exp(1i*k3UpperFluid*Half);
    e_ = exp(-1i*k3LowerFluid*Half);
    Y = NaN(size(Wavenumber));
    for l = 1:width(Wavenumber)
        for j = 1:height(Wavenumber)
            if  ~isnan(Wavenumber(j,l))
                M(1,1) = a1(j,l)*SinxH(j,l);
                M(1,2) = a1(j,l)*CosxH(j,l);
                M(1,3) = -a4(j,l)*CosyH(j,l);
                M(1,4) = a4(j,l)*SinyH(j,l);
                M(1,5) = a6*e(j,l);
                M(2,1) = a3(j,l)*CosxH(j,l);
                M(2,2) = -a3(j,l)*SinxH(j,l);
                M(2,3) = a2(j,l)*SinyH(j,l);
                M(2,4) = a2(j,l)*CosyH(j,l);
                M(3,1) = x(j,l)*CosxH(j,l);
                M(3,2) = -x(j,l)*SinxH(j,l);
                M(3,3) = a5(j,l)*SinyH(j,l);
                M(3,4) = a5(j,l)*CosyH(j,l);
                M(3,5) = a8(j,l)*e(j,l);
                M(4,1) = a1(j,l)*Sin_xH(j,l);
                M(4,2) = a1(j,l)*Cos_xH(j,l);
                M(4,3) = -a4(j,l)*Cos_yH(j,l);
                M(4,4) = a4(j,l)*Sin_yH(j,l);
                M(4,6) = a7*e_(j,l);
                M(5,1) = a3(j,l)*Cos_xH(j,l);
                M(5,2) = -a3(j,l)*Sin_xH(j,l);
                M(5,3) = a2(j,l)*Sin_yH(j,l);
                M(5,4) = a2(j,l)*Cos_yH(j,l);
                M(6,1) = x(j,l)*Cos_xH(j,l);
                M(6,2) = -x(j,l)*Sin_xH(j,l);
                M(6,3) = a5(j,l)*Sin_yH(j,l);
                M(6,4) = a5(j,l)*Cos_yH(j,l);
                M(6,6) = a9(j,l)*e_(j,l);
                Y(j,l) = abs(det(M));
            end
        end
    end
end