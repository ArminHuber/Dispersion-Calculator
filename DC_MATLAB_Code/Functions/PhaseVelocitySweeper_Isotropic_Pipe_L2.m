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
function HL = PhaseVelocitySweeper_Isotropic_Pipe_L2(Material,Symmetric,Viscoelastic,FluidLoading,Ro,Ri,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Frequency)
% Note: The lower limit of SweepRangeImag is taken as -1e-10 instead
% of 0 because when there is only fluid inside the pipe and no outer fluid 
% and there is a sink inside, the dispersion equation amplitude has a sharp 
% drop at SweepRangeImag = 0. This seems to happen only below the fluid
% velocity.

Resolution = 1e-6; % (m/s)
VelocityStep = 1; % (m/s); for coarse sweep
VelocityStepFine = .05; % (m/s); for first iteration in the fine sweeps
EnlargementReal = 3; % SweepRangeReal = XRough(p,1)+EnlargementReal*VelocityStep*[-1 1];
EnlargementImag = 2; % SweepRangeImag = XRough(p,2)+EnlargementImag*VelocityStep*[-1 1];
ImagRealThreshold1 = .1; % exclude wavenumbers with ratios of imaginary part/real part greater than this threshold
ImagRealThreshold2 = .5; % for above slow fluid velocity
DensityThreshold = .01; % rho_fluid/rho_solid
% disp(1+2*VelocityStep/VelocityStepFine*[EnlargementReal EnlargementImag]) % array size of first iteration in the fine sweeps

%#ok<*AGROW>
if  FluidLoading
    if  ~Symmetric && ToggleOuterFluid && ToggleInnerFluid
        if  OuterFluid.Velocity > InnerFluid.Velocity
            FastFluidVelocity = OuterFluid.Velocity;
            SlowFluidVelocity = InnerFluid.Velocity;
            Density_ = InnerFluid.Density/Material.Density;
        else
            FastFluidVelocity = InnerFluid.Velocity;
            SlowFluidVelocity = OuterFluid.Velocity;
            Density_ = OuterFluid.Density/Material.Density;
        end
    elseif Symmetric || (ToggleOuterFluid && ~ToggleInnerFluid)
        SlowFluidVelocity = OuterFluid.Velocity;
        Density_ = OuterFluid.Density/Material.Density;
    elseif ~ToggleOuterFluid && ToggleInnerFluid
        SlowFluidVelocity = InnerFluid.Velocity;
        Density_ = InnerFluid.Density/Material.Density;
    end  
    if  Density_ > DensityThreshold
        RangeImag = .2;
    else
        RangeImag = .04;
    end
else
    RangeImag = .04;
end
H = cell(1,2);
Ro2 = Ro^2;
Ri2 = Ri^2;
Densityi = InnerFluid.Density/Material.Density;
Densityo = OuterFluid.Density/Material.Density;
AngularFrequency = 2*pi*Frequency*1e3;
kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
kInnerFluid2 = (AngularFrequency/InnerFluid.Velocity)^2;
kOuterFluid2 = (AngularFrequency/OuterFluid.Velocity)^2;
for r = 1:2
    if  r == 2 && (~FluidLoading || ~(~Symmetric && ToggleOuterFluid && ToggleInnerFluid)) 
        break
    end
    XRough = [];
    if  r == 1
        if  FluidLoading
            SweepRangeReal = 50:VelocityStep:SlowFluidVelocity+VelocityStep;
            if  Viscoelastic
                SweepRangeImag = -1e-10:-VelocityStep:-RangeImag*SlowFluidVelocity;
            else
                SweepRangeImag = -1e-10;
            end
        else
            SweepRangeReal = 50:VelocityStep:.99*Material.TransverseVelocity;
            if  Viscoelastic
                SweepRangeImag = -1e-10:-VelocityStep:-RangeImag*Material.TransverseVelocity;
            else
                SweepRangeImag = -1e-10;
            end
        end
        SweepRangeImag(abs(SweepRangeImag/SweepRangeReal(end)) > ImagRealThreshold1) = [];
    elseif r == 2
        SweepRangeReal = SlowFluidVelocity-VelocityStep:VelocityStep:FastFluidVelocity+VelocityStep;
        if  ~Viscoelastic && ~Sink && OuterFluid.Velocity > InnerFluid.Velocity
            SweepRangeImag = -1e-10;
        else
            SweepRangeImag = -1e-10:-VelocityStep:-RangeImag*FastFluidVelocity;
            SweepRangeImag(abs(SweepRangeImag/SweepRangeReal(end)) > ImagRealThreshold2) = [];
        end
    end
    Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
    if  r == 1
        Wavenumber(imag(Wavenumber)./real(Wavenumber) > ImagRealThreshold1) = NaN;
    elseif r == 2
        Wavenumber(imag(Wavenumber)./real(Wavenumber) > ImagRealThreshold2) = NaN;
    end
    Y = Computer(r,Wavenumber,kL2,kT2,Ro,Ri,Ro2,Ri2,Densityo,Densityi,FluidLoading,ToggleInnerFluid,ToggleOuterFluid,Sink,InnerFluid,OuterFluid,kInnerFluid2,kOuterFluid2);
    Y(isinf(Y)) = NaN;
    if  ~Viscoelastic && (r == 1 || (r == 2 && ~Sink && OuterFluid.Velocity > InnerFluid.Velocity))
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
        if  r == 1 && FluidLoading
            XRough(end+1,:) = [SlowFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        elseif r == 2
            XRough(end+1,:) = [FastFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        end
    else
        if  r == 1 && FluidLoading
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
                    if  FluidLoading
                        SweepRangeReal(SweepRangeReal > SlowFluidVelocity) = [];
                    else
                        SweepRangeReal(SweepRangeReal > .99*Material.TransverseVelocity) = [];
                    end
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
            Y = Computer(r,Wavenumber,kL2,kT2,Ro,Ri,Ro2,Ri2,Densityo,Densityi,FluidLoading,ToggleInnerFluid,ToggleOuterFluid,Sink,InnerFluid,OuterFluid,kInnerFluid2,kOuterFluid2);
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
HL = [H{1};H{2}];
if  isempty(HL)
    disp('No higher order modes found!')
    return
end
String = ['cp @ ',num2str(Frequency),' kHz:',newline,'Mode     cp(m/ms)  alpha(Np/m)'];
for p = 1:height(HL)
    if  p < 10
        String = append(String,newline,'L(0,',num2str(p),')    ',num2str(HL(p,1)/1e3,'%.5f'),'   ',num2str(HL(p,2)));
    else
        String = append(String,newline,'L(0,',num2str(p),')   ',num2str(HL(p,1)/1e3,'%.5f'),'   ',num2str(HL(p,2)));
    end
end
disp(append(String,newline,newline,'L: ',num2str(height(HL)),newline,'----------------'));
end
function Y = Computer(r,Wavenumber,kL2,kT2,Ro,Ri,Ro2,Ri2,Densityo,Densityi,FluidLoading,ToggleInnerFluid,ToggleOuterFluid,Sink,InnerFluid,OuterFluid,kInnerFluid2,kOuterFluid2)
    k2 = Wavenumber.^2;
    x = sqrt(k2-kL2);
    y = sqrt(k2-kT2);
    y2 = kT2-k2;
    xRo = x*Ro;
    yRo = y*Ro;
    xRi = x*Ri;
    yRi = y*Ri;
    Z0xi = besseli(0,xRi);
    Z0xo = besseli(0,xRo);
    Z0yi = besseli(0,yRi);
    Z0yo = besseli(0,yRo);
    W0xi = besselk(0,xRi);
    W0xo = besselk(0,xRo);
    W0yi = besselk(0,yRi);
    W0yo = besselk(0,yRo);
    Z1xiRi = -xRi.*besseli(1,xRi);
    Z1xoRo = -xRo.*besseli(1,xRo);
    Z1yiRi = -yRi.*besseli(1,yRi);
    Z1yoRo = -yRo.*besseli(1,yRo);
    W1xiRi = xRi.*besselk(1,xRi);
    W1xoRo = xRo.*besselk(1,xRo);
    W1yiRi = yRi.*besselk(1,yRi);
    W1yoRo = yRo.*besselk(1,yRo);
    if  FluidLoading
        if  r == 1
            zRi = -sqrt(kInnerFluid2-k2)*Ri;
            zRo = -sqrt(kOuterFluid2-k2)*Ro;
        else
            if  OuterFluid.Velocity > InnerFluid.Velocity
                zRi = sqrt(kInnerFluid2-k2)*Ri;
                zRo = -sqrt(kOuterFluid2-k2)*Ro;
            else
                zRi = -sqrt(kInnerFluid2-k2)*Ri;
                zRo = sqrt(kOuterFluid2-k2)*Ro;
            end
        end
    end
    if  ToggleInnerFluid
        if  Sink
            Z0zi = besselh(0,2,zRi);
            Z1ziRi = -zRi.*besselh(1,2,zRi);
        else
            Z0zi = besselj(0,zRi);
            Z1ziRi = -zRi.*besselj(1,zRi);
        end
    end
    if  ToggleOuterFluid
        H0zo = besselh(0,zRo);
        H1zoRo = -zRo.*besselh(1,zRo);
    end
    Y = NaN(size(Wavenumber));
    for l = 1:width(Wavenumber)
        for j = 1:height(Wavenumber)
            if  ~isnan(Wavenumber(j,l))
                if  ToggleInnerFluid && ToggleOuterFluid
                    M(1,1) = Ri2*(.5*kT2-k2(j,l))*Z0xi(j,l)-Z1xiRi(j,l);
                    M(1,2) = Ri2*(.5*kT2-k2(j,l))*W0xi(j,l)-W1xiRi(j,l);
                    M(1,3) = Wavenumber(j,l)*(Z1yiRi(j,l)-y2(j,l)*Ri2*Z0yi(j,l));
                    M(1,4) = Wavenumber(j,l)*(W1yiRi(j,l)-y2(j,l)*Ri2*W0yi(j,l));
                    M(1,5) = .5*kT2*Ri2*Densityi*Z0zi(j,l);
                    M(2,1) = -2*Wavenumber(j,l)*Z1xiRi(j,l);
                    M(2,2) = -2*Wavenumber(j,l)*W1xiRi(j,l);
                    M(2,3) = (k2(j,l)-y2(j,l))*Z1yiRi(j,l);
                    M(2,4) = (k2(j,l)-y2(j,l))*W1yiRi(j,l);
                    M(3,1) = -Z1xiRi(j,l);
                    M(3,2) = -W1xiRi(j,l);
                    M(3,3) = Wavenumber(j,l)*Z1yiRi(j,l);
                    M(3,4) = Wavenumber(j,l)*W1yiRi(j,l);
                    M(3,5) = Z1ziRi(j,l);
                    M(4,1) = Ro2*(.5*kT2-k2(j,l))*Z0xo(j,l)-Z1xoRo(j,l);
                    M(4,2) = Ro2*(.5*kT2-k2(j,l))*W0xo(j,l)-W1xoRo(j,l);
                    M(4,3) = Wavenumber(j,l)*(Z1yoRo(j,l)-y2(j,l)*Ro2*Z0yo(j,l));
                    M(4,4) = Wavenumber(j,l)*(W1yoRo(j,l)-y2(j,l)*Ro2*W0yo(j,l));
                    M(4,6) = .5*kT2*Ro2*Densityo*H0zo(j,l);
                    M(5,1) = -2*Wavenumber(j,l)*Z1xoRo(j,l);
                    M(5,2) = -2*Wavenumber(j,l)*W1xoRo(j,l);
                    M(5,3) = (k2(j,l)-y2(j,l))*Z1yoRo(j,l);
                    M(5,4) = (k2(j,l)-y2(j,l))*W1yoRo(j,l);
                    M(6,1) = -Z1xoRo(j,l);
                    M(6,2) = -W1xoRo(j,l);
                    M(6,3) = Wavenumber(j,l)*Z1yoRo(j,l);
                    M(6,4) = Wavenumber(j,l)*W1yoRo(j,l);
                    M(6,6) = H1zoRo(j,l);
                elseif ToggleInnerFluid && ~ToggleOuterFluid
                    M(1,1) = Ri2*(.5*kT2-k2(j,l))*Z0xi(j,l)-Z1xiRi(j,l);
                    M(1,2) = Ri2*(.5*kT2-k2(j,l))*W0xi(j,l)-W1xiRi(j,l);
                    M(1,3) = Wavenumber(j,l)*(Z1yiRi(j,l)-y2(j,l)*Ri2*Z0yi(j,l));
                    M(1,4) = Wavenumber(j,l)*(W1yiRi(j,l)-y2(j,l)*Ri2*W0yi(j,l));
                    M(1,5) = .5*kT2*Ri2*Densityi*Z0zi(j,l);
                    M(2,1) = -2*Wavenumber(j,l)*Z1xiRi(j,l);
                    M(2,2) = -2*Wavenumber(j,l)*W1xiRi(j,l);
                    M(2,3) = (k2(j,l)-y2(j,l))*Z1yiRi(j,l);
                    M(2,4) = (k2(j,l)-y2(j,l))*W1yiRi(j,l);
                    M(3,1) = -Z1xiRi(j,l);
                    M(3,2) = -W1xiRi(j,l);
                    M(3,3) = Wavenumber(j,l)*Z1yiRi(j,l);
                    M(3,4) = Wavenumber(j,l)*W1yiRi(j,l);
                    M(3,5) = Z1ziRi(j,l);
                    M(4,1) = Ro2*(.5*kT2-k2(j,l))*Z0xo(j,l)-Z1xoRo(j,l);
                    M(4,2) = Ro2*(.5*kT2-k2(j,l))*W0xo(j,l)-W1xoRo(j,l);
                    M(4,3) = Wavenumber(j,l)*(Z1yoRo(j,l)-y2(j,l)*Ro2*Z0yo(j,l));
                    M(4,4) = Wavenumber(j,l)*(W1yoRo(j,l)-y2(j,l)*Ro2*W0yo(j,l));
                    M(5,1) = -2*Wavenumber(j,l)*Z1xoRo(j,l);
                    M(5,2) = -2*Wavenumber(j,l)*W1xoRo(j,l);
                    M(5,3) = (k2(j,l)-y2(j,l))*Z1yoRo(j,l);
                    M(5,4) = (k2(j,l)-y2(j,l))*W1yoRo(j,l);
                elseif ~ToggleInnerFluid && ToggleOuterFluid
                    M(1,1) = Ri2*(.5*kT2-k2(j,l))*Z0xi(j,l)-Z1xiRi(j,l);
                    M(1,2) = Ri2*(.5*kT2-k2(j,l))*W0xi(j,l)-W1xiRi(j,l);
                    M(1,3) = Wavenumber(j,l)*(Z1yiRi(j,l)-y2(j,l)*Ri2*Z0yi(j,l));
                    M(1,4) = Wavenumber(j,l)*(W1yiRi(j,l)-y2(j,l)*Ri2*W0yi(j,l));
                    M(2,1) = -2*Wavenumber(j,l)*Z1xiRi(j,l);
                    M(2,2) = -2*Wavenumber(j,l)*W1xiRi(j,l);
                    M(2,3) = (k2(j,l)-y2(j,l))*Z1yiRi(j,l);
                    M(2,4) = (k2(j,l)-y2(j,l))*W1yiRi(j,l);
                    M(3,1) = Ro2*(.5*kT2-k2(j,l))*Z0xo(j,l)-Z1xoRo(j,l);
                    M(3,2) = Ro2*(.5*kT2-k2(j,l))*W0xo(j,l)-W1xoRo(j,l);
                    M(3,3) = Wavenumber(j,l)*(Z1yoRo(j,l)-y2(j,l)*Ro2*Z0yo(j,l));
                    M(3,4) = Wavenumber(j,l)*(W1yoRo(j,l)-y2(j,l)*Ro2*W0yo(j,l));
                    M(3,5) = .5*kT2*Ro2*Densityo*H0zo(j,l);
                    M(4,1) = -2*Wavenumber(j,l)*Z1xoRo(j,l);
                    M(4,2) = -2*Wavenumber(j,l)*W1xoRo(j,l);
                    M(4,3) = (k2(j,l)-y2(j,l))*Z1yoRo(j,l);
                    M(4,4) = (k2(j,l)-y2(j,l))*W1yoRo(j,l);
                    M(5,1) = -Z1xoRo(j,l);
                    M(5,2) = -W1xoRo(j,l);
                    M(5,3) = Wavenumber(j,l)*Z1yoRo(j,l);
                    M(5,4) = Wavenumber(j,l)*W1yoRo(j,l);
                    M(5,5) = H1zoRo(j,l);
                elseif ~ToggleInnerFluid && ~ToggleOuterFluid
                    M(1,1) = Ri2*(.5*kT2-k2(j,l))*Z0xi(j,l)-Z1xiRi(j,l);
                    M(1,2) = Ri2*(.5*kT2-k2(j,l))*W0xi(j,l)-W1xiRi(j,l);
                    M(1,3) = Wavenumber(j,l)*(Z1yiRi(j,l)-y2(j,l)*Ri2*Z0yi(j,l));
                    M(1,4) = Wavenumber(j,l)*(W1yiRi(j,l)-y2(j,l)*Ri2*W0yi(j,l));
                    M(2,1) = -2*Wavenumber(j,l)*Z1xiRi(j,l);
                    M(2,2) = -2*Wavenumber(j,l)*W1xiRi(j,l);
                    M(2,3) = (k2(j,l)-y2(j,l))*Z1yiRi(j,l);
                    M(2,4) = (k2(j,l)-y2(j,l))*W1yiRi(j,l);
                    M(3,1) = Ro2*(.5*kT2-k2(j,l))*Z0xo(j,l)-Z1xoRo(j,l);
                    M(3,2) = Ro2*(.5*kT2-k2(j,l))*W0xo(j,l)-W1xoRo(j,l);
                    M(3,3) = Wavenumber(j,l)*(Z1yoRo(j,l)-y2(j,l)*Ro2*Z0yo(j,l));
                    M(3,4) = Wavenumber(j,l)*(W1yoRo(j,l)-y2(j,l)*Ro2*W0yo(j,l));
                    M(4,1) = -2*Wavenumber(j,l)*Z1xoRo(j,l);
                    M(4,2) = -2*Wavenumber(j,l)*W1xoRo(j,l);
                    M(4,3) = (k2(j,l)-y2(j,l))*Z1yoRo(j,l);
                    M(4,4) = (k2(j,l)-y2(j,l))*W1yoRo(j,l);
                end
                Y(j,l) = abs(det(M));
            end
        end
    end
end