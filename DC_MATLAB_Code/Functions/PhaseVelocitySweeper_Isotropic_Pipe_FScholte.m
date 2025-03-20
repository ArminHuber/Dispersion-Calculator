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
function HF = PhaseVelocitySweeper_Isotropic_Pipe_FScholte(Material,Symmetric,Viscoelastic,FluidLoading,Ro,Ri,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Frequency,FlexuralModeOrders)
% Note: The lower limit of SweepRangeImag is taken as -1e-10 instead
% of 0 because when there is only fluid inside the pipe and no outer fluid 
% and there is a sink inside, the dispersion equation amplitude has a sharp 
% drop at SweepRangeImag = 0. This seems to happen only below the fluid
% velocity.

Resolution = 1e-6; % (m/s); CAUTION: Computer_Isotropic_Pipe_FnScholte_D2 ==> Hdiff(abs(H(:,1)-Material.TransverseVelocity) < 1e-6) = []; H(abs(H(:,1)-Material.TransverseVelocity) < 1e-6,:) = [];
VelocityStep = 1; % (m/s); for coarse sweep
VelocityStepFine = .05; % (m/s); for first iteration in the fine sweeps
VelocityStepSuperFine = .0025; % (m/s); for the additional search for the solution very close to transverse velocity
EnlargementReal = 3; % SweepRangeReal = XRough(p,1)+EnlargementReal*VelocityStep*[-1 1];
EnlargementImag = 2; % SweepRangeImag = XRough(p,2)+EnlargementImag*VelocityStep*[-1 1];
EnlargementRealFine = 10; % SweepRangeReal = XRough(p,1)+EnlargementRealFine*VelocityStepFine*[-1 1]; for the additional search for the solution very close to transverse velocity
EnlargementImagFine = 2; % SweepRangeImag = XRough(p,2)+EnlargementImagFine*VelocityStepFine*[-1 1];
ImagRealThreshold1 = .1; % exclude wavenumbers with ratios of imaginary part/real part greater than this threshold
ImagRealThreshold2 = .5; % for above slow fluid velocity
DensityThreshold = .01; % rho_fluid/rho_solid
% disp(1+2*VelocityStep/VelocityStepFine*[EnlargementReal EnlargementImag]) % array size of first iteration in the fine sweeps
% disp(1+2*VelocityStepFine/VelocityStepSuperFine*[EnlargementRealFine EnlargementImagFine]) % for the additional search for the solution very close to transverse velocity

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
Ro2 = Ro^2;
Ri2 = Ri^2;
Densityi = InnerFluid.Density/Material.Density;
Densityo = OuterFluid.Density/Material.Density;
AngularFrequency = 2*pi*Frequency*1e3;
kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
kInnerFluid2 = (AngularFrequency/InnerFluid.Velocity)^2;
kOuterFluid2 = (AngularFrequency/OuterFluid.Velocity)^2;
TotalModeNumber = 0;
disp(['cp @ ',num2str(Frequency),' kHz:',newline,'Mode     cp(m/ms)  alpha(Np/m)']);
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',FlexuralModeOrders),'Name',['F-mode orders at ',num2str(Frequency),' kHz...']);
for n = 1:FlexuralModeOrders
    H = cell(1,2);
    n2 = n^2;
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
                SweepRangeReal = 50:VelocityStep:.999*Material.TransverseVelocity;
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
        Y = Computer(n,n2,r,Wavenumber,kL2,kT2,Ro,Ri,Ro2,Ri2,Densityo,Densityi,FluidLoading,ToggleInnerFluid,ToggleOuterFluid,Sink,InnerFluid,OuterFluid,kInnerFluid2,kOuterFluid2);
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
            if  abs(XRough(p,1)-Material.TransverseVelocity) <= VelocityStep && (XRough(p,2) == -1e-10 || (XRough(p,2) < -1e-10 && mod(XRough(p,1),VelocityStep) > 0))
                SweepRangeReal = XRough(p,1)+EnlargementRealFine*VelocityStepFine*[-1 1];
                SweepRangeImag = XRough(p,2)+EnlargementImagFine*VelocityStepFine*[-1 1];
            else
                SweepRangeReal = XRough(p,1)+EnlargementReal*VelocityStep*[-1 1];
                SweepRangeImag = XRough(p,2)+EnlargementImag*VelocityStep*[-1 1];
            end
            for k = 1:1e2
                if  k == 1
                    if  abs(XRough(p,1)-Material.TransverseVelocity) <= VelocityStep && (XRough(p,2) == -1e-10 || (XRough(p,2) < -1e-10 && mod(XRough(p,1),VelocityStep) > 0))
                        SweepRangeReal = SweepRangeReal(1):VelocityStepSuperFine:SweepRangeReal(end);
                        SweepRangeImag = SweepRangeImag(1):VelocityStepSuperFine:SweepRangeImag(end);
                    else
                        SweepRangeReal = SweepRangeReal(1):VelocityStepFine:SweepRangeReal(end);
                        SweepRangeImag = SweepRangeImag(1):VelocityStepFine:SweepRangeImag(end);
                    end
                    if  r == 1
                        if  FluidLoading
                            SweepRangeReal(SweepRangeReal > SlowFluidVelocity) = [];
                        else
                            SweepRangeReal(SweepRangeReal > .999*Material.TransverseVelocity) = [];
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
                if  ~isempty(H{r})
                    for j = 2:height(Wavenumber)-1
                        for m = 1:height(H{r})
                            if  SweepRangeReal(j-1) < H{r}(m,3) && SweepRangeReal(j+1) > H{r}(m,3)
                                Wavenumber(j,:) = NaN;
                            end
                        end
                    end
                end
                Y = Computer(n,n2,r,Wavenumber,kL2,kT2,Ro,Ri,Ro2,Ri2,Densityo,Densityi,FluidLoading,ToggleInnerFluid,ToggleOuterFluid,Sink,InnerFluid,OuterFluid,kInnerFluid2,kOuterFluid2);
                Y(:,end+1) = max(Y,[],'all'); % add wall below zero attenuation to allow finding minimum at zero attenuation
% if  p==3&&k <= 1
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
                if  abs(H{r}(end,1)-Material.TransverseVelocity) < Resolution % additional search for solution very close to transverse velocity
                    XRough(end+1,:) = [SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))]; % real velocity (m/s), imaginary velocity (m/s)
                end
            end
            p = p+1;
        end
        % H{r}(abs(H{r}(:,1)-Material.TransverseVelocity) < Resolution,:) = [];
        if  ~isempty(H{r})
            H{r} = sortrows(H{r},1);
        end
    end
    HF{n} = [H{1};H{2}];
    if  isempty(HF{n})
        disp('No higher order modes found!')
        close(h)
        return
    end
    String = '';
    Index = 1;
    for p = 1:height(HF{n})
        if  abs(HF{n}(p,1)-Material.TransverseVelocity) > Resolution
            if  Index < 10
                String = append(String,newline,'F(',num2str(n),',',num2str(Index),')    ',num2str(HF{n}(p,1)/1e3,'%.5f'),'   ',num2str(HF{n}(p,2)));
            else
                String = append(String,newline,'F(',num2str(n),',',num2str(Index),')   ',num2str(HF{n}(p,1)/1e3,'%.5f'),'   ',num2str(HF{n}(p,2)));  
            end
            Index = Index+1;
        end
    end
    disp(String)
    waitbar(n/FlexuralModeOrders,h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',n,FlexuralModeOrders,100*n/FlexuralModeOrders,toc))
    TotalModeNumber = TotalModeNumber+Index-1;
end
close(h)
disp([newline,'F: ',num2str(TotalModeNumber),newline,'----------------']);
end
function Y = Computer(n,n2,r,Wavenumber,kL2,kT2,Ro,Ri,Ro2,Ri2,Densityo,Densityi,FluidLoading,ToggleInnerFluid,ToggleOuterFluid,Sink,InnerFluid,OuterFluid,kInnerFluid2,kOuterFluid2)
    k2 = Wavenumber.^2;
    x = sqrt(k2-kL2);
    y = sqrt(k2-kT2);
    y2 = kT2-k2;
    xRo = x*Ro;
    yRo = y*Ro;
    xRi = x*Ri;
    yRi = y*Ri;
    Znxi = besseli(n,xRi);
    Znxo = besseli(n,xRo);
    Znyi = besseli(n,yRi);
    Znyo = besseli(n,yRo);
    Wnxi = besselk(n,xRi);
    Wnxo = besselk(n,xRo);
    Wnyi = besselk(n,yRi);
    Wnyo = besselk(n,yRo);
    dZnxiRi = n*Znxi+xRi.*besseli(n+1,xRi);
    dZnxoRo = n*Znxo+xRo.*besseli(n+1,xRo);
    dZnyiRi = n*Znyi+yRi.*besseli(n+1,yRi);
    dZnyoRo = n*Znyo+yRo.*besseli(n+1,yRo);
    dWnxiRi = n*Wnxi-xRi.*besselk(n+1,xRi);
    dWnxoRo = n*Wnxo-xRo.*besselk(n+1,xRo);
    dWnyiRi = n*Wnyi-yRi.*besselk(n+1,yRi);
    dWnyoRo = n*Wnyo-yRo.*besselk(n+1,yRo);
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
            Znzi = besselh(n,2,zRi);
            dZnziRi = n*Znzi-zRi.*besselh(n+1,2,zRi);
        else
            Znzi = besselj(n,zRi);
            dZnziRi = n*Znzi-zRi.*besselj(n+1,zRi);
        end
    end
    if  ToggleOuterFluid
        Hnzo = besselh(n,zRo);
        dHnzoRo = n*Hnzo-zRo.*besselh(n+1,zRo);
    end
    Y = NaN(size(Wavenumber));
    for l = 1:width(Wavenumber)
        for j = 1:height(Wavenumber)
            if  ~isnan(Wavenumber(j,l))
                if  ToggleInnerFluid && ToggleOuterFluid
                    M(1,1) = dZnxiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-n2)*Znxi(j,l);
                    M(1,2) = dWnxiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-n2)*Wnxi(j,l);
                    M(1,3) = -Wavenumber(j,l)*(dZnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Znyi(j,l));
                    M(1,4) = -Wavenumber(j,l)*(dWnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Wnyi(j,l));
                    M(1,5) = n*(dZnyiRi(j,l)-Znyi(j,l));
                    M(1,6) = n*(dWnyiRi(j,l)-Wnyi(j,l));
                    M(1,7) = .5*kT2*Ri2*Densityi*Znzi(j,l);
                    M(2,1) = 2*n*(dZnxiRi(j,l)-Znxi(j,l));
                    M(2,2) = 2*n*(dWnxiRi(j,l)-Wnxi(j,l));
                    M(2,3) = 2*Wavenumber(j,l)*n*(Znyi(j,l)-dZnyiRi(j,l));
                    M(2,4) = 2*Wavenumber(j,l)*n*(Wnyi(j,l)-dWnyiRi(j,l));
                    M(2,5) = 2*dZnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Znyi(j,l);
                    M(2,6) = 2*dWnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Wnyi(j,l);
                    M(3,1) = 2*Wavenumber(j,l)*dZnxiRi(j,l);
                    M(3,2) = 2*Wavenumber(j,l)*dWnxiRi(j,l);
                    M(3,3) = (y2(j,l)-k2(j,l))*dZnyiRi(j,l);
                    M(3,4) = (y2(j,l)-k2(j,l))*dWnyiRi(j,l);
                    M(3,5) = -Wavenumber(j,l)*n*Znyi(j,l);
                    M(3,6) = -Wavenumber(j,l)*n*Wnyi(j,l);
                    M(4,1) = dZnxiRi(j,l);
                    M(4,2) = dWnxiRi(j,l);
                    M(4,3) = -Wavenumber(j,l)*dZnyiRi(j,l);
                    M(4,4) = -Wavenumber(j,l)*dWnyiRi(j,l);
                    M(4,5) = -n*Znyi(j,l);
                    M(4,6) = -n*Wnyi(j,l);
                    M(4,7) = dZnziRi(j,l);
                    M(5,1) = dZnxoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-n2)*Znxo(j,l);
                    M(5,2) = dWnxoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-n2)*Wnxo(j,l);
                    M(5,3) = -Wavenumber(j,l)*(dZnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Znyo(j,l));
                    M(5,4) = -Wavenumber(j,l)*(dWnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Wnyo(j,l));
                    M(5,5) = n*(dZnyoRo(j,l)-Znyo(j,l));
                    M(5,6) = n*(dWnyoRo(j,l)-Wnyo(j,l));
                    M(5,8) = .5*kT2*Ro2*Densityo*Hnzo(j,l);
                    M(6,1) = 2*n*(dZnxoRo(j,l)-Znxo(j,l));
                    M(6,2) = 2*n*(dWnxoRo(j,l)-Wnxo(j,l));
                    M(6,3) = 2*Wavenumber(j,l)*n*(Znyo(j,l)-dZnyoRo(j,l));
                    M(6,4) = 2*Wavenumber(j,l)*n*(Wnyo(j,l)-dWnyoRo(j,l));
                    M(6,5) = 2*dZnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Znyo(j,l);
                    M(6,6) = 2*dWnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Wnyo(j,l);
                    M(7,1) = 2*Wavenumber(j,l)*dZnxoRo(j,l);
                    M(7,2) = 2*Wavenumber(j,l)*dWnxoRo(j,l);
                    M(7,3) = (y2(j,l)-k2(j,l))*dZnyoRo(j,l);
                    M(7,4) = (y2(j,l)-k2(j,l))*dWnyoRo(j,l);
                    M(7,5) = -Wavenumber(j,l)*n*Znyo(j,l);
                    M(7,6) = -Wavenumber(j,l)*n*Wnyo(j,l);
                    M(8,1) = dZnxoRo(j,l);
                    M(8,2) = dWnxoRo(j,l);
                    M(8,3) = -Wavenumber(j,l)*dZnyoRo(j,l);
                    M(8,4) = -Wavenumber(j,l)*dWnyoRo(j,l);
                    M(8,5) = -n*Znyo(j,l);
                    M(8,6) = -n*Wnyo(j,l);
                    M(8,8) = dHnzoRo(j,l);
                elseif ToggleInnerFluid && ~ToggleOuterFluid
                    M(1,1) = dZnxiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-n2)*Znxi(j,l);
                    M(1,2) = dWnxiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-n2)*Wnxi(j,l);
                    M(1,3) = -Wavenumber(j,l)*(dZnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Znyi(j,l));
                    M(1,4) = -Wavenumber(j,l)*(dWnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Wnyi(j,l));
                    M(1,5) = n*(dZnyiRi(j,l)-Znyi(j,l));
                    M(1,6) = n*(dWnyiRi(j,l)-Wnyi(j,l));
                    M(1,7) = .5*kT2*Ri2*Densityi*Znzi(j,l);
                    M(2,1) = 2*n*(dZnxiRi(j,l)-Znxi(j,l));
                    M(2,2) = 2*n*(dWnxiRi(j,l)-Wnxi(j,l));
                    M(2,3) = 2*Wavenumber(j,l)*n*(Znyi(j,l)-dZnyiRi(j,l));
                    M(2,4) = 2*Wavenumber(j,l)*n*(Wnyi(j,l)-dWnyiRi(j,l));
                    M(2,5) = 2*dZnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Znyi(j,l);
                    M(2,6) = 2*dWnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Wnyi(j,l);
                    M(3,1) = 2*Wavenumber(j,l)*dZnxiRi(j,l);
                    M(3,2) = 2*Wavenumber(j,l)*dWnxiRi(j,l);
                    M(3,3) = (y2(j,l)-k2(j,l))*dZnyiRi(j,l);
                    M(3,4) = (y2(j,l)-k2(j,l))*dWnyiRi(j,l);
                    M(3,5) = -Wavenumber(j,l)*n*Znyi(j,l);
                    M(3,6) = -Wavenumber(j,l)*n*Wnyi(j,l);
                    M(4,1) = dZnxiRi(j,l);
                    M(4,2) = dWnxiRi(j,l);
                    M(4,3) = -Wavenumber(j,l)*dZnyiRi(j,l);
                    M(4,4) = -Wavenumber(j,l)*dWnyiRi(j,l);
                    M(4,5) = -n*Znyi(j,l);
                    M(4,6) = -n*Wnyi(j,l);
                    M(4,7) = dZnziRi(j,l);
                    M(5,1) = dZnxoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-n2)*Znxo(j,l);
                    M(5,2) = dWnxoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-n2)*Wnxo(j,l);
                    M(5,3) = -Wavenumber(j,l)*(dZnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Znyo(j,l));
                    M(5,4) = -Wavenumber(j,l)*(dWnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Wnyo(j,l));
                    M(5,5) = n*(dZnyoRo(j,l)-Znyo(j,l));
                    M(5,6) = n*(dWnyoRo(j,l)-Wnyo(j,l));
                    M(6,1) = 2*n*(dZnxoRo(j,l)-Znxo(j,l));
                    M(6,2) = 2*n*(dWnxoRo(j,l)-Wnxo(j,l));
                    M(6,3) = 2*Wavenumber(j,l)*n*(Znyo(j,l)-dZnyoRo(j,l));
                    M(6,4) = 2*Wavenumber(j,l)*n*(Wnyo(j,l)-dWnyoRo(j,l));
                    M(6,5) = 2*dZnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Znyo(j,l);
                    M(6,6) = 2*dWnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Wnyo(j,l);
                    M(7,1) = 2*Wavenumber(j,l)*dZnxoRo(j,l);
                    M(7,2) = 2*Wavenumber(j,l)*dWnxoRo(j,l);
                    M(7,3) = (y2(j,l)-k2(j,l))*dZnyoRo(j,l);
                    M(7,4) = (y2(j,l)-k2(j,l))*dWnyoRo(j,l);
                    M(7,5) = -Wavenumber(j,l)*n*Znyo(j,l);
                    M(7,6) = -Wavenumber(j,l)*n*Wnyo(j,l);
                elseif ~ToggleInnerFluid && ToggleOuterFluid
                    M(1,1) = dZnxiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-n2)*Znxi(j,l);
                    M(1,2) = dWnxiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-n2)*Wnxi(j,l);
                    M(1,3) = -Wavenumber(j,l)*(dZnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Znyi(j,l));
                    M(1,4) = -Wavenumber(j,l)*(dWnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Wnyi(j,l));
                    M(1,5) = n*(dZnyiRi(j,l)-Znyi(j,l));
                    M(1,6) = n*(dWnyiRi(j,l)-Wnyi(j,l));
                    M(2,1) = 2*n*(dZnxiRi(j,l)-Znxi(j,l));
                    M(2,2) = 2*n*(dWnxiRi(j,l)-Wnxi(j,l));
                    M(2,3) = 2*Wavenumber(j,l)*n*(Znyi(j,l)-dZnyiRi(j,l));
                    M(2,4) = 2*Wavenumber(j,l)*n*(Wnyi(j,l)-dWnyiRi(j,l));
                    M(2,5) = 2*dZnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Znyi(j,l);
                    M(2,6) = 2*dWnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Wnyi(j,l);
                    M(3,1) = 2*Wavenumber(j,l)*dZnxiRi(j,l);
                    M(3,2) = 2*Wavenumber(j,l)*dWnxiRi(j,l);
                    M(3,3) = (y2(j,l)-k2(j,l))*dZnyiRi(j,l);
                    M(3,4) = (y2(j,l)-k2(j,l))*dWnyiRi(j,l);
                    M(3,5) = -Wavenumber(j,l)*n*Znyi(j,l);
                    M(3,6) = -Wavenumber(j,l)*n*Wnyi(j,l);
                    M(4,1) = dZnxoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-n2)*Znxo(j,l);
                    M(4,2) = dWnxoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-n2)*Wnxo(j,l);
                    M(4,3) = -Wavenumber(j,l)*(dZnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Znyo(j,l));
                    M(4,4) = -Wavenumber(j,l)*(dWnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Wnyo(j,l));
                    M(4,5) = n*(dZnyoRo(j,l)-Znyo(j,l));
                    M(4,6) = n*(dWnyoRo(j,l)-Wnyo(j,l));
                    M(4,7) = .5*kT2*Ro2*Densityo*Hnzo(j,l);
                    M(5,1) = 2*n*(dZnxoRo(j,l)-Znxo(j,l));
                    M(5,2) = 2*n*(dWnxoRo(j,l)-Wnxo(j,l));
                    M(5,3) = 2*Wavenumber(j,l)*n*(Znyo(j,l)-dZnyoRo(j,l));
                    M(5,4) = 2*Wavenumber(j,l)*n*(Wnyo(j,l)-dWnyoRo(j,l));
                    M(5,5) = 2*dZnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Znyo(j,l);
                    M(5,6) = 2*dWnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Wnyo(j,l);
                    M(6,1) = 2*Wavenumber(j,l)*dZnxoRo(j,l);
                    M(6,2) = 2*Wavenumber(j,l)*dWnxoRo(j,l);
                    M(6,3) = (y2(j,l)-k2(j,l))*dZnyoRo(j,l);
                    M(6,4) = (y2(j,l)-k2(j,l))*dWnyoRo(j,l);
                    M(6,5) = -Wavenumber(j,l)*n*Znyo(j,l);
                    M(6,6) = -Wavenumber(j,l)*n*Wnyo(j,l);
                    M(7,1) = dZnxoRo(j,l);
                    M(7,2) = dWnxoRo(j,l);
                    M(7,3) = -Wavenumber(j,l)*dZnyoRo(j,l);
                    M(7,4) = -Wavenumber(j,l)*dWnyoRo(j,l);
                    M(7,5) = -n*Znyo(j,l);
                    M(7,6) = -n*Wnyo(j,l);
                    M(7,7) = dHnzoRo(j,l);
                elseif ~ToggleInnerFluid && ~ToggleOuterFluid
                    M(1,1) = dZnxiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-n2)*Znxi(j,l);
                    M(1,2) = dWnxiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-n2)*Wnxi(j,l);
                    M(1,3) = -Wavenumber(j,l)*(dZnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Znyi(j,l));
                    M(1,4) = -Wavenumber(j,l)*(dWnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Wnyi(j,l));
                    M(1,5) = n*(dZnyiRi(j,l)-Znyi(j,l));
                    M(1,6) = n*(dWnyiRi(j,l)-Wnyi(j,l));
                    M(2,1) = 2*n*(dZnxiRi(j,l)-Znxi(j,l));
                    M(2,2) = 2*n*(dWnxiRi(j,l)-Wnxi(j,l));
                    M(2,3) = 2*Wavenumber(j,l)*n*(Znyi(j,l)-dZnyiRi(j,l));
                    M(2,4) = 2*Wavenumber(j,l)*n*(Wnyi(j,l)-dWnyiRi(j,l));
                    M(2,5) = 2*dZnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Znyi(j,l);
                    M(2,6) = 2*dWnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Wnyi(j,l);
                    M(3,1) = 2*Wavenumber(j,l)*dZnxiRi(j,l);
                    M(3,2) = 2*Wavenumber(j,l)*dWnxiRi(j,l);
                    M(3,3) = (y2(j,l)-k2(j,l))*dZnyiRi(j,l);
                    M(3,4) = (y2(j,l)-k2(j,l))*dWnyiRi(j,l);
                    M(3,5) = -Wavenumber(j,l)*n*Znyi(j,l);
                    M(3,6) = -Wavenumber(j,l)*n*Wnyi(j,l);
                    M(4,1) = dZnxoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-n2)*Znxo(j,l);
                    M(4,2) = dWnxoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-n2)*Wnxo(j,l);
                    M(4,3) = -Wavenumber(j,l)*(dZnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Znyo(j,l));
                    M(4,4) = -Wavenumber(j,l)*(dWnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Wnyo(j,l));
                    M(4,5) = n*(dZnyoRo(j,l)-Znyo(j,l));
                    M(4,6) = n*(dWnyoRo(j,l)-Wnyo(j,l));
                    M(5,1) = 2*n*(dZnxoRo(j,l)-Znxo(j,l));
                    M(5,2) = 2*n*(dWnxoRo(j,l)-Wnxo(j,l));
                    M(5,3) = 2*Wavenumber(j,l)*n*(Znyo(j,l)-dZnyoRo(j,l));
                    M(5,4) = 2*Wavenumber(j,l)*n*(Wnyo(j,l)-dWnyoRo(j,l));
                    M(5,5) = 2*dZnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Znyo(j,l);
                    M(5,6) = 2*dWnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Wnyo(j,l);
                    M(6,1) = 2*Wavenumber(j,l)*dZnxoRo(j,l);
                    M(6,2) = 2*Wavenumber(j,l)*dWnxoRo(j,l);
                    M(6,3) = (y2(j,l)-k2(j,l))*dZnyoRo(j,l);
                    M(6,4) = (y2(j,l)-k2(j,l))*dWnyoRo(j,l);
                    M(6,5) = -Wavenumber(j,l)*n*Znyo(j,l);
                    M(6,6) = -Wavenumber(j,l)*n*Wnyo(j,l);
                end
                Y(j,l) = abs(det(M));
            end
        end
    end
end