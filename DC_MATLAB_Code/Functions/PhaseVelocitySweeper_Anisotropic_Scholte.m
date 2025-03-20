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
function [HSScholte,HAScholte,HBScholte] = PhaseVelocitySweeper_Anisotropic_Scholte(c,Delta,Material,Viscoelastic,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Frequency,I,I1,Pattern,SuperLayerSize,LayerThicknesses,Symmetric,SymmetricSystem,Decoupled)
Resolution = 1e-6; % (m/s)
VelocityStep = 2; % (m/s); for coarse sweep
VelocityStepFine = .5; % (m/s); for first iteration in the fine sweeps
EnlargementReal = 3; % SweepRangeReal = XRough(p,1)+EnlargementReal*VelocityStep*[-1 1];
EnlargementImag = 2; % SweepRangeImag = XRough(p,2)+EnlargementImag*VelocityStep*[-1 1];
ImagRealThreshold1 = .1; % exclude wavenumbers with ratios of imaginary part/real part greater than this threshold
ImagRealThreshold2 = .5; % for above slow fluid velocity
DensityThreshold = .01; % rho_fluid/rho_solid
% disp(1+2*VelocityStep/VelocityStepFine*[EnlargementReal EnlargementImag]) % array size of first iteration in the fine sweeps

%#ok<*AGROW>
HSScholte=[];HAScholte=[];HBScholte=[];
if  ToggleUpperFluid && ToggleLowerFluid && ~strcmp(UpperFluid.Name,LowerFluid.Name)
    if  UpperFluid.Velocity > LowerFluid.Velocity
        FastFluidVelocity = UpperFluid.Velocity;
        SlowFluidVelocity = LowerFluid.Velocity;
        Density_ = LowerFluid.Density/Material{end}.Density;
    else
        FastFluidVelocity = LowerFluid.Velocity;
        SlowFluidVelocity = UpperFluid.Velocity;
        Density_ = UpperFluid.Density/Material{1}.Density;
    end
elseif ToggleUpperFluid && (~ToggleLowerFluid || (ToggleLowerFluid && strcmp(UpperFluid.Name,LowerFluid.Name)))
    SlowFluidVelocity = UpperFluid.Velocity;
    Density_ = UpperFluid.Density/Material{1}.Density;
    FastFluidVelocity = 0;
elseif ~ToggleUpperFluid && ToggleLowerFluid
    SlowFluidVelocity = LowerFluid.Velocity;
    Density_ = LowerFluid.Density/Material{end}.Density;
    FastFluidVelocity = 0;
end
if  Density_ > DensityThreshold
    RangeImag = .2;
else
    RangeImag = .04;
end
XS = cell(1,2);
XA = cell(1,2);
AngularFrequency = 2*pi*Frequency*1e3;
AngularFrequency2 = AngularFrequency^2;
kUpperFluid2 = AngularFrequency2/UpperFluid.Velocity^2;
kLowerFluid2 = AngularFrequency2/LowerFluid.Velocity^2;
gUpperFluid = 1i*UpperFluid.Density*AngularFrequency2;
gLowerFluid = 1i*LowerFluid.Density*AngularFrequency2;
for m = 1:SuperLayerSize
    rw2(m) = Material{m}.Density*AngularFrequency2;
    if  ~Decoupled
        r2w4 = rw2(m)^2;
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        b12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))*rw2(m)/Delta(m);
        b22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))*rw2(m)/Delta(m);
        b23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))*r2w4/Delta(m);
        b32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))*rw2(m)/Delta(m);
        b33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))*r2w4/Delta(m);
        b34(m) = -rw2(m)^3/Delta(m);
        A1=0;
    else
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a31(m) = c{m}(1,1)*c{m}(5,5);
        b22(m) = -(c{m}(3,3)+c{m}(5,5))*rw2(m);
        b32(m) = -(c{m}(1,1)+c{m}(5,5))*rw2(m);
        b33(m) = rw2(m)^2;
        a11=0;b12=0;b23=0;b34=0;
    end
end
for r = 1:2
    XSRough = [];
    XARough = [];
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
    if  Symmetric
        if  ~Decoupled
            [YS,YA] = Computer_Coupled(0,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
        else
            [YS,YA] = Computer_Decoupled(0,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33);
        end
        YA(isinf(YA)) = NaN;
    else
        if  ~Decoupled
            YS = Computer_Coupled(0,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
        else
            YS = Computer_Decoupled(0,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33);
        end
    end
    YS(isinf(YS)) = NaN;
    if  ~Viscoelastic && r == 1
        for j = 2:size(YS,1)-1
            if  abs(YS(j)) < abs(YS(j-1)) && abs(YS(j)) < abs(YS(j+1)) && abs(diff(angle([YS(j-1) YS(j+1)]))) > pi/2
                XSRough(end+1,:) = [SweepRangeReal(j) SweepRangeImag]; % real velocity (m/s), imaginary velocity (m/s)
            end
            if  Symmetric && abs(YA(j)) < abs(YA(j-1)) && abs(YA(j)) < abs(YA(j+1)) && abs(diff(angle([YA(j-1) YA(j+1)]))) > pi/2
                XARough(end+1,:) = [SweepRangeReal(j) SweepRangeImag];
            end
        end
% figure,plot(SweepRangeReal,20*log10(abs(YS)))
% figure,plot(SweepRangeReal,20*log10(abs(YA)))
    else
        YS = [max(YS,[],'all')*ones(height(YS),1) YS]; % CAUTION: XRough is at SweepRangeImag(l-1) instead of SweepRangeImag(l) if wall is added at the START of SweepRangeImag instead of at the END!
        if  Symmetric
            YA = [max(YA,[],'all')*ones(height(YA),1) YA];
        end
        for l = 2:size(YS,2)-1
            for j = 2:size(YS,1)-1
                if  abs(YS(j,l)) < abs(YS(j-1,l)) && abs(YS(j,l)) < abs(YS(j+1,l)) && abs(YS(j,l)) < abs(YS(j,l-1)) && abs(YS(j,l)) < abs(YS(j,l+1))
                    if  l == 2
                        YangleDiff = abs(diff(angle([YS(j-1,l) YS(j+1,l)])));
                    else
                        YangleDiff = abs(diff(angle([YS(j-1,l) YS(j,l-1);YS(j+1,l) YS(j,l+1)])));
                    end
                    if  any(YangleDiff > pi/2)
                        XSRough(end+1,:) = [SweepRangeReal(j) SweepRangeImag(l-1)]; % real velocity (m/s), imaginary velocity (m/s)
                    end
                end
                if  Symmetric && abs(YA(j,l)) < abs(YA(j-1,l)) && abs(YA(j,l)) < abs(YA(j+1,l)) && abs(YA(j,l)) < abs(YA(j,l-1)) && abs(YA(j,l)) < abs(YA(j,l+1))
                    if  l == 2
                        YangleDiff = abs(diff(angle([YA(j-1,l) YA(j+1,l)])));
                    else
                        YangleDiff = abs(diff(angle([YA(j-1,l) YA(j,l-1);YA(j+1,l) YA(j,l+1)])));
                    end
                    if  any(YangleDiff > pi/2)
                        XARough(end+1,:) = [SweepRangeReal(j) SweepRangeImag(l-1)]; % real velocity (m/s), imaginary velocity (m/s)
                    end
                end
            end
        end
% figure;surf([-SweepRangeImag(2) SweepRangeImag],SweepRangeReal,20*log10(abs(YS)),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
% figure;surf([-SweepRangeImag(2) SweepRangeImag],SweepRangeReal,20*log10(abs(YA)),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
    end
% XSRough
% XARough
    if  ~isempty(XSRough)
        if  r == 1
            XSRough(end+1,:) = [SlowFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        elseif r == 2
            XSRough(end+1,:) = [FastFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        end
    else
        if  r == 1
            XSRough = [SlowFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        elseif r == 2
            XSRough = [FastFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        end
    end
    if  Symmetric
        if  ~isempty(XARough)
            if  r == 1
                XARough(end+1,:) = [SlowFluidVelocity-EnlargementReal*VelocityStep -1e-10];
            elseif r == 2
                XARough(end+1,:) = [FastFluidVelocity-EnlargementReal*VelocityStep -1e-10];
            end
        else
            if  r == 1
                XARough = [SlowFluidVelocity-EnlargementReal*VelocityStep -1e-10];
            elseif r == 2
                XARough = [FastFluidVelocity-EnlargementReal*VelocityStep -1e-10];
            end
        end
    end
    XS{r} = Converger(1,r,XS{r},XSRough,FastFluidVelocity,SlowFluidVelocity,Resolution,EnlargementReal,EnlargementImag,VelocityStep,VelocityStepFine,Decoupled,Wavenumber,AngularFrequency,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a11,a21,a31,b12,b22,b23,b32,b33,b34);
    if  Symmetric
        XA{r} = Converger(2,r,XA{r},XARough,FastFluidVelocity,SlowFluidVelocity,Resolution,EnlargementReal,EnlargementImag,VelocityStep,VelocityStepFine,Decoupled,Wavenumber,AngularFrequency,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a11,a21,a31,b12,b22,b23,b32,b33,b34);
    end
    if  ~(ToggleUpperFluid && ToggleLowerFluid && ~strcmp(UpperFluid.Name,LowerFluid.Name))
        break
    end
end
XS = [XS{1};XS{2}];
XA = [XA{1};XA{2}];
if  isempty(XS) && isempty(XA)
    disp('No higher order Scholte modes found!')
    return
end
String = ['cp @ ',num2str(Frequency),' kHz:',newline,'Mode        cp(m/ms)  alpha(Np/m)'];
if  Symmetric
    HSScholte = XS(:,1:4);
    HAScholte = XA(:,1:4);
    if  any(HAScholte)
        for i = 1:height(HAScholte)
            if  i < 11
                String = append(String,newline,'AScholte',num2str(i-1),'    ',num2str(HAScholte(i,1)/1e3,'%.5f'),'   ',num2str(HAScholte(i,2)));
            else
                String = append(String,newline,'AScholte',num2str(i-1),'   ',num2str(HAScholte(i,1)/1e3,'%.5f'),'   ',num2str(HAScholte(i,2)));
            end
        end
    end
    String = append(String,newline);
    if  any(HSScholte)
        for i = 1:height(HSScholte)
            if  i < 11
                String = append(String,newline,'SScholte',num2str(i-1),'    ',num2str(HSScholte(i,1)/1e3,'%.5f'),'   ',num2str(HSScholte(i,2)));
            else
                String = append(String,newline,'SScholte',num2str(i-1),'   ',num2str(HSScholte(i,1)/1e3,'%.5f'),'   ',num2str(HSScholte(i,2)));
            end
        end
    end
    String = append(String,newline,newline);
    if  any(HAScholte)
        String = append(String,'AScholte: ',num2str(height(HAScholte)));
    end
    String = append(String,newline);
    if  any(HSScholte)
        String = append(String,'SScholte: ',num2str(height(HSScholte)));
    end
else
    HBScholte = XS(:,1:4);
    if  any(HBScholte)
        for i = 1:height(HBScholte)
            if  i < 11
                String = append(String,newline,'BScholte',num2str(i-1),'    ',num2str(HBScholte(i,1)/1e3,'%.5f'),'   ',num2str(HBScholte(i,2)));
            else
                String = append(String,newline,'BScholte',num2str(i-1),'   ',num2str(HBScholte(i,1)/1e3,'%.5f'),'   ',num2str(HBScholte(i,2)));
            end
        end
    end
    String = append(String,newline,newline);
    if  any(HBScholte)
        String = append(String,'BScholte: ',num2str(height(HBScholte)));
    end
end
disp([String,newline,'----------------'])
end
function X = Converger(ModeFamily,r,X,XRough,FastFluidVelocity,SlowFluidVelocity,Resolution,EnlargementReal,EnlargementImag,VelocityStep,VelocityStepFine,Decoupled,Wavenumber,AngularFrequency,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a11,a21,a31,b12,b22,b23,b32,b33,b34)
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
                SweepRangeReal = SweepRangeReal(1):(SweepRangeReal(end)-SweepRangeReal(1))/4:SweepRangeReal(end);
                SweepRangeImag = SweepRangeImag(1):(SweepRangeImag(end)-SweepRangeImag(1))/4:SweepRangeImag(end);
            end
            SweepRangeImag(SweepRangeImag > -1e-10) = [];
            Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
            if  Symmetric
                if  ModeFamily == 1
                    if  ~Decoupled
                        Y = Computer_Coupled(ModeFamily,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
                    else
                        Y = Computer_Decoupled(ModeFamily,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33);
                    end
                elseif ModeFamily == 2
                    if  ~Decoupled
                        [~,Y] = Computer_Coupled(ModeFamily,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
                    else
                        [~,Y] = Computer_Decoupled(ModeFamily,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33);
                    end
                end
            else
                if  ~Decoupled
                    Y = Computer_Coupled(ModeFamily,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
                else
                    Y = Computer_Decoupled(ModeFamily,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33);
                end
            end
            if  ~isempty(X)
                for j = 2:height(Wavenumber)-1
                    for m = 1:height(X)
                        if  SweepRangeReal(j-1) < X(m,3) && SweepRangeReal(j+1) > X(m,3)
                            Y(j,:) = NaN;
                        end
                    end
                end
            end
            Y(:,end+1) = max(Y,[],'all'); % add wall below zero attenuation to allow finding minimum at zero attenuation
% if  p >= 17 && k <= 1
% f = figure;surf([SweepRangeImag SweepRangeImag(end)-(SweepRangeImag(1)-SweepRangeImag(2))],SweepRangeReal,20*log10(abs(Y)))
% close(f)
% end
            Min = zeros(size(Y));
            for l = 2:size(Y,2)-1
                for j = 2:size(Y,1)-1
                    if  abs(Y(j,l)) < abs(Y(j-1,l)) && abs(Y(j,l)) < abs(Y(j+1,l)) && abs(Y(j,l)) < abs(Y(j,l-1)) && abs(Y(j,l)) < abs(Y(j,l+1))
                        Min(j,l) = 1;
                    end
                end
            end
            [b1,b2] = find(Min);
            if  ~isempty(b1)
                for l = 1:length(b1)
                    if  l == 1
                        MIN = [b1(1) b2(1)];
                    elseif k == 1 && l > 1 && ~any(all([SweepRangeReal(b1(l)) SweepRangeImag(b2(l))] == XRough,2)) % add minima only at k == 1 because there can be instability minima when we converge deeply into the mimima; add only unique mimima
                        XRough(end+1,:) = [SweepRangeReal(b1(l)) SweepRangeImag(b2(l))]; % real velocity (m/s), imaginary velocity (m/s)
                    end
                end
            else
                MIN = 0;
                break
            end
            if  k == 1
                if  MIN(2) == length(SweepRangeImag)
                    YangleDiff = abs(diff(angle([Y(MIN(1)-1,MIN(2)) Y(MIN(1)+1,MIN(2))])));
                else
                    YangleDiff = abs(diff(angle([Y(MIN(1)-1,MIN(2)) Y(MIN(1),MIN(2)-1);Y(MIN(1)+1,MIN(2)) Y(MIN(1),MIN(2)+1)])));
                end
                if  all(YangleDiff < pi/2)
                    MIN = 0;
                    break
                end
            elseif k == 100 || Resolution > abs(SweepRangeReal(1)-SweepRangeReal(2))
                break
            end
            SweepRangeReal = SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))*[-1 1];
            SweepRangeImag = SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))*[-1 1];
        end
        if  any(MIN)
            X(end+1,:) = [AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) imag(Wavenumber(MIN(1),MIN(2))) SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2)) Wavenumber(MIN(1),MIN(2))]; % phase velocity (m/s), attenuation (Np/m), real velocity (m/s), imaginary velocity (m/s), wavenumber (1/m)
        end
        p = p+1;
    end
    if  ~isempty(X)
        X = sortrows(X,1);
    end
end
% function [YS,YA] = Computer_Coupled(ModeFamily,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34)
%     Size = size(Wavenumber);
%     Wavenumber = reshape(Wavenumber,1,1,[]);
%     Wavenumber2 = Wavenumber.^2;
%     Wavenumber4 = Wavenumber2.^2;
%     Wavenumber6 = Wavenumber2.^3;
%     Length = length(Wavenumber);
%     YS = NaN(size(Wavenumber)); % to be discarded once pagedet exists
%     YA = YS;
%     for m = 1:SuperLayerSize
%         A1 = a11(m)*Wavenumber2+b12(m);
%         A2 = a21(m)*Wavenumber4+b22(m)*Wavenumber2+b23(m);
%         A3 = a31(m)*Wavenumber6+b32(m)*Wavenumber4+b33(m)*Wavenumber2+b34(m);
%         d1 = A1/3;
%         d2 = A2/3-d1.^2;
%         d3 = d1.^3-d1.*A2/2+A3/2;
%         d4 = (sqrt(d2.^3+d3.^2)-d3).^(1/3);
%         d5 = d2./d4;
%         d6 = (d5-d4)/2-d1;
%         d7 = (d5+d4)/2i*sqrt(3);
%         k3(1,1,:) = sqrt(d6+d7);
%         k3(1,2,:) = sqrt(d6-d7);
%         k3(1,3,:) = sqrt(d4-d5-d1);
%         k32 = k3.^2;
%         k3k = k3.*Wavenumber;
%         m11 = c{m}(1,1)*Wavenumber2+c{m}(5,5)*k32-rw2(m);
%         m22 = c{m}(6,6)*Wavenumber2+c{m}(4,4)*k32-rw2(m);
%         m12 = c{m}(1,6)*Wavenumber2+c{m}(4,5)*k32;
%         m13 = (c{m}(1,3)+c{m}(5,5))*k3k;
%         m23 = (c{m}(3,6)+c{m}(4,5))*k3k;
%         m1 = m13.*m22-m12.*m23;
%         V = (m11.*m23-m13.*m12)./m1;
%         W = (m11.*m22-m12.^2)./-m1;
%         e1 = Wavenumber.*W+k3;
%         e2 = k3.*V;
%         D3 = 1i*((c{m}(1,3)+c{m}(3,6)*V).*Wavenumber+c{m}(3,3)*k3.*W);
%         D4 = 1i*(c{m}(4,5)*e1+c{m}(4,4)*e2);
%         D5 = 1i*(c{m}(5,5)*e1+c{m}(4,5)*e2);
%         if  Symmetric && SuperLayerSize == 1
%             E = exp(.5i*k3*LayerThicknesses);
%         else
%             E = exp(1i*k3*LayerThicknesses(m));
%         end
%         L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
%         L2 = [ones(1,3,Length) E;V V.*E;W -W.*E;E ones(1,3,Length);V.*E V;W.*E -W];
%         L{m} = pagemrdivide(L1,L2);
%     end
%     M{1} = L{1};
%     for m = 2:SuperLayerSize
%         M0 = L{m}(1:3,1:3,:)-M{1}(4:6,4:6,:);
%         M1 = pagemrdivide(M{1}(1:3,4:6,:),M0);
%         M2 = pagemrdivide(L{m}(4:6,1:3,:),M0);
%         M{1} = [M{1}(1:3,1:3,:)+pagemtimes(M1,M{1}(4:6,1:3,:)) -pagemtimes(M1,L{m}(1:3,4:6,:));pagemtimes(M2,M{1}(4:6,1:3,:)) L{m}(4:6,4:6,:)-pagemtimes(M2,L{m}(1:3,4:6,:))];
%     end
%     for m = 1:length(Pattern)
%         M0 = M{Pattern(m)}(1:3,1:3,:)-M{m}(4:6,4:6,:);
%         M1 = pagemrdivide(M{m}(1:3,4:6,:),M0);
%         M2 = pagemrdivide(M{Pattern(m)}(4:6,1:3,:),M0);
%         M{m+1} = [M{m}(1:3,1:3,:)+pagemtimes(M1,M{m}(4:6,1:3,:)) -pagemtimes(M1,M{Pattern(m)}(1:3,4:6,:));pagemtimes(M2,M{m}(4:6,1:3,:)) M{Pattern(m)}(4:6,4:6,:)-pagemtimes(M2,M{Pattern(m)}(1:3,4:6,:))];
%     end
%     if  Symmetric
%         M{end} = pageinv(M{end});
%         if  ModeFamily <= 1
%             YS = reshape((-sqrt(kUpperFluid2-Wavenumber2)/gUpperFluid+M{end}(3,1,:)-M{end}(3,4,:).*M{end}(6,1,:)./M{end}(6,4,:)),Size);
%         end
%         if  ModeFamily ~= 1
%             YA = reshape((-sqrt(kUpperFluid2-Wavenumber2)/gUpperFluid+M{end}(3,1,:)-M{end}(3,5,:).*M{end}(5,1,:)./M{end}(5,5,:)-(M{end}(3,5,:).*M{end}(5,6,:)./M{end}(5,5,:)-M{end}(3,6,:)).*(M{end}(4,1,:).*M{end}(5,5,:)-M{end}(4,5,:).*M{end}(5,1,:))./(M{end}(4,5,:).*M{end}(5,6,:)-M{end}(4,6,:).*M{end}(5,5,:))),Size);
%         end
%     else
%         if  SymmetricSystem
%             M0 = M{end}(4:6,4:6,:).*I-M{end}(4:6,4:6,:);
%             M1 = pagemrdivide(M{end}(1:3,4:6,:),M0);
%             M2 = pagemrdivide(M{end}(1:3,4:6,:).*I,M0);
%             M{end} = [M{end}(1:3,1:3,:)+pagemtimes(M1,M{end}(4:6,1:3,:)) -pagemtimes(M1,M{end}(4:6,1:3,:).*I);pagemtimes(M2,M{end}(4:6,1:3,:)) M{end}(1:3,1:3,:).*I-pagemtimes(M2,M{end}(4:6,1:3,:).*I)];
%         end
%         M{end} = pageinv(M{end});
%         if  ToggleUpperFluid && ToggleLowerFluid
%             if  r == 1
%                 k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2);
%                 k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2);
%             else
%                 if  UpperFluid.Velocity > LowerFluid.Velocity
%                     k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2);
%                     k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
%                 else
%                     k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
%                     k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2);
%                 end
%             end
%             WUpperFluid = k3UpperFluid./Wavenumber;
%             WLowerFluid = k3LowerFluid./Wavenumber;
%             DUpperFluid = gUpperFluid./Wavenumber;
%             DLowerFluid = gLowerFluid./Wavenumber;
%             YS = reshape(((WLowerFluid-M{end}(6,4,:).*DLowerFluid).*(WUpperFluid+M{end}(3,1,:).*DUpperFluid)+M{end}(3,4,:).*M{end}(6,1,:).*DUpperFluid.*DLowerFluid),Size); 
%         elseif ToggleUpperFluid && ~ToggleLowerFluid
%             YS = reshape((-sqrt(kUpperFluid2-Wavenumber2)/gUpperFluid+M{end}(3,1,:)),Size);
%         elseif ~ToggleUpperFluid && ToggleLowerFluid
%             YS = reshape((-sqrt(kLowerFluid2-Wavenumber2)/gLowerFluid-M{end}(6,4,:)),Size);
%         end
%     end
% end
function [YS,YA] = Computer_Coupled(ModeFamily,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34)
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Wavenumber2 = Wavenumber.^2;
    Wavenumber4 = Wavenumber2.^2;
    Wavenumber6 = Wavenumber2.^3;
    Length = length(Wavenumber);
    YS = NaN(size(Wavenumber)); % to be discarded once pagedet exists
    YA = YS;
    for m = 1:SuperLayerSize
        A1 = a11(m)*Wavenumber2+b12(m);
        A2 = a21(m)*Wavenumber4+b22(m)*Wavenumber2+b23(m);
        A3 = a31(m)*Wavenumber6+b32(m)*Wavenumber4+b33(m)*Wavenumber2+b34(m);
        d1 = A1/3;
        d2 = A2/3-d1.^2;
        d3 = d1.^3-d1.*A2/2+A3/2;
        d4 = (sqrt(d2.^3+d3.^2)-d3).^(1/3);
        d5 = d2./d4;
        d6 = (d5-d4)/2-d1;
        d7 = (d5+d4)/2i*sqrt(3);
        k3(1,1,:) = sqrt(d6+d7);
        k3(1,2,:) = sqrt(d6-d7);
        k3(1,3,:) = sqrt(d4-d5-d1);
        k32 = k3.^2;
        k3k = k3.*Wavenumber;
        m11 = c{m}(1,1)*Wavenumber2+c{m}(5,5)*k32-rw2(m);
        m22 = c{m}(6,6)*Wavenumber2+c{m}(4,4)*k32-rw2(m);
        m12 = c{m}(1,6)*Wavenumber2+c{m}(4,5)*k32;
        m13 = (c{m}(1,3)+c{m}(5,5))*k3k;
        m23 = (c{m}(3,6)+c{m}(4,5))*k3k;
        m1 = m13.*m22-m12.*m23;
        V = (m11.*m23-m13.*m12)./m1;
        W = (m11.*m22-m12.^2)./-m1;
        e1 = Wavenumber.*W+k3;
        e2 = k3.*V;
        D3 = 1i*((c{m}(1,3)+c{m}(3,6)*V).*Wavenumber+c{m}(3,3)*k3.*W);
        D4 = 1i*(c{m}(4,5)*e1+c{m}(4,4)*e2);
        D5 = 1i*(c{m}(5,5)*e1+c{m}(4,5)*e2);
        if  Symmetric && SuperLayerSize == 1
            E = exp(.5i*k3*LayerThicknesses);
        else
            E = exp(1i*k3*LayerThicknesses(m));
        end
        L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
        L2 = [ones(1,3,Length) E;V V.*E;W -W.*E;E ones(1,3,Length);V.*E V;W.*E -W];
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    for m = 2:SuperLayerSize
        M0 = L{m}(1:3,1:3,:)-M{1}(4:6,4:6,:);
        M1 = pagemrdivide(M{1}(1:3,4:6,:),M0);
        M2 = pagemrdivide(L{m}(4:6,1:3,:),M0);
        M{1} = [M{1}(1:3,1:3,:)+pagemtimes(M1,M{1}(4:6,1:3,:)) -pagemtimes(M1,L{m}(1:3,4:6,:));pagemtimes(M2,M{1}(4:6,1:3,:)) L{m}(4:6,4:6,:)-pagemtimes(M2,L{m}(1:3,4:6,:))];
    end
    for m = 1:length(Pattern)
        M0 = M{Pattern(m)}(1:3,1:3,:)-M{m}(4:6,4:6,:);
        M1 = pagemrdivide(M{m}(1:3,4:6,:),M0);
        M2 = pagemrdivide(M{Pattern(m)}(4:6,1:3,:),M0);
        M{m+1} = [M{m}(1:3,1:3,:)+pagemtimes(M1,M{m}(4:6,1:3,:)) -pagemtimes(M1,M{Pattern(m)}(1:3,4:6,:));pagemtimes(M2,M{m}(4:6,1:3,:)) M{Pattern(m)}(4:6,4:6,:)-pagemtimes(M2,M{Pattern(m)}(1:3,4:6,:))];
    end
    if  Symmetric
        WFluid = -sqrt(kUpperFluid2-Wavenumber2)./Wavenumber;
        DFluid = gUpperFluid./Wavenumber;
        if  ModeFamily <= 1
            MS = [M{end}([1:3 5:6],[1:2 4:5],:) [-M{end}(1,3,:).*WFluid-DFluid;-M{end}([2:3 5:6],3,:).*WFluid]];
            for j = 1:length(Wavenumber)
                YS(j) = det(MS(:,:,j));
            end
        end
        if  ModeFamily ~= 1
            MA = [M{end}(1:4,[1:2 6],:) [-M{end}(1,3,:).*WFluid-DFluid;-M{end}(2:4,3,:).*WFluid]];
            for j = 1:length(Wavenumber)
                YA(j) = det(MA(:,:,j));
            end
            YA = reshape(YA,Size);
        end
    else
        if  SymmetricSystem
            M0 = M{end}(4:6,4:6,:).*I-M{end}(4:6,4:6,:);
            M1 = pagemrdivide(M{end}(1:3,4:6,:),M0);
            M2 = pagemrdivide(M{end}(1:3,4:6,:).*I,M0);
            M{end} = [M{end}(1:3,1:3,:)+pagemtimes(M1,M{end}(4:6,1:3,:)) -pagemtimes(M1,M{end}(4:6,1:3,:).*I);pagemtimes(M2,M{end}(4:6,1:3,:)) M{end}(1:3,1:3,:).*I-pagemtimes(M2,M{end}(4:6,1:3,:).*I)];
        end
        if  ToggleUpperFluid && ToggleLowerFluid
            if  r == 1
                k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2);
                k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2);
            else
                if  UpperFluid.Velocity > LowerFluid.Velocity
                    k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                else
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2);
                end
            end
            WUpperFluid = k3UpperFluid./Wavenumber;
            WLowerFluid = k3LowerFluid./Wavenumber;
            DUpperFluid = gUpperFluid./Wavenumber;
            DLowerFluid = gLowerFluid./Wavenumber;
            M{end} = [M{end}(:,[1:2 4:5],:) [-M{end}(1,3,:).*WUpperFluid-DUpperFluid;-M{end}(2:6,3,:).*WUpperFluid] [M{end}(1:3,6,:).*WLowerFluid;M{end}(4,6,:).*WLowerFluid-DLowerFluid;M{end}(5:6,6,:).*WLowerFluid]];
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            WUpperFluid = -sqrt(kUpperFluid2-Wavenumber2)./Wavenumber;
            DUpperFluid = gUpperFluid./Wavenumber;
            M{end} = [M{end}(:,[1:2 4:6],:) [-M{end}(1,3,:).*WUpperFluid-DUpperFluid;-M{end}(2:6,3,:).*WUpperFluid]];
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            WLowerFluid = -sqrt(kLowerFluid2-Wavenumber2)./Wavenumber;
            DLowerFluid = gLowerFluid./Wavenumber;
            M{end} = [M{end}(:,[1:2 4:5],:) -M{end}(:,3,:) [M{end}(1:3,6,:).*WLowerFluid;M{end}(4,6,:).*WLowerFluid-DLowerFluid;M{end}(5:6,6,:).*WLowerFluid]];
        end
        for j = 1:length(Wavenumber)
            YS(j) = det(M{end}(:,:,j));
        end
    end
    YS = reshape(YS,Size);
end
% function [YS,YA] = Computer_Decoupled(ModeFamily,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33)
%     Size = size(Wavenumber);
%     Wavenumber = reshape(Wavenumber,1,1,[]);
%     Wavenumber2 = Wavenumber.^2;
%     Wavenumber4 = Wavenumber2.^2;
%     Length = length(Wavenumber);
%     YS = NaN(size(Wavenumber)); % to be discarded once pagedet exists
%     YA = YS;
%     for m = 1:SuperLayerSize
%         A2 = a21(m)*Wavenumber2+b22(m);
%         A3 = a31(m)*Wavenumber4+b32(m)*Wavenumber2+b33(m);
%         d1 = sqrt(A2.^2-2*A1(m)*A3);
%         k3(1,1,:) = sqrt((-A2+d1)/A1(m));
%         k3(1,2,:) = sqrt((-A2-d1)/A1(m));
%         W = (rw2(m)-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k3);
%         D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*k3.*W);
%         D5 = 1i*(c{m}(5,5)*(k3+Wavenumber.*W));
%         if  Symmetric && SuperLayerSize == 1
%             E = exp(.5i*k3*LayerThicknesses);
%         else
%             E = exp(1i*k3*LayerThicknesses(m));
%         end
%         L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
%         L2 = [ones(1,2,Length) E;W -W.*E;E ones(1,2,Length);W.*E -W];
%         L{m} = pagemrdivide(L1,L2);
%     end
%     M{1} = L{1};
%     for m = 2:SuperLayerSize
%         M0 = L{m}(1:2,1:2,:)-M{1}(3:4,3:4,:);
%         M1 = pagemrdivide(M{1}(1:2,3:4,:),M0);
%         M2 = pagemrdivide(L{m}(3:4,1:2,:),M0);
%         M{1} = [M{1}(1:2,1:2,:)+pagemtimes(M1,M{1}(3:4,1:2,:)) -pagemtimes(M1,L{m}(1:2,3:4,:));pagemtimes(M2,M{1}(3:4,1:2,:)) L{m}(3:4,3:4,:)-pagemtimes(M2,L{m}(1:2,3:4,:))];
%     end
%     for m = 1:length(Pattern)
%         M0 = M{Pattern(m)}(1:2,1:2,:)-M{m}(3:4,3:4,:);
%         M1 = pagemrdivide(M{m}(1:2,3:4,:),M0);
%         M2 = pagemrdivide(M{Pattern(m)}(3:4,1:2,:),M0);
%         M{m+1} = [M{m}(1:2,1:2,:)+pagemtimes(M1,M{m}(3:4,1:2,:)) -pagemtimes(M1,M{Pattern(m)}(1:2,3:4,:));pagemtimes(M2,M{m}(3:4,1:2,:)) M{Pattern(m)}(3:4,3:4,:)-pagemtimes(M2,M{Pattern(m)}(1:2,3:4,:))];
%     end
%     if  Symmetric
%         M{end} = pageinv(M{end});
%         if  ModeFamily <= 1
%             YS = reshape((-sqrt(kUpperFluid2-Wavenumber2)/gUpperFluid+M{end}(2,1,:)-M{end}(2,3,:).*M{end}(4,1,:)./M{end}(4,3,:)),Size);
%         end
%         if  ModeFamily ~= 1
%             YA = reshape((-sqrt(kUpperFluid2-Wavenumber2)/gUpperFluid+M{end}(2,1,:)-M{end}(2,4,:).*M{end}(3,1,:)./M{end}(3,4,:)),Size);
%         end
%     else
%         if  SymmetricSystem
%             M0 = M{end}(3:4,3:4,:).*I1-M{end}(3:4,3:4,:);
%             M1 = pagemrdivide(M{end}(1:2,3:4,:),M0);
%             M2 = pagemrdivide(M{end}(1:2,3:4,:).*I1,M0);
%             M{end} = [M{end}(1:2,1:2,:)+pagemtimes(M1,M{end}(3:4,1:2,:)) -pagemtimes(M1,M{end}(3:4,1:2,:).*I1);pagemtimes(M2,M{end}(3:4,1:2,:)) M{end}(1:2,1:2,:).*I1-pagemtimes(M2,M{end}(3:4,1:2,:).*I1)];
%         end
%         M{end} = pageinv(M{end});
%         if  ToggleUpperFluid && ToggleLowerFluid
%             if  r == 1
%                 k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2);
%                 k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2);
%             else
%                 if  UpperFluid.Velocity > LowerFluid.Velocity
%                     k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2);
%                     k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
%                 else
%                     k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
%                     k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2);
%                 end
%             end
%             WUpperFluid = k3UpperFluid./Wavenumber;
%             WLowerFluid = k3LowerFluid./Wavenumber;
%             DUpperFluid = gUpperFluid./Wavenumber;
%             DLowerFluid = gLowerFluid./Wavenumber;
%             YS = reshape(((WLowerFluid-M{end}(4,3,:).*DLowerFluid).*(WUpperFluid+M{end}(2,1,:).*DUpperFluid)+M{end}(2,3,:).*M{end}(4,1,:).*DUpperFluid.*DLowerFluid),Size);
%         elseif ToggleUpperFluid && ~ToggleLowerFluid
%             YS = reshape((-sqrt(kUpperFluid2-Wavenumber2)/gUpperFluid+M{end}(2,1,:)),Size);
%         elseif ~ToggleUpperFluid && ToggleLowerFluid
%             YS = reshape((-sqrt(kLowerFluid2-Wavenumber2)/gLowerFluid-M{end}(4,3,:)),Size);
%         end
%     end
% end
function [YS,YA] = Computer_Decoupled(ModeFamily,r,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33)
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Wavenumber2 = Wavenumber.^2;
    Wavenumber4 = Wavenumber2.^2;
    Length = length(Wavenumber);
    YS = NaN(size(Wavenumber)); % to be discarded once pagedet exists
    YA = YS;
    for m = 1:SuperLayerSize
        A2 = a21(m)*Wavenumber2+b22(m);
        A3 = a31(m)*Wavenumber4+b32(m)*Wavenumber2+b33(m);
        d1 = sqrt(A2.^2-2*A1(m)*A3);
        k3(1,1,:) = sqrt((-A2+d1)/A1(m));
        k3(1,2,:) = sqrt((-A2-d1)/A1(m));
        W = (rw2(m)-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k3);
        D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*k3.*W);
        D5 = 1i*(c{m}(5,5)*(k3+Wavenumber.*W));
        if  Symmetric && SuperLayerSize == 1
            E = exp(.5i*k3*LayerThicknesses);
        else
            E = exp(1i*k3*LayerThicknesses(m));
        end
        L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
        L2 = [ones(1,2,Length) E;W -W.*E;E ones(1,2,Length);W.*E -W];
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    for m = 2:SuperLayerSize
        M0 = L{m}(1:2,1:2,:)-M{1}(3:4,3:4,:);
        M1 = pagemrdivide(M{1}(1:2,3:4,:),M0);
        M2 = pagemrdivide(L{m}(3:4,1:2,:),M0);
        M{1} = [M{1}(1:2,1:2,:)+pagemtimes(M1,M{1}(3:4,1:2,:)) -pagemtimes(M1,L{m}(1:2,3:4,:));pagemtimes(M2,M{1}(3:4,1:2,:)) L{m}(3:4,3:4,:)-pagemtimes(M2,L{m}(1:2,3:4,:))];
    end
    for m = 1:length(Pattern)
        M0 = M{Pattern(m)}(1:2,1:2,:)-M{m}(3:4,3:4,:);
        M1 = pagemrdivide(M{m}(1:2,3:4,:),M0);
        M2 = pagemrdivide(M{Pattern(m)}(3:4,1:2,:),M0);
        M{m+1} = [M{m}(1:2,1:2,:)+pagemtimes(M1,M{m}(3:4,1:2,:)) -pagemtimes(M1,M{Pattern(m)}(1:2,3:4,:));pagemtimes(M2,M{m}(3:4,1:2,:)) M{Pattern(m)}(3:4,3:4,:)-pagemtimes(M2,M{Pattern(m)}(1:2,3:4,:))];
    end
    if  Symmetric
        WFluid = -sqrt(kUpperFluid2-Wavenumber2)./Wavenumber;
        DFluid = gUpperFluid./Wavenumber;
        if  ModeFamily <= 1
            MS = [M{end}([1:2 4],[1 3],:) [-M{end}(1,2,:).*WFluid-DFluid;-M{end}([2 4],2,:).*WFluid]];
            for j = 1:length(Wavenumber)
                YS(j) = det(MS(:,:,j));
            end
        end
        if  ModeFamily ~= 1
            MA = [M{end}(1:3,[1 4],:) [-M{end}(1,2,:).*WFluid-DFluid;-M{end}(2:3,2,:).*WFluid]];
            for j = 1:length(Wavenumber)
                YA(j) = det(MA(:,:,j));
            end
            YA = reshape(YA,Size);
        end
    else
        if  SymmetricSystem
            M0 = M{end}(3:4,3:4,:).*I1-M{end}(3:4,3:4,:);
            M1 = pagemrdivide(M{end}(1:2,3:4,:),M0);
            M2 = pagemrdivide(M{end}(1:2,3:4,:).*I1,M0);
            M{end} = [M{end}(1:2,1:2,:)+pagemtimes(M1,M{end}(3:4,1:2,:)) -pagemtimes(M1,M{end}(3:4,1:2,:).*I1);pagemtimes(M2,M{end}(3:4,1:2,:)) M{end}(1:2,1:2,:).*I1-pagemtimes(M2,M{end}(3:4,1:2,:).*I1)];
        end
        if  ToggleUpperFluid && ToggleLowerFluid
            if  r == 1
                k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2);
                k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2);
            else
                if  UpperFluid.Velocity > LowerFluid.Velocity
                    k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                else
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2);
                end
            end
            WUpperFluid = k3UpperFluid./Wavenumber;
            WLowerFluid = k3LowerFluid./Wavenumber;
            DUpperFluid = gUpperFluid./Wavenumber;
            DLowerFluid = gLowerFluid./Wavenumber;
            M{end} = [M{end}(:,[1 3],:) [-M{end}(1,2,:).*WUpperFluid-DUpperFluid;-M{end}(2:4,2,:).*WUpperFluid] [M{end}(1:2,4,:).*WLowerFluid;M{end}(3,4,:).*WLowerFluid-DLowerFluid;M{end}(4,4,:).*WLowerFluid]];
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            WUpperFluid = -sqrt(kUpperFluid2-Wavenumber2)./Wavenumber;
            DUpperFluid = gUpperFluid./Wavenumber;
            M{end} = [M{end}(:,[1 3:4],:) [-M{end}(1,2,:).*WUpperFluid-DUpperFluid;-M{end}(2:4,2,:).*WUpperFluid]];
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            WLowerFluid = -sqrt(kLowerFluid2-Wavenumber2)./Wavenumber;
            DLowerFluid = gLowerFluid./Wavenumber;
            M{end} = [M{end}(:,[1 3],:) -M{end}(:,2,:) [M{end}(1:2,4,:).*WLowerFluid;M{end}(3,4,:).*WLowerFluid-DLowerFluid;M{end}(4,4,:).*WLowerFluid]];
        end
        for j = 1:length(Wavenumber)
            YS(j) = det(M{end}(:,:,j));
        end
    end
    YS = reshape(YS,Size);
end