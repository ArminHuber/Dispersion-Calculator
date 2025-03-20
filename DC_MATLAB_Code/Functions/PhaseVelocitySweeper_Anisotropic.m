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
function [FSLamb,FALamb,FAScholte] = PhaseVelocitySweeper_Anisotropic(c,Delta,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Frequency,I,I1,Pattern,SuperLayerSize,LayerThicknesses,Symmetric,SymmetricSystem,Decoupled,XS0)
Resolution = 1e-6; % (m/s)
DensityThreshold = .01; % rho_fluid/rho_solid

%#ok<*AGROW>
AngularFrequency = 2*pi*Frequency*1e3;
AngularFrequency2 = AngularFrequency^2;
kUpperFluid2 = AngularFrequency2/UpperFluid.Velocity^2;
kLowerFluid2 = AngularFrequency2/LowerFluid.Velocity^2;
gUpperFluid = 1i*UpperFluid.Density*AngularFrequency2;
gLowerFluid = 1i*LowerFluid.Density*AngularFrequency2;
if  ToggleUpperFluid && ToggleLowerFluid
    Density_ = .5*(UpperFluid.Density/Material{1}.Density+LowerFluid.Density/Material{end}.Density);
elseif ToggleUpperFluid && ~ToggleLowerFluid
    Density_ = .5*UpperFluid.Density/Material{1}.Density;
elseif ~ToggleUpperFluid && ToggleLowerFluid
    Density_ = .5*LowerFluid.Density/Material{end}.Density;
end

Range1 = .02; % for the weakly damped solutions near the real axis A0Scholte/B0Scholte
if  FluidLoading && Density_ > DensityThreshold % for the strongly damped solutions A0/B0 
    Range2 = 1.8;
else
    Range2 = .1;
end
Range3 = .05; % for the weakly damped solutions near the real axis S0/B1 and S1/B2
Steps1 = 50;
Steps2 = 50;
Steps3 = 50;

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

% rough search for undamped A0Scholte/B0Scholte mode on the real axis
SweepRange = 1:1000;
Wavenumber = AngularFrequency./SweepRange;
if  ~Decoupled
    Y = Computer_Coupled('A','Solid','SMM',Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34); 
else
    Y = Computer_Decoupled('A','Solid','SMM',Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33);
end
for i = 2:length(SweepRange)-1
    if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
        XRough = SweepRange(i);
        break
    end
end
% figure,plot(SweepRange,20*log10(Y))

% fine search for undamped A0Scholte/B0Scholte mode on the real axis
Bisections = ceil(log2(Resolution/(abs(SweepRange(1)-SweepRange(2))))/log2(2/4));
PhaseVelocity = XRough+(SweepRange(2)-SweepRange(1))*[-1 1];
for o = 1:Bisections
    PhaseVelocity = PhaseVelocity(1):(PhaseVelocity(end)-PhaseVelocity(1))/4:PhaseVelocity(end);
    Wavenumber = AngularFrequency./PhaseVelocity;
    if  ~Decoupled
        Y = Computer_Coupled('A','Solid','SMM',Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
    else
        Y = Computer_Decoupled('A','Solid','SMM',Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33);
    end
    for i = 2:length(PhaseVelocity)-1
        if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
            if  o == Bisections
                FALamb_UD = PhaseVelocity(i);
            end
            PhaseVelocity = [PhaseVelocity(i-1) PhaseVelocity(i+1)];
            break
        end
    end
end

% search for the weakly damped Scholte mode near the real axis
if  Viscoelastic

    % rough search for the weakly damped Scholte mode near the real axis
    XRough = [];
    SweepRangeReal = (1-.5*Range1)*FALamb_UD:Range1*FALamb_UD/Steps1:(1+.5*Range1)*FALamb_UD;
    SweepRangeImag = -Range1*FALamb_UD:Range1*FALamb_UD/Steps1:0;
    Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
    if  ~Decoupled
        Y = Computer_Coupled('A','Fluid2','TMM',Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
    else
        Y = Computer_Decoupled('A','Fluid2','TMM',Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33);
    end
    if  ~FluidLoading
        Y(:,end+1) = max(Y,[],'all');
    end
% if  ~FluidLoading
% figure,surf([SweepRangeImag SweepRangeImag(end)-(SweepRangeImag(1)-SweepRangeImag(2))],SweepRangeReal,20*log10(Y))
% else
% figure,surf(SweepRangeImag,SweepRangeReal,20*log10(Y))
% end
    for l = 2:size(Y,2)-1
        for j = 2:size(Y,1)-1
            if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                XRough = [SweepRangeReal(j) SweepRangeImag(l)]; % real velocity (m/s), imaginary velocity (m/s)
                break
            end
        end
    end

    % fine search for the weakly damped Scholte mode near the real axis
    if  ~isempty(XRough)
        XFine = Converger('A','Fluid2','TMM',Decoupled,XRough,Resolution,AngularFrequency,SweepRangeReal,SweepRangeImag,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a11,a21,a31,b12,b22,b23,b32,b33,b34);
    else
        XFine = [FALamb_UD 0 FALamb_UD 0]; 
    end
else
    XFine = [FALamb_UD 0 FALamb_UD 0];
end

% search for the strongly damped A0/B0 when fluid-loading is present
if  FluidLoading

    % rough search for the strongly damped A0/B0 when fluid-loading is present
    XRough = [];
    for o = [1 2 5]
        SweepRangeReal = (1-.5*Range2)*FALamb_UD:Range2*FALamb_UD/Steps2/o:(1+.5*Range2)*FALamb_UD;
        SweepRangeImag = -Range2*FALamb_UD:Range2*FALamb_UD/Steps2/o:0;
        Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
        if  ~Decoupled
            Y = Computer_Coupled('A','Solid','TMM',Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
        else
            Y = Computer_Decoupled('A','Solid','TMM',Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33);
        end
% figure,surf(SweepRangeImag,SweepRangeReal,20*log10(Y))
        for l = 2:size(Y,2)-1
            for j = 2:size(Y,1)-1
                if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                    XRough = [SweepRangeReal(j) SweepRangeImag(l)]; % real velocity (m/s), imaginary velocity (m/s)
                    break
                end
            end
        end
        if  ~isempty(XRough)
            break
        end
    end

    % fine search for the strongly damped A0/B0 when fluid-loading is present
    if  ~isempty(XRough)
        FALamb = Converger('A','Solid','TMM',Decoupled,XRough,Resolution,AngularFrequency,SweepRangeReal,SweepRangeImag,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a11,a21,a31,b12,b22,b23,b32,b33,b34);
    else
        FALamb = XFine; 
    end
    FAScholte = XFine;
else
    FALamb = XFine;
    FAScholte = zeros(1,4);
end

% search for the weakly damped S0/B1 and S1/B2 modes near the real axis
for p = 1:length(XS0)

    % rough search for the weakly damped S0/B1 and S1/B2 modes near the real axis
    XRough = [];
    SweepRangeReal = (1-.5*Range3/2)*XS0(p):Range3/2*XS0(p)/Steps3:(1+.5*Range3/2)*XS0(p);
    SweepRangeImag = -Range3*XS0(p):Range3*XS0(p)/Steps3:0;
    Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
    if  ~Decoupled
        Y = Computer_Coupled('S','Solid','TMM',Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
    else
        Y = Computer_Decoupled('S','Solid','TMM',Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33);
    end
    Y(:,end+1) = max(Y,[],'all');
% figure,surf([SweepRangeImag SweepRangeImag(end)-(SweepRangeImag(1)-SweepRangeImag(2))],SweepRangeReal,20*log10(Y))
    for l = 2:size(Y,2)-1
        for j = 2:size(Y,1)-1
            if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                XRough = [SweepRangeReal(j) SweepRangeImag(l)]; % real velocity (m/s), imaginary velocity (m/s)
                break
            end
        end
    end
    
    % fine search for the weakly damped S0/B1 and S1/B2 modes near the real axis
    if  ~isempty(XRough)
        FSLamb(p,:) = Converger('S','Solid','TMM',Decoupled,XRough,Resolution,AngularFrequency,SweepRangeReal,SweepRangeImag,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a11,a21,a31,b12,b22,b23,b32,b33,b34);
    else
        FSLamb(p,:) = [XS0(p) 0 XS0(p) 0];
    end
end
% if  Symmetric
%     ModeNames = {'A0Scholte ','A0        ','S0        ','S1        '};
% else
%     ModeNames = {'B0Scholte ','B0        ','B1        ','B2        '};
% end
% String = ['Modes @ ',num2str(Frequency),' kHz:',newline...
%           ModeNames{1},num2str(FAScholte(1)),', ',num2str(FAScholte(2)),', ',num2str(FAScholte(3)),', ',num2str(FAScholte(4)),newline...
%           ModeNames{2},num2str(FALamb(1)),', ',num2str(FALamb(2)),', ',num2str(FALamb(3)),', ',num2str(FALamb(4)),newline...
%           ModeNames{3},num2str(FSLamb(1,1)),', ',num2str(FSLamb(1,2)),', ',num2str(FSLamb(1,3)),', ',num2str(FSLamb(1,4))];
% if  ~Decoupled
%     String = append(String,newline,ModeNames{4},num2str(FSLamb(2,1)),', ',num2str(FSLamb(2,2)),', ',num2str(FSLamb(2,3)),', ',num2str(FSLamb(2,4)));
% end
% disp([String,newline,'-----------------------'])
end
function X = Converger(ModeFamily,ModeType,MatrixMethod,Decoupled,XRough,Resolution,AngularFrequency,SweepRangeReal,SweepRangeImag,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a11,a21,a31,b12,b22,b23,b32,b33,b34)
    Bisections = ceil(log2(Resolution/(abs(SweepRangeReal(1)-SweepRangeReal(2))))/log2(2/4));
    RangeReal = XRough(1)+(SweepRangeReal(2)-SweepRangeReal(1))*[-2 2];
    RangeImag = XRough(2)+(SweepRangeImag(2)-SweepRangeImag(1))*[-2 2];
    for k = 1:Bisections
        if  k == 1
            RangeReal = RangeReal(1):(RangeReal(end)-RangeReal(1))/8:RangeReal(end);
            RangeImag = RangeImag(1):(RangeImag(end)-RangeImag(1))/8:RangeImag(end);
        else
            RangeReal = RangeReal(1):(RangeReal(end)-RangeReal(1))/4:RangeReal(end);
            RangeImag = RangeImag(1):(RangeImag(end)-RangeImag(1))/4:RangeImag(end);
        end
        RangeImag(RangeImag > 0) = [];
        Wavenumber = AngularFrequency./(RangeReal'+RangeImag*1i);
        if  ~Decoupled
            Y = Computer_Coupled(ModeFamily,ModeType,MatrixMethod,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
        else
            Y = Computer_Decoupled(ModeFamily,ModeType,MatrixMethod,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33);
        end
        Y(:,end+1) = max(Y,[],'all');
% f = figure;surf([RangeImag RangeImag(end)-(RangeImag(1)-RangeImag(2))],RangeReal,20*log10(Y))
% close(f)
        for l = 2:size(Y,2)-1
            for j = 2:size(Y,1)-1
                if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                    MIN = [j l];
                    break
                end
            end
        end
        if  k < Bisections
            RangeReal = RangeReal(MIN(1))+(RangeReal(2)-RangeReal(1))*[-1 1];
            RangeImag = RangeImag(MIN(2))+(RangeImag(2)-RangeImag(1))*[-1 1];
        end
    end
    X = [AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) imag(Wavenumber(MIN(1),MIN(2))) RangeReal(MIN(1)) RangeImag(MIN(2))]; % phase velocity (m/s), attenuation (Np/m), real velocity (m/s), imaginary velocity (m/s)
end
function Y = Computer_Coupled(ModeFamily,ModeType,MatrixMethod,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34)
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Wavenumber2 = Wavenumber.^2;
    Wavenumber4 = Wavenumber2.^2;
    Wavenumber6 = Wavenumber2.^3;
    Length = length(Wavenumber);
    Y = NaN(size(Wavenumber)); % to be discarded once pagedet exists
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
            Phi = .5i*k3*LayerThicknesses;
        else
            Phi = 1i*k3*LayerThicknesses(m);
        end
        E = exp(Phi);
        if  strcmp(MatrixMethod,'TMM')
            E_ = exp(-Phi);
            L1 = [E E_;V.*E V.*E_;W.*E -W.*E_;D3.*E D3.*E_;D5.*E -D5.*E_;D4.*E -D4.*E_];
            L2 = [ones(1,6,Length);V V;W -W;D3 D3;D5 -D5;D4 -D4];
        elseif strcmp(MatrixMethod,'SMM')
            L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
            L2 = [ones(1,3,Length) E;V V.*E;W -W.*E;E ones(1,3,Length);V.*E V;W.*E -W];
        end
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    if  strcmp(MatrixMethod,'TMM')
        for m = 2:SuperLayerSize
            M{1} = pagemtimes(M{1},L{m});
        end
        for m = 1:length(Pattern)
            M{m+1} = pagemtimes(M{m},M{Pattern(m)});
        end
        if  Symmetric
            if  FluidLoading
                k3Fluid = sqrt(kUpperFluid2-Wavenumber2);
                if  strcmp(ModeType,'Fluid1') || strcmp(ModeType,'Fluid2')
                    k3Fluid = -k3Fluid;
                end
            end
            if  strcmp(ModeFamily,'S')
                if  FluidLoading
                    M{end} = [M{end}(3:6,[1:2 4],:) [k3Fluid./Wavenumber;gUpperFluid./Wavenumber;zeros(2,1,Length)]];
                    for j = 1:length(Wavenumber)
                        Y(j) = abs(det(M{end}(:,:,j)));
                    end
                else
                    for j = 1:length(Wavenumber)
                        Y(j) = abs(det(M{end}(4:6,[1:2 4],j)));
                    end
                end
            elseif strcmp(ModeFamily,'A')
                if  FluidLoading
                    M{end} = [M{end}(3:6,[3 5:6],:) [k3Fluid./Wavenumber;gUpperFluid./Wavenumber;zeros(2,1,Length)]];
                    for j = 1:length(Wavenumber)
                        Y(j) = abs(det(M{end}(:,:,j)));
                    end
                else
                    for j = 1:length(Wavenumber)
                        Y(j) = abs(det(M{end}(4:6,[3 5:6],j)));
                    end
                end
            end
        else
            if  SymmetricSystem
                M2{1} = L{end};
                for m = SuperLayerSize-1:-1:1
                    M2{1} = pagemtimes(M2{1},L{m});
                end
                for m = 1:length(Pattern)
                    M2{m+1} = pagemtimes(M2{m},M2{Pattern(m)});
                end
                M{end} = pagemtimes(M{end},M2{end});
            end
            if  FluidLoading
                if  ToggleUpperFluid && ToggleLowerFluid
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                    if  strcmp(ModeType,'Fluid1') && UpperFluid.Velocity > LowerFluid.Velocity
                        k3UpperFluid = -k3UpperFluid;
                    elseif strcmp(ModeType,'Fluid1') && UpperFluid.Velocity < LowerFluid.Velocity
                        k3LowerFluid = -k3LowerFluid;
                    elseif strcmp(ModeType,'Fluid2')
                        k3UpperFluid = -k3UpperFluid;
                        k3LowerFluid = -k3LowerFluid;
                    end
                    QUpperFluid = gUpperFluid/k3UpperFluid;
                    QLowerFluid = gLowerFluid/k3LowerFluid;
                    for j = 1:length(Wavenumber)
                        M11(1,1,j) = det(M{end}([3 5:6],1:3,j));
                        M12(1,1,j) = det(M{end}([3 5:6],[1:2 4],j));
                        M21(1,1,j) = det(M{end}(4:6,1:3,j));
                        M22(1,1,j) = det(M{end}(4:6,[1:2 4],j));
                    end
                    Y = abs(M21-QLowerFluid.*M22-QUpperFluid.*(M11-QLowerFluid.*M12));
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    if  strcmp(ModeType,'Fluid1') || strcmp(ModeType,'Fluid2')
                        k3UpperFluid = -k3UpperFluid;
                    end
                    for j = 1:length(Wavenumber)
                        M11(1,1,j) = det(M{end}([3 5:6],1:3,j));
                        M21(1,1,j) = det(M{end}(4:6,1:3,j));
                    end
                    Y = abs(M21-gUpperFluid/k3UpperFluid.*M11);
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                    if  strcmp(ModeType,'Fluid1') || strcmp(ModeType,'Fluid2')
                        k3LowerFluid = -k3LowerFluid;
                    end
                    for j = 1:length(Wavenumber)
                        M21(1,1,j) = det(M{end}(4:6,1:3,j));
                        M22(1,1,j) = det(M{end}(4:6,[1:2 4],j));
                    end
                    Y = abs(M21-gLowerFluid/k3LowerFluid.*M22);
                end
            else
                for j = 1:length(Wavenumber)
                    Y(j) = abs(det(M{end}(4:6,1:3,j)));
                end
            end
        end
    elseif strcmp(MatrixMethod,'SMM')
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
            if  FluidLoading
                M{end} = pageinv(M{end});
                k3Fluid = sqrt(kUpperFluid2-Wavenumber2);
                if  strcmp(ModeType,'Fluid1') || strcmp(ModeType,'Fluid2')
                    k3Fluid = -k3Fluid;
                end
            end
            if  strcmp(ModeFamily,'S')
                if  FluidLoading
                    Y = abs(k3Fluid/gUpperFluid+M{end}(3,1,:)-M{end}(3,4,:).*M{end}(6,1,:)./M{end}(6,4,:));
                else
                    for j = 1:length(Wavenumber)
                        Y(j) = abs(det(M{end}([1:3 5:6],1:5,j)));
                    end
                end
            elseif strcmp(ModeFamily,'A')
                if  FluidLoading
                    Y = abs(k3Fluid/gUpperFluid+M{end}(3,1,:)-M{end}(3,5,:).*M{end}(5,1,:)./M{end}(5,5,:)-(M{end}(3,5,:).*M{end}(5,6,:)./M{end}(5,5,:)-M{end}(3,6,:)).*(M{end}(4,1,:).*M{end}(5,5,:)-M{end}(4,5,:).*M{end}(5,1,:))./(M{end}(4,5,:).*M{end}(5,6,:)-M{end}(4,6,:).*M{end}(5,5,:)));
                else
                    for j = 1:length(Wavenumber)
                        Y(j) = abs(det(M{end}(1:4,[1:3 6],j)));
                    end
                end
            end
        else
            if  SymmetricSystem
                M0 = M{end}(4:6,4:6,:).*I-M{end}(4:6,4:6,:);
                M1 = pagemrdivide(M{end}(1:3,4:6,:),M0);
                M2 = pagemrdivide(M{end}(1:3,4:6,:).*I,M0);
                M{end} = [M{end}(1:3,1:3,:)+pagemtimes(M1,M{end}(4:6,1:3,:)) -pagemtimes(M1,M{end}(4:6,1:3,:).*I);pagemtimes(M2,M{end}(4:6,1:3,:)) M{end}(1:3,1:3,:).*I-pagemtimes(M2,M{end}(4:6,1:3,:).*I)];
            end
            if  FluidLoading
                M{end} = pageinv(M{end});
                if  ToggleUpperFluid && ToggleLowerFluid
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                    if  strcmp(ModeType,'Fluid1') && UpperFluid.Velocity > LowerFluid.Velocity
                        k3UpperFluid = -k3UpperFluid;
                    elseif strcmp(ModeType,'Fluid1') && UpperFluid.Velocity < LowerFluid.Velocity
                        k3LowerFluid = -k3LowerFluid;
                    elseif strcmp(ModeType,'Fluid2')
                        k3UpperFluid = -k3UpperFluid;
                        k3LowerFluid = -k3LowerFluid;
                    end
                    WUpperFluid = k3UpperFluid./Wavenumber;
                    WLowerFluid = k3LowerFluid./Wavenumber;
                    DUpperFluid = gUpperFluid./Wavenumber;
                    DLowerFluid = gLowerFluid./Wavenumber;
                    Y = abs((WLowerFluid-M{end}(6,4,:).*DLowerFluid).*(WUpperFluid+M{end}(3,1,:).*DUpperFluid)+M{end}(3,4,:).*M{end}(6,1,:).*DUpperFluid.*DLowerFluid); 
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    if  strcmp(ModeType,'Fluid1') || strcmp(ModeType,'Fluid2')
                        k3UpperFluid = -k3UpperFluid;
                    end
                    Y = abs(k3UpperFluid/gUpperFluid+M{end}(3,1,:));
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                    if  strcmp(ModeType,'Fluid1') || strcmp(ModeType,'Fluid2')
                        k3LowerFluid = -k3LowerFluid;
                    end
                    Y = abs(k3LowerFluid/gLowerFluid-M{end}(6,4,:));
                end
            else
                for j = 1:length(Wavenumber)
                    Y(j) = abs(det(M{end}(:,:,j)));
                end
            end
        end
    end
    Y = reshape(Y,Size);
end
function Y = Computer_Decoupled(ModeFamily,ModeType,MatrixMethod,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33)
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Wavenumber2 = Wavenumber.^2;
    Wavenumber4 = Wavenumber2.^2;
    Length = length(Wavenumber);
    Y = NaN(size(Wavenumber)); % to be discarded once pagedet exists
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
            Phi = .5i*k3*LayerThicknesses;
        else
            Phi = 1i*k3*LayerThicknesses(m);
        end
        E = exp(Phi);
        if  strcmp(MatrixMethod,'TMM')
            E_ = exp(-Phi);
            L1 = [E E_;W.*E -W.*E_;D3.*E D3.*E_;D5.*E -D5.*E_];
            L2 = [ones(1,4,Length);W -W;D3 D3;D5 -D5];
        elseif strcmp(MatrixMethod,'SMM')
            L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
            L2 = [ones(1,2,Length) E;W -W.*E;E ones(1,2,Length);W.*E -W];
        end
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    if  strcmp(MatrixMethod,'TMM')
        for m = 2:SuperLayerSize
            M{1} = pagemtimes(M{1},L{m});
        end
        for m = 1:length(Pattern)
            M{m+1} = pagemtimes(M{m},M{Pattern(m)});
        end
        if  Symmetric
            if  FluidLoading
                k3Fluid = sqrt(kUpperFluid2-Wavenumber2);
                if  strcmp(ModeType,'Fluid1') || strcmp(ModeType,'Fluid2')
                    k3Fluid = -k3Fluid;
                end
            end
            if  strcmp(ModeFamily,'S')
                if  FluidLoading
                    M{end} = [M{end}(2:4,[1 3],:) [k3Fluid./Wavenumber;gUpperFluid./Wavenumber;zeros(1,1,Length)]];
                    for j = 1:length(Wavenumber)
                        Y(j) = abs(det(M{end}(:,:,j)));
                    end
                else
                    for j = 1:length(Wavenumber)
                        Y(j) = abs(det(M{end}(3:4,[1 3],j)));
                    end
                end
            elseif strcmp(ModeFamily,'A')
                if  FluidLoading
                    M{end} = [M{end}(2:4,[2 4],:) [k3Fluid./Wavenumber;gUpperFluid./Wavenumber;zeros(1,1,Length)]];
                    for j = 1:length(Wavenumber)
                        Y(j) = abs(det(M{end}(:,:,j)));
                    end
                else
                    for j = 1:length(Wavenumber)
                        Y(j) = abs(det(M{end}(3:4,[2 4],j)));
                    end
                end
            end
        else
            if  SymmetricSystem
                M2{1} = L{end};
                for m = SuperLayerSize-1:-1:1
                    M2{1} = pagemtimes(M2{1},L{m});
                end
                for m = 1:length(Pattern)
                    M2{m+1} = pagemtimes(M2{m},M2{Pattern(m)});
                end 
                M{end} = pagemtimes(M{end},M2{end});
            end
            if  FluidLoading
                if  ToggleUpperFluid && ToggleLowerFluid
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                    if  strcmp(ModeType,'Fluid1') && UpperFluid.Velocity > LowerFluid.Velocity
                        k3UpperFluid = -k3UpperFluid;
                    elseif strcmp(ModeType,'Fluid1') && UpperFluid.Velocity < LowerFluid.Velocity
                        k3LowerFluid = -k3LowerFluid;
                    elseif strcmp(ModeType,'Fluid2')
                        k3UpperFluid = -k3UpperFluid;
                        k3LowerFluid = -k3LowerFluid;
                    end
                    QUpperFluid = gUpperFluid/k3UpperFluid;
                    QLowerFluid = gLowerFluid/k3LowerFluid;
                    for j = 1:length(Wavenumber)
                        M11(1,1,j) = det(M{end}([2 4],1:2,j));
                        M12(1,1,j) = det(M{end}([2 4],[1 3],j));
                        M21(1,1,j) = det(M{end}(3:4,1:2,j));
                        M22(1,1,j) = det(M{end}(3:4,[1 3],j));
                    end
                    Y = abs(M21-QLowerFluid.*M22-QUpperFluid.*(M11-QLowerFluid.*M12));
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    if  strcmp(ModeType,'Fluid1') || strcmp(ModeType,'Fluid2')
                        k3UpperFluid = -k3UpperFluid;
                    end
                    for j = 1:length(Wavenumber)
                        M11(1,1,j) = det(M{end}([2 4],1:2,j));
                        M21(1,1,j) = det(M{end}(3:4,1:2,j));
                    end
                    Y = abs(M21-gUpperFluid/k3UpperFluid.*M11);
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                    if  strcmp(ModeType,'Fluid1') || strcmp(ModeType,'Fluid2')
                        k3LowerFluid = -k3LowerFluid;
                    end
                    for j = 1:length(Wavenumber)
                        M21(1,1,j) = det(M{end}(3:4,1:2,j));
                        M22(1,1,j) = det(M{end}(3:4,[1 3],j));
                    end
                    Y = abs(M21-gLowerFluid/k3LowerFluid.*M22);
                end
            else
                for j = 1:length(Wavenumber)
                    Y(j) = abs(det(M{end}(3:4,1:2,j)));
                end
            end
        end
    elseif strcmp(MatrixMethod,'SMM')
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
            if  FluidLoading
                M{end} = pageinv(M{end});
                k3Fluid = sqrt(kUpperFluid2-Wavenumber2);
                if  strcmp(ModeType,'Fluid1') || strcmp(ModeType,'Fluid2')
                    k3Fluid = -k3Fluid;
                end
            end
            if  strcmp(ModeFamily,'S')
                if  FluidLoading
                    Y = abs(k3Fluid/gUpperFluid+M{end}(2,1,:)-M{end}(2,3,:).*M{end}(4,1,:)./M{end}(4,3,:));
                else
                    for j = 1:length(Wavenumber)
                        Y(j) = abs(det(M{end}([1:2 4],1:3,j)));
                    end
                end
            elseif strcmp(ModeFamily,'A')
                if  FluidLoading
                    Y = abs(k3Fluid/gUpperFluid+M{end}(2,1,:)-M{end}(2,4,:).*M{end}(3,1,:)./M{end}(3,4,:));
                else
                    for j = 1:length(Wavenumber)
                        Y(j) = abs(det(M{end}(1:3,[1:2 4],j)));
                    end
                end
            end
        else
            if  SymmetricSystem
                M0 = M{end}(3:4,3:4,:).*I1-M{end}(3:4,3:4,:);
                M1 = pagemrdivide(M{end}(1:2,3:4,:),M0);
                M2 = pagemrdivide(M{end}(1:2,3:4,:).*I1,M0);
                M{end} = [M{end}(1:2,1:2,:)+pagemtimes(M1,M{end}(3:4,1:2,:)) -pagemtimes(M1,M{end}(3:4,1:2,:).*I1);pagemtimes(M2,M{end}(3:4,1:2,:)) M{end}(1:2,1:2,:).*I1-pagemtimes(M2,M{end}(3:4,1:2,:).*I1)];
            end
            if  FluidLoading
                M{end} = pageinv(M{end});
                if  ToggleUpperFluid && ToggleLowerFluid
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                    if  strcmp(ModeType,'Fluid1') && UpperFluid.Velocity > LowerFluid.Velocity
                        k3UpperFluid = -k3UpperFluid;
                    elseif strcmp(ModeType,'Fluid1') && UpperFluid.Velocity < LowerFluid.Velocity
                        k3LowerFluid = -k3LowerFluid;
                    elseif strcmp(ModeType,'Fluid2')
                        k3UpperFluid = -k3UpperFluid;
                        k3LowerFluid = -k3LowerFluid;
                    end
                    WUpperFluid = k3UpperFluid./Wavenumber;
                    WLowerFluid = k3LowerFluid./Wavenumber;
                    DUpperFluid = gUpperFluid./Wavenumber;
                    DLowerFluid = gLowerFluid./Wavenumber;
                    Y = abs((WLowerFluid-M{end}(4,3,:).*DLowerFluid).*(WUpperFluid+M{end}(2,1,:).*DUpperFluid)+M{end}(2,3,:).*M{end}(4,1,:).*DUpperFluid.*DLowerFluid);
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    if  strcmp(ModeType,'Fluid1') || strcmp(ModeType,'Fluid2')
                        k3UpperFluid = -k3UpperFluid;
                    end
                    Y = abs(k3UpperFluid/gUpperFluid+M{end}(2,1,:));
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                    if  strcmp(ModeType,'Fluid1') || strcmp(ModeType,'Fluid2')
                        k3LowerFluid = -k3LowerFluid;
                    end
                    Y = abs(k3LowerFluid/gLowerFluid-M{end}(4,3,:));
                end
            else
                for j = 1:length(Wavenumber)
                    Y(j) = abs(det(M{end}(:,:,j)));
                end
            end
        end
    end
    Y = reshape(Y,Size);
end