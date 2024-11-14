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
function FLambF = PhaseVelocitySweeper_Isotropic_Pipe_F(FLamb,Material,Viscoelastic,Ro,Ri,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Frequency)
Resolution = 1e-6; % (m/s)
Range1 = .02; % for the weakly damped solution near the real axis
Steps1 = 2e2;
Range2 = .3; % for the strongly damped solution when fluid-loading is present
Steps2 = 2e2;

%#ok<*AGROW>
Ro2 = Ro^2;
Ri2 = Ri^2;
Densityi = InnerFluid.Density/Material.Density;
Densityo = OuterFluid.Density/Material.Density;
AngularFrequency = 2*pi*Frequency*1e3;
kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
kInnerFluid2 = (AngularFrequency/InnerFluid.Velocity)^2;
kOuterFluid2 = (AngularFrequency/OuterFluid.Velocity)^2;

% search for damped solution near the real axis
if  Viscoelastic

    % rough search for damped solution near the real axis
    XRough = [];
    SweepRangeReal = (1-.5*Range1)*FLamb:Range1*FLamb/Steps1:(1+.5*Range1)*FLamb;
    SweepRangeImag = -Range1*FLamb:Range1*FLamb/Steps1:0;
    Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
    Y = Computer(-1,Wavenumber,kL2,kT2,Ro,Ri,Ro2,Ri2,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink,kInnerFluid2,kOuterFluid2);
% figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
    for l = 2:size(Y,2)-1
        for j = 2:size(Y,1)-1
            if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                XRough = [SweepRangeReal(j) SweepRangeImag(l)]; % real velocity (m/s), imaginary velocity (m/s)
            end
        end
    end
    
    % fine search for damped solution near the real axis
    if  ~isempty(XRough)
        Bisections = ceil(log2(Resolution/(2*(SweepRangeReal(2)-SweepRangeReal(1))))/log2(2*.25));
        RangeReal = [XRough(1)-2*(SweepRangeReal(2)-SweepRangeReal(1)) XRough(1)+2*(SweepRangeReal(2)-SweepRangeReal(1))];
        RangeImag = [XRough(2)-2*(SweepRangeImag(2)-SweepRangeImag(1)) XRough(2)+2*(SweepRangeImag(2)-SweepRangeImag(1))];
        for k = 1:Bisections % search minimum in characteristic equation and converge upon it
            if  k == 1
                RangeReal = RangeReal(1):.125*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                RangeImag = RangeImag(1):.125*(RangeImag(end)-RangeImag(1)):RangeImag(end);
            else
                RangeReal = RangeReal(1):.25*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                RangeImag = RangeImag(1):.25*(RangeImag(end)-RangeImag(1)):RangeImag(end);
            end
            Wavenumber = AngularFrequency./(RangeReal'+RangeImag*1i);
            Y = Computer(-1,Wavenumber,kL2,kT2,Ro,Ri,Ro2,Ri2,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink,kInnerFluid2,kOuterFluid2);
            if  abs(RangeImag(end)) < 1e-3
                Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
            end
% if  abs(RangeImag(end)) < 1e-3
% f = figure;surf(horzcat(RangeImag,-RangeImag(end-1)),RangeReal,20*log10(Y))    
% else
% f = figure;surf(RangeImag,RangeReal,20*log10(Y))
% end
% close(f)
            for l = 2:size(Y,2)-1
                for j = 2:size(Y,1)-1
                    if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                        MIN = [j l];
                    end
                end
            end
            if  k < Bisections % set the new search area around the found minimum
                RangeReal = [RangeReal(MIN(1))-(RangeReal(2)-RangeReal(1)) RangeReal(MIN(1))+(RangeReal(2)-RangeReal(1))];
                RangeImag = [RangeImag(MIN(2))-(RangeImag(2)-RangeImag(1)) RangeImag(MIN(2))+(RangeImag(2)-RangeImag(1))];
            end
        end
        FLambF(1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
        FLambF(2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
        FLambF(3) = RangeReal(MIN(1)); % real velocity (m/s)
        FLambF(4) = RangeImag(MIN(2)); % imaginary velocity (m/s)
    else
        FLambF = [FLamb 0 FLamb 0];
    end
else
    FLambF = [FLamb 0 FLamb 0];
end
if  FLambF(1) < Material.TransverseVelocity+2 && FLambF(1) > Material.TransverseVelocity-2
    FLambF = [0 0 0 0];
end

% search for strongly damped solutions when fluid-loading is present
for p = 1%:2
    if  (p == 1 && ~ToggleOuterFluid && FLambF(1) > 0) || (p == 2 && (~ToggleInnerFluid || (ToggleInnerFluid && ~Sink)))
        % if  p == 2
            FLambF(3,:) = 0;
        % end
        continue
    end

    % rough search for strongly damped solutions when fluid-loading is present
    XRough = [];
    for o = [1 2 5]
        SweepRangeReal = 0:Range2*FLamb/Steps2/o:Range2*FLamb;
        if  p == 1
            SweepRangeImag = -SweepRangeReal;    
        elseif p == 2
            SweepRangeImag = SweepRangeReal;
        end
        Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
        Y = Computer(1,Wavenumber,kL2,kT2,Ro,Ri,Ro2,Ri2,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink,kInnerFluid2,kOuterFluid2); 
% figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
        for l = 2:size(Y,2)-1
            for j = 2:size(Y,1)-1
                if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                    XRough(end+1,:) = [SweepRangeReal(j) SweepRangeImag(l) imag(Wavenumber(j,l))]; % real velocity (m/s), imaginary velocity (m/s)
                end
            end
        end
        if  ~isempty(XRough)
            break
        end
    end
    if  height(XRough) > 1
        [~,z] = min(abs(XRough(:,3)));
        XRough = XRough(z,:);
    end

    % fine search for strongly damped solutions when fluid-loading is present
    if  ~isempty(XRough)
        Bisections = ceil(log2(Resolution/(2*(SweepRangeReal(2)-SweepRangeReal(1))))/log2(2*.25));
        RangeReal = [XRough(1)-2*(SweepRangeReal(2)-SweepRangeReal(1)) XRough(1)+2*(SweepRangeReal(2)-SweepRangeReal(1))];
        RangeImag = [XRough(2)-2*(SweepRangeImag(2)-SweepRangeImag(1)) XRough(2)+2*(SweepRangeImag(2)-SweepRangeImag(1))];
        for k = 1:Bisections % search minimum in characteristic equation and converge upon it
            if  k == 1
                RangeReal = RangeReal(1):.125*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                RangeImag = RangeImag(1):.125*(RangeImag(end)-RangeImag(1)):RangeImag(end);
            else
                RangeReal = RangeReal(1):.25*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                RangeImag = RangeImag(1):.25*(RangeImag(end)-RangeImag(1)):RangeImag(end);
            end
            Wavenumber = AngularFrequency./(RangeReal'+RangeImag*1i);
            Y = Computer(1,Wavenumber,kL2,kT2,Ro,Ri,Ro2,Ri2,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink,kInnerFluid2,kOuterFluid2);
            if  abs(RangeImag(end)) < 1e-3
                Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
            end
% if  abs(RangeImag(end)) < 1e-3
% f = figure;surf(horzcat(RangeImag,-RangeImag(end-1)),RangeReal,20*log10(Y))    
% else
% f = figure;surf(RangeImag,RangeReal,20*log10(Y))
% end
% close(f)
            for l = 2:size(Y,2)-1
                for j = 2:size(Y,1)-1
                    if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                        MIN = [j l];
                    end
                end
            end
            if  k < Bisections % set the new search area around the found minimum
                RangeReal = [RangeReal(MIN(1))-(RangeReal(2)-RangeReal(1)) RangeReal(MIN(1))+(RangeReal(2)-RangeReal(1))];
                RangeImag = [RangeImag(MIN(2))-(RangeImag(2)-RangeImag(1)) RangeImag(MIN(2))+(RangeImag(2)-RangeImag(1))];
            end
        end
        FLambF(p+1,1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
        FLambF(p+1,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
        FLambF(p+1,3) = RangeReal(MIN(1)); % real phase velocity (m/s)
        FLambF(p+1,4) = RangeImag(MIN(2)); % imaginary phase velocity (m/s)
    end
end
% if  FLambF(1) == 0
%     String = ['Undamped --',newline];
% else
%     String = ['Undamped ',num2str(FLamb),newline];
% end
% for p = 1:height(FLambF)
%     String = append(String,'F(1,',num2str(p),')   ',num2str(FLambF(p,1)),', ',num2str(FLambF(p,2)),', ',num2str(FLambF(p,3)),', ',num2str(FLambF(p,4)),newline);
% end
% disp([String,'-----------------------------------------'])
end
function Y = Computer(Sign,Wavenumber,kL2,kT2,Ro,Ri,Ro2,Ri2,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink,kInnerFluid2,kOuterFluid2)
    k2 = Wavenumber.^2;
    x = sqrt(k2-kL2);
    y = sqrt(k2-kT2);
    y2 = kT2-k2;
    xRo = x*Ro;
    yRo = y*Ro;
    xRi = x*Ri;
    yRi = y*Ri;
    Z1xi = besseli(1,xRi);
    Z1xo = besseli(1,xRo);
    Z1yi = besseli(1,yRi);
    Z1yo = besseli(1,yRo);
    W1xi = besselk(1,xRi);
    W1xo = besselk(1,xRo);
    W1yi = besselk(1,yRi);
    W1yo = besselk(1,yRo);
    dZ1xiRi = Z1xi+xRi.*besseli(2,xRi);
    dZ1xoRo = Z1xo+xRo.*besseli(2,xRo);
    dZ1yiRi = Z1yi+yRi.*besseli(2,yRi);
    dZ1yoRo = Z1yo+yRo.*besseli(2,yRo);
    dW1xiRi = W1xi-xRi.*besselk(2,xRi);
    dW1xoRo = W1xo-xRo.*besselk(2,xRo);
    dW1yiRi = W1yi-yRi.*besselk(2,yRi);
    dW1yoRo = W1yo-yRo.*besselk(2,yRo);
    if  ToggleInnerFluid
        zRi = Sign*sqrt(kInnerFluid2-k2)*Ri;
        if  Sink
            Z1zi = besselh(1,2,zRi);
            dZ1ziRi = Z1zi-zRi.*besselh(2,2,zRi);
        else
            Z1zi = besselj(1,zRi);
            dZ1ziRi = Z1zi-zRi.*besselj(2,zRi);
        end
    end
    if  ToggleOuterFluid
        zRo = Sign*sqrt(kOuterFluid2-k2)*Ro;
        H1zo = besselh(1,zRo);
        dH1zoRo = H1zo-zRo.*besselh(2,zRo);
    end
    Y = NaN(size(Wavenumber));
    for l = 1:width(Wavenumber)
        for j = 1:height(Wavenumber)
            if  ToggleInnerFluid && ToggleOuterFluid
                M(1,1) = dZ1xiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-1)*Z1xi(j,l);
                M(1,2) = dW1xiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-1)*W1xi(j,l);
                M(1,3) = -Wavenumber(j,l)*(dZ1yiRi(j,l)+(y2(j,l)*Ri2-1)*Z1yi(j,l));
                M(1,4) = -Wavenumber(j,l)*(dW1yiRi(j,l)+(y2(j,l)*Ri2-1)*W1yi(j,l));
                M(1,5) = dZ1yiRi(j,l)-Z1yi(j,l);
                M(1,6) = dW1yiRi(j,l)-W1yi(j,l);
                M(1,7) = .5*kT2*Ri2*Densityi*Z1zi(j,l);
                M(2,1) = 2*(dZ1xiRi(j,l)-Z1xi(j,l));
                M(2,2) = 2*(dW1xiRi(j,l)-W1xi(j,l));
                M(2,3) = 2*Wavenumber(j,l)*(Z1yi(j,l)-dZ1yiRi(j,l));
                M(2,4) = 2*Wavenumber(j,l)*(W1yi(j,l)-dW1yiRi(j,l));
                M(2,5) = 2*dZ1yiRi(j,l)+(y2(j,l)*Ri2-2)*Z1yi(j,l);
                M(2,6) = 2*dW1yiRi(j,l)+(y2(j,l)*Ri2-2)*W1yi(j,l);
                M(3,1) = 2*Wavenumber(j,l)*dZ1xiRi(j,l);
                M(3,2) = 2*Wavenumber(j,l)*dW1xiRi(j,l);
                M(3,3) = (y2(j,l)-k2(j,l))*dZ1yiRi(j,l);
                M(3,4) = (y2(j,l)-k2(j,l))*dW1yiRi(j,l);
                M(3,5) = -Wavenumber(j,l)*Z1yi(j,l);
                M(3,6) = -Wavenumber(j,l)*W1yi(j,l);
                M(4,1) = dZ1xiRi(j,l);
                M(4,2) = dW1xiRi(j,l);
                M(4,3) = -Wavenumber(j,l)*dZ1yiRi(j,l);
                M(4,4) = -Wavenumber(j,l)*dW1yiRi(j,l);
                M(4,5) = -Z1yi(j,l);
                M(4,6) = -W1yi(j,l);
                M(4,7) = dZ1ziRi(j,l);
                M(5,1) = dZ1xoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-1)*Z1xo(j,l);
                M(5,2) = dW1xoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-1)*W1xo(j,l);
                M(5,3) = -Wavenumber(j,l)*(dZ1yoRo(j,l)+(y2(j,l)*Ro2-1)*Z1yo(j,l));
                M(5,4) = -Wavenumber(j,l)*(dW1yoRo(j,l)+(y2(j,l)*Ro2-1)*W1yo(j,l));
                M(5,5) = dZ1yoRo(j,l)-Z1yo(j,l);
                M(5,6) = dW1yoRo(j,l)-W1yo(j,l);
                M(5,8) = .5*kT2*Ro2*Densityo*H1zo(j,l);
                M(6,1) = 2*(dZ1xoRo(j,l)-Z1xo(j,l));
                M(6,2) = 2*(dW1xoRo(j,l)-W1xo(j,l));
                M(6,3) = 2*Wavenumber(j,l)*(Z1yo(j,l)-dZ1yoRo(j,l));
                M(6,4) = 2*Wavenumber(j,l)*(W1yo(j,l)-dW1yoRo(j,l));
                M(6,5) = 2*dZ1yoRo(j,l)+(y2(j,l)*Ro2-2)*Z1yo(j,l);
                M(6,6) = 2*dW1yoRo(j,l)+(y2(j,l)*Ro2-2)*W1yo(j,l);
                M(7,1) = 2*Wavenumber(j,l)*dZ1xoRo(j,l);
                M(7,2) = 2*Wavenumber(j,l)*dW1xoRo(j,l);
                M(7,3) = (y2(j,l)-k2(j,l))*dZ1yoRo(j,l);
                M(7,4) = (y2(j,l)-k2(j,l))*dW1yoRo(j,l);
                M(7,5) = -Wavenumber(j,l)*Z1yo(j,l);
                M(7,6) = -Wavenumber(j,l)*W1yo(j,l);
                M(8,1) = dZ1xoRo(j,l);
                M(8,2) = dW1xoRo(j,l);
                M(8,3) = -Wavenumber(j,l)*dZ1yoRo(j,l);
                M(8,4) = -Wavenumber(j,l)*dW1yoRo(j,l);
                M(8,5) = -Z1yo(j,l);
                M(8,6) = -W1yo(j,l);
                M(8,8) = dH1zoRo(j,l);
            elseif ToggleInnerFluid && ~ToggleOuterFluid
                M(1,1) = dZ1xiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-1)*Z1xi(j,l);
                M(1,2) = dW1xiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-1)*W1xi(j,l);
                M(1,3) = -Wavenumber(j,l)*(dZ1yiRi(j,l)+(y2(j,l)*Ri2-1)*Z1yi(j,l));
                M(1,4) = -Wavenumber(j,l)*(dW1yiRi(j,l)+(y2(j,l)*Ri2-1)*W1yi(j,l));
                M(1,5) = dZ1yiRi(j,l)-Z1yi(j,l);
                M(1,6) = dW1yiRi(j,l)-W1yi(j,l);
                M(1,7) = .5*kT2*Ri2*Densityi*Z1zi(j,l);
                M(2,1) = 2*(dZ1xiRi(j,l)-Z1xi(j,l));
                M(2,2) = 2*(dW1xiRi(j,l)-W1xi(j,l));
                M(2,3) = 2*Wavenumber(j,l)*(Z1yi(j,l)-dZ1yiRi(j,l));
                M(2,4) = 2*Wavenumber(j,l)*(W1yi(j,l)-dW1yiRi(j,l));
                M(2,5) = 2*dZ1yiRi(j,l)+(y2(j,l)*Ri2-2)*Z1yi(j,l);
                M(2,6) = 2*dW1yiRi(j,l)+(y2(j,l)*Ri2-2)*W1yi(j,l);
                M(3,1) = 2*Wavenumber(j,l)*dZ1xiRi(j,l);
                M(3,2) = 2*Wavenumber(j,l)*dW1xiRi(j,l);
                M(3,3) = (y2(j,l)-k2(j,l))*dZ1yiRi(j,l);
                M(3,4) = (y2(j,l)-k2(j,l))*dW1yiRi(j,l);
                M(3,5) = -Wavenumber(j,l)*Z1yi(j,l);
                M(3,6) = -Wavenumber(j,l)*W1yi(j,l);
                M(4,1) = dZ1xiRi(j,l);
                M(4,2) = dW1xiRi(j,l);
                M(4,3) = -Wavenumber(j,l)*dZ1yiRi(j,l);
                M(4,4) = -Wavenumber(j,l)*dW1yiRi(j,l);
                M(4,5) = -Z1yi(j,l);
                M(4,6) = -W1yi(j,l);
                M(4,7) = dZ1ziRi(j,l);
                M(5,1) = dZ1xoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-1)*Z1xo(j,l);
                M(5,2) = dW1xoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-1)*W1xo(j,l);
                M(5,3) = -Wavenumber(j,l)*(dZ1yoRo(j,l)+(y2(j,l)*Ro2-1)*Z1yo(j,l));
                M(5,4) = -Wavenumber(j,l)*(dW1yoRo(j,l)+(y2(j,l)*Ro2-1)*W1yo(j,l));
                M(5,5) = dZ1yoRo(j,l)-Z1yo(j,l);
                M(5,6) = dW1yoRo(j,l)-W1yo(j,l);
                M(6,1) = 2*(dZ1xoRo(j,l)-Z1xo(j,l));
                M(6,2) = 2*(dW1xoRo(j,l)-W1xo(j,l));
                M(6,3) = 2*Wavenumber(j,l)*(Z1yo(j,l)-dZ1yoRo(j,l));
                M(6,4) = 2*Wavenumber(j,l)*(W1yo(j,l)-dW1yoRo(j,l));
                M(6,5) = 2*dZ1yoRo(j,l)+(y2(j,l)*Ro2-2)*Z1yo(j,l);
                M(6,6) = 2*dW1yoRo(j,l)+(y2(j,l)*Ro2-2)*W1yo(j,l);
                M(7,1) = 2*Wavenumber(j,l)*dZ1xoRo(j,l);
                M(7,2) = 2*Wavenumber(j,l)*dW1xoRo(j,l);
                M(7,3) = (y2(j,l)-k2(j,l))*dZ1yoRo(j,l);
                M(7,4) = (y2(j,l)-k2(j,l))*dW1yoRo(j,l);
                M(7,5) = -Wavenumber(j,l)*Z1yo(j,l);
                M(7,6) = -Wavenumber(j,l)*W1yo(j,l);
            elseif ~ToggleInnerFluid && ToggleOuterFluid
                M(1,1) = dZ1xiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-1)*Z1xi(j,l);
                M(1,2) = dW1xiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-1)*W1xi(j,l);
                M(1,3) = -Wavenumber(j,l)*(dZ1yiRi(j,l)+(y2(j,l)*Ri2-1)*Z1yi(j,l));
                M(1,4) = -Wavenumber(j,l)*(dW1yiRi(j,l)+(y2(j,l)*Ri2-1)*W1yi(j,l));
                M(1,5) = dZ1yiRi(j,l)-Z1yi(j,l);
                M(1,6) = dW1yiRi(j,l)-W1yi(j,l);
                M(2,1) = 2*(dZ1xiRi(j,l)-Z1xi(j,l));
                M(2,2) = 2*(dW1xiRi(j,l)-W1xi(j,l));
                M(2,3) = 2*Wavenumber(j,l)*(Z1yi(j,l)-dZ1yiRi(j,l));
                M(2,4) = 2*Wavenumber(j,l)*(W1yi(j,l)-dW1yiRi(j,l));
                M(2,5) = 2*dZ1yiRi(j,l)+(y2(j,l)*Ri2-2)*Z1yi(j,l);
                M(2,6) = 2*dW1yiRi(j,l)+(y2(j,l)*Ri2-2)*W1yi(j,l);
                M(3,1) = 2*Wavenumber(j,l)*dZ1xiRi(j,l);
                M(3,2) = 2*Wavenumber(j,l)*dW1xiRi(j,l);
                M(3,3) = (y2(j,l)-k2(j,l))*dZ1yiRi(j,l);
                M(3,4) = (y2(j,l)-k2(j,l))*dW1yiRi(j,l);
                M(3,5) = -Wavenumber(j,l)*Z1yi(j,l);
                M(3,6) = -Wavenumber(j,l)*W1yi(j,l);
                M(4,1) = dZ1xoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-1)*Z1xo(j,l);
                M(4,2) = dW1xoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-1)*W1xo(j,l);
                M(4,3) = -Wavenumber(j,l)*(dZ1yoRo(j,l)+(y2(j,l)*Ro2-1)*Z1yo(j,l));
                M(4,4) = -Wavenumber(j,l)*(dW1yoRo(j,l)+(y2(j,l)*Ro2-1)*W1yo(j,l));
                M(4,5) = dZ1yoRo(j,l)-Z1yo(j,l);
                M(4,6) = dW1yoRo(j,l)-W1yo(j,l);
                M(4,7) = .5*kT2*Ro2*Densityo*H1zo(j,l);
                M(5,1) = 2*(dZ1xoRo(j,l)-Z1xo(j,l));
                M(5,2) = 2*(dW1xoRo(j,l)-W1xo(j,l));
                M(5,3) = 2*Wavenumber(j,l)*(Z1yo(j,l)-dZ1yoRo(j,l));
                M(5,4) = 2*Wavenumber(j,l)*(W1yo(j,l)-dW1yoRo(j,l));
                M(5,5) = 2*dZ1yoRo(j,l)+(y2(j,l)*Ro2-2)*Z1yo(j,l);
                M(5,6) = 2*dW1yoRo(j,l)+(y2(j,l)*Ro2-2)*W1yo(j,l);
                M(6,1) = 2*Wavenumber(j,l)*dZ1xoRo(j,l);
                M(6,2) = 2*Wavenumber(j,l)*dW1xoRo(j,l);
                M(6,3) = (y2(j,l)-k2(j,l))*dZ1yoRo(j,l);
                M(6,4) = (y2(j,l)-k2(j,l))*dW1yoRo(j,l);
                M(6,5) = -Wavenumber(j,l)*Z1yo(j,l);
                M(6,6) = -Wavenumber(j,l)*W1yo(j,l);
                M(7,1) = dZ1xoRo(j,l);
                M(7,2) = dW1xoRo(j,l);
                M(7,3) = -Wavenumber(j,l)*dZ1yoRo(j,l);
                M(7,4) = -Wavenumber(j,l)*dW1yoRo(j,l);
                M(7,5) = -Z1yo(j,l);
                M(7,6) = -W1yo(j,l);
                M(7,7) = dH1zoRo(j,l);
            elseif ~ToggleInnerFluid && ~ToggleOuterFluid
                M(1,1) = dZ1xiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-1)*Z1xi(j,l);
                M(1,2) = dW1xiRi(j,l)+(.5*kT2*Ri2-k2(j,l)*Ri2-1)*W1xi(j,l);
                M(1,3) = -Wavenumber(j,l)*(dZ1yiRi(j,l)+(y2(j,l)*Ri2-1)*Z1yi(j,l));
                M(1,4) = -Wavenumber(j,l)*(dW1yiRi(j,l)+(y2(j,l)*Ri2-1)*W1yi(j,l));
                M(1,5) = dZ1yiRi(j,l)-Z1yi(j,l);
                M(1,6) = dW1yiRi(j,l)-W1yi(j,l);
                M(2,1) = 2*(dZ1xiRi(j,l)-Z1xi(j,l));
                M(2,2) = 2*(dW1xiRi(j,l)-W1xi(j,l));
                M(2,3) = 2*Wavenumber(j,l)*(Z1yi(j,l)-dZ1yiRi(j,l));
                M(2,4) = 2*Wavenumber(j,l)*(W1yi(j,l)-dW1yiRi(j,l));
                M(2,5) = 2*dZ1yiRi(j,l)+(y2(j,l)*Ri2-2)*Z1yi(j,l);
                M(2,6) = 2*dW1yiRi(j,l)+(y2(j,l)*Ri2-2)*W1yi(j,l);
                M(3,1) = 2*Wavenumber(j,l)*dZ1xiRi(j,l);
                M(3,2) = 2*Wavenumber(j,l)*dW1xiRi(j,l);
                M(3,3) = (y2(j,l)-k2(j,l))*dZ1yiRi(j,l);
                M(3,4) = (y2(j,l)-k2(j,l))*dW1yiRi(j,l);
                M(3,5) = -Wavenumber(j,l)*Z1yi(j,l);
                M(3,6) = -Wavenumber(j,l)*W1yi(j,l);
                M(4,1) = dZ1xoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-1)*Z1xo(j,l);
                M(4,2) = dW1xoRo(j,l)+(.5*kT2*Ro2-k2(j,l)*Ro2-1)*W1xo(j,l);
                M(4,3) = -Wavenumber(j,l)*(dZ1yoRo(j,l)+(y2(j,l)*Ro2-1)*Z1yo(j,l));
                M(4,4) = -Wavenumber(j,l)*(dW1yoRo(j,l)+(y2(j,l)*Ro2-1)*W1yo(j,l));
                M(4,5) = dZ1yoRo(j,l)-Z1yo(j,l);
                M(4,6) = dW1yoRo(j,l)-W1yo(j,l);
                M(5,1) = 2*(dZ1xoRo(j,l)-Z1xo(j,l));
                M(5,2) = 2*(dW1xoRo(j,l)-W1xo(j,l));
                M(5,3) = 2*Wavenumber(j,l)*(Z1yo(j,l)-dZ1yoRo(j,l));
                M(5,4) = 2*Wavenumber(j,l)*(W1yo(j,l)-dW1yoRo(j,l));
                M(5,5) = 2*dZ1yoRo(j,l)+(y2(j,l)*Ro2-2)*Z1yo(j,l);
                M(5,6) = 2*dW1yoRo(j,l)+(y2(j,l)*Ro2-2)*W1yo(j,l);
                M(6,1) = 2*Wavenumber(j,l)*dZ1xoRo(j,l);
                M(6,2) = 2*Wavenumber(j,l)*dW1xoRo(j,l);
                M(6,3) = (y2(j,l)-k2(j,l))*dZ1yoRo(j,l);
                M(6,4) = (y2(j,l)-k2(j,l))*dW1yoRo(j,l);
                M(6,5) = -Wavenumber(j,l)*Z1yo(j,l);
                M(6,6) = -Wavenumber(j,l)*W1yo(j,l);
            end
            Y(j,l) = abs(det(M));
        end
    end
end