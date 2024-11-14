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
function LScholte = Computer_Isotropic_Pipe_LScholte_D2(Multithreading,Q,ax,Symmetric,Viscoelastic,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRange,Ro,Ri,H,FrequencyResolution,Resolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples)        
% Note: The lower limit of SweepRangeImag is taken as -1e-10 instead
% of 0 because when there is only fluid inside the pipe and no outer fluid 
% and there is a sink inside, the dispersion equation amplitude has a sharp 
% drop at SweepRangeImag = 0. This seems to happen only below the fluid
% velocity.

PhaseVelocityOffset = 1; % (m/s)
% PhaseVelocityOffset = 0;
if  ToggleInnerFluid && ~Sink
    FrequencyRange_Fundamental = FrequencyRange;
else
    x = 100/(Ro-Ri)/1e3; % limit up to which a 5 times smaller frequency step is used (kHz/mm)
    FrequencyRange_Fundamental = [FrequencyRange(1) 2*FrequencyRange(1):FrequencyResolution/5:x x+FrequencyResolution:FrequencyResolution:FrequencyRange(end)];
    FrequencyRange_Fundamental(FrequencyRange_Fundamental > FrequencyRange(end)) = [];
    if  FrequencyRange_Fundamental(end) ~= FrequencyRange(end)
        FrequencyRange_Fundamental(end) = FrequencyRange(end);
    end
end
FrequencyRange_Fundamental = fliplr(FrequencyRange_Fundamental);

%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS>
global Stop
Stop = 0;
LScholte{1} = [];
if  ~Multithreading
    for p = 1:height(H)
        g(p) = animatedline(ax,'LineStyle','-.','color','r');
    end
end
if  ~Symmetric && ToggleOuterFluid && ToggleInnerFluid
    if  OuterFluid.Velocity > InnerFluid.Velocity
        FastFluidVelocity = OuterFluid.Velocity;
        SlowFluidVelocity = InnerFluid.Velocity;
    else
        FastFluidVelocity = InnerFluid.Velocity;
        SlowFluidVelocity = OuterFluid.Velocity;
    end
elseif Symmetric || (ToggleOuterFluid && ~ToggleInnerFluid)
    FastFluidVelocity = OuterFluid.Velocity;
    SlowFluidVelocity = OuterFluid.Velocity;
elseif ~ToggleOuterFluid && ToggleInnerFluid
    FastFluidVelocity = InnerFluid.Velocity;
    SlowFluidVelocity = InnerFluid.Velocity;
end
Ro2 = Ro^2;
Ri2 = Ri^2;
Densityo = OuterFluid.Density/Material.Density;
Densityi = InnerFluid.Density/Material.Density;
MissingModes(height(H)) = 0;
Diff = diff(H(:,3));
Hdiff = min([[Diff;Inf] [Inf;Diff]],[],2);
% Hdiff(abs(H(:,1)-Material.TransverseVelocity) < 1e-6) = [];
% H(abs(H(:,1)-Material.TransverseVelocity) < 1e-6,:) = [];
for p = 1:height(H)
    X = H(p,:);
    if  Multithreading
        send(Q,[FrequencyRange_Fundamental(1),X(1)/1e3,p])
    else
        addpoints(g(p),FrequencyRange_Fundamental(1),X(1)/1e3);
        drawnow limitrate
    end
    Misses = 0;
    for i = 2:length(FrequencyRange_Fundamental)
        if  Stop
            return
        end
        AngularFrequency = 2*pi*FrequencyRange_Fundamental(i)*1e3;
        kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
        kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
        if  ToggleOuterFluid
            kOuterFluid2 = (AngularFrequency/OuterFluid.Velocity)^2;
        end
        if  ToggleInnerFluid
            kInnerFluid2 = (AngularFrequency/InnerFluid.Velocity)^2;
        end
        X(i,1) = 0;
        Neighbors = [];
        if  p > 1
            for j = 1:length(LScholte)
                if  i <= height(LScholte{j})
                    Neighbors(j,:) = [LScholte{j}(i,7) LScholte{j}(i,8)];
                end
            end
        end
        NeighborsNumber = height(Neighbors);
        for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal x SweepRangeImag)
            if  numel(find(X(:,1) ~= 0)) < 4
                if  Hdiff(p) > .001*X(end-1,3)
                    SweepRangeReal = [1.001*X(end-1,3) .999*X(end-1,3)];
                else
                    SweepRangeReal = [X(end-1,3)+Hdiff(p) X(end-1,3)-Hdiff(p)];
                end
            else
                SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
            end
            if  SweepRangeReal(1) == SweepRangeReal(2)
                SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
            end
            if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
            end
            if  SweepRangeReal(2) < 0
                SweepRangeReal(2) = 0;
            end
            if  numel(find(X(:,1) ~= 0)) < 4
                SweepRangeImag= [1.1*X(end-1,4) .9*X(end-1,4)];
            else
                SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary velocity corresponds to the attenuation
                if  -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                    SweepRangeImag(2) = -1e-10;
                end
            end
            if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
            end
            if  SweepRangeImag(2) > -1e-10
                SweepRangeImag(2) = -1e-10;
            end
            for o = 1:SearchAreaSections % increase search resolution
                if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                    break
                end 
                for k = 1:1e2 % search minimum in characteristic equation and converge upon it 
                    if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                        if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                            SweepRangeReal = SweepRangeReal(1):.25^o/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                            SweepRangeImag = SweepRangeImag(1):.25/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                        elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            SweepRangeReal = SweepRangeReal(1):.25/q/o*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                            SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                        elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            SweepRangeReal = SweepRangeReal(1):.25/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                            SweepRangeImag = SweepRangeImag(1):.25^o/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                        end
                    else
                        if  length(SweepRangeReal) == 2
                            SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                        end
                        if  length(SweepRangeImag) == 2
                            SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                        end
                    end
                    if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                        break
                    end
                    Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
                    cp = AngularFrequency./real(Wavenumber);
                    if  ~isempty(Neighbors) % remove solutions of previously found lower modes
                        if  Viscoelastic
                            for l = 2:width(Wavenumber)-1
                                for j = 2:height(Wavenumber)-1
                                    for m = 1:NeighborsNumber
                                        if  SweepRangeReal(j-1) > Neighbors(m,1) && SweepRangeReal(j+1) < Neighbors(m,1) && SweepRangeImag(l-1) < Neighbors(m,2) && SweepRangeImag(l+1) > Neighbors(m,2)
                                            Wavenumber(j,l) = NaN;
                                        end
                                    end
                                end
                            end
                        else
                            for j = 2:height(Wavenumber)-1
                                for m = 1:NeighborsNumber
                                    if  SweepRangeReal(j-1) > Neighbors(m,1) && SweepRangeReal(j+1) < Neighbors(m,1)
                                        Wavenumber(j,:) = NaN;
                                    end
                                end
                            end
                        end
                    end
                    for l = 1:width(Wavenumber)
                        for j = 2:height(Wavenumber)-1
                            if  cp(j-1,l) > Material.TransverseVelocity && cp(j+1,l) < Material.TransverseVelocity
                                Wavenumber(j,l) = NaN;
                            end
                        end
                    end
                    k2 = Wavenumber.^2;
                    if  (numel(find(X(:,1) ~= 0)) <= 3 && H(p,1) > Material.LongitudinalVelocity) || (numel(find(X(:,1) ~= 0)) > 3 && X(end-1,1) > Material.LongitudinalVelocity)
                        x = sqrt(kL2-k2);
                        y = sqrt(kT2-k2);
                    elseif (numel(find(X(:,1) ~= 0)) <= 3 && H(p,1) < Material.LongitudinalVelocity && H(p,1) > Material.TransverseVelocity) || (numel(find(X(:,1) ~= 0)) > 3 && X(end-1,1) < Material.LongitudinalVelocity && X(end-1,1) > Material.TransverseVelocity)
                        x = sqrt(k2-kL2);
                        y = sqrt(kT2-k2);
                    elseif (numel(find(X(:,1) ~= 0)) <= 3 && H(p,1) < Material.TransverseVelocity) || (numel(find(X(:,1) ~= 0)) > 3 && X(end-1,1) < Material.TransverseVelocity)
                        x = sqrt(k2-kL2);
                        y = sqrt(k2-kT2);
                    end
                    y2 = kT2-k2;
                    xRo = x*Ro;
                    yRo = y*Ro;
                    xRi = x*Ri;
                    yRi = y*Ri;
                    if  (numel(find(X(:,1) ~= 0)) <= 3 && H(p,1) > Material.LongitudinalVelocity) || (numel(find(X(:,1) ~= 0)) > 3 && X(end-1,1) > Material.LongitudinalVelocity)
                        Z0xi = besselj(0,xRi);
                        Z0xo = besselj(0,xRo);
                        Z0yi = besselj(0,yRi);
                        Z0yo = besselj(0,yRo);
                        W0xi = bessely(0,xRi);
                        W0xo = bessely(0,xRo);
                        W0yi = bessely(0,yRi);
                        W0yo = bessely(0,yRo);
                        Z1xiRi = xRi.*besselj(1,xRi);
                        Z1xoRo = xRo.*besselj(1,xRo);
                        Z1yiRi = yRi.*besselj(1,yRi);
                        Z1yoRo = yRo.*besselj(1,yRo);
                        W1xiRi = xRi.*bessely(1,xRi);
                        W1xoRo = xRo.*bessely(1,xRo);
                        W1yiRi = yRi.*bessely(1,yRi);
                        W1yoRo = yRo.*bessely(1,yRo);
                    elseif (numel(find(X(:,1) ~= 0)) <= 3 && H(p,1) < Material.LongitudinalVelocity && H(p,1) > Material.TransverseVelocity) || (numel(find(X(:,1) ~= 0)) > 3 && X(end-1,1) < Material.LongitudinalVelocity && X(end-1,1) > Material.TransverseVelocity)
                        Z0xi = besseli(0,xRi);
                        Z0xo = besseli(0,xRo);
                        Z0yi = besselj(0,yRi);
                        Z0yo = besselj(0,yRo);
                        W0xi = besselk(0,xRi);
                        W0xo = besselk(0,xRo);
                        W0yi = bessely(0,yRi);
                        W0yo = bessely(0,yRo);
                        Z1xiRi = -xRi.*besseli(1,xRi);
                        Z1xoRo = -xRo.*besseli(1,xRo);
                        Z1yiRi = yRi.*besselj(1,yRi);
                        Z1yoRo = yRo.*besselj(1,yRo);
                        W1xiRi = xRi.*besselk(1,xRi);
                        W1xoRo = xRo.*besselk(1,xRo);
                        W1yiRi = yRi.*bessely(1,yRi);
                        W1yoRo = yRo.*bessely(1,yRo);
                    elseif (numel(find(X(:,1) ~= 0)) <= 3 && H(p,1) < Material.TransverseVelocity) || (numel(find(X(:,1) ~= 0)) > 3 && X(end-1,1) < Material.TransverseVelocity)
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
                    end
                    if  ToggleInnerFluid && ToggleOuterFluid
                        if  Symmetric
                            zRi = -sqrt(kInnerFluid2-k2)*Ri;
                            zRo = -sqrt(kOuterFluid2-k2)*Ro;
                        else
                            if  X(end-1,1) < SlowFluidVelocity
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
                    elseif ToggleInnerFluid && ~ToggleOuterFluid
                        zRi = -sqrt(kInnerFluid2-k2)*Ri;
                    elseif ~ToggleInnerFluid && ToggleOuterFluid
                        zRo = -sqrt(kOuterFluid2-k2)*Ro;
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
                                end
                                Y(j,l) = abs(det(M));
                            end
                        end
                    end
                    if  abs(SweepRangeImag(end)) < 1e-3
                        Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                    end
% if  p==1&&i>=267
% if  abs(SweepRangeImag(end)) < 1e-3
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y))
% else
% f = figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y))
% end
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
                    if  ~isempty(b1) % one or multiple minima are found
                        if  isscalar(b1)
                            MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                        else
                            delta = [];
                            for l = 1:length(b1) % calculate which minimum lies closest to previous solution
                                cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
                                delta(l) = abs(cp(l)-X(end-1,1));
                            end
                            [~,l] = min(delta);
                            MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                        end
                    else
                        Min = zeros(size(Y)); % find border minima
                        for j = 2:size(Y,1)-1
                            if  Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)
                                Min(j,1) = 1;
                            end
                            if  Y(j,size(Y,2)) < Y(j-1,size(Y,2)) && Y(j,size(Y,2)) < Y(j+1,size(Y,2)) && Y(j,size(Y,2)) < Y(j,size(Y,2)-1)
                                Min(j,size(Y,2)) = 1;
                            end
                        end
                        for j = 2:size(Y,2)-1
                            if  Y(1,j) < Y(1,j-1) && Y(1,j) < Y(1,j+1) && Y(1,j) < Y(2,j)
                                Min(1,j) = 1;
                            end
                            if  Y(size(Y,1),j) < Y(size(Y,1),j-1) && Y(size(Y,1),j) < Y(size(Y,1),j+1) && Y(size(Y,1),j) < Y(size(Y,1)-1,j)
                                Min(size(Y,1),j) = 1;
                            end
                        end
                        [b1,b2] = find(Min);
                        if  ~isempty(b1) % one or multiple BORDER minima are found
                            if  isscalar(b1)
                                MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                            else
                                Value = [];
                                for l = 1:length(b1)
                                    Value(l) = Y(b1(l),b2(l));
                                end
                                [~,l] = min(Value);
                                MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                            end
                        else
                            if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal x SweepRangeImag)
                                SweepRangeReal = [SweepRangeReal(1)+(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                if  SweepRangeReal(2) < 0
                                    SweepRangeReal(2) = 0;
                                end
                                if  SweepRangeImag(2) > -1e-10
                                    SweepRangeImag(2) = -1e-10;
                                end
                            end
                            MIN = 0;
                            break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                        end
                    end
                    if  k == 100 || (Resolution > abs(SweepRangeReal(1)-SweepRangeReal(2)) && Resolution > abs(SweepRangeImag(1)-SweepRangeImag(2)))
                        break
                    end
                    if  MIN(1) == 1 || MIN(2) == 1 || MIN(1) == size(Y,1) || MIN(2) == size(Y,2) % set the new search area around the found minimum
                        if  MIN(1) == 1
                            if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                                SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                            end
                            SweepRangeReal = [SweepRangeReal(1)+4*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)];
                        elseif MIN(2) == 1
                            if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                            end
                            SweepRangeImag = [SweepRangeImag(1)-4*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
                        elseif MIN(1) == size(Y,1)
                            if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                                SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                            end
                            SweepRangeReal = [SweepRangeReal(1) SweepRangeReal(end)-4*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                        elseif MIN(2) == size(Y,2)
                            if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                            end
                            SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                        end
                    else
                        if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                            if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                                SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                            end
                        elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                                SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                            end
                            if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                            end
                        elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                            end
                        end
                    end
                    if  SweepRangeReal(2) < 0
                        SweepRangeReal(2) = 0;
                    end
                    if  SweepRangeImag(2) > -1e-10
                        SweepRangeImag(2) = -1e-10;
                    end
                end
                if  any(MIN)
                    if  numel(find(X(:,1) ~= 0)) < 4
                        Outlier = 0;
                    else
                        z1 = isoutlier(vertcat(X(1:end-1,3),SweepRangeReal(MIN(1))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                        z2 = isoutlier(vertcat(X(1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                        if  (z1(end) && abs(X(end-1,3)-SweepRangeReal(MIN(1))) > 1) || (z2(end) && abs(X(end-1,4)-SweepRangeImag(MIN(2))) > 1)
                            Outlier = 1;
                        else
                            Outlier = 0;
                        end
                    end
                    if  ~Outlier || all(X(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                        X(i,1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
                        X(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
                        X(i,3) = SweepRangeReal(MIN(1)); % real velocity (m/s)
                        X(i,4) = SweepRangeImag(MIN(2)); % imaginary velocity (m/s)
                        Misses(i) = 0;
                        break
                    end
                    if  numel(find(X(:,1) ~= 0)) < 4
                        if  Hdiff(p) > .001*X(end-1,3)
                            SweepRangeReal = [1.001*X(end-1,3) .999*X(end-1,3)];
                        else
                            SweepRangeReal = [X(end-1,3)+Hdiff(p) X(end-1,3)-Hdiff(p)];
                        end
                    else
                        SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                    end
                    if  SweepRangeReal(1) == SweepRangeReal(2)
                        SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
                    end
                    if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                        SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                        SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                    end
                    if  SweepRangeReal(2) < 0
                        SweepRangeReal(2) = 0;
                    end
                    if  numel(find(X(:,1) ~= 0)) < 4
                        SweepRangeImag= [1.1*X(end-1,4) .9*X(end-1,4)];
                    else
                        SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary velocity corresponds to the attenuation
                        if  -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                            SweepRangeImag(2) = -1e-10;
                        end
                    end
                    if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                        SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                        SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                    end
                    if  SweepRangeImag(2) > -1e-10
                        SweepRangeImag(2) = -1e-10;
                    end
                end
            end
            if  X(i,1) > 0 % stop q-loop if minimum has been found
                break
            end
        end
        if  X(i,1) == 0 % fit phase velocity, attenuation, and imaginary velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
            if  isscalar(find(X(:,1) > 0))
                MissingModes(p) = 1;
                break
            else
                Smooth1 = filloutliers(X(1:i-1,1),'spline','movmedian',5,'ThresholdFactor',1);
                Smooth2 = filloutliers(X(1:i-1,2),'spline','movmedian',5,'ThresholdFactor',1);
                Smooth3 = filloutliers(X(1:i-1,3),'spline','movmedian',5,'ThresholdFactor',1);
                Smooth4 = filloutliers(X(1:i-1,4),'spline','movmedian',5,'ThresholdFactor',1);
                Fit1 = fit(FrequencyRange_Fundamental(1:i-1)',Smooth1,'cubicspline');
                Fit2 = fit(FrequencyRange_Fundamental(1:i-1)',Smooth2,'cubicspline');
                Fit3 = fit(FrequencyRange_Fundamental(1:i-1)',Smooth3,'cubicspline');
                Fit4 = fit(FrequencyRange_Fundamental(1:i-1)',Smooth4,'cubicspline');
                X(i,1) = Fit1(FrequencyRange_Fundamental(i));
                X(i,2) = Fit2(FrequencyRange_Fundamental(i));
                X(i,3) = Fit3(FrequencyRange_Fundamental(i));
                X(i,4) = Fit4(FrequencyRange_Fundamental(i));
                if  X(i,2) < 0 % exclude negative attenuation
                    X(i,[2 4]) = 0;
                end
                Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
            end
        end
        if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % stop tracing the current dispersion curve if more than 'MissingSamples' have been missed
            X(end-MissingSamples:end,:) = [];
            Misses(end-MissingSamples:end) = 0;
            break
        end
        if  X(end,1) > FastFluidVelocity-PhaseVelocityOffset || (~(~Sink && OuterFluid.Velocity > InnerFluid.Velocity) && ((X(1) < SlowFluidVelocity && X(end,1) > SlowFluidVelocity-PhaseVelocityOffset) || (X(1) > SlowFluidVelocity && X(end,1) < SlowFluidVelocity+PhaseVelocityOffset)))
            X(end,:) = [];
            if  height(X) <= 10
                MissingModes(p) = 1;
            end
            break
        end
        if  Multithreading
            send(Q,[FrequencyRange_Fundamental(i),X(i)/1e3,p])
        else
            addpoints(g(p),FrequencyRange_Fundamental(i),X(i)/1e3);
            drawnow limitrate
        end
% if  p == 1
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange_Fundamental(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
% end
    end
    if  ~MissingModes(p)
        X(Misses(1:height(X)) == 1,:) = NaN;
        LScholte{p}(:,1) = FrequencyRange_Fundamental(1:height(X));
        LScholte{p}(:,2) = FrequencyRange_Fundamental(1:height(X))/1e3;
        LScholte{p}(:,3) = FrequencyRange_Fundamental(1:height(X))*(Ro-Ri);
        LScholte{p}(:,4) = fillmissing(X(:,1),'spline')/1e3;
        LScholte{p}(:,6) = fillmissing(X(:,2),'spline');
        LScholte{p}(:,7) = fillmissing(X(:,3),'spline'); % real velocity (m/s)
        LScholte{p}(:,8) = fillmissing(X(:,4),'spline'); % imaginary velocity (m/s)
        LScholte{p}(LScholte{p}(:,6) < 0,6) = 0; % exclude negative attenuation
% LScholte{p}(:,5) = smooth(((LScholte{p}(:,4)).^2)./(LScholte{p}(:,4)-LScholte{p}(:,1).*differentiate(fit(LScholte{p}(:,1),LScholte{p}(:,4),'cubicspline'),LScholte{p}(:,1))));
    else
        if  p == 1
            LScholte{1}(1,1) = 0;
        else
            LScholte{p} = LScholte{p-1};
        end
    end
end
LScholte(MissingModes == 1) = [];
for p = 1:length(LScholte)
    LScholte{p}(:,7:8) = [];
    LScholte{p} = flipud(LScholte{p});
end
if  isempty(LScholte)
    LScholte{1} = [];
end