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
function FScholte_ = Computer_Isotropic_Pipe_FnScholte_D2(Multithreading,Q,ax,Symmetric,Viscoelastic,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRange,Ro,Ri,n,H,LineColor,FrequencyResolution,Resolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples)        
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
FScholte_{1} = [];
if  ~Multithreading
    for p = 1:height(H)
        if  ToggleOuterFluid || ToggleInnerFluid
            g(p) = animatedline(ax,'LineStyle','-.','color',LineColor);
        else
            g(p) = animatedline(ax,'color',LineColor);
        end
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
n2 = n^2;
Ro2 = Ro^2;
Ri2 = Ri^2;
Densityo = OuterFluid.Density/Material.Density;
Densityi = InnerFluid.Density/Material.Density;
MissingModes(height(H)) = 0;
Diff = diff(H(:,3));
Hdiff = min([[Diff;Inf] [Inf;Diff]],[],2);
Hdiff(abs(H(:,1)-Material.TransverseVelocity) < 1e-6) = [];
H(abs(H(:,1)-Material.TransverseVelocity) < 1e-6,:) = [];
H(abs(H(:,1)-SlowFluidVelocity) < 1e-6,:) = [];
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
            for j = 1:length(FScholte_)
                if  i <= height(FScholte_{j})
                    Neighbors(j,:) = [FScholte_{j}(i,7) FScholte_{j}(i,8)];
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
                    if  p > 1 % remove solutions of previously found lower modes
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
                    for j = 2:height(Wavenumber)-1
                        if  SweepRangeReal(j-1) > SlowFluidVelocity && SweepRangeReal(j+1) < SlowFluidVelocity
                            Wavenumber(j,end) = NaN;
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
                        Znxi = besselj(n,xRi);
                        Znxo = besselj(n,xRo);
                        Znyi = besselj(n,yRi);
                        Znyo = besselj(n,yRo);
                        Wnxi = bessely(n,xRi);
                        Wnxo = bessely(n,xRo);
                        Wnyi = bessely(n,yRi);
                        Wnyo = bessely(n,yRo);
                        dZnxiRi = n*Znxi-xRi.*besselj(n+1,xRi);
                        dZnxoRo = n*Znxo-xRo.*besselj(n+1,xRo);
                        dZnyiRi = n*Znyi-yRi.*besselj(n+1,yRi);
                        dZnyoRo = n*Znyo-yRo.*besselj(n+1,yRo);
                        dWnxiRi = n*Wnxi-xRi.*bessely(n+1,xRi);
                        dWnxoRo = n*Wnxo-xRo.*bessely(n+1,xRo);
                        dWnyiRi = n*Wnyi-yRi.*bessely(n+1,yRi);
                        dWnyoRo = n*Wnyo-yRo.*bessely(n+1,yRo);
                    elseif (numel(find(X(:,1) ~= 0)) <= 3 && H(p,1) < Material.LongitudinalVelocity && H(p,1) > Material.TransverseVelocity) || (numel(find(X(:,1) ~= 0)) > 3 && X(end-1,1) < Material.LongitudinalVelocity && X(end-1,1) > Material.TransverseVelocity)
                        Znxi = besseli(n,xRi);
                        Znxo = besseli(n,xRo);
                        Znyi = besselj(n,yRi);
                        Znyo = besselj(n,yRo);
                        Wnxi = besselk(n,xRi);
                        Wnxo = besselk(n,xRo);
                        Wnyi = bessely(n,yRi);
                        Wnyo = bessely(n,yRo);
                        dZnxiRi = n*Znxi+xRi.*besseli(n+1,xRi);
                        dZnxoRo = n*Znxo+xRo.*besseli(n+1,xRo);
                        dZnyiRi = n*Znyi-yRi.*besselj(n+1,yRi);
                        dZnyoRo = n*Znyo-yRo.*besselj(n+1,yRo);
                        dWnxiRi = n*Wnxi-xRi.*besselk(n+1,xRi);
                        dWnxoRo = n*Wnxo-xRo.*besselk(n+1,xRo);
                        dWnyiRi = n*Wnyi-yRi.*bessely(n+1,yRi);
                        dWnyoRo = n*Wnyo-yRo.*bessely(n+1,yRo);
                    elseif (numel(find(X(:,1) ~= 0)) <= 3 && H(p,1) < Material.TransverseVelocity) || (numel(find(X(:,1) ~= 0)) > 3 && X(end-1,1) < Material.TransverseVelocity)
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
                    if  abs(SweepRangeImag(end)) < 1e-3
                        Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                    end
% if  p==4&&i>=2
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
        if  ((ToggleOuterFluid || ToggleInnerFluid) && (X(end,1) > FastFluidVelocity-PhaseVelocityOffset || (~(~Sink && OuterFluid.Velocity > InnerFluid.Velocity) && ((X(1) < SlowFluidVelocity && X(end,1) > SlowFluidVelocity-PhaseVelocityOffset) || (X(1) > SlowFluidVelocity && X(end,1) < SlowFluidVelocity+PhaseVelocityOffset))))) ||...
            (~ToggleOuterFluid && ~ToggleInnerFluid && X(end,1) > Material.TransverseVelocity)
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
        FScholte_{p}(:,1) = FrequencyRange_Fundamental(1:height(X));
        FScholte_{p}(:,2) = FrequencyRange_Fundamental(1:height(X))/1e3;
        FScholte_{p}(:,3) = FrequencyRange_Fundamental(1:height(X))*(Ro-Ri);
        FScholte_{p}(:,4) = fillmissing(X(:,1),'spline')/1e3;
        FScholte_{p}(:,6) = fillmissing(X(:,2),'spline');
        FScholte_{p}(:,7) = fillmissing(X(:,3),'spline'); % real velocity (m/s)
        FScholte_{p}(:,8) = fillmissing(X(:,4),'spline'); % imaginary velocity (m/s)
        FScholte_{p}(FScholte_{p}(:,6) < 0,6) = 0; % exclude negative attenuation
    else
        if  p == 1
            FScholte_{1}(1,1) = 0;
        else
            FScholte_{p} = FScholte_{p-1};
        end
    end
end
FScholte_(MissingModes == 1) = [];
for p = 1:length(FScholte_)
    FScholte_{p}(:,7:8) = [];
    FScholte_{p} = flipud(FScholte_{p});
end
if  isempty(FScholte_)
    FScholte_{1} = [];
end