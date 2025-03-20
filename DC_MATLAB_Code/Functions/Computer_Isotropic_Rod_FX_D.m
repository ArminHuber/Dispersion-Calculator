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
function [FScholte_,F] = Computer_Isotropic_Rod_FX_D(Multithreading,Q,ax,FluidLoading,Fluid,Viscoelastic,Material,FrequencyRange,R,FScholte,FLambF,Resolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS>
global Stop
Stop = 0;
F{1}=[];FScholte_{1}=[];
if  ~Multithreading
    g(1) = animatedline(ax,'color','b');
    g(2) = animatedline(ax,'color','b');
end
R2 = R^2;
Density = Fluid.Density/Material.Density;
if  FluidLoading && Viscoelastic
    T = Fluid.Density*Fluid.Velocity/(Fluid.Density*Fluid.Velocity+Material.Density*Material.PlateVelocity)+Material.LongitudinalAttenuation+Material.TransverseAttenuation;
elseif FluidLoading && ~Viscoelastic
    T = Fluid.Density*Fluid.Velocity/(Fluid.Density*Fluid.Velocity+Material.Density*Material.PlateVelocity);
elseif ~FluidLoading && Viscoelastic
    T = Material.LongitudinalAttenuation+Material.TransverseAttenuation;
end
for p = 1:2
    if  p == 1
        if  FScholte(1,1) == 0
            FX{1}(length(FrequencyRange),8) = 0;
            continue
        end
        X = FScholte(1,:);
    elseif p == 2
        X = FLambF(1,:);
        Misses = 0;
    end
    if  Multithreading
        send(Q,[FrequencyRange(1),X(1)/1e3,1])
    else
        addpoints(g(p),FrequencyRange(1),X(1)/1e3);
        drawnow limitrate
    end
    for i = 2:length(FrequencyRange)
        if  Stop
            return
        end
        AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
        kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
        kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
        if  FluidLoading
            kF2 = (AngularFrequency/Fluid.Velocity)^2;
        end
        X(i,1) = 0;
        Neighbors = [];
        if  p == 2
            if  i <= height(FX{1})
                Neighbors = [FX{1}(i,7) FX{1}(i,8)];
            end
            NeighborsNumber = 1;
        end
        for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
            if  p == 1
                if  i <= 5
                    SweepRangeReal = [5*X(end-1,3) X(end-1,3)];
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
                if  i <= 5
                    SweepRangeImag = [5*X(end-1,4) X(end-1,4)];
                else
                    SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                    if  all(SweepRangeImag == [0 0])
                        SweepRangeImag = [-2*T 0];
                    elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                        SweepRangeImag(2) = 0;
                    end
                end
                if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                    SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                    SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                end
                if  SweepRangeImag(2) > 0
                    SweepRangeImag(2) = 0;
                end
            elseif p == 2
                if  i == 2
                    if  FluidLoading
                        SweepRangeReal = [15*X(end-1,3) 5*X(end-1,3)];
                    else
                        SweepRangeReal = [5*X(end-1,3) X(end-1,3)];
                    end
                elseif i == 3
                    if  FluidLoading
                        SweepRangeReal = [2.5*X(end-1,3) 1.5*X(end-1,3)];
                    else
                        SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
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
                if  i == 2
                    if  FluidLoading
                        SweepRangeImag = [15*X(end-1,4) 5*X(end-1,4)];
                    else
                        SweepRangeImag = [5*X(end-1,4) X(end-1,4)];
                    end
                elseif i == 3
                    if  FluidLoading
                        SweepRangeImag = [2.5*X(end-1,4) 1.5*X(end-1,4)];
                    else
                        SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)];
                    end
                else
                    SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                    if  all(SweepRangeImag == [0 0])
                        SweepRangeImag = [-20*T 0];
                    elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                        SweepRangeImag(2) = 0;
                    end
                end
                if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                    SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                    SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                end
                if  SweepRangeImag(2) > 0
                    SweepRangeImag(2) = 0;
                end
                if  i > 2 && SweepRangeImag(1) > min(X(:,4)) && Scholte(i-1)
                    SweepRangeImag(1) = min(X(:,4));
                end
            end
            for o = 1:SearchAreaSections % increase search resolution
                if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                    break
                end
                for k = 1:1e2
                    if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                        if  p == 1
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
                        elseif p == 2
                            SweepRangeReal = SweepRangeReal(1):.25/q/o*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                            SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
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
                    if  ~isempty(Neighbors)
                        for l = 2:width(Wavenumber)-1 % remove solutions of previously found lower modes
                            for j = 2:height(Wavenumber)-1
                                for n = 1:NeighborsNumber
                                    if  SweepRangeReal(j-1) > Neighbors(n,1) && SweepRangeReal(j+1) < Neighbors(n,1) && SweepRangeImag(l-1) < Neighbors(n,2) && SweepRangeImag(l+1) > Neighbors(n,2)
                                        Wavenumber(j,l) = NaN;
                                    end
                                end
                            end
                        end
                    end
                    k2 = Wavenumber.^2;
                    y2 = kT2-k2;
                    xR = sqrt(kL2-k2)*R;
                    yR = sqrt(kT2-k2)*R;
                    if  FluidLoading
                        if  p == 2
                            for j = 2:height(Wavenumber)-1 % avoid jumping to the Scholte mode
                                if  SweepRangeReal(j-1) > Fluid.Velocity && SweepRangeReal(j+1) < Fluid.Velocity && SweepRangeImag(end) == 0
                                    Wavenumber(j,end) = NaN;
                                end
                            end
                        end
                        J1x = besselj(1,xR);
                        J1y = besselj(1,yR);
                        dJ1xR = J1x-xR.*besselj(2,xR);
                        dJ1yR = J1y-yR.*besselj(2,yR);
                        zR = sqrt(kF2-k2)*R;
                        if  p == 1 && X(end-1,1) < Fluid.Velocity
                            zR = -zR;
                        end
                        H1z = besselh(1,zR);
                        dH1zR = H1z-zR.*besselh(2,zR);
                        Y = NaN(size(Wavenumber));
                        for l = 1:width(Wavenumber)
                            for j = 1:height(Wavenumber)
                                if  ~isnan(Wavenumber(j,l))
                                    M(1,1) = dJ1xR(j,l)+(.5*kT2*R2-k2(j,l)*R2-1)*J1x(j,l);
                                    M(1,2) = -Wavenumber(j,l)*(dJ1yR(j,l)+(y2(j,l)*R2-1)*J1y(j,l));
                                    M(1,3) = dJ1yR(j,l)-J1y(j,l);
                                    M(1,4) = .5*kT2*R2*Density*H1z(j,l);
                                    M(2,1) = 2*(dJ1xR(j,l)-J1x(j,l));
                                    M(2,2) = 2*Wavenumber(j,l)*(J1y(j,l)-dJ1yR(j,l));
                                    M(2,3) = 2*dJ1yR(j,l)+(y2(j,l)*R2-2)*J1y(j,l);
                                    M(3,1) = 2*Wavenumber(j,l)*dJ1xR(j,l);
                                    M(3,2) = (y2(j,l)-k2(j,l))*dJ1yR(j,l);
                                    M(3,3) = -Wavenumber(j,l)*J1y(j,l);
                                    M(4,1) = dJ1xR(j,l);
                                    M(4,2) = -Wavenumber(j,l)*dJ1yR(j,l);
                                    M(4,3) = -J1y(j,l);
                                    M(4,4) = dH1zR(j,l);
                                    Y(j,l) = abs(det(M));
                                end
                            end
                        end
                    else
                        kR2 = k2*R2;
                        kR4 = kR2.^2;
                        yR2 = y2*R2;
                        yR4 = yR2.^2;
                        ZxR = xR.*besselj(0,xR)./besselj(1,xR);
                        ZyR = yR.*besselj(0,yR)./besselj(1,yR);
                        ZyR2 = ZyR.^2;
                        F1 = 2*(yR2-kR2).^2;
                        F2 = 2*yR2.*(yR2+5*kR2);
                        F3 = yR2.^3-10*yR4-2*yR4.*kR2+2*yR2.*kR2+yR2.*kR4-4*kR4;
                        F4 = 2*yR2.*(2*yR2.*kR2-yR2-9*kR2);
                        F5 = yR2.*(-yR4+8*yR2-2*yR2.*kR2+8*kR2-kR4);
                        Y = abs(F1+F2.*ZxR./ZyR+F3./ZyR+F4.*ZxR./ZyR2+F5./ZyR2); 
                    end
                    for l = 1:width(Wavenumber)
                        for j = 2:height(Wavenumber)-1
                            if  cp(j-1,l) > Material.TransverseVelocity && cp(j+1,l) < Material.TransverseVelocity
                                Y(j,l) = NaN;
                            end
                        end
                    end
                    if  abs(SweepRangeImag(end)) < 1e-3
                        Y(:,end+1) = max(Y,[],'all'); % add wall below zero attenuation to allow finding minimum at zero attenuation
                    end
% if i==3
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
                            cp = [];
                            for l = 1:length(b1) % calculate which minimum lies closest to last solution
                                cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
                                delta(l) = abs(cp(l)-X(end-1,1));
                            end
                            if  p == 1 && i > 2 && any(cp > Fluid.Velocity) && any(cp < Fluid.Velocity) && X(end-1,1)+(X(end-1,1)-X(end-2,1)) > Fluid.Velocity
                                l = find(cp > Fluid.Velocity);
                                MIN = [b1(l(1)) b2(l(1))];
                            else
                                [~,l] = min(delta);
                                MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                            end
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
                            if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                                SweepRangeReal = [SweepRangeReal(1)+(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                SweepRangeImag = [SweepRangeImag(1)-(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                if  SweepRangeReal(2) < 0
                                    SweepRangeReal(2) = 0;
                                end
                                if  SweepRangeImag(2) > 0
                                    SweepRangeImag(2) = 0;
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
                    if  SweepRangeImag(2) > 0
                        SweepRangeImag(2) = 0;
                    end
                end
                if  any(MIN)
                    if  i < 4
                        if  (p == 1 && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) > Fluid.Velocity) || (p == 2 && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) > Material.PlateVelocity)
                            Outlier = 1;
                        else
                            Outlier = 0;
                        end
                    else
                        z1 = isoutlier(vertcat(X(1:end-1,3),SweepRangeReal(MIN(1))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                        if  p == 1
                            if  Viscoelastic % to overcome the point where the mode jumps above the fluid velocity in a perfectly elastic case; damping rise significantly from zero, which is normally filtered out as an outlier
                                z2 = isoutlier(vertcat(X(1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                            end
                            if  (z1(end) && abs(X(end-1,3)-SweepRangeReal(MIN(1))) > 1) || (Viscoelastic && (z2(end) && abs(X(end-1,4)-SweepRangeImag(MIN(2))) > 1))
                                Outlier = 1;
                            else
                                Outlier = 0;
                            end
                        elseif p == 2
                            z2 = isoutlier(vertcat(X(1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
                            z3 = isoutlier(vertcat(X(1:end-1,1),AngularFrequency/real(Wavenumber(MIN(1),MIN(2)))),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
                            if  (z1(end) && abs(X(end-1,3)-SweepRangeReal(MIN(1))) > 1) || (z2(end) && abs(X(end-1,4)-SweepRangeImag(MIN(2))) > 1) || (z3(end) && abs(X(end-1,1)-AngularFrequency/real(Wavenumber(MIN(1),MIN(2)))) > 1)
                                Outlier = 1;
                            else
                                Outlier = 0;
                            end
                        end
                    end
                    if  ~Outlier || all(X(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                        X(i,1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
                        X(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
                        X(i,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
                        X(i,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                        Misses(i) = 0;
                        if  p == 2
                            if  SweepRangeImag(end) == 0
                                Scholte(i) = 1;
                            else
                                Scholte(i) = 0;
                            end
                        end
                        break
                    end
                    if  p == 1
                        if  i <= 5
                            SweepRangeReal = [5*X(end-1,3) X(end-1,3)];
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
                        if  i <= 5
                            SweepRangeImag = [5*X(end-1,4) X(end-1,4)];
                        else
                            SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                            if  all(SweepRangeImag == [0 0])
                                SweepRangeImag = [-2*T 0];
                            elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                                SweepRangeImag(2) = 0;
                            end
                        end
                        if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                            SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                            SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                        end
                        if  SweepRangeImag(2) > 0
                            SweepRangeImag(2) = 0;
                        end
                    elseif p == 2
                        if  i == 2
                            if  FluidLoading
                                SweepRangeReal = [15*X(end-1,3) 5*X(end-1,3)];
                            else
                                SweepRangeReal = [5*X(end-1,3) X(end-1,3)];
                            end
                        elseif i == 3
                            if  FluidLoading
                                SweepRangeReal = [2.5*X(end-1,3) 1.5*X(end-1,3)];
                            else
                                SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
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
                        if  i == 2
                            if  FluidLoading
                                SweepRangeImag = [15*X(end-1,4) 5*X(end-1,4)];
                            else
                                SweepRangeImag = [5*X(end-1,4) X(end-1,4)];
                            end
                        elseif i == 3
                            if  FluidLoading
                                SweepRangeImag = [2.5*X(end-1,4) 1.5*X(end-1,4)];
                            else
                                SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)];
                            end
                        else
                            SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                            if  all(SweepRangeImag == [0 0])
                                SweepRangeImag = [-20*T 0];
                            elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                                SweepRangeImag(2) = 0;
                            end
                        end
                        if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                            SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                            SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                        end
                        if  SweepRangeImag(2) > 0
                            SweepRangeImag(2) = 0;
                        end
                        if  i > 2 && SweepRangeImag(1) > min(X(:,4)) && Scholte(i-1)
                            SweepRangeImag(1) = min(X(:,4));
                        end
                    end
                end
            end
            if  X(i,1) > 0 % stop q-loop if minimum has been found
                break
            end
        end
        if  X(i,1) == 0 % fit phase velocity, attenuation, and imaginary phase velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
            Smooth1 = filloutliers(X(1:i-1,1),'spline','movmedian',5,'ThresholdFactor',1);
            Smooth2 = filloutliers(X(1:i-1,2),'spline','movmedian',5,'ThresholdFactor',1);
            Smooth3 = filloutliers(X(1:i-1,3),'spline','movmedian',5,'ThresholdFactor',1);
            Smooth4 = filloutliers(X(1:i-1,4),'spline','movmedian',5,'ThresholdFactor',1);
            Fit1 = fit(FrequencyRange(1:i-1)',Smooth1,'cubicspline');
            Fit2 = fit(FrequencyRange(1:i-1)',Smooth2,'cubicspline');
            Fit3 = fit(FrequencyRange(1:i-1)',Smooth3,'cubicspline');
            Fit4 = fit(FrequencyRange(1:i-1)',Smooth4,'cubicspline');
            X(i,1) = Fit1(FrequencyRange(i));
            X(i,2) = Fit2(FrequencyRange(i));
            X(i,3) = Fit3(FrequencyRange(i));
            X(i,4) = Fit4(FrequencyRange(i));
            if  X(i,2) < 0 % negative attenuation is impossible
                X(i,[2 4]) = 0;
            end
            Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
            if  p == 2
                Scholte(i) = 0;
            end
        end
        if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % if more than 'MissingSamples' sample points are missing, the damped tracing stops, and only the nondamped tracing continues to have that data for the next mode
            X(end-MissingSamples:end,:) = [];
            Misses(end-MissingSamples:end) = 0;
            break
        end
        if  Multithreading
            send(Q,[FrequencyRange(i),X(i)/1e3,1])
        else
            addpoints(g(p),FrequencyRange(i),X(i)/1e3);
            drawnow limitrate
        end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
% if  Misses(i) == 1
%     disp(['p = ',num2str(p),' f = ',num2str(FrequencyRangeF(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' Miss'])
% end
    end
    X(Misses(1:height(X)) == 1,:) = NaN;
    FX{p}(:,1) = FrequencyRange(1:height(X));
    FX{p}(:,2) = FrequencyRange(1:height(X))/1e3;
    FX{p}(:,3) = FrequencyRange(1:height(X))*2*R;
    FX{p}(:,4) = fillmissing(X(:,1),'spline')/1e3;
    FX{p}(:,6) = fillmissing(X(:,2),'spline');
    FX{p}(:,7) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
    FX{p}(:,8) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
    FX{p}(FX{p}(:,6) < 0,6) = 0; % negative attenuation is impossible
end
if  FScholte(1,1) == 0
    F{1} = FX{2};
    FScholte_ = [];
else
    if  max(FX{1}(:,4)) > max(FX{2}(:,4))
        F{1} = FX{1};
        FScholte_{1} = FX{2};
    else
        F{1} = FX{2};
        FScholte_{1} = FX{1};
    end
end