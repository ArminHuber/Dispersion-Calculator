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
function L = Computer_Isotropic_Rod_L_D(Multithreading,Q1,Q2,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRange,R,HigherOrderModes,FLambL,H,FrequencyResolution,PhaseVelocityLimit,Resolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS> 
global Stop 
Stop = 0;
L{1} = [];
if  ~Multithreading && HigherOrderModes && any(H)
    for p = 1:length(H)+2
        g(p) = animatedline(ax,'color','r');
        g1(p) = animatedline(ax,'color','r');
    end
else
    g = animatedline(ax,'color','r');
    g(2) = animatedline(ax,'color','r');
end
R2 = R^2;
Density = Fluid.Density/Material.Density;
X0 = Material.CylinderVelocity;
if  FluidLoading && Viscoelastic
    T = Fluid.Density*Fluid.Velocity/(Fluid.Density*Fluid.Velocity+Material.Density*X0)+Material.LongitudinalAttenuation+Material.TransverseAttenuation;
elseif FluidLoading && ~Viscoelastic
    T = Fluid.Density*Fluid.Velocity/(Fluid.Density*Fluid.Velocity+Material.Density*X0);
elseif ~FluidLoading && Viscoelastic
    T = Material.LongitudinalAttenuation+Material.TransverseAttenuation;
end
% if  isempty(FLambL)
    L{1}(length(FrequencyRange),8) = 0;
    MissingModes(1) = 1;
% else
%     MissingModes(1) = 0;
%     X = FLambL(1,:);
%     if  Multithreading
%         send(Q1,[FrequencyRange(1),X(1)/1e3,1])
%     else
%         addpoints(g(1),FrequencyRange(1),X(1)/1e3);
%         drawnow limitrate
%     end
% %     Misses = 0;
% %     Scholte = 0;
%     for i = 2:length(FrequencyRange)
%         if  Stop
%             return
%         end
%         AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
%         kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
%         kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
%         kF2 = (AngularFrequency/Fluid.Velocity)^2;
%         X(i,1) = 0;
% %         Neighbors = [];
% %         if  i <= height(L{1})
% %             Neighbors = [L{1}(i,7) L{1}(i,8)];
% %         end
% %         NeighborsNumber = 1;
%         for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
%             if  i == 2
%                 SweepRangeReal = [15*X(end-1,3) 5*X(end-1,3)];
%             elseif i == 3
%                 SweepRangeReal = [2.5*X(end-1,3) 1.5*X(end-1,3)];
%             else
%                 SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
%             end
%             if  SweepRangeReal(1) == SweepRangeReal(2)
%                 SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
%             end
%             if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
%                 SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
%                 SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
%             end
%             if  SweepRangeReal(2) < 0
%                 SweepRangeReal(2) = 0;
%             end
%             if  i == 2
%                 SweepRangeImag = [15*X(end-1,4) 5*X(end-1,4)];
%             elseif i == 3
%                 SweepRangeImag = [2.5*X(end-1,4) 1.5*X(end-1,4)];
%             else
%                 SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
%                 if  all(SweepRangeImag == [0 0])
%                     SweepRangeImag = [20*T 0];
%                 elseif SweepRangeImag(1)-SweepRangeImag(2) > SweepRangeImag(2)
%                     SweepRangeImag(2) = 0;
%                 end
%             end
%             if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
%                 SweepRangeImag(1) = SweepRangeImag(1)+Resolution;
%                 SweepRangeImag(2) = SweepRangeImag(2)-Resolution;
%             end
%             if  SweepRangeImag(2) < 0
%                 SweepRangeImag(2) = 0;
%             end
%             if  i > 2 && SweepRangeImag(1) < max(X(:,4)) && Scholte(i-1)
%                 SweepRangeImag(1) = max(X(:,4));
%             end
%             for o = 1:SearchAreaSections % increase search resolution
%                 if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
%                     break
%                 end
%                 for k = 1:1e2
%                     if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
%                         SweepRangeReal = SweepRangeReal(1):.25/q/o*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
%                         SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
%                     else
%                         if  length(SweepRangeReal) == 2
%                             SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
%                         end
%                         if  length(SweepRangeImag) == 2
%                             SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
%                         end
%                     end
%                     if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
%                         break
%                     end
%                     Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
                    % cp = AngularFrequency./real(Wavenumber);
%                     for j = 2:height(Wavenumber)-1 % avoid jumping to the Scholte mode
%                         if  SweepRangeReal(j-1) > Fluid.Velocity && SweepRangeReal(j+1) < Fluid.Velocity && SweepRangeImag(end) == 0
%                             Wavenumber(j,end) = NaN;
%                         end
%                     end
                    % for l = 1:width(Wavenumber)
                    %     for j = 2:height(Wavenumber)-1
                    %         if  cp(j-1,l) > Material.TransverseVelocity && cp(j+1,l) < Material.TransverseVelocity
                    %             Wavenumber(j,l) = NaN;
                    %         end
                    %     end
                    % end
% %                     if  ~isempty(Neighbors)
% %                         for l = 2:width(Wavenumber)-1 % remove solutions of previously found lower modes
% %                             for j = 2:height(Wavenumber)-1
% %                                 for m = 1:NeighborsNumber
% %                                     if  Neighbors(m,2) ~= 0
% %                                         if  SweepRangeReal(j-1) > Neighbors(m,1) && SweepRangeReal(j+1) < Neighbors(m,1) && SweepRangeImag(l-1) > Neighbors(m,2) && SweepRangeImag(l+1) < Neighbors(m,2)
% %                                             Wavenumber(j,l) = NaN;
% %                                         end
% %                                     else
% %                                         if  SweepRangeReal(j-1) > Neighbors(m,1) && SweepRangeReal(j+1) < Neighbors(m,1)
% %                                             Wavenumber(j,l) = NaN;
% %                                         end
% %                                     end
% %                                 end
% %                             end
% %                         end
% %                     end
%                     k2 = Wavenumber.^2;
%                     y2 = kT2-k2;
%                     xR = sqrt(kL2-k2)*R;
%                     yR = sqrt(kT2-k2)*R;
%                     J0x = besselj(0,xR);
%                     J0y = besselj(0,yR);
%                     J1xR = xR.*besselj(1,xR);
%                     J1yR = yR.*besselj(1,yR);
%                     zR = sqrt(kF2-k2)*R;
%                     if  X(end-1,1) < Fluid.Velocity
%                         zR = -zR;
%                     end
%                     H0z = besselh(0,zR);
%                     H1zR = -zR.*besselh(1,zR);
%                     Y = NaN(size(Wavenumber));
%                     for l = 1:width(Wavenumber)
%                         for j = 1:height(Wavenumber)
%                             if  ~isnan(Wavenumber(j,l))
%                                 M(1,1) = R2*(.5*kT2-k2(j,l))*J0x(j,l)-J1xR(j,l);
%                                 M(1,2) = Wavenumber(j,l)*(J1yR(j,l)-y2(j,l)*R2*J0y(j,l));
%                                 M(1,3) = .5*kT2*R2*Density*H0z(j,l);
%                                 M(2,1) = -2*Wavenumber(j,l)*J1xR(j,l);
%                                 M(2,2) = (k2(j,l)-y2(j,l))*J1yR(j,l);
%                                 M(3,1) = -J1xR(j,l);
%                                 M(3,2) = Wavenumber(j,l)*J1yR(j,l);
%                                 M(3,3) = H1zR(j,l);
%                                 Y(j,l) = abs(det(M));
%                             end
%                         end
%                     end
%                     if  abs(SweepRangeImag(end)) < 1e-3
%                         Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
%                     end
% % if  i<=3&& k == 1
% % if  abs(SweepRangeImag(end)) < 1e-3
% % f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y))    
% % else
% % f = figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y))
% % end
% % close(f)
% % end
%                     Min = zeros(size(Y));
%                     for l = 2:size(Y,2)-1
%                         for j = 2:size(Y,1)-1
%                             if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
%                                 Min(j,l) = 1;
%                             end
%                         end
%                     end
%                     [b1,b2] = find(Min);
%                     if  ~isempty(b1) % one or multiple minima are found
%                         if  isscalar(b1)
%                             MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
%                         else
%                             delta = [];
%                             for l = 1:length(b1) % calculate which minimum lies closest to last solution
%                                 cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
%                                 if  i == 1
%                                     delta(l) = abs(cp(l)-X0);
%                                 else
%                                     delta(l) = abs(cp(l)-X(end-1,1));
%                                 end
%                             end
%                             [~,l] = min(delta);
%                             MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
%                         end
%                     else
% %                         if  i <= 3 && k == 1
% %                             if  q < SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
% %                                 SweepRangeReal = [SweepRangeReal(1)+q*o*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)];
% %                                 SweepRangeImag = [SweepRangeImag(1)+q*o*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
% %                             else
% %                                 SweepRangeReal = [SweepRangeReal(1)+(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end))];
% %                                 SweepRangeImag = [SweepRangeImag(1)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
% %                                 if  SweepRangeReal(2) < 0
% %                                     SweepRangeReal(2) = 0;
% %                                 end
% %                                 if  SweepRangeImag(2) < 0
% %                                     SweepRangeImag(2) = 0;
% %                                 end
% %                             end
% %                             MIN = 0;
% %                             break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
% %                         else
%                             Min = zeros(size(Y)); % find border minima
%                             for j = 2:size(Y,1)-1
%                                 if  Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)
%                                     Min(j,1) = 1;
%                                 end
%                                 if  Y(j,size(Y,2)) < Y(j-1,size(Y,2)) && Y(j,size(Y,2)) < Y(j+1,size(Y,2)) && Y(j,size(Y,2)) < Y(j,size(Y,2)-1)
%                                     Min(j,size(Y,2)) = 1;
%                                 end
%                             end
%                             for j = 2:size(Y,2)-1
%                                 if  Y(1,j) < Y(1,j-1) && Y(1,j) < Y(1,j+1) && Y(1,j) < Y(2,j)
%                                     Min(1,j) = 1;
%                                 end
%                                 if  Y(size(Y,1),j) < Y(size(Y,1),j-1) && Y(size(Y,1),j) < Y(size(Y,1),j+1) && Y(size(Y,1),j) < Y(size(Y,1)-1,j)
%                                     Min(size(Y,1),j) = 1;
%                                 end
%                             end
%                             [b1,b2] = find(Min);
%                             if  ~isempty(b1) % one or multiple BORDER minima are found
%                                 if  isscalar(b1)
%                                     MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
%                                 else
%                                     Value = [];
%                                     for l = 1:length(b1)
%                                         Value(l) = Y(b1(l),b2(l));
%                                     end
%                                     [~,l] = min(Value);
%                                     MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
%                                 end
%                             else
%                                 if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
%                                     SweepRangeReal = [SweepRangeReal(1)+(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end))];
%                                     SweepRangeImag = [SweepRangeImag(1)+(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)-(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end))];
%                                     if  SweepRangeReal(2) < 0
%                                         SweepRangeReal(2) = 0;
%                                     end
%                                     if  SweepRangeImag(2) < 0
%                                         SweepRangeImag(2) = 0;
%                                     end
%                                 end
%                                 MIN = 0;
%                                 break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
%                             end
% %                         end
%                     end
%                     if  k == 100 || (Resolution > abs(SweepRangeReal(1)-SweepRangeReal(2)) && Resolution > abs(SweepRangeImag(1)-SweepRangeImag(2)))
%                         break
%                     end
%                     if  MIN(1) == 1 || MIN(2) == 1 || MIN(1) == size(Y,1) || MIN(2) == size(Y,2) % set the new search area around the found minimum
%                         if  MIN(1) == 1
%                             if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
%                                 SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
%                             end
%                             SweepRangeReal = [SweepRangeReal(1)+4*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)];
%                         elseif MIN(2) == 1
%                             if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
%                                 SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
%                             end
%                             SweepRangeImag = [SweepRangeImag(1)-4*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
%                         elseif MIN(1) == size(Y,1)
%                             if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
%                                 SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
%                             end
%                             SweepRangeReal = [SweepRangeReal(1) SweepRangeReal(end)-4*abs(SweepRangeReal(1)-SweepRangeReal(end))];
%                         elseif MIN(2) == size(Y,2)
%                             if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
%                                 SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
%                             end
%                             SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*abs(SweepRangeImag(1)-SweepRangeImag(end))];
%                         end
%                     else
%                         if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
%                             if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
%                                 SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
%                             end
%                         elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
%                             if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
%                                 SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
%                             end
%                             if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
%                                 SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
%                             end
%                         elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
%                             if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
%                                 SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
%                             end
%                         end
%                     end
%                     if  SweepRangeReal(2) < 0
%                         SweepRangeReal(2) = 0;
%                     end
%                     if  SweepRangeImag(2) < 0
%                         SweepRangeImag(2) = 0;
%                     end
%                 end
%                 if  any(MIN)
%                     if  i < 4
%                         if  AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) > Material.PlateVelocity
%                             Outlier = 1;
%                         else
%                             Outlier = 0;
%                         end
%                     else
%                         z1 = isoutlier(vertcat(X(1:end-1,3),SweepRangeReal(MIN(1))),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
%                         z2 = isoutlier(vertcat(X(1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
%                         z3 = isoutlier(vertcat(X(1:end-1,1),AngularFrequency/real(Wavenumber(MIN(1),MIN(2)))),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
%                         if  (z1(end) && abs(X(end-1,3)-SweepRangeReal(MIN(1))) > 1) || (z2(end) && abs(X(end-1,4)-SweepRangeImag(MIN(2))) > 1) || (z3(end) && abs(X(end-1,1)-AngularFrequency/real(Wavenumber(MIN(1),MIN(2)))) > 1)
%                             Outlier = 1;
%                         else
%                             Outlier = 0;
%                         end
%                     end
%                     if  ~Outlier || all(X(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
%                         X(i,1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
%                         X(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
%                         X(i,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
%                         X(i,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
%                         Misses(i) = 0;
%                         if  SweepRangeImag(end) == 0
%                             Scholte(i) = 1;
%                         else
%                             Scholte(i) = 0;
%                         end
%                         break
%                     end
%                     if  i == 2
%                         SweepRangeReal = [15*X(end-1,3) 5*X(end-1,3)];
%                     elseif i == 3
%                         SweepRangeReal = [2.5*X(end-1,3) 1.5*X(end-1,3)];
%                     else
%                         SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
%                     end
%                     if  SweepRangeReal(1) == SweepRangeReal(2)
%                         SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
%                     end
%                     if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
%                         SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
%                         SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
%                     end
%                     if  SweepRangeReal(2) < 0
%                         SweepRangeReal(2) = 0;
%                     end
%                     if  i == 2
%                         SweepRangeImag = [15*X(end-1,4) 5*X(end-1,4)];
%                     elseif i == 3
%                         SweepRangeImag = [2.5*X(end-1,4) 1.5*X(end-1,4)];
%                     else
%                         SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
%                         if  all(SweepRangeImag == [0 0])
%                             SweepRangeImag = [20*T 0];
%                         elseif SweepRangeImag(1)-SweepRangeImag(2) > SweepRangeImag(2)
%                             SweepRangeImag(2) = 0;
%                         end
%                     end
%                     if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
%                         SweepRangeImag(1) = SweepRangeImag(1)+Resolution;
%                         SweepRangeImag(2) = SweepRangeImag(2)-Resolution;
%                     end
%                     if  SweepRangeImag(2) < 0
%                         SweepRangeImag(2) = 0;
%                     end
%                     if  i > 2 && SweepRangeImag(1) < max(X(:,4)) && Scholte(i-1)
%                         SweepRangeImag(1) = max(X(:,4));
%                     end
%                 end
%             end
%             if  X(i,1) > 0 % stop q-loop if minimum has been found
%                 break
%             end
%         end
%         if  X(i,1) == 0 % fit phase velocity, attenuation, and imaginary phase velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
%             Smooth1 = filloutliers(X(1:i-1,1),'spline','movmedian',5,'ThresholdFactor',1);
%             Smooth2 = filloutliers(X(1:i-1,2),'spline','movmedian',5,'ThresholdFactor',1);
%             Smooth3 = filloutliers(X(1:i-1,3),'spline','movmedian',5,'ThresholdFactor',1);
%             Smooth4 = filloutliers(X(1:i-1,4),'spline','movmedian',5,'ThresholdFactor',1);
%             Fit1 = fit(FrequencyRange(1:i-1)',Smooth1,'cubicspline');
%             Fit2 = fit(FrequencyRange(1:i-1)',Smooth2,'cubicspline');
%             Fit3 = fit(FrequencyRange(1:i-1)',Smooth3,'cubicspline');
%             Fit4 = fit(FrequencyRange(1:i-1)',Smooth4,'cubicspline');
%             X(i,1) = Fit1(FrequencyRange(i));
%             X(i,2) = Fit2(FrequencyRange(i));
%             X(i,3) = Fit3(FrequencyRange(i));
%             X(i,4) = Fit4(FrequencyRange(i));
% %             if  X(i,2) < 0 % negative attenuation is impossible
% %                 X(i,[2 4]) = 0;
% %             end
%             Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
%             Scholte(i) = 0;
%         end
%         if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % if more than 'MissingSamples' sample points are missing, the damped tracing stops, and only the nondamped tracing continues to have that data for the next mode
%             X(end-MissingSamples:end,:) = [];
%             Misses(end-MissingSamples:end) = 0;
%             break
%         end
%         if  Multithreading
%             send(Q1,[FrequencyRange(i),X(i)/1e3,1])
%         else
%             addpoints(g(1),FrequencyRange(i),X(i)/1e3);
%             drawnow limitrate
%         end   
% % String = ['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% % if  Misses(i) == 1 
% %     String = append(String,' Miss');
% % end
% % disp(String)
% % if  Misses(i) == 1
% %     disp(['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' Miss'])   
% % end
%     end
%     X(Misses(1:height(X)) == 1,:) = NaN;
%     L{1}(:,1) = FrequencyRange(1:height(X));
%     L{1}(:,2) = FrequencyRange(1:height(X))/1e3;
%     L{1}(:,3) = FrequencyRange(1:height(X))*2*R;
%     L{1}(:,4) = fillmissing(X(:,1),'spline')/1e3;
%     L{1}(:,6) = fillmissing(X(:,2),'spline');
%     L{1}(:,7) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
%     L{1}(:,8) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
% %     L{1}(L{1}(:,6) < 0,6) = 0; % negative attenuation is impossible
% end
X = [];
Misses = 0;
for i = 1:length(FrequencyRange)
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
    for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
        if  i < 4
            SweepRangeReal = [1.01*X0 .99*X0];
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
        if  i == 1
            SweepRangeImag = [-10*T 0]; % the imaginary phase velocity corresponds to the attenuation
        elseif i == 2 || i == 3
            SweepRangeImag = [-20*T 0]; % the imaginary phase velocity corresponds to the attenuation
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
        for o = 1:SearchAreaSections % increase search resolution
            if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                break
            end
            for k = 1:1e2
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
                k2 = Wavenumber.^2;
                x = sqrt(kL2-k2);
                y = sqrt(kT2-k2);
                y2 = kT2-k2;
                xR = x*R;
                yR = y*R;
                if  FluidLoading
                    for j = 2:height(Wavenumber)-1 % avoid jumping to the Scholte mode
                        if  SweepRangeReal(j-1) > Fluid.Velocity && SweepRangeReal(j+1) < Fluid.Velocity && SweepRangeImag(end) == 0
                            Wavenumber(j,end) = NaN;
                        end
                    end
                    J0x = besselj(0,xR);
                    J0y = besselj(0,yR);
                    J1xR = xR.*besselj(1,xR);
                    J1yR = yR.*besselj(1,yR);
                    zR = sqrt(kF2-k2)*R;
                    H0z = besselh(0,zR);
                    H1zR = -zR.*besselh(1,zR);
                    Y = NaN(size(Wavenumber));
                    for l = 1:width(Wavenumber)
                        for j = 1:height(Wavenumber)
                            if  ~isnan(Wavenumber(j,l))
                                M(1,1) = R2*(.5*kT2-k2(j,l))*J0x(j,l)-J1xR(j,l);
                                M(1,2) = Wavenumber(j,l)*(J1yR(j,l)-y2(j,l)*R2*J0y(j,l));
                                M(1,3) = .5*kT2*R2*Density*H0z(j,l);
                                M(2,1) = -2*Wavenumber(j,l)*J1xR(j,l);
                                M(2,2) = (k2(j,l)-y2(j,l))*J1yR(j,l);
                                M(3,1) = -J1xR(j,l);
                                M(3,2) = Wavenumber(j,l)*J1yR(j,l);
                                M(3,3) = H1zR(j,l);
                                Y(j,l) = abs(det(M));
                            end
                        end
                    end
                else
                    Y = abs(2*x/R.*(y2+k2)-(y2-k2).^2.*besselj(0,xR)./besselj(1,xR)-4*k2.*x.*y.*besselj(0,yR)./besselj(1,yR));
                end
                for l = 1:width(Wavenumber)
                    for j = 2:height(Wavenumber)-1
                        if  cp(j-1,l) > Material.TransverseVelocity && cp(j+1,l) < Material.TransverseVelocity
                            Y(j,l) = NaN;
                        end
                    end
                end
                if  abs(SweepRangeImag(end)) < 1e-3
                    Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                end
% if  i==10 && k == 8
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
                        for l = 1:length(b1) % calculate which minimum lies closest to last solution
                            cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
                            if  i == 1
                                delta(l) = abs(cp(l)-X0);
                            else
                                delta(l) = abs(cp(l)-X(end-1,1));
                            end
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
                    X(i,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
                    X(i,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                    Misses(i) = 0;
                    break
                end
                if  i < 4
                    SweepRangeReal = [1.01*X0 .99*X0];
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
                if  i == 1
                    SweepRangeImag = [-10*T 0]; % the imaginary phase velocity corresponds to the attenuation
                elseif i == 2 || i == 3
                    SweepRangeImag = [-20*T 0]; % the imaginary phase velocity corresponds to the attenuation
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
    end
    if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % if more than 'MissingSamples' sample points are missing, the damped tracing stops, and only the nondamped tracing continues to have that data for the next mode
        X(end-MissingSamples:end,:) = [];
        Misses(end-MissingSamples:end) = 0;
        break
    end
    if  FluidLoading && i > 2 && X(end,1) < Fluid.Velocity % avoid jumping to Scholte modes
        X(end,:) = [];
        break
    end
    if  Multithreading
        send(Q1,[FrequencyRange(i),X(i)/1e3,1])
    else
        addpoints(g(2),FrequencyRange(i),X(i)/1e3);
        drawnow limitrate
    end
% String = ['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
% if  Misses(i) == 1
%     disp(['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' Miss'])   
% end
end
X(Misses(1:height(X)) == 1,:) = NaN;
L{2}(:,1) = FrequencyRange(1:height(X));
L{2}(:,2) = FrequencyRange(1:height(X))/1e3;
L{2}(:,3) = FrequencyRange(1:height(X))*2*R;
L{2}(:,4) = fillmissing(X(:,1),'spline')/1e3;
L{2}(:,6) = fillmissing(X(:,2),'spline');
L{2}(:,7) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
L{2}(:,8) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
L{2}(L{2}(:,6) < 0,6) = 0; % negative attenuation is impossible
if  HigherOrderModes && any(H)
    MissingModes(length(H)+2) = 0;
    for p = 1:length(H)
        X = [];
        Misses = 0;
        BelowCutoff = 0;
        for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange)
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
            for j = 1:length(L)
                if  FrequencyRange(i) <= L{j}(end,1)
                    z = find(abs(FrequencyRange(i)-L{j}(:,1)) == min(abs(FrequencyRange(i)-L{j}(:,1))));
                    Neighbors(j,:) = [L{j}(z(1),7) L{j}(z(1),8)];
                end
            end
            NeighborsNumber = height(Neighbors);
            for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                if  all(X(:,1) == 0) 
                    SweepRangeReal = [1.1*PhaseVelocityLimit X0];
                elseif isscalar(find(X(:,1) ~= 0))
                    SweepRangeReal = [PhaseVelocityLimit X0];
                elseif numel(find(X(:,1) ~= 0)) == 2
                    SweepRangeReal = [1.1*X(end-1,3) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
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
                if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < X0
                    SweepRangeReal(2) = X0;
                elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                    SweepRangeReal(2) = 0;
                end
                if  all(X(:,1) == 0)
                    SweepRangeImag = [-1000*T 0];
                elseif isscalar(find(X(:,1) ~= 0))
                    SweepRangeImag = [4*SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; 
                else
                    SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                end
                if  all(SweepRangeImag == [0 0])
                    SweepRangeImag = [-20*T 0];
                elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                    SweepRangeImag(2) = 0;
                end
                if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                    SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                    SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                end
                if  SweepRangeImag(2) > 0
                    SweepRangeImag(2) = 0;
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
                        % cp = AngularFrequency./real(Wavenumber);
                        k2 = Wavenumber.^2;
                        x = sqrt(kL2-k2);
                        y = sqrt(kT2-k2);
                        y2 = kT2-k2;
                        xR = x*R;
                        yR = y*R;
                        if  FluidLoading
                            J0x = besselj(0,xR);
                            J0y = besselj(0,yR);
                            J1xR = xR.*besselj(1,xR);
                            J1yR = yR.*besselj(1,yR);
                            zR = sqrt(kF2-k2)*R;
                            H0z = besselh(0,zR);
                            H1zR = -zR.*besselh(1,zR);
                            Y = NaN(size(Wavenumber));
                            for l = 1:width(Wavenumber)
                                for j = 1:height(Wavenumber)
                                    M(1,1) = R2*(.5*kT2-k2(j,l))*J0x(j,l)-J1xR(j,l);
                                    M(1,2) = Wavenumber(j,l)*(J1yR(j,l)-y2(j,l)*R2*J0y(j,l));
                                    M(1,3) = .5*kT2*R2*Density*H0z(j,l);
                                    M(2,1) = -2*Wavenumber(j,l)*J1xR(j,l);
                                    M(2,2) = (k2(j,l)-y2(j,l))*J1yR(j,l);
                                    M(3,1) = -J1xR(j,l);
                                    M(3,2) = Wavenumber(j,l)*J1yR(j,l);
                                    M(3,3) = H1zR(j,l);
                                    Y(j,l) = abs(det(M));
                                end
                            end
                        else
                            Y = abs(2*x/R.*(y2+k2)-(y2-k2).^2.*besselj(0,xR)./besselj(1,xR)-4*k2.*x.*y.*besselj(0,yR)./besselj(1,yR));
                        end
                        % for l = 1:width(Wavenumber)
                        %     for j = 2:height(Wavenumber)-1
                        %         if  cp(j-1,l) > Material.TransverseVelocity && cp(j+1,l) < Material.TransverseVelocity
                        %             Y(j,l) = NaN;
                        %         end
                        %     end
                        % end
                        if  numel(find(X(:,1) ~= 0)) <= 3
                            for l = 2:width(Wavenumber)-1 % remove solutions of previously found lower modes
                                for j = 2:height(Wavenumber)-1
                                    for n = 1:NeighborsNumber
                                        if  SweepRangeReal(j-1) > Neighbors(n,1) && SweepRangeReal(j+1) < Neighbors(n,1)    
                                            Y(j,l) = NaN;
                                        end
                                    end
                                end
                            end
                        else
                            for l = 2:width(Wavenumber)-1 % remove solutions of previously found lower modes
                                for j = 2:height(Wavenumber)-1
                                    for n = 1:NeighborsNumber
                                        if  SweepRangeReal(j-1) > Neighbors(n,1) && SweepRangeReal(j+1) < Neighbors(n,1) && SweepRangeImag(l-1) < Neighbors(n,2) && SweepRangeImag(l+1) > Neighbors(n,2)
                                            Y(j,l) = NaN;
                                        end
                                    end
                                end
                            end
                        end
                        if  abs(SweepRangeImag(end)) < 1e-3 && numel(find(X(:,1) ~= 0)) > 3
                            Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                        end
% if  p==1&&i>=609
% if  abs(SweepRangeImag(end)) < 1e-3 && numel(find(X(:,1) ~= 0)) > 3
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
                                for l = 1:length(b1) % calculate which minimum lies closest to last solution
                                    cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
                                    delta(l) = abs(cp(l)-X(end-1,1));
                                end
                                if  all(X(:,1) == 0)
                                    [~,l] = max(delta);
                                else
                                    [~,l] = min(delta);
                                end
                                MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                            end
                        else
                            if  numel(find(X(:,1) ~= 0)) <= 3 && k <= 2
                                if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                                    SweepRangeReal = [SweepRangeReal(1)+(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                    SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                    if  SweepRangeReal(1) > 1.1*PhaseVelocityLimit
                                        SweepRangeReal(1) = 1.1*PhaseVelocityLimit;
                                    end
                                    if  SweepRangeReal(2) < X0
                                        SweepRangeReal(2) = X0;
                                    end
                                    if  SweepRangeImag(2) > 0
                                        SweepRangeImag(2) = 0;
                                    end
                                end
                                MIN = 0;
                                break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
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
                                        SweepRangeReal = [SweepRangeReal(1)+(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                        SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
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
                        if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(1) > 1.1*PhaseVelocityLimit
                            SweepRangeReal(1) = 1.1*PhaseVelocityLimit;
                        end
                        if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < X0
                            SweepRangeReal(2) = X0;
                        elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                            SweepRangeReal(2) = 0;
                        end
                        if  SweepRangeImag(2) > 0
                            SweepRangeImag(2) = 0;
                        end
                    end
                    if  any(MIN)
                        if  (numel(find(X(:,1) ~= 0)) <= 3 && i <= height(L{p+1}) && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) < L{p+1}(i,4)*1e3) ||...
                            (numel(find(X(:,1) ~= 0)) >= 1 && numel(find(X(:,1) ~= 0)) <= 3 && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) > X(end-1,1)) ||...
                            SweepRangeReal(MIN(1)) == 0
                            Neighbors(end+1,:) = [SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))];
                            NeighborsNumber = height(Neighbors);
                        else
                            if  numel(find(X(:,1) ~= 0)) <= 5
                                Outlier = 0;
                            else
                                z1 = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:end-1,3),SweepRangeReal(MIN(1))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                                z2 = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                                if  (z1(end) && abs(X(end-1,3)-SweepRangeReal(MIN(1))) > 1) || (z2(end) && abs(X(end-1,4)-SweepRangeImag(MIN(2))) > 1)
                                    Outlier = 1;
                                else
                                    Outlier = 0;
                                end
                            end
                            if  ~Outlier || all(X(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                                X(i,1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
                                X(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
                                X(i,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
                                X(i,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                                Misses(i) = 0;
                                BelowCutoff(i) = 0;
                                break
                            end
                        end
                        if  all(X(:,1) == 0) 
                            SweepRangeReal = [1.1*PhaseVelocityLimit X0];
                        elseif isscalar(find(X(:,1) ~= 0))
                            SweepRangeReal = [PhaseVelocityLimit X0];
                        elseif numel(find(X(:,1) ~= 0)) == 2
                            SweepRangeReal = [1.1*X(end-1,3) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
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
                        if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < X0
                            SweepRangeReal(2) = X0;
                        elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                            SweepRangeReal(2) = 0;
                        end
                        if  all(X(:,1) == 0)
                            SweepRangeImag = [-1000*T 0];
                        elseif isscalar(find(X(:,1) ~= 0))
                            SweepRangeImag = [4*SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; 
                        else
                            SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                        end
                        if  all(SweepRangeImag == [0 0])
                            SweepRangeImag = [-20*T 0];
                        elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                            SweepRangeImag(2) = 0;
                        end
                        if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                            SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                            SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                        end
                        if  SweepRangeImag(2) > 0
                            SweepRangeImag(2) = 0;
                        end
                    end
                end
                if  X(i,1) > 0 % stop q-loop if minimum has been found
                    break
                end
            end
            if  X(i,1) == 0 && any(X(:,1)) % fit phase velocity, attenuation, and imaginary phase velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
                if  isscalar(find(X(:,1) > 0))
                    X(i-1,:) = 0;
                    Misses(i) = 0;
                    BelowCutoff(i) = 1;
                else
                    Smooth1 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,1),'spline','movmedian',5,'ThresholdFactor',1);
                    Smooth2 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,2),'spline','movmedian',5,'ThresholdFactor',1);
                    Smooth3 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,3),'spline','movmedian',5,'ThresholdFactor',1);
                    Smooth4 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,4),'spline','movmedian',5,'ThresholdFactor',1);
                    Fit1 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth1,'cubicspline');
                    Fit2 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth2,'cubicspline');
                    Fit3 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth3,'cubicspline');
                    Fit4 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth4,'cubicspline');
                    X(i,1) = Fit1(FrequencyRange(i));
                    X(i,2) = Fit2(FrequencyRange(i));
                    X(i,3) = Fit3(FrequencyRange(i));
                    X(i,4) = Fit4(FrequencyRange(i));
                    if  X(i,2) < 0 % negative attenuation is impossible
                        X(i,[2 4]) = 0;
                    end
                    Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
                    BelowCutoff(i) = 0;
                end
            elseif all(X(:,1) == 0) % if we are scanning below the cut-off frequency of the damped mode
                Misses(i) = 0;
                BelowCutoff(i) = 1;
            end
            if  FrequencyResolution*length(BelowCutoff(BelowCutoff == 1)) > BelowCutoffWidth/R/2e3 ||... % if we exceed the allowed scanning width (in kHz*mm) below the cut-off frequency of a damped mode without finding it; it can happen that the damped cut-off is higher than the non-damped from which we start the search, or there exists no damped mode for every non-damped one
                length(Misses) == ceil(H(p)/FrequencyResolution)+6 && numel(find(Misses == 1)) >= 2
                MissingModes(p+2) = 1;
                break
            end
            if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % if more than 'MissingSamples' sample points are missing, the damped tracing stops, and only the nondamped tracing continues to have that data for the next mode
                X(end-MissingSamples:end,:) = [];
                Misses(end-MissingSamples:end) = 0;
                break
            end
            if  FluidLoading && any(X(:,1)) && X(end,1) < Fluid.Velocity % avoid jumping to Scholte modes
                X(end,:) = [];
                break
            end
            if  Multithreading && X(i,1) > 0 && BelowCutoff(i) == 0
                send(Q1,[FrequencyRange(i),X(i)/1e3,p+2])
            elseif ~Multithreading && X(i,1) > 0 && BelowCutoff(i) == 0
                addpoints(g(p+2),FrequencyRange(i),X(i)/1e3);
                drawnow limitrate
            end
% if  p == 2
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRangeF(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
% end
% if  Misses(i) == 1
%     disp(['p = ',num2str(p),' f = ',num2str(FrequencyRangeF(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' Miss'])
% end
        end
        if  all(X(:,1) == 0)
            MissingModes(p+2) = 1;
        end
        if  ~MissingModes(p+2)
            [~,z] = max(X(:,1)); % find non-zero data below the cut-off frequency
            X(1:z-1,:) = 0; % remove them
            Misses(1:z-1) = 0; % remove also misses
            if  X(z,2) == 0 % if frequency cut-off is at zero damping, remove it also
                X(z,:) = 0;
            end
            X(Misses(1:height(X)) == 1,:) = NaN;
            L{p+2}(:,1) = FrequencyRange(1:height(X));
            L{p+2}(:,2) = FrequencyRange(1:height(X))/1e3;
            L{p+2}(:,3) = FrequencyRange(1:height(X))*2*R;
            L{p+2}(:,4) = fillmissing(X(:,1),'spline')/1e3;
            L{p+2}(:,6) = fillmissing(X(:,2),'spline');
            L{p+2}(:,7) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
            L{p+2}(:,8) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
            L{p+2}(L{p+2}(:,6) < 0,6) = 0; % negative attenuation is impossible
        else
            L{p+2} = L{p+1};
            X1{p+2}(1,1) = 0;
            continue
        end      
        if  max(L{p+2}(:,4))*1e3 > PhaseVelocityLimit
            X1{p+2}(1,1) = 0;
        else
            [Max,MaxInd] = max(L{p+2}(:,7));
            PhaseVelocityRange = Max+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
            for i = 1:length(PhaseVelocityRange)
                if  Stop
                    return
                end
                X1{p+2}(i,1) = 0;
                for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeFrq X SweepRangeImag)
                    if  i == 1
                        SweepRangeFrq = [L{p+2}(MaxInd,1)-FrequencyOffset/R/2e3 L{p+2}(MaxInd,1)+4/R/2e3];
                        SweepRangeImag = [SearchWidthImag(1)*L{p+2}(MaxInd,8) SearchWidthImag(2)*L{p+2}(MaxInd,8)];
                    else
                        SweepRangeFrq = [X1{p+2}(end-1,1)-FrequencyOffset/R/2e3 X1{p+2}(end-1,1)+4/R/2e3];
                        SweepRangeImag = [SearchWidthImag(1)*X1{p+2}(end-1,3) SearchWidthImag(2)*X1{p+2}(end-1,3)];
                    end
                    if  all(SweepRangeImag == [0 0])
                        SweepRangeImag = [-20*T 0];
                    end
                    if  SweepRangeImag(2) > 0
                        SweepRangeImag(2) = 0;
                    end
                    for o = 1:SearchAreaSections % increase search resolution
                        if  isempty(SweepRangeFrq) || isempty(SweepRangeImag)
                            break
                        end
                        for k = 1:1e2 % search minimum in characteristic equation and converge upon it 
                            if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                                if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                                    SweepRangeFrq = SweepRangeFrq(1):.25^o/q*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                                    SweepRangeImag = SweepRangeImag(1):.25/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) <= 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeFrq(1)-SweepRangeFrq(end)) >= .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    SweepRangeFrq = SweepRangeFrq(1):.25/q/o*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                                    SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) < .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    SweepRangeFrq = SweepRangeFrq(1):.25/q*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                                    SweepRangeImag = SweepRangeImag(1):.25^o/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                end
                            else
                                if  length(SweepRangeFrq) == 2
                                    SweepRangeFrq = SweepRangeFrq(1):.25*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                                end
                                if  length(SweepRangeImag) == 2
                                    SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                end
                            end
                            if  isempty(SweepRangeFrq) || isempty(SweepRangeImag)
                                break
                            end
                            AngularFrequency = 2*pi*SweepRangeFrq*1e3;
                            Wavenumber = AngularFrequency'./(PhaseVelocityRange(i)+SweepRangeImag*1i);
                            k2 = Wavenumber.^2;
                            kT2 = (AngularFrequency'/Material.TransverseVelocity_complex).^2;
                            x = sqrt((AngularFrequency'/Material.LongitudinalVelocity_complex).^2-k2);
                            y = sqrt(kT2-k2);
                            y2 = kT2-k2;
                            xR = x*R;
                            yR = y*R;
                            if  FluidLoading
                                J0x = besselj(0,xR);
                                J0y = besselj(0,yR);
                                J1xR = xR.*besselj(1,xR);
                                J1yR = yR.*besselj(1,yR);
                                zR = sqrt((AngularFrequency'/Fluid.Velocity).^2-k2)*R;
                                H0z = besselh(0,zR);
                                H1zR = -zR.*besselh(1,zR);
                                Y = NaN(size(Wavenumber));
                                for l = 1:width(Wavenumber)
                                    for j = 1:height(Wavenumber)
                                        M(1,1) = R2*(.5*kT2(j)-k2(j,l))*J0x(j,l)-J1xR(j,l);
                                        M(1,2) = Wavenumber(j,l)*(J1yR(j,l)-y2(j,l)*R2*J0y(j,l));
                                        M(1,3) = .5*kT2(j)*R2*Density*H0z(j,l);
                                        M(2,1) = -2*Wavenumber(j,l)*J1xR(j,l);
                                        M(2,2) = (k2(j,l)-y2(j,l))*J1yR(j,l);
                                        M(3,1) = -J1xR(j,l);
                                        M(3,2) = Wavenumber(j,l)*J1yR(j,l);
                                        M(3,3) = H1zR(j,l);
                                        Y(j,l) = abs(det(M));
                                    end
                                end
                            else
                                Y = abs(2*x/R.*(y2+k2)-(y2-k2).^2.*besselj(0,xR)./besselj(1,xR)-4*k2.*x.*y.*besselj(0,yR)./besselj(1,yR));
                            end
                            if  abs(SweepRangeImag(end)) < 1e-3 % in case the fluid has low density, this allows finding solutions at zero attenuation for A0
                                Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                            end
% if  p==4
% if  abs(SweepRangeImag(end)) < 1e-3
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeFrq,20*log10(Y))
% else
% f = figure;surf(SweepRangeImag,SweepRangeFrq,20*log10(Y))
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
                                    MIN = [b1 b2]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                                else
                                    delta = [];
                                    for l = 1:length(b1) % calculate which minimum lies closest to last solution
                                        frq(l) = AngularFrequency(b1(l))/(2*pi)/1e3;
                                        if  i == 1
                                            delta(l) = abs(frq(l)-L{p+2}(MaxInd,1));
                                        else
                                            delta(l) = abs(frq(l)-X1{p+2}(end-1,1));
                                        end
                                    end
                                    [~,l] = min(delta);
                                    MIN = [b1(l) b2(l)]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
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
                                        MIN = [b1 b2]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                                    else
                                        Value = [];
                                        for l = 1:length(b1)
                                            Value(l) = Y(b1(l),b2(l));
                                        end
                                        [~,l] = min(Value);
                                        MIN = [b1(l) b2(l)]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                                    end
                                else
                                    if  q > 1% && o == 1 % q > 1 runs are with extended search area (SweepRangeFrq X SweepRangeImag)
                                        SweepRangeFrq = [SweepRangeFrq(1)+(q-1)*abs(SweepRangeFrq(1)-SweepRangeFrq(end)) SweepRangeFrq(end)-(q-1)*abs(SweepRangeFrq(1)-SweepRangeFrq(end))];
                                        SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                        if  SweepRangeImag(2) > 0
                                            SweepRangeImag(2) = 0;
                                        end
                                    end
                                    MIN = 0;
                                    break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                                end
                            end 
                            if  k == 100 || (Resolution > abs(SweepRangeFrq(1)-SweepRangeFrq(2)) && Resolution > abs(SweepRangeImag(1)-SweepRangeImag(2)))
                                break
                            end
                            if  MIN(1) == 1 || MIN(2) == 1 || MIN(1) == size(Y,1) || MIN(2) == size(Y,2) % set the new search area around the found minimum
                                if  MIN(1) == 1
                                    if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeFrq(1)-SweepRangeFrq(end))
                                        SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                    end
                                    SweepRangeFrq = [SweepRangeFrq(1)+4*abs(SweepRangeFrq(1)-SweepRangeFrq(end)) SweepRangeFrq(end)];
                                elseif MIN(2) == 1
                                    if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                        SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                                    end
                                    SweepRangeImag = [SweepRangeImag(1)-4*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
                                elseif MIN(1) == size(Y,1)
                                    if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeFrq(1)-SweepRangeFrq(end))
                                        SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                    end
                                    SweepRangeFrq = [SweepRangeFrq(1) SweepRangeFrq(end)-4*abs(SweepRangeFrq(1)-SweepRangeFrq(end))];
                                elseif MIN(2) == size(Y,2)
                                    if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                        SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                                    end
                                    SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                end
                            else
                                if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                                    if  Resolution < abs(SweepRangeFrq(1)-SweepRangeFrq(2))
                                        SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                                    end
                                elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) <= 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeFrq(1)-SweepRangeFrq(end)) >= .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    if  Resolution < abs(SweepRangeFrq(1)-SweepRangeFrq(2))
                                        SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                                    end
                                    if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                        SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                    end
                                elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) < .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                        SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                    end
                                end
                            end
                            if  SweepRangeImag(2) > 0
                                SweepRangeImag(2) = 0;
                            end
                        end
                        if  any(MIN)
                            z = isoutlier(vertcat(X1{p+2}(1:end-1,1),AngularFrequency(MIN(1))/(2*pi)/1e3),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
                            if  ~z(end) || all(X1{p+2}(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                                X1{p+2}(i,4) = AngularFrequency(MIN(1))/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
                                X1{p+2}(i,1) = AngularFrequency(MIN(1))/(2*pi)/1e3; % frequency (kHz)
                                X1{p+2}(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
                                X1{p+2}(i,3) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                                break
                            end
                            if  i == 1
                                SweepRangeFrq = [L{p+2}(MaxInd,1)-FrequencyOffset/R/2e3 L{p+2}(MaxInd,1)+4/R/2e3];
                                SweepRangeImag = [SearchWidthImag(1)*L{p+2}(MaxInd,8) SearchWidthImag(2)*L{p+2}(MaxInd,8)];
                            else
                                SweepRangeFrq = [X1{p+2}(end-1,1)-FrequencyOffset/R/2e3 X1{p+2}(end-1,1)+4/R/2e3];
                                SweepRangeImag = [SearchWidthImag(1)*X1{p+2}(end-1,3) SearchWidthImag(2)*X1{p+2}(end-1,3)];
                            end
                            if  all(SweepRangeImag == [0 0])
                                SweepRangeImag = [-20*T 0];
                            end
                            if  SweepRangeImag(2) > 0
                                SweepRangeImag(2) = 0;
                            end
                        end
                    end
                    if  X1{p+2}(i) > 0 % stop q-loop if minimum has been found
                        break
                    end
                end
                if  X1{p+2}(i) == 0 || X1{p+2}(i,4) > PhaseVelocityLimit
                    break
                end            
                if  X1{p+2}(i,2) < 0 % negative attenuation is impossible
                    X1{p+2}(i,2:3) = 0;
                end
                if  Multithreading
                    send(Q2,[X1{p+2}(i,1),X1{p+2}(i,4)/1e3,p+2])
                else
                    addpoints(g1(p+2),X1{p+2}(i,1),X1{p+2}(i,4)/1e3);
                    drawnow limitrate
                end
% disp(['c = ',num2str(PhaseVelocityRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)]);
            end
            if  length(X1) < p+2
                X1{p+2}(1,1) = 0;
            end
            if  X1{p+2}(1,1) > 0 
                X1{p+2}(:,6) = X1{p+2}(:,2);
                X1{p+2}(:,2) = X1{p+2}(:,1)/1e3;
                X1{p+2}(:,3) = X1{p+2}(:,1)*2*R;
                X1{p+2}(:,4) = X1{p+2}(:,4)/1e3;
                X1{p+2}(X1{p+2}(:,1) == 0,:) = []; 
                X1{p+2} = flipud(X1{p+2});
            end
        end
    end
    L{1}(:,7:8) = [];
    L{2}(:,7:8) = [];
    for p = 3:length(L)
        L{p}(L{p}(:,4) == 0,:) = [];
        L{p}(:,7:8) = [];
        if  X1{p}(1,1) > 0
            L{p} = vertcat(X1{p},L{p});
        end
    end    
end
L(MissingModes == 1) = []; 