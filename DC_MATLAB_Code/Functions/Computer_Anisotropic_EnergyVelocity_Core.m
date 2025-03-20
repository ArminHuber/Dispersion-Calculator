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
function [X,Counter] = Computer_Anisotropic_EnergyVelocity_Core(X,ModeType,ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34)
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
for p = 1:length(X)
    PhaseVelocity = reshape(X{p}(:,4),1,1,[])*1e3;
    Attenuation = PhaseVelocity.*reshape(X{p}(:,7)./X{p}(:,1),1,1,[])/1e3; % Np/wavelength
    AngularFrequency = reshape(X{p}(:,1),1,1,[])*pi*2e3;
    AngularFrequency2 = AngularFrequency.^2;
    Wavenumber = AngularFrequency./PhaseVelocity.*(1+1i*Attenuation/2/pi);
    Wavenumber2 = Wavenumber.^2;
    Length = length(Wavenumber);
    KineticEnergyDensity = zeros(Height,1,Length);
    A = cell(0);
    if  ~Decoupled
        Wavenumber4 = Wavenumber2.^2;
        Wavenumber6 = Wavenumber2.^3;
        k3 = zeros(1,3,Length);
        v = zeros(Height,3,Length);
        epsilon = zeros(Height,6,Length);
        sigma = epsilon;
        for m = 1:SuperLayerSize
            rw2 = Material{m}.Density*AngularFrequency2;
            r2w4 = rw2.^2;
            A1 = a11(m)*Wavenumber2+a12(m)*rw2;
            A2 = a21(m)*Wavenumber4+a22(m)*rw2.*Wavenumber2+a23(m)*r2w4;
            A3 = a31(m)*Wavenumber6+a32(m)*rw2.*Wavenumber4+a33(m)*r2w4.*Wavenumber2+a34(m)*rw2.^3;
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
            m11 = c{m}(1,1)*Wavenumber2+c{m}(5,5)*k32-rw2;
            m22 = c{m}(6,6)*Wavenumber2+c{m}(4,4)*k32-rw2;
            m12 = c{m}(1,6)*Wavenumber2+c{m}(4,5)*k32;
            m13 = (c{m}(1,3)+c{m}(5,5))*k3k;
            m23 = (c{m}(3,6)+c{m}(4,5))*k3k;
            m1 = m13.*m22-m12.*m23;
            V = (m11.*m23-m13.*m12)./m1;
            W = (m11.*m22-m12.^2)./-m1;
            e1 = Wavenumber.*W+k3;
            e2 = k3.*V;
            e3 = k3.*W;
            D3 = 1i*((c{m}(1,3)+c{m}(3,6)*V).*Wavenumber+c{m}(3,3)*e3);
            D4 = 1i*(c{m}(4,5)*e1+c{m}(4,4)*e2);
            D5 = 1i*(c{m}(5,5)*e1+c{m}(4,5)*e2);
            E = exp(1i*k3*LayerThicknesses(m));
            L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
            A{m,3} = [ones(1,3,Length) E;V V.*E;W -W.*E;E ones(1,3,Length);V.*E V;W.*E -W]; % L2
            A{m,1} = pagemrdivide(L1,A{m,3}); % L
            A{m,4} = LayerThicknesses(m);
            A{m,5} = Material{m}.Density;
            A{m,6} = k3;
            A{m,7} = V;
            A{m,8} = W;
            A{m,9} = 1i*e3; % epsilon33
            A{m,10} = 1i*e2; % epsilon23
            A{m,11} = 1i*e1; % epsilon13
            A{m,12} = 1i*Wavenumber.*V; % epsilon12
            A{m,13} = 1i*((c{m}(1,1)+c{m}(1,6)*V).*Wavenumber+c{m}(1,3)*e3); % sigma11
            A{m,14} = 1i*((c{m}(1,2)+c{m}(2,6)*V).*Wavenumber+c{m}(2,3)*e3); % sigma22
            A{m,15} = D3; % sigma33
            A{m,16} = D4; % sigma23
            A{m,17} = D5; % sigma13
            A{m,18} = 1i*((c{m}(1,6)+c{m}(6,6)*V).*Wavenumber+c{m}(3,6)*e3); % sigma12
        end
        A = repmat(A,Repetitions,1);
        M{1} = A{1};
        for m = 2:SuperLayerSize
            M0 = A{m,1}(1:3,1:3,:)-M{1}(4:6,4:6,:);
            M1 = pagemrdivide(M{1}(1:3,4:6,:),M0);
            M2 = pagemrdivide(A{m,1}(4:6,1:3,:),M0);
            M{1} = [M{1}(1:3,1:3,:)+pagemtimes(M1,M{1}(4:6,1:3,:)) -pagemtimes(M1,A{m,1}(1:3,4:6,:));pagemtimes(M2,M{1}(4:6,1:3,:)) A{m,1}(4:6,4:6,:)-pagemtimes(M2,A{m,1}(1:3,4:6,:))];
        end
        for m = 1:length(Pattern)
            M0 = M{Pattern(m)}(1:3,1:3,:)-M{m}(4:6,4:6,:);
            M1 = pagemrdivide(M{m}(1:3,4:6,:),M0);
            M2 = pagemrdivide(M{Pattern(m)}(4:6,1:3,:),M0);
            M{m+1} = [M{m}(1:3,1:3,:)+pagemtimes(M1,M{m}(4:6,1:3,:)) -pagemtimes(M1,M{Pattern(m)}(1:3,4:6,:));pagemtimes(M2,M{m}(4:6,1:3,:)) M{Pattern(m)}(4:6,4:6,:)-pagemtimes(M2,M{Pattern(m)}(1:3,4:6,:))];
        end
        if  SymmetricSystem
            A = [A;flipud(A)];
            M0 = M{end}(4:6,4:6,:).*I-M{end}(4:6,4:6,:);
            M1 = pagemrdivide(M{end}(1:3,4:6,:),M0);
            M2 = pagemrdivide(M{end}(1:3,4:6,:).*I,M0);
            M{end} = [M{end}(1:3,1:3,:)+pagemtimes(M1,M{end}(4:6,1:3,:)) -pagemtimes(M1,M{end}(4:6,1:3,:).*I);pagemtimes(M2,M{end}(4:6,1:3,:)) M{end}(1:3,1:3,:).*I-pagemtimes(M2,M{end}(4:6,1:3,:).*I)];
        end
        if  FluidLoading
            M{end} = pageinv(M{end});
            k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-Wavenumber2);
            if  strcmp(ModeType,'Scholte') && max(PhaseVelocity) < UpperFluid.Velocity
                k3UpperFluid = -k3UpperFluid;
            end
            DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2./Wavenumber; % sigma11, sigma22, sigma33 in the upper fluid
            DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2./Wavenumber; % in the lower fluid
            ULowerFluid = -(k3UpperFluid./Wavenumber+M{end}(3,1,:).*DUpperFluid)./(M{end}(3,4,:).*DLowerFluid);
            uInterfaces = M{end}(1:3,1,:).*DUpperFluid+M{end}(1:3,4,:).*DLowerFluid.*ULowerFluid;
            uInterfaces(:,Layers+1,:) = M{end}(4:6,1,:).*DUpperFluid+M{end}(4:6,4,:).*DLowerFluid.*ULowerFluid;
        else
            Z1 = [pagemtimes(-M{end}(1:3,1:3,:),[ones(1,3,Length);A{1,7};-A{1,8}]) pagemtimes(-M{end}(1:3,4:6,:),[ones(1,3,Length);A{end,7};A{end,8}]);pagemtimes(M{end}(4:6,1:3,:),[ones(1,3,Length);A{1,7};-A{1,8}]) pagemtimes(M{end}(4:6,4:6,:),[ones(1,3,Length);A{end,7};A{end,8}])];
            Z2 = [pagemtimes(M{end}(1:3,1:3,:),[ones(1,3,Length);A{1,7};A{1,8}]);pagemtimes(-M{end}(4:6,1:3,:),[ones(1,3,Length);A{1,7};A{1,8}])];
            RT = pagemldivide(Z1,Z2(:,1,:));
            uInterfaces = [ones(1,1,Length);A{1,7}(1,1,:);A{1,8}(1,1,:)]+pagemtimes([ones(1,3,Length);A{1,7};-A{1,8}],RT(1:3,1,:));
            uInterfaces(:,Layers+1,:) = pagemtimes([ones(1,3,Length);A{end,7};A{end,8}],RT(4:6,1,:));
        end
        A{1,2} = A{1};
        for m = 2:Layers
            M0 = A{m,1}(1:3,1:3,:)-A{m-1,2}(4:6,4:6,:);
            M1 = pagemrdivide(A{m-1,2}(1:3,4:6,:),M0);
            M2 = pagemrdivide(A{m,1}(4:6,1:3,:),M0);
            A{m,2} = [A{m-1,2}(1:3,1:3,:)+pagemtimes(M1,A{m-1,2}(4:6,1:3,:)) -pagemtimes(M1,A{m,1}(1:3,4:6,:));pagemtimes(M2,A{m-1,2}(4:6,1:3,:)) A{m,1}(4:6,4:6,:)-pagemtimes(M2,A{m,1}(1:3,4:6,:))];
        end
        for m = Layers:-1:2
            M0 = A{m,1}(1:3,1:3,:)-A{m-1,2}(4:6,4:6,:);
            uInterfaces(:,m,:) = pagemtimes(pagemldivide(M0,A{m-1,2}(4:6,1:3,:)),uInterfaces(:,1,:))-pagemtimes(pagemldivide(M0,A{m,1}(1:3,4:6,:)),uInterfaces(:,m+1,:));
        end
        for m = 1:Layers
            x3 = (0:A{m,4}/SamplesX3:A{m,4})';
            if  m == 1
                x3Total = x3;
            else
                x3Total = [x3Total;x3Total(end)+x3];
            end
            r = (m-1)*length(x3)+1:m*length(x3);
            E = exp(1i*[A{m,6}.*x3 A{m,6}.*(A{m,4}-x3)]);
            U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]);
            v(r,1,:) = -1i*AngularFrequency.*pagemtimes(E,U); % v1
            v(r,2,:) = -1i*AngularFrequency.*pagemtimes([A{m,7} A{m,7}].*E,U); % v2
            v(r,3,:) = -1i*AngularFrequency.*pagemtimes([A{m,8} -A{m,8}].*E,U); % v3
            epsilon(r,1,:) = 1i*Wavenumber.*pagemtimes(E,U); % epsilon11
            epsilon(r,3,:) = pagemtimes([A{m,9} A{m,9}].*E,U); % epsilon33
            epsilon(r,4,:) = pagemtimes([A{m,10} -A{m,10}].*E,U); % epsilon23
            epsilon(r,5,:) = pagemtimes([A{m,11} -A{m,11}].*E,U); % epsilon13
            epsilon(r,6,:) = pagemtimes([A{m,12} A{m,12}].*E,U); % epsilon12
            sigma(r,1,:) = pagemtimes([A{m,13} A{m,13}].*E,U); % sigma11
            sigma(r,2,:) = pagemtimes([A{m,14} A{m,14}].*E,U); % sigma22
            sigma(r,3,:) = pagemtimes([A{m,15} A{m,15}].*E,U); % sigma33
            sigma(r,4,:) = pagemtimes([A{m,16} -A{m,16}].*E,U); % sigma23
            sigma(r,5,:) = pagemtimes([A{m,17} -A{m,17}].*E,U); % sigma13
            sigma(r,6,:) = pagemtimes([A{m,18} A{m,18}].*E,U); % sigma12
            KineticEnergyDensity(r,1,:) = A{m,5}*sum(abs(v(r,:,:)).^2,2)/2;
        end
        PowerFlowDensity = -real(sigma(:,1,:).*conj(v(:,1,:))+sigma(:,6,:).*conj(v(:,2,:))+sigma(:,5,:).*conj(v(:,3,:)))/2;
        PowerFlowDensity(:,2,:) = -real(sigma(:,6,:).*conj(v(:,1,:))+sigma(:,2,:).*conj(v(:,2,:))+sigma(:,4,:).*conj(v(:,3,:)))/2;
    elseif Decoupled && ~strcmp(ModeType,'Shear')
        Wavenumber4 = Wavenumber2.^2;
        k3 = zeros(1,2,Length);
        v = zeros(Height,2,Length);
        epsilon = zeros(Height,3,Length);
        sigma = epsilon;
        for m = 1:SuperLayerSize
            rw2 = Material{m}.Density*AngularFrequency2;
            A2 = a21(m)*Wavenumber2+a22(m)*rw2;
            A3 = a31(m)*Wavenumber4+a32(m)*rw2.*Wavenumber2+rw2.^2;
            d1 = sqrt(A2.^2-2*A1(m)*A3);
            k3(1,1,:) = sqrt((-A2+d1)/A1(m));
            k3(1,2,:) = sqrt((-A2-d1)/A1(m));
            W = (rw2-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k3);
            e1 = Wavenumber.*W+k3;
            e3 = k3.*W;
            D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*e3);
            D5 = 1i*c{m}(5,5)*e1;
            E = exp(1i*k3*LayerThicknesses(m));
            L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
            A{m,3} = [ones(1,2,Length) E;W -W.*E;E ones(1,2,Length);W.*E -W]; % L2
            A{m,1} = pagemrdivide(L1,A{m,3}); % L
            A{m,4} = LayerThicknesses(m);
            A{m,5} = Material{m}.Density;
            A{m,6} = k3;
            A{m,7} = W;
            A{m,8} = 1i*e3; % epsilon33
            A{m,9} = 1i*e1; % epsilon13
            A{m,10} = 1i*(c{m}(1,1)*Wavenumber+c{m}(1,3)*e3); % sigma11
            A{m,11} = D3; % sigma33
            A{m,12} = D5; % sigma13
        end
        A = repmat(A,Repetitions,1);
        M{1} = A{1};
        for m = 2:SuperLayerSize
            M0 = A{m,1}(1:2,1:2,:)-M{1}(3:4,3:4,:);
            M1 = pagemrdivide(M{1}(1:2,3:4,:),M0);
            M2 = pagemrdivide(A{m,1}(3:4,1:2,:),M0);
            M{1} = [M{1}(1:2,1:2,:)+pagemtimes(M1,M{1}(3:4,1:2,:)) -pagemtimes(M1,A{m,1}(1:2,3:4,:));pagemtimes(M2,M{1}(3:4,1:2,:)) A{m,1}(3:4,3:4,:)-pagemtimes(M2,A{m,1}(1:2,3:4,:))];
        end
        for m = 1:length(Pattern)
            M0 = M{Pattern(m)}(1:2,1:2,:)-M{m}(3:4,3:4,:);
            M1 = pagemrdivide(M{m}(1:2,3:4,:),M0);
            M2 = pagemrdivide(M{Pattern(m)}(3:4,1:2,:),M0);
            M{m+1} = [M{m}(1:2,1:2,:)+pagemtimes(M1,M{m}(3:4,1:2,:)) -pagemtimes(M1,M{Pattern(m)}(1:2,3:4,:));pagemtimes(M2,M{m}(3:4,1:2,:)) M{Pattern(m)}(3:4,3:4,:)-pagemtimes(M2,M{Pattern(m)}(1:2,3:4,:))];
        end
        if  SymmetricSystem
            A = [A;flipud(A)];
            M0 = M{end}(3:4,3:4,:).*I1-M{end}(3:4,3:4,:);
            M1 = pagemrdivide(M{end}(1:2,3:4,:),M0);
            M2 = pagemrdivide(M{end}(1:2,3:4,:).*I1,M0);
            M{end} = [M{end}(1:2,1:2,:)+pagemtimes(M1,M{end}(3:4,1:2,:)) -pagemtimes(M1,M{end}(3:4,1:2,:).*I1);pagemtimes(M2,M{end}(3:4,1:2,:)) M{end}(1:2,1:2,:).*I1-pagemtimes(M2,M{end}(3:4,1:2,:).*I1)];
        end
        if  FluidLoading
            M{end} = pageinv(M{end});
            k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-Wavenumber2);
            if  strcmp(ModeType,'Scholte') && max(PhaseVelocity) < UpperFluid.Velocity
                k3UpperFluid = -k3UpperFluid;
            end
            DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2./Wavenumber; % sigma11, sigma22, sigma33 in the upper fluid
            DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2./Wavenumber; % in the lower fluid
            ULowerFluid = -(k3UpperFluid./Wavenumber+M{end}(2,1,:).*DUpperFluid)./(M{end}(2,3,:).*DLowerFluid);
            uInterfaces = M{end}(1:2,1,:).*DUpperFluid+M{end}(1:2,3,:).*DLowerFluid.*ULowerFluid;
            uInterfaces(:,Layers+1,:) = M{end}(3:4,1,:).*DUpperFluid+M{end}(3:4,3,:).*DLowerFluid.*ULowerFluid;
        else
            Z1 = [pagemtimes(-M{end}(1:2,1:2,:),[ones(1,2,Length);-A{1,7}]) pagemtimes(-M{end}(1:2,3:4,:),[ones(1,2,Length);A{end,7}]);pagemtimes(M{end}(3:4,1:2,:),[ones(1,2,Length);-A{1,7}]) pagemtimes(M{end}(3:4,3:4,:),[ones(1,2,Length);A{end,7}])];
            Z2 = [pagemtimes(M{end}(1:2,1:2,:),[ones(1,2,Length);A{1,7}]);pagemtimes(-M{end}(3:4,1:2,:),[ones(1,2,Length);A{1,7}])];
            RT = pagemldivide(Z1,Z2(:,1,:));
            uInterfaces = [ones(1,1,Length);A{1,7}(1,1,:)]+pagemtimes([ones(1,2,Length);-A{1,7}],RT(1:2,1,:));
            uInterfaces(:,Layers+1,:) = pagemtimes([ones(1,2,Length);A{end,7}],RT(3:4,1,:));
        end
        A{1,2} = A{1};
        for m = 2:Layers
            M0 = A{m,1}(1:2,1:2,:)-A{m-1,2}(3:4,3:4,:);
            M1 = pagemrdivide(A{m-1,2}(1:2,3:4,:),M0);
            M2 = pagemrdivide(A{m,1}(3:4,1:2,:),M0);
            A{m,2} = [A{m-1,2}(1:2,1:2,:)+pagemtimes(M1,A{m-1,2}(3:4,1:2,:)) -pagemtimes(M1,A{m,1}(1:2,3:4,:));pagemtimes(M2,A{m-1,2}(3:4,1:2,:)) A{m,1}(3:4,3:4,:)-pagemtimes(M2,A{m,1}(1:2,3:4,:))];
        end
        for m = Layers:-1:2
            M0 = A{m,1}(1:2,1:2,:)-A{m-1,2}(3:4,3:4,:);
            uInterfaces(:,m,:) = pagemtimes(pagemldivide(M0,A{m-1,2}(3:4,1:2,:)),uInterfaces(:,1,:))-pagemtimes(pagemldivide(M0,A{m,1}(1:2,3:4,:)),uInterfaces(:,m+1,:));
        end
        for m = 1:Layers
            x3 = (0:A{m,4}/SamplesX3:A{m,4})';
            if  m == 1
                x3Total = x3;
            else
                x3Total = [x3Total;x3Total(end)+x3];
            end
            r = (m-1)*length(x3)+1:m*length(x3);
            E = exp(1i*[A{m,6}.*x3 A{m,6}.*(A{m,4}-x3)]);
            U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]);
            v(r,1,:) = -1i*AngularFrequency.*pagemtimes(E,U); % v1
            v(r,2,:) = -1i*AngularFrequency.*pagemtimes([A{m,7} -A{m,7}].*E,U); % v3
            epsilon(r,1,:) = 1i*Wavenumber.*pagemtimes(E,U); % epsilon11
            epsilon(r,2,:) = pagemtimes([A{m,8} A{m,8}].*E,U); % epsilon33
            epsilon(r,3,:) = pagemtimes([A{m,9} -A{m,9}].*E,U); % epsilon13
            sigma(r,1,:) = pagemtimes([A{m,10} A{m,10}].*E,U); % sigma11
            sigma(r,2,:) = pagemtimes([A{m,11} A{m,11}].*E,U); % sigma33
            sigma(r,3,:) = pagemtimes([A{m,12} -A{m,12}].*E,U); % sigma13
            KineticEnergyDensity(r,1,:) = A{m,5}*sum(abs(v(r,:,:)).^2,2)/2;
        end
        PowerFlowDensity = -real(sigma(:,1,:).*conj(v(:,1,:))+sigma(:,3,:).*conj(v(:,2,:)))/2;
    elseif Decoupled && strcmp(ModeType,'Shear')
        k3 = zeros(1,1,Length);
        v = zeros(Height,1,Length);
        epsilon = zeros(Height,2,Length);
        sigma = epsilon;
        for m = 1:SuperLayerSize
            k3(1,1,:) = sqrt((Material{m}.Density*AngularFrequency2-Wavenumber2*c{m}(6,6))/c{m}(4,4));
            if  k3(1) == 0
                k3(1,1,:) = 1e-10;
            end
            D4 = 1i*k3*c{m}(4,4);
            E = exp(1i*k3*LayerThicknesses(m));
            E2 = E.^2;
            A{m,1} = D4./(E2-ones(1,1,Length)).*[-ones(1,1,Length)-E2 2*E;-2*E ones(1,1,Length)+E2]; % L
            A{m,3} = [ones(1,1,Length) E;E ones(1,1,Length)]; % L2
            A{m,4} = LayerThicknesses(m);
            A{m,5} = Material{m}.Density;
            A{m,6} = 1i*k3; % epsilon23
            A{m,7} = D4; % sigma23
            A{m,8} = 1i*Wavenumber*c{m}(6,6); % sigma12
        end
        A = repmat(A,Repetitions,1);
        M{1} = A{1};
        for m = 2:SuperLayerSize
            M0 = A{m,1}(1,1,:)-M{1}(2,2,:);
            M1 = M{1}(1,2,:)./M0;
            M2 = A{m,1}(2,1,:)./M0;
            M{1} = [M{1}(1,1,:)+M1.*M{1}(2,1,:) -M1.*A{m,1}(1,2,:);M2.*M{1}(2,1,:) A{m,1}(2,2,:)-M2.*A{m,1}(1,2,:)];
        end
        for m = 1:length(Pattern)
            M0 = M{Pattern(m)}(1,1,:)-M{m}(2,2,:);
            M1 = M{m}(1,2,:)./M0;
            M2 = M{Pattern(m)}(2,1,:)./M0;
            M{m+1} = [M{m}(1,1,:)+M1.*M{m}(2,1,:) -M1.*M{Pattern(m)}(1,2,:);M2.*M{m}(2,1,:) M{Pattern(m)}(2,2,:)-M2.*M{Pattern(m)}(1,2,:)];
        end
        if  SymmetricSystem
            A = [A;flipud(A)];
            M1 = -M{end}(1,2,:)./(2*M{end}(2,2,:));
            M{end} = [M{end}(1,1,:)+M1.*M{end}(2,1,:) M1.*M{end}(2,1,:);-M1.*M{end}(2,1,:) -M{end}(1,1,:)-M1.*M{end}(2,1,:)];
        end
        Z1 = [M{end}(1,1,:) -M{end}(1,2,:);M{end}(2,:,:)]; % changed sign in Z1(1,1)!
        Z2 = [M{end}(1,1,:);-M{end}(2,1,:)];
        RT = pagemldivide(Z1,Z2(:,1,:));
        uInterfaces = ones(1,1,Length)+RT(1,1,:);
        uInterfaces(:,Layers+1,:) = RT(2,1,:);
        A{1,2} = A{1};
        for m = 2:Layers
            M0 = A{m,1}(1,1,:)-A{m-1,2}(2,2,:);
            M1 = A{m-1,2}(1,2,:)./M0;
            M2 = A{m,1}(2,1,:)./M0;
            A{m,2} = [A{m-1,2}(1,1,:)+M1.*A{m-1,2}(2,1,:) -M1.*A{m,1}(1,2,:);M2.*A{m-1,2}(2,1,:) A{m,1}(2,2,:)-M2.*A{m,1}(1,2,:)];
        end
        for m = Layers:-1:2
            M0 = A{m,1}(1,1,:)-A{m-1,2}(2,2,:);
            uInterfaces(:,m,:) = M0.\A{m-1,2}(2,1,:).*uInterfaces(:,1,:)-M0.\A{m,1}(1,2,:).*uInterfaces(:,m+1,:);
        end
        for m = 1:Layers
            x3 = (0:A{m,4}/SamplesX3:A{m,4})';
            if  m == 1
                x3Total = x3;
            else
                x3Total = [x3Total;x3Total(end)+x3];
            end
            r = (m-1)*length(x3)+1:m*length(x3);
            E = exp(A{m,6}.*[x3 A{m,4}-x3]);
            U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]);
            v(r,1,:) = -1i*AngularFrequency.*pagemtimes(E,U); % v2
            epsilon(r,1,:) = pagemtimes([A{m,6} -A{m,6}].*E,U); % epsilon23
            epsilon(r,2,:) = 1i*Wavenumber.*pagemtimes(E,U); % epsilon12
            sigma(r,1,:) = pagemtimes([A{m,7} -A{m,7}].*E,U); % sigma23
            sigma(r,2,:) = pagemtimes([A{m,8} A{m,8}].*E,U); % sigma12
            KineticEnergyDensity(r,1,:) = A{m,5}*abs(v(r,1,:)).^2/2;
        end
        PowerFlowDensity = -real(sigma(:,2,:).*conj(v(:,1,:)))/2;
    end
    StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
    PowerFlow = trapz(x3Total,reshape(PowerFlowDensity(:,1,:),Height,Length));
    if  ~Decoupled
        PowerFlow(2,:) = trapz(x3Total,reshape(PowerFlowDensity(:,2,:),Height,Length));
    end
    TotalEnergy = trapz(x3Total,reshape(StrainEnergyDensity,Height,Length)+reshape(KineticEnergyDensity,Height,Length))/2;
    X{p}(:,5:4+size(PowerFlowDensity,2)) = fillmissing(filloutliers((PowerFlow./TotalEnergy)'/1e3,'spline','movmedian',5,'ThresholdFactor',1),'spline'); % ce1,ce2 (m/ms)
    Counter = Counter+1;
    waitbar(Counter/ModeTotal,h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
    if  Stop
        return
    end
end