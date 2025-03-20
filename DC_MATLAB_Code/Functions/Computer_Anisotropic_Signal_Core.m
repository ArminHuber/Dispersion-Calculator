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
function [uSum,Counter,ExcitationSpectrumRange] = Computer_Anisotropic_Signal_Core(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,CoherenceTime,Counter,Distance,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,ModeType,ModeTotal,PhaseVelocity,Attenuation,Time,TimeLimit,z1,z2,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory)
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
uSum{length(PhaseVelocity)} = [];
ExcitationSpectrumRange{length(PhaseVelocity)} = [];
for p = 1:length(PhaseVelocity)
    if  any(PhaseVelocity{p})
        ExcitationSpectrumRange{p} = ExcitationSpectrum(:,z1(p):z2(p)); % extract from the frequency spectrum the range for which we have phase velocities
        if  CoherenceTime < TimeLimit % if the coherence time is smaller than the calculated time range, we have to expect seeing unwanted twins of the wave packet
            Fit1 = fit(ExcitationSpectrumRange{p}(1,:)',ExcitationSpectrumRange{p}(2,:)','cubicspline'); % fit the spectral amplitudes
            Fit2 = fit(ExcitationSpectrumRange{p}(1,:)',PhaseVelocity{p},'cubicspline'); % fit the phase velocity
            Fit3 = fit(ExcitationSpectrumRange{p}(1,:)',Attenuation{p},'cubicspline'); % fit the attenuation
            ExcitationSpectrumRange{p} = ExcitationSpectrumRange{p}(1,1):FrequencyResolution2:ExcitationSpectrumRange{p}(1,end); % generate the new frequency range with smaller steps
            ExcitationSpectrumRange{p}(2,:) = Fit1(ExcitationSpectrumRange{p}(1,:))/FrequencyResolution*FrequencyResolution2; % interpolate the spectral amplitudes
            PhaseVelocity{p} = Fit2(ExcitationSpectrumRange{p}(1,:)); % interpolate the phase velocity
            Attenuation{p} = Fit3(ExcitationSpectrumRange{p}(1,:)); % interpolate the attenuation
        end
        Length = length(ExcitationSpectrumRange{p});
        if  Length*length(Time)*8 > .8*Memory.MaxPossibleArrayBytes % Bytes needed for u = real(exp(1i*(reshape(Wavenumber,[],1)*Distance-reshape(AngularFrequency,[],1).*Time)));
            errordlg('The data size exceeds your physical RAM! Decrease the parameter settings.','Error');
            return
        end
        PhaseVelocity{p} = reshape(PhaseVelocity{p},1,1,[]);
        Attenuation{p} = reshape(Attenuation{p},1,1,[]);
        AngularFrequency = reshape(ExcitationSpectrumRange{p}(1,:),1,1,[])*pi*2e3;
        AngularFrequency2 = AngularFrequency.^2;
        Wavenumber = AngularFrequency./PhaseVelocity{p}.*(1+1i*Attenuation{p}/2/pi);
        Wavenumber2 = Wavenumber.^2;
        A = cell(0);
        if  ~Decoupled
            Wavenumber4 = Wavenumber2.^2;
            Wavenumber6 = Wavenumber2.^3;
            k3 = zeros(1,3,Length);
            v = zeros(Height,3,Length);
            sigma = v;
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
                A{m,5} = k3;
                A{m,6} = V;
                A{m,7} = W;
                A{m,8} = 1i*((c{m}(1,1)+c{m}(1,6)*V).*Wavenumber+c{m}(1,3)*e3); % sigma11
                A{m,9} = D5; % sigma13
                A{m,10} = 1i*((c{m}(1,6)+c{m}(6,6)*V).*Wavenumber+c{m}(3,6)*e3); % sigma12
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
                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2./Wavenumber; % sigma11, sigma22, sigma33 in the upper fluid
                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2./Wavenumber; % in the lower fluid
                ULowerFluid = -(k3UpperFluid./Wavenumber+M{end}(3,1,:).*DUpperFluid)./(M{end}(3,4,:).*DLowerFluid);
                uInterfaces = M{end}(1:3,1,:).*DUpperFluid+M{end}(1:3,4,:).*DLowerFluid.*ULowerFluid;
                uInterfaces(:,Layers+1,:) = M{end}(4:6,1,:).*DUpperFluid+M{end}(4:6,4,:).*DLowerFluid.*ULowerFluid;
            else
                Z1 = [pagemtimes(-M{end}(1:3,1:3,:),[ones(1,3,Length);A{1,6};-A{1,7}]) pagemtimes(-M{end}(1:3,4:6,:),[ones(1,3,Length);A{end,6};A{end,7}]);pagemtimes(M{end}(4:6,1:3,:),[ones(1,3,Length);A{1,6};-A{1,7}]) pagemtimes(M{end}(4:6,4:6,:),[ones(1,3,Length);A{end,6};A{end,7}])];
                Z2 = [pagemtimes(M{end}(1:3,1:3,:),[ones(1,3,Length);A{1,6};A{1,7}]);pagemtimes(-M{end}(4:6,1:3,:),[ones(1,3,Length);A{1,6};A{1,7}])];
                RT = pagemldivide(Z1,Z2(:,1,:));
                uInterfaces = [ones(1,1,Length);A{1,6}(1,1,:);A{1,7}(1,1,:)]+pagemtimes([ones(1,3,Length);A{1,6};-A{1,7}],RT(1:3,1,:));
                uInterfaces(:,Layers+1,:) = pagemtimes([ones(1,3,Length);A{end,6};A{end,7}],RT(4:6,1,:));
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
                U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]);
                if  m == 1
                    U0 = U;
                    x30 = x3;
                    x3Total = x3;
                else
                    x3Total = [x3Total;x3Total(end)+x3];
                end
                r = (m-1)*length(x3)+1:m*length(x3);
                E = exp(1i*[A{m,5}.*x3 A{m,5}.*(A{m,4}-x3)]);
                v(r,1,:) = -1i*AngularFrequency.*pagemtimes(E,U); % v1
                v(r,2,:) = -1i*AngularFrequency.*pagemtimes([A{m,6} A{m,6}].*E,U); % v2
                v(r,3,:) = -1i*AngularFrequency.*pagemtimes([A{m,7} -A{m,7}].*E,U); % v3
                sigma(r,1,:) = pagemtimes([A{m,8} A{m,8}].*E,U); % sigma11
                sigma(r,3,:) = pagemtimes([A{m,9} -A{m,9}].*E,U); % sigma13
                sigma(r,2,:) = pagemtimes([A{m,10} A{m,10}].*E,U); % sigma12
            end
        elseif Decoupled && ~strcmp(ModeType,'Shear')
            Wavenumber4 = Wavenumber2.^2;
            k3 = zeros(1,2,Length);
            v = zeros(Height,2,Length);
            sigma = v;
            for m = 1:SuperLayerSize
                rw2 = Material{m}.Density*AngularFrequency2;
                A2 = a21(m)*Wavenumber2+a22(m)*rw2;
                A3 = a31(m)*Wavenumber4+a32(m)*rw2.*Wavenumber2+rw2.^2;
                d1 = sqrt(A2.^2-2*A1(m)*A3);
                k3(1,1,:) = sqrt((-A2+d1)/A1(m));
                k3(1,2,:) = sqrt((-A2-d1)/A1(m));
                W = (rw2-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k3);
                e3 = k3.*W;
                D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*e3);
                D5 = 1i*c{m}(5,5)*(Wavenumber.*W+k3);
                E = exp(1i*k3*LayerThicknesses(m));
                L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                A{m,3} = [ones(1,2,Length) E;W -W.*E;E ones(1,2,Length);W.*E -W]; % L2
                A{m,1} = pagemrdivide(L1,A{m,3}); % L
                A{m,4} = LayerThicknesses(m);
                A{m,5} = k3;
                A{m,7} = W;
                A{m,6} = 1i*(c{m}(1,1)*Wavenumber+c{m}(1,3)*e3); % sigma11
                A{m,8} = D5; % sigma13
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
                U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]);
                if  m == 1
                    U0 = U;
                    x30 = x3;
                    x3Total = x3;
                else
                    x3Total = [x3Total;x3Total(end)+x3];
                end
                r = (m-1)*length(x3)+1:m*length(x3);
                E = exp(1i*[A{m,5}.*x3 A{m,5}.*(A{m,4}-x3)]);
                v(r,1,:) = -1i*AngularFrequency.*pagemtimes(E,U); % v1
                v(r,2,:) = -1i*AngularFrequency.*pagemtimes([A{m,7} -A{m,7}].*E,U); % v3
                sigma(r,1,:) = pagemtimes([A{m,6} A{m,6}].*E,U); % sigma11
                sigma(r,2,:) = pagemtimes([A{m,8} -A{m,8}].*E,U); % sigma13
            end
        elseif Decoupled && strcmp(ModeType,'Shear') && DisplacementComponent == 2
            k3 = zeros(1,1,Length);
            v = zeros(Height,1,Length);
            sigma = v;
            for m = 1:SuperLayerSize
                k3(1,1,:) = sqrt((Material{m}.Density*AngularFrequency2-Wavenumber2*c{m}(6,6))/c{m}(4,4));
                if  k3(1) == 0
                    k3(1,1,:) = 1e-10;
                end
                E = exp(1i*k3*LayerThicknesses(m));
                E2 = E.^2;
                A{m,1} = 1i*k3*c{m}(4,4)./(E2-ones(1,1,Length)).*[-ones(1,1,Length)-E2 2*E;-2*E ones(1,1,Length)+E2]; % L
                A{m,3} = [ones(1,1,Length) E;E ones(1,1,Length)]; % L2
                A{m,4} = LayerThicknesses(m);
                A{m,5} = k3;
                A{m,6} = 1i*Wavenumber*c{m}(6,6); % sigma12
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
                U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]); 
                if  m == 1
                    U0 = U;
                    x30 = x3;
                    x3Total = x3;
                else
                    x3Total = [x3Total;x3Total(end)+x3];
                end
                r = (m-1)*length(x3)+1:m*length(x3);
                E = exp(1i*A{m,5}.*[x3 (A{m,4}-x3)]);
                v(r,1,:) = -1i*AngularFrequency.*pagemtimes(E,U); % v2
                sigma(r,1,:) = pagemtimes([A{m,6} A{m,6}].*E,U); % sigma12
            end
        end
        if  strcmp(ModeType,'Shear') && DisplacementComponent == 1
            uSum{p} = zeros(1,length(Time));
        else
            PowerFlow = trapz(x3Total,-real(sum(sigma.*conj(v),2))/2);
            E = exp(1i*(Wavenumber*Distance+[A{1,5}.*x30(1:2) A{1,5}.*(A{1,4}-x30(1:2))]));
            if  DisplacementComponent == 1 % u3
                unorm = pagemtimes([A{1,7} -A{1,7}].*E,U0)./sqrt(PowerFlow);
            elseif DisplacementComponent == 2 % u1,u2
                unorm = pagemtimes(E,U0)./sqrt(PowerFlow);
            end
            unorm = abs(real(unorm(1,1,:).*exp(-1i*angle(unorm(2,1,:)))));
            uNorm = fillmissing(filloutliers(unorm,'spline','movmedian',5,'ThresholdFactor',1),'spline');
            u = real(exp(1i*(reshape(Wavenumber,[],1)*Distance-reshape(AngularFrequency,[],1).*Time)));
            u = u./max(u,[],2).*reshape(uNorm,[],1).*ExcitationSpectrumRange{p}(2,:)';
            uSum{p} = sum(u);
% figure
% hold on
% plot(ExcitationSpectrumRange{p}(1,:),reshape(unorm,1,Length),'linewidth',4,'color','r')
% plot(ExcitationSpectrumRange{p}(1,:),reshape(uNorm,1,Length),'linewidth',1.5,'color','g')
% figure
% hold on
% for i = 1:length(PhaseVelocity{p})
% plot(u(i,:))
% end
        end
        Counter = Counter+1;
        if  ModeTotal > 0
            waitbar(Counter/ModeTotal,h1,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
        end
    end
    if  Stop
        return
    end
end