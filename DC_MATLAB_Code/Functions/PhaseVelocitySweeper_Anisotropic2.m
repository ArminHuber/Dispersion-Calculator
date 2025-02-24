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
function HBScholte = PhaseVelocitySweeper_Anisotropic2(c,Delta,Material,Viscoelastic,UpperFluid,LowerFluid,Frequency,I,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled)
Resolution = 1e-6; % (m/s)
VelocityStep = 2; % (m/s); for coarse sweep
VelocityStepFine = .5; % (m/s); for first iteration in the fine sweeps
EnlargementReal = 3; % SweepRangeReal = XRough(p,1)+EnlargementReal*VelocityStep*[-1 1];
EnlargementImag = 2; % SweepRangeImag = XRough(p,2)+EnlargementImag*VelocityStep*[-1 1];
ImagRealThreshold1 = .1; % exclude wavenumbers with ratios of imaginary part/real part greater than this threshold
ImagRealThreshold2 = .5; % for above slow fluid velocity
DensityThreshold = .01; % rho_fluid/rho_solid
% disp(1+2*VelocityStep/VelocityStepFine*[EnlargementReal EnlargementImag]) % array size of first iteration in the fine sweeps

%#ok<*MINV>
%#ok<*AGROW>
h = msgbox(['Searching Scholte modes at ',num2str(Frequency),' kHz...']);
if  UpperFluid.Velocity > LowerFluid.Velocity
    FastFluidVelocity = UpperFluid.Velocity;
    SlowFluidVelocity = LowerFluid.Velocity;
    Density_ = LowerFluid.Density/Material{end}.Density;
else
    FastFluidVelocity = LowerFluid.Velocity;
    SlowFluidVelocity = UpperFluid.Velocity;
    Density_ = UpperFluid.Density/Material{1}.Density;
end
if  Density_ > DensityThreshold
    RangeImag = .2;
else
    RangeImag = .04;
end
H = cell(1,2);
AngularFrequency = 2*pi*Frequency*1e3;
AngularFrequency2 = AngularFrequency^2;
kUpperFluid2 = AngularFrequency2/UpperFluid.Velocity^2;
kLowerFluid2 = AngularFrequency2/LowerFluid.Velocity^2;
for r = 1:2
    XRough = [];
    for m = 1:SuperLayerSize
        if  ~Decoupled
            rw2(m) = Material{m}.Density*AngularFrequency2;
            r2w4(m) = rw2(m)^2;
            r3w6(m) = rw2(m)^3;
            a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
            a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
            a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
            a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
            a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
            a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
            a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
            a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
            a34(m) = -1/Delta(m);
            b12(m) = a12(m)*rw2(m);
            b22(m) = a22(m)*rw2(m);
            b23(m) = a23(m)*r2w4(m);
            b32(m) = a32(m)*rw2(m);
            b33(m) = a33(m)*r2w4(m);
            b34(m) = a34(m)*r3w6(m);
        else
            rw2(m) = Material{m}.Density*AngularFrequency2;
            r2w4(m) = rw2(m)^2;
            A1(m) = 2*c{m}(3,3)*c{m}(5,5);
            a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
            a22(m) = -c{m}(3,3)-c{m}(5,5);
            a31(m) = c{m}(1,1)*c{m}(5,5);
            a32(m) = -c{m}(1,1)-c{m}(5,5);
            b22(m) = a22(m)*rw2(m);
            b32(m) = a32(m)*rw2(m);
            b33(m) = r2w4(m);
        end
    end
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
    if  ~Decoupled
        Y = Computer_Coupled(r,AngularFrequency2,Wavenumber,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34,I,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,LowerFluid,UpperFluid,kLowerFluid2,kUpperFluid2);
    else
        Y = Computer_Decoupled(r,AngularFrequency2,Wavenumber,c,A1,rw2,a21,a31,b22,b32,b33,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,LowerFluid,UpperFluid,kLowerFluid2,kUpperFluid2);
    end
    Y(isinf(Y)) = NaN;
    if  ~Viscoelastic && r == 1
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
        if  r == 1
            XRough(end+1,:) = [SlowFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        elseif r == 2
            XRough(end+1,:) = [FastFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        end
    else
        if  r == 1
            XRough = [SlowFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        elseif r == 2
            XRough = [FastFluidVelocity-EnlargementReal*VelocityStep -1e-10];
        end
    end
% disp(XRough)
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
            if  ~Decoupled
                Y = Computer_Coupled(r,AngularFrequency2,Wavenumber,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34,I,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,LowerFluid,UpperFluid,kLowerFluid2,kUpperFluid2);
            else
                Y = Computer_Decoupled(r,AngularFrequency2,Wavenumber,c,A1,rw2,a21,a31,b22,b32,b33,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,LowerFluid,UpperFluid,kLowerFluid2,kUpperFluid2);
            end
            Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
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
        end
        p = p+1;
    end
    if  ~isempty(H{r})
        H{r} = sortrows(H{r},1);
    end
end
HBScholte = [H{1};H{2}];
if  isempty(HBScholte)
    disp('No higher order Scholte modes found!')
    close(h)
    return
end
String = ['cp @ ',num2str(Frequency),' kHz:',newline,'Mode        cp(m/ms)  alpha(Np/m)'];
for p = 1:height(HBScholte)
    if  p < 11
        String = append(String,newline,'BScholte',num2str(p-1),'    ',num2str(HBScholte(p,1)/1e3,'%.5f'),'   ',num2str(HBScholte(p,2)));
    else
        String = append(String,newline,'BScholte',num2str(p-1),'   ',num2str(HBScholte(p,1)/1e3,'%.5f'),'   ',num2str(HBScholte(p,2)));
    end
end
disp(append(String,newline,newline,'BScholte: ',num2str(height(HBScholte)),newline,'----------------'));
close(h)
end
function Y = Computer_Coupled(r,AngularFrequency2,Wavenumber,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34,I,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,LowerFluid,UpperFluid,kLowerFluid2,kUpperFluid2)
    Wavenumber2 = Wavenumber.^2;
    Wavenumber4 = Wavenumber2.^2;
    Wavenumber6 = Wavenumber2.^3;
    Y = NaN(size(Wavenumber));
    n1 = 1:height(Wavenumber);
    n2 = 1:width(Wavenumber);
    if  SuperLayerSize > 1
        for m = 1:SuperLayerSize 
            A1 = a11(m)*Wavenumber2+b12(m);
            A2 = a21(m)*Wavenumber4+b22(m)*Wavenumber2+b23(m);
            A3 = a31(m)*Wavenumber6+b32(m)*Wavenumber4+b33(m)*Wavenumber2+b34(m);
            d1 = A2/3-A1.^2/9;
            d2 = A1.^3/27-A1.*A2/6+A3/2;
            d3 = (sqrt(d2.^2+d1.^3)-d2).^(1/3);
            d4 = d1./(2*d3)-d3/2;
            d5 = d1./d3;
            d6 = (sqrt(3)*(d3+d5)*1i)/2;
            k31 = sqrt(d4-d6-A1/3);
            k32 = sqrt(d4+d6-A1/3);
            k33 = -sqrt(d3-d5-A1/3);
            k312 = k31.^2;
            k322 = k32.^2;
            k332 = k33.^2;
            V1(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31));
            V2(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32));
            V3(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33));
            W1(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k312).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31)));
            W2(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k322).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32)));
            W3(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k332).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33)));
            D31(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V1(n1,n2,m)+c{m}(3,3)*k31.*W1(n1,n2,m)); % sigma33
            D32(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V2(n1,n2,m)+c{m}(3,3)*k32.*W2(n1,n2,m)); % sigma33
            D33(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V3(n1,n2,m)+c{m}(3,3)*k33.*W3(n1,n2,m)); % sigma33
            D41(n1,n2,m) = 1i*(c{m}(4,5)*(k31+Wavenumber.*W1(n1,n2,m))+c{m}(4,4)*k31.*V1(n1,n2,m)); % sigma23
            D42(n1,n2,m) = 1i*(c{m}(4,5)*(k32+Wavenumber.*W2(n1,n2,m))+c{m}(4,4)*k32.*V2(n1,n2,m)); % sigma23
            D43(n1,n2,m) = 1i*(c{m}(4,5)*(k33+Wavenumber.*W3(n1,n2,m))+c{m}(4,4)*k33.*V3(n1,n2,m)); % sigma23
            D51(n1,n2,m) = 1i*(c{m}(5,5)*(k31+Wavenumber.*W1(n1,n2,m))+c{m}(4,5)*k31.*V1(n1,n2,m)); % sigma13
            D52(n1,n2,m) = 1i*(c{m}(5,5)*(k32+Wavenumber.*W2(n1,n2,m))+c{m}(4,5)*k32.*V2(n1,n2,m)); % sigma13
            D53(n1,n2,m) = 1i*(c{m}(5,5)*(k33+Wavenumber.*W3(n1,n2,m))+c{m}(4,5)*k33.*V3(n1,n2,m)); % sigma13
            E1(n1,n2,m) = exp(1i*k31*LayerThicknesses(m));
            E2(n1,n2,m) = exp(1i*k32*LayerThicknesses(m));
            E3(n1,n2,m) = exp(1i*k33*LayerThicknesses(m));
        end
    end
    for l = n2
        for j = n1
            if  ~isnan(Wavenumber(j,l))
                if  SuperLayerSize == 1
                    A1 = a11*Wavenumber2(j,l)+b12;
                    A2 = a21*Wavenumber4(j,l)+b22*Wavenumber2(j,l)+b23;
                    A3 = a31*Wavenumber6(j,l)+b32*Wavenumber4(j,l)+b33*Wavenumber2(j,l)+b34;
                    d1 = A2/3-A1^2/9;
                    d2 = A1^3/27-A1*A2/6+A3/2;
                    d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                    d4 = d1/(2*d3)-d3/2;
                    d5 = d1/d3;
                    d6 = (sqrt(3)*(d3+d5)*1i)/2;
                    k3(1) = sqrt(d4-d6-A1/3);
                    k3(2) = sqrt(d4+d6-A1/3);
                    k3(3) = -sqrt(d3-d5-A1/3);
                    k32 = k3.^2;
                    m11 = c{1}(1,1)*Wavenumber2(j,l)-rw2+c{1}(5,5)*k32;
                    m12 = c{1}(1,6)*Wavenumber2(j,l)+c{1}(4,5)*k32;
                    m13 = (c{1}(1,3)+c{1}(5,5))*Wavenumber(j,l)*k3;
                    m22 = c{1}(6,6)*Wavenumber2(j,l)-rw2+c{1}(4,4)*k32;
                    m23 = (c{1}(3,6)+c{1}(4,5))*Wavenumber(j,l)*k3;
                    V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                    W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                    D3 = 1i*(c{1}(1,3)*Wavenumber(j,l)+c{1}(3,6)*Wavenumber(j,l)*V+c{1}(3,3)*k3.*W); % sigma33
                    D4 = 1i*(c{1}(4,5)*(k3+Wavenumber(j,l)*W)+c{1}(4,4)*k3.*V); % sigma23
                    D5 = 1i*(c{1}(5,5)*(k3+Wavenumber(j,l)*W)+c{1}(4,5)*k3.*V); % sigma13
                    E = exp(1i*k3*LayerThicknesses);
%                     L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
%                     L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                    L1(1,1)=D3(1);L1(1,2)=D3(2);L1(1,3)=D3(3);L1(1,4)=D3(1)*E(1);L1(1,5)=D3(2)*E(2);L1(1,6)=D3(3)*E(3);L1(2,1)=D5(1);L1(2,2)=D5(2);L1(2,3)=D5(3);L1(2,4)=-D5(1)*E(1);L1(2,5)=-D5(2)*E(2);L1(2,6)=-D5(3)*E(3);L1(3,1)=D4(1);L1(3,2)=D4(2);L1(3,3)=D4(3);L1(3,4)=-D4(1)*E(1);L1(3,5)=-D4(2)*E(2);L1(3,6)=-D4(3)*E(3);L1(4,1)=D3(1)*E(1);L1(4,2)=D3(2)*E(2);L1(4,3)=D3(3)*E(3);L1(4,4)=D3(1);L1(4,5)=D3(2);L1(4,6)=D3(3);L1(5,1)=D5(1)*E(1);L1(5,2)=D5(2)*E(2);L1(5,3)=D5(3)*E(3);L1(5,4)=-D5(1);L1(5,5)=-D5(2);L1(5,6)=-D5(3);L1(6,1)=D4(1)*E(1);L1(6,2)=D4(2)*E(2);L1(6,3)=D4(3)*E(3);L1(6,4)=-D4(1);L1(6,5)=-D4(2);L1(6,6)=-D4(3);
                    L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(1);L2(2,2)=V(2);L2(2,3)=V(3);L2(2,4)=V(1)*E(1);L2(2,5)=V(2)*E(2);L2(2,6)=V(3)*E(3);L2(3,1)=W(1);L2(3,2)=W(2);L2(3,3)=W(3);L2(3,4)=-W(1)*E(1);L2(3,5)=-W(2)*E(2);L2(3,6)=-W(3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(1)*E(1);L2(5,2)=V(2)*E(2);L2(5,3)=V(3)*E(3);L2(5,4)=V(1);L2(5,5)=V(2);L2(5,6)=V(3);L2(6,1)=W(1)*E(1);L2(6,2)=W(2)*E(2);L2(6,3)=W(3)*E(3);L2(6,4)=-W(1);L2(6,5)=-W(2);L2(6,6)=-W(3);
                    L{1} = L1/L2;
                else
                    for m = 1:SuperLayerSize
%                         L1 = [D31(j,l,m) D32(j,l,m) D33(j,l,m) D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m);D51(j,l,m) D52(j,l,m) D53(j,l,m) -D51(j,l,m)*E1(j,l,m) -D52(j,l,m)*E2(j,l,m) -D53(j,l,m)*E3(j,l,m);D41(j,l,m) D42(j,l,m) D43(j,l,m) -D41(j,l,m)*E1(j,l,m) -D42(j,l,m)*E2(j,l,m) -D43(j,l,m)*E3(j,l,m);D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m) D31(j,l,m) D32(j,l,m) D33(j,l,m);D51(j,l,m)*E1(j,l,m) D52(j,l,m)*E2(j,l,m) D53(j,l,m)*E3(j,l,m) -D51(j,l,m) -D52(j,l,m) -D53(j,l,m);D41(j,l,m)*E1(j,l,m) D42(j,l,m)*E2(j,l,m) D43(j,l,m)*E3(j,l,m) -D41(j,l,m) -D42(j,l,m) -D43(j,l,m)];
%                         L2 = [1 1 1 E1(j,l,m) E2(j,l,m) E3(j,l,m);V1(j,l,m) V2(j,l,m) V3(j,l,m) V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m);W1(j,l,m) W2(j,l,m) W3(j,l,m) -W1(j,l,m)*E1(j,l,m) -W2(j,l,m)*E2(j,l,m) -W3(j,l,m)*E3(j,l,m);E1(j,l,m) E2(j,l,m) E3(j,l,m) 1 1 1;V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m) V1(j,l,m) V2(j,l,m) V3(j,l,m);W1(j,l,m)*E1(j,l,m) W2(j,l,m)*E2(j,l,m) W3(j,l,m)*E3(j,l,m) -W1(j,l,m) -W2(j,l,m) -W3(j,l,m)];
                        L1(1,1)=D31(j,l,m);L1(1,2)=D32(j,l,m);L1(1,3)=D33(j,l,m);L1(1,4)=D31(j,l,m)*E1(j,l,m);L1(1,5)=D32(j,l,m)*E2(j,l,m);L1(1,6)=D33(j,l,m)*E3(j,l,m);L1(2,1)=D51(j,l,m);L1(2,2)=D52(j,l,m);L1(2,3)=D53(j,l,m);L1(2,4)=-D51(j,l,m)*E1(j,l,m);L1(2,5)=-D52(j,l,m)*E2(j,l,m);L1(2,6)=-D53(j,l,m)*E3(j,l,m);L1(3,1)=D41(j,l,m);L1(3,2)=D42(j,l,m);L1(3,3)=D43(j,l,m);L1(3,4)=-D41(j,l,m)*E1(j,l,m);L1(3,5)=-D42(j,l,m)*E2(j,l,m);L1(3,6)=-D43(j,l,m)*E3(j,l,m);L1(4,1)=D31(j,l,m)*E1(j,l,m);L1(4,2)=D32(j,l,m)*E2(j,l,m);L1(4,3)=D33(j,l,m)*E3(j,l,m);L1(4,4)=D31(j,l,m);L1(4,5)=D32(j,l,m);L1(4,6)=D33(j,l,m);L1(5,1)=D51(j,l,m)*E1(j,l,m);L1(5,2)=D52(j,l,m)*E2(j,l,m);L1(5,3)=D53(j,l,m)*E3(j,l,m);L1(5,4)=-D51(j,l,m);L1(5,5)=-D52(j,l,m);L1(5,6)=-D53(j,l,m);L1(6,1)=D41(j,l,m)*E1(j,l,m);L1(6,2)=D42(j,l,m)*E2(j,l,m);L1(6,3)=D43(j,l,m)*E3(j,l,m);L1(6,4)=-D41(j,l,m);L1(6,5)=-D42(j,l,m);L1(6,6)=-D43(j,l,m);
                        L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E1(j,l,m);L2(1,5)=E2(j,l,m);L2(1,6)=E3(j,l,m);L2(2,1)=V1(j,l,m);L2(2,2)=V2(j,l,m);L2(2,3)=V3(j,l,m);L2(2,4)=V1(j,l,m)*E1(j,l,m);L2(2,5)=V2(j,l,m)*E2(j,l,m);L2(2,6)=V3(j,l,m)*E3(j,l,m);L2(3,1)=W1(j,l,m);L2(3,2)=W2(j,l,m);L2(3,3)=W3(j,l,m);L2(3,4)=-W1(j,l,m)*E1(j,l,m);L2(3,5)=-W2(j,l,m)*E2(j,l,m);L2(3,6)=-W3(j,l,m)*E3(j,l,m);L2(4,1)=E1(j,l,m);L2(4,2)=E2(j,l,m);L2(4,3)=E3(j,l,m);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V1(j,l,m)*E1(j,l,m);L2(5,2)=V2(j,l,m)*E2(j,l,m);L2(5,3)=V3(j,l,m)*E3(j,l,m);L2(5,4)=V1(j,l,m);L2(5,5)=V2(j,l,m);L2(5,6)=V3(j,l,m);L2(6,1)=W1(j,l,m)*E1(j,l,m);L2(6,2)=W2(j,l,m)*E2(j,l,m);L2(6,3)=W3(j,l,m)*E3(j,l,m);L2(6,4)=-W1(j,l,m);L2(6,5)=-W2(j,l,m);L2(6,6)=-W3(j,l,m);                    
                        L{m} = L1/L2;
                    end
                end
                M = L{1};
                for m = 2:SuperLayerSize
                    N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                    M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)]; 
                end
                MM{1} = M;
                for m = 2:log2(Repetitions)+1
                    N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                    MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                end
                for m = m+1:length(Pattern)
                    N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                    MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                end
                if  SymmetricSystem
                    N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                    MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                end
                G = inv(MM{end});
                if  r == 1
                    k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2(j,l));
                    k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2(j,l));
                else
                    if  UpperFluid.Velocity > LowerFluid.Velocity
                        k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2(j,l));
                        k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2(j,l));
                    else
                        k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2(j,l));
                        k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2(j,l));
                    end
                end
                WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2/Wavenumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2/Wavenumber(j,l); % in the lower fluid
                Y(j,l) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
            end
        end
    end
end
function Y = Computer_Decoupled(r,AngularFrequency2,Wavenumber,c,A1,rw2,a21,a31,b22,b32,b33,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,LowerFluid,UpperFluid,kLowerFluid2,kUpperFluid2)
    Wavenumber2 = Wavenumber.^2;
    Wavenumber4 = Wavenumber2.^2;
    Y = NaN(size(Wavenumber));
    for l = 1:width(Wavenumber)
        for j = 1:height(Wavenumber)
            if  ~isnan(Wavenumber(j,l))
                for m = 1:SuperLayerSize
                    A2 = a21(m)*Wavenumber2(j,l)+b22(m);
                    A3 = a31(m)*Wavenumber4(j,l)+b32(m)*Wavenumber2(j,l)+b33(m);
                    k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                    k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                    W = (rw2(m)-c{m}(1,1)*Wavenumber2(j,l)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber(j,l)*k3);
                    D3 = 1i*(c{m}(1,3)*Wavenumber(j,l)+c{m}(3,3)*k3.*W); % sigma33
                    D5 = 1i*(c{m}(5,5)*(k3+Wavenumber(j,l)*W)); % sigma13
                    E = exp(1i*k3*LayerThicknesses(m));
%                     L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
%                     L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
                    L1(1,1)=D3(1);L1(1,2)=D3(2);L1(1,3)=D3(1)*E(1);L1(1,4)=D3(2)*E(2);L1(2,1)=D5(1);L1(2,2)=D5(2);L1(2,3)=-D5(1)*E(1);L1(2,4)=-D5(2)*E(2);L1(3,1)=D3(1)*E(1);L1(3,2)=D3(2)*E(2);L1(3,3)=D3(1);L1(3,4)=D3(2);L1(4,1)=D5(1)*E(1);L1(4,2)=D5(2)*E(2);L1(4,3)=-D5(1);L1(4,4)=-D5(2);
                    L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(1);L2(2,2)=W(2);L2(2,3)=-W(1)*E(1);L2(2,4)=-W(2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(1)*E(1);L2(4,2)=W(2)*E(2);L2(4,3)=-W(1);L2(4,4)=-W(2);
                    L{m} = L1/L2;
                end
                M = L{1};
                for m = 2:SuperLayerSize
                    N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                    M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                end
                MM{1} = M;
                for m = 2:log2(Repetitions)+1
                    N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                    MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
                end
                for m = m+1:length(Pattern)
                    N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
                    MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
                end
                if  SymmetricSystem
                    N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
                    MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
                end
                G = inv(MM{end});
                if  r == 1
                    k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2(j,l));
                    k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2(j,l));
                else
                    if  UpperFluid.Velocity > LowerFluid.Velocity
                        k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2(j,l));
                        k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2(j,l));
                    else
                        k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2(j,l));
                        k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2(j,l));
                    end
                end
                WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2/Wavenumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2/Wavenumber(j,l); % in the lower fluid
                Y(j,l) = abs((WLowerFluid-G(4,3)*DLowerFluid)*(WUpperFluid+G(2,1)*DUpperFluid)+G(2,3)*G(4,1)*DUpperFluid*DLowerFluid);
            end
        end
    end
end