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
function [MatrixMethodLimit,EvanescenceLimit,XS0,XA0,XSH0] = MatrixMethodLimitComputer(c,Delta,Material,Hybrid,Pattern,SuperLayerSize,LayerThicknesses,Symmetric,SymmetricSystem,Decoupled,FrequencyRange)
Resolution = 1e-6; % (m/s)
InstabilityLimitFactor = 4; % multiples of the phase velocity limit at which TMM becomes unstable; below this multiple, SMM is used instead of TMM

%#ok<*AGROW>
XSH0 = [];
for m = 1:SuperLayerSize
    c{m} = real(c{m});
    if  ~Decoupled
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34(m) = -1/Delta(m);
        A1=0;
    else
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22(m) = -c{m}(3,3)-c{m}(5,5);
        a31(m) = c{m}(1,1)*c{m}(5,5);
        a32(m) = -c{m}(1,1)-c{m}(5,5);
        a11=0;a12=0;a23=0;a33=0;a34=0;
    end
end

% find S0/B1 and S1/B2 at low frequency
AngularFrequency = 2*pi*FrequencyRange(1)*1e3;
SweepRange = 200:10:25e3;
Wavenumber = AngularFrequency./SweepRange;
if  ~Decoupled
    [Y,EvanescenceLimit] = Computer_Coupled('S',1,SweepRange,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,a11,a21,a31,a12,a22,a23,a32,a33,a34);
else
    [Y,EvanescenceLimit] = Computer_Decoupled('S',1,SweepRange,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,A1,a21,a31,a22,a32);
end
XRough = [];
for i = 2:length(SweepRange)-1 % remove strange spikes crossing through zero and causing wrong solutions
    if  sign(Y(i)) ~= sign(Y(i-1)) && sign(Y(i)) ~= sign(Y(i+1))
        Y(i) = NaN;
    end
end
for i = 2:length(SweepRange)-1
    if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1)) && sign(Y(i-1)) ~= sign(Y(i+1))
        XRough(end+1) = SweepRange(i);
        if  (~Decoupled && length(XRough) == 2) || (Decoupled && isscalar(XRough))
            break
        end
    end
end
% figure,plot(SweepRange,20*log10(abs(Y)))
% figure,plot(SweepRange,real(Y)),yline(0),ax = gca;ax.YLim = 1e19*[-1 1];
XS0 = Converger('S',0,0,Decoupled,XRough,AngularFrequency,SweepRange,Resolution,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,A1,a11,a21,a31,a12,a22,a23,a32,a33,a34);

% find A0/B0 at low frequency
SweepRange = .1:200;
Wavenumber = AngularFrequency./SweepRange;
if  ~Decoupled
    [Y,~] = Computer_Coupled('A',0,SweepRange,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,a11,a21,a31,a12,a22,a23,a32,a33,a34);
else
    [Y,~] = Computer_Decoupled('A',0,SweepRange,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,A1,a21,a31,a22,a32);
end
for i = 2:length(SweepRange)-1
    if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1)) && sign(Y(i-1)) ~= sign(Y(i+1))
        XRough = SweepRange(i);
        break
    end
end
XA0 = Converger('A',0,0,Decoupled,XRough,AngularFrequency,SweepRange,Resolution,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,A1,a11,a21,a31,a12,a22,a23,a32,a33,a34);

% find SSH0/BSH0 at low frequency
if  Decoupled
    if  Hybrid
        SweepRange = XA0:10:XS0;
        Wavenumber = AngularFrequency./SweepRange;
        Y = Computer_Decoupled_SH(SweepRange,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,c);
        for i = 2:length(SweepRange)-1
            if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1)) && sign(Y(i-1)) ~= sign(Y(i+1))
                XRough = SweepRange(i);
                break
            end
        end
        XSH0 = Converger('',0,1,Decoupled,XRough,AngularFrequency,SweepRange,Resolution,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,A1,a11,a21,a31,a12,a22,a23,a32,a33,a34);
    else
        XSH0 = sqrt(c{1}(6,6)/Material{1}.Density);
    end
end

% find instability limit at high frequency
AngularFrequency = 2*pi*FrequencyRange(end)*1e3;
SweepRange = 0:XS0(end);
Wavenumber = AngularFrequency./SweepRange;
if  ~Decoupled
    [Y,~] = Computer_Coupled('S',0,SweepRange,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,a11,a21,a31,a12,a22,a23,a32,a33,a34);
else
    [Y,~] = Computer_Decoupled('S',0,SweepRange,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,A1,a21,a31,a22,a32);
end
% figure,plot(SweepRange,20*log10(abs(Y)))
XRough = [];
for i = length(SweepRange)-1:-1:2
    if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1))
        XRough(end+1) = SweepRange(i);
    end
end
Diff = abs(diff(XRough));
InstabilityLimit = XS0(end)/4;
for i = 3:length(Diff)
    if  all(Diff(i-2:i) <= 10*(SweepRange(2)-SweepRange(1)))
        InstabilityLimit = XRough(i-2);
        break
    end
end
MatrixMethodLimit = InstabilityLimitFactor*InstabilityLimit*FrequencyRange/FrequencyRange(end);
MatrixMethodLimit(MatrixMethodLimit > EvanescenceLimit) = EvanescenceLimit;
end
function X = Converger(ModeFamily,GetLimit,SHMode,Decoupled,XRough,AngularFrequency,SweepRange,Resolution,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,A1,a11,a21,a31,a12,a22,a23,a32,a33,a34)
    X = [];
    Bisections = ceil(log2(abs(SweepRange(1)-SweepRange(2))/Resolution));
    for j = 1:length(XRough)
        PhaseVelocity = XRough(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            PhaseVelocity = PhaseVelocity(1):(PhaseVelocity(end)-PhaseVelocity(1))/4:PhaseVelocity(end);
            Wavenumber = AngularFrequency./PhaseVelocity;
            if  ~Decoupled
                [Y,~] = Computer_Coupled(ModeFamily,GetLimit,PhaseVelocity,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,a11,a21,a31,a12,a22,a23,a32,a33,a34);
            else
                if  ~SHMode
                    [Y,~] = Computer_Decoupled(ModeFamily,GetLimit,PhaseVelocity,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,A1,a21,a31,a22,a32);
                else
                    Y = Computer_Decoupled_SH(PhaseVelocity,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,c);
                end
            end
            for i = 2:length(PhaseVelocity)-1
                if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1))
                    if  o == Bisections
                        X(end+1) = PhaseVelocity(i);
                    end
                    PhaseVelocity = [PhaseVelocity(i-1) PhaseVelocity(i+1)];
                    break
                end
            end
        end
    end
end
function [Y,EvanescenceLimit] = Computer_Coupled(ModeFamily,GetLimit,PhaseVelocity,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,a11,a21,a31,a12,a22,a23,a32,a33,a34)
    Size = size(Wavenumber);
    PhaseVelocity2 = reshape(PhaseVelocity.^2,1,1,[]);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Length = length(Wavenumber);
    Y = NaN(size(Wavenumber)); % to be discarded once pagedet exists
    for m = 1:SuperLayerSize
        rc2 = Material{m}.Density*PhaseVelocity2;
        r2c4 = rc2.^2;
        A1 = a11(m)+a12(m)*rc2;
        A2 = a21(m)+a22(m)*rc2+a23(m)*r2c4;
        A3 = a31(m)+a32(m)*rc2+a33(m)*r2c4+a34(m)*rc2.^3;
        d1 = A1/3;
        d2 = A2/3-d1.^2;
        d3 = d1.^3-d1.*A2/2+A3/2;
        d4 = (sqrt(d2.^3+d3.^2)-d3).^(1/3);
        d5 = d2./d4;
        d6 = (d5-d4)/2-d1;
        d7 = (d5+d4)/2i*sqrt(3);
        Alpha(1,1,:) = sqrt(d6+d7);
        Alpha(1,2,:) = sqrt(d6-d7);
        Alpha(1,3,:) = sqrt(d4-d5-d1);
        Alpha2 = Alpha.^2;
        m11 = c{m}(1,1)+c{m}(5,5)*Alpha2-rc2;
        m22 = c{m}(6,6)+c{m}(4,4)*Alpha2-rc2;
        m12 = c{m}(1,6)+c{m}(4,5)*Alpha2;
        m13 = (c{m}(1,3)+c{m}(5,5))*Alpha;
        m23 = (c{m}(3,6)+c{m}(4,5))*Alpha;
        m1 = m13.*m22-m12.*m23;
        V = (m11.*m23-m13.*m12)./m1;
        W = (m11.*m22-m12.^2)./-m1;
        e1 = Alpha+W;
        e2 = Alpha.*V;
        D3 = c{m}(1,3)+c{m}(3,6)*V+c{m}(3,3)*Alpha.*W;
        D4 = c{m}(4,5)*e1+c{m}(4,4)*e2;
        D5 = c{m}(5,5)*e1+c{m}(4,5)*e2;
        if  Symmetric && SuperLayerSize == 1
            Phi = .5i*Wavenumber.*Alpha*LayerThicknesses;
        else
            Phi = 1i*Wavenumber.*Alpha*LayerThicknesses(m);
        end
        E = exp(Phi);
        E_ = exp(-Phi);
        L1 = [E E_;V.*E V.*E_;W.*E -W.*E_;D3.*E D3.*E_;D5.*E -D5.*E_;D4.*E -D4.*E_];
        L2 = [ones(1,6,Length);V V;W -W;D3 D3;D5 -D5;D4 -D4];
        L{m} = pagemrdivide(L1,L2);
        if  GetLimit
            Propagating(3*m-2:3*m,:) = abs(imag(Alpha)) < 1e-10;
        end
    end
    M{1} = L{1};
    for m = 2:SuperLayerSize
        M{1} = pagemtimes(M{1},L{m});
    end
    for m = 1:length(Pattern)
        M{m+1} = pagemtimes(M{m},M{Pattern(m)});
    end
    if  Symmetric
        if  strcmp(ModeFamily,'S')
            for j = 1:length(Wavenumber)
                Y(j) = real(det(M{end}(4:6,[1:2 4],j)));
            end
        elseif strcmp(ModeFamily,'A')
            for j = 1:length(Wavenumber)
                Y(j) = imag(det(M{end}(4:6,[3 5:6],j)));
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
        for j = 1:length(Wavenumber)
            Y(j) = imag(det(M{end}(4:6,1:3,j)));
        end
    end
    Y = reshape(Y,Size);
    if  GetLimit
        for i = 1:width(Propagating)
            if  all(Propagating(:,i))
                EvanescenceLimit = PhaseVelocity(i);
                break
            end
        end
    else
        EvanescenceLimit = 0;
    end
end
function [Y,EvanescenceLimit] = Computer_Decoupled(ModeFamily,GetLimit,PhaseVelocity,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,c,A1,a21,a31,a22,a32)
    Size = size(Wavenumber);
    PhaseVelocity2 = reshape(PhaseVelocity.^2,1,1,[]);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Length = length(Wavenumber);
    Y = NaN(size(Wavenumber)); % to be discarded once pagedet exists
    for m = 1:SuperLayerSize
        rc2 = Material{m}.Density*PhaseVelocity2;
        A2 = a21(m)+a22(m)*rc2;
        A3 = a31(m)+a32(m)*rc2+rc2.^2;
        d1 = sqrt(A2.^2-2*A1(m)*A3);
        Alpha(1,1,:) = sqrt((-A2+d1)/A1(m));
        Alpha(1,2,:) = sqrt((-A2-d1)/A1(m));
        W = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha.^2)./((c{m}(1,3)+c{m}(5,5))*Alpha);
        D3 = c{m}(1,3)+c{m}(3,3)*Alpha.*W;
        D5 = c{m}(5,5)*(Alpha+W);
        if  Symmetric && SuperLayerSize == 1
            Phi = .5i*Wavenumber.*Alpha*LayerThicknesses;
        else
            Phi = 1i*Wavenumber.*Alpha*LayerThicknesses(m);
        end
        E = exp(Phi);
        E_ = exp(-Phi);
        L1 = [E E_;W.*E -W.*E_;D3.*E D3.*E_;D5.*E -D5.*E_];
        L2 = [ones(1,4,Length);W -W;D3 D3;D5 -D5];
        L{m} = pagemrdivide(L1,L2);
        if  GetLimit
            Propagating(2*m-1:2*m,:) = abs(imag(Alpha)) < 1e-10;
        end
    end
    M{1} = L{1};
    for m = 2:SuperLayerSize
        M{1} = pagemtimes(M{1},L{m});
    end
    for m = 1:length(Pattern)
        M{m+1} = pagemtimes(M{m},M{Pattern(m)});
    end
    if  Symmetric
        if  strcmp(ModeFamily,'S')
            for j = 1:length(Wavenumber)
                Y(j) = imag(det(M{end}(3:4,[1 3],j)));
            end
        elseif strcmp(ModeFamily,'A')
            for j = 1:length(Wavenumber)
                Y(j) = imag(det(M{end}(3:4,[2 4],j)));
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
        for j = 1:length(Wavenumber)
            Y(j) = real(det(M{end}(3:4,1:2,j)));
        end
    end
    Y = reshape(Y,Size);
    if  GetLimit
        for i = 1:width(Propagating)
            if  all(Propagating(:,i))
                EvanescenceLimit = PhaseVelocity(i);
                break
            end
        end
    else
        EvanescenceLimit = 0;
    end
end
function Y = Computer_Decoupled_SH(PhaseVelocity,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,c)
    Size = size(Wavenumber);
    PhaseVelocity2 = reshape(PhaseVelocity.^2,1,1,[]);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    for m = 1:SuperLayerSize
        Alpha = sqrt((Material{m}.Density*PhaseVelocity2-c{m}(6,6))/c{m}(4,4));
        D = Alpha*c{m}(4,4);
        G = Wavenumber.*Alpha*LayerThicknesses(m);
        CosG = cos(G);
        SinG = sin(G);
        L{m} = [CosG 1i*SinG./D;1i*SinG.*D CosG];
    end
    M{1} = L{1};
    for m = 2:SuperLayerSize
        M{1} = pagemtimes(M{1},L{m});
    end
    for m = 1:length(Pattern)
        M{m+1} = pagemtimes(M{m},M{Pattern(m)});
    end
    Y = reshape(M{end}(2,1,:),Size);
end