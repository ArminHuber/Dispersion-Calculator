% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2026 DLR
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
function [HSLamb,HALamb,HBLamb,HSShear,HAShear] = FrequencySweeper_Isotropic(SweepRange,Material,PhaseVelocity,Half,Symmetric,OutputWindow1aUI1,OutputWindow1bUI1,OutputWindow2aUI1,OutputWindow2bUI1)
Bisections = 17; % ceil(log2(1e5))

%#ok<*AGROW>
HSLamb=[];HALamb=[];HBLamb=[];XRoughS=[];XRoughA=[];
cL2 = Material.LongitudinalVelocity^2;
cT2 = Material.TransverseVelocity^2;
AngularFrequency2 = (pi*SweepRange).^2*4e6;
PhaseVelocity2 = PhaseVelocity^2;
k2 = AngularFrequency2/PhaseVelocity2;
y2 = AngularFrequency2/cT2-k2;
x = sqrt(AngularFrequency2/cL2-k2);
y = sqrt(y2);
a1 = (y2-k2).^2./y;
a2 = 4*k2.*x;
a3 = tan(x*Half);
a4 = tan(y*Half);
YS = abs(a1./a3+a2./a4);
YA = abs(a1.*a3+a2.*a4);
for i = 2:length(SweepRange)-1
    if  YS(i) < YS(i-1) && YS(i) < YS(i+1)
        XRoughS(end+1) = SweepRange(i);
    end
    if  YA(i) < YA(i-1) && YA(i) < YA(i+1)
        XRoughA(end+1)  = SweepRange(i);
    end
end
% figure,plot(SweepRange,20*log10(YS))
a = PhaseVelocity*Material.TransverseVelocity/sqrt(PhaseVelocity2-cT2)/Half/4e3;
n = floor(SweepRange(end)/a);
HAShear = (1:2:n)*a;
HSShear = (2:2:n)*a;
if  isempty(XRoughS) && isempty(XRoughA) && isempty(HSShear) && isempty(HAShear)
    String = 'No higher order modes found!';
    OutputWindow1aUI1.String = String;
    OutputWindow1bUI1.String = '';
    OutputWindow2aUI1.String = '';
    OutputWindow2bUI1.String = '';
    disp(String)
    return
end    
HSLamb = Converger('S',XRoughS,Bisections,PhaseVelocity2,SweepRange,Half,cL2,cT2);
HALamb = Converger('A',XRoughA,Bisections,PhaseVelocity2,SweepRange,Half,cL2,cT2);
ModeNameLength = 1+5;
String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode     Frq.(kHz)'];
if  Symmetric
    if  any(HALamb)
        for i = 1:length(HALamb)
            String = append(String,newline,pad(['A',num2str(i)],ModeNameLength),pad(num2str(HALamb(i),'%.3f'),11,'left'));
        end
    end
else
    HBLamb = sort(horzcat(HALamb,HSLamb));
    if  any(HBLamb)
        for i = 1:length(HBLamb)
            String = append(String,newline,pad(['B',num2str(i+1)],ModeNameLength),pad(num2str(HBLamb(i),'%.3f'),11,'left'));
        end
    end
end
String = append(String,newline,' ');
if  any(HAShear)
    for i = 1:length(HAShear)
        String = append(String,newline,pad(['ASH',num2str(i)],ModeNameLength),pad(num2str(HAShear(i),'%.3f'),11,'left'));
    end
end
OutputWindow1aUI1.String = String;
disp(String)
String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode     Frq.(kHz)'];
if  any(HSLamb) && Symmetric
    for i = 1:length(HSLamb)
        String = append(String,newline,pad(['S',num2str(i)],ModeNameLength),pad(num2str(HSLamb(i),'%.3f'),11,'left'));
    end
end
String = append(String,newline,' ');
if  any(HSShear)
    for i = 1:length(HSShear)
        String = append(String,newline,pad(['SSH',num2str(i)],ModeNameLength),pad(num2str(HSShear(i),'%.3f'),11,'left'));
    end
end
OutputWindow1bUI1.String = String;
disp(extractAfter(String,')'))
disp(' ')
String = '';
if  Symmetric
    if  any(HALamb)
        String = append(String,'A:   ',num2str(length(HALamb)));
    end
else
    if  any(HBLamb)
        String = append(String,'B:   ',num2str(length(HBLamb)));
    end
end
if  any(HAShear)
    String = append(String,newline,'ASH: ',num2str(length(HAShear)));
end
OutputWindow2aUI1.String = String;
disp(String)
String = '';
if  any(HSLamb) && Symmetric
    String = append(String,'S:   ',num2str(length(HSLamb)));
end
if  any(HSShear)
    String = append(String,newline,'SSH: ',num2str(length(HSShear)));
end
OutputWindow2bUI1.String = String;
disp([String,newline,'----------------'])
end
function X = Converger(ModeFamily,XRough,Bisections,PhaseVelocity2,SweepRange,Half,cL2,cT2)
    X = [];
    for j = 1:length(XRough)
        Frequency = XRough(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            Frequency = Frequency(1):(Frequency(end)-Frequency(1))/4:Frequency(end);
            AngularFrequency2 = (pi*Frequency).^2*4e6;
            k2 = AngularFrequency2/PhaseVelocity2;
            y2 = AngularFrequency2/cT2-k2;
            x = sqrt(AngularFrequency2/cL2-k2);
            y = sqrt(y2);
            if  strcmp(ModeFamily,'S')
                Yc = (y2-k2).^2./y./tan(x*Half)+4*k2.*x./tan(y*Half);
            elseif strcmp(ModeFamily,'A')
                Yc = (y2-k2).^2./y.*tan(x*Half)+4*k2.*x.*tan(y*Half);
            end
            Y = abs(Yc);
            for i = 2:length(Frequency)-1
                if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
                    if  o == Bisections && sign(Yc(i-1)) ~= sign(Yc(i+1))
                        X(end+1) = Frequency(i);
                    end
                    Frequency = [Frequency(i-1) Frequency(i+1)];
                    break
                end
            end
        end
    end
end