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
function [AModes,ALambModes,AShearModes,BModes,BLambModes,BShearModes,SModes,SLambModes,SShearModes,Frequency] = ModeFinder_Anisotropic(A,ALamb,AShear,B,BLamb,BShear,S,SLamb,SShear,Frequency)
%#ok<*EFIND>
AModes = false(1);
ALambModes = false(1);
AShearModes = false(1);
BModes = false(1);
BLambModes = false(1);
BShearModes = false(1);
SModes = false(1);
SLambModes = false(1);
SShearModes = false(1);
String = 'MODES:';
if  ~isempty(A{1})
    for p = 1:length(A)
        q = find(A{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(A{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(A{p}(end,1)) || (isempty(z) && Frequency < A{p}(1,1)) || (~isempty(z) && Frequency < A{p}(z,1))
                AModes(p) = 0;
            else
                [~,q] = min(abs(A{p}(:,1)-Frequency));
                Frequency = A{p}(q,1);
                AModes(p) = 1;
            end
        else
            AModes(p) = 1;
        end
        if  p == 10
            break
        end
    end
    if  any(AModes)
        String = append(String,newline,'A: ',num2str(numel(find(AModes == 1))),'  [',num2str(AModes),']');
    end
end
if  ~isempty(S{1})
    for p = 1:length(S)
        q = find(S{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(S{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(S{p}(end,1)) || (isempty(z) && Frequency < S{p}(1,1)) || (~isempty(z) && Frequency < S{p}(z,1))
                SModes(p) = 0;
            else
                [~,q] = min(abs(S{p}(:,1)-Frequency));
                Frequency = S{p}(q,1);
                SModes(p) = 1;
            end
        else
            SModes(p) = 1;
        end
        if  p == 10
            break
        end
    end
    if  any(SModes)
        String = append(String,newline,'S: ',num2str(numel(find(SModes == 1))),'  [',num2str(SModes),']');
    end
end
if  ~isempty(B{1})
    for p = 1:length(B)
        q = find(B{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(B{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(B{p}(end,1)) || (isempty(z) && Frequency < B{p}(1,1)) || (~isempty(z) && Frequency < B{p}(z,1))
                BModes(p) = 0;
            else
                [~,q] = min(abs(B{p}(:,1)-Frequency));
                Frequency = B{p}(q,1);
                BModes(p) = 1;
            end
        else
            BModes(p) = 1;
        end
        if  p == 10
            break
        end
    end
    if  any(BModes)
        String = append(String,newline,'B: ',num2str(numel(find(BModes == 1))),'  [',num2str(BModes),']');
    end
end
if  ~isempty(ALamb{1})
    for p = 1:length(ALamb)
        q = find(ALamb{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(ALamb{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(ALamb{p}(end,1)) || (isempty(z) && Frequency < ALamb{p}(1,1)) || (~isempty(z) && Frequency < ALamb{p}(z,1))
                ALambModes(p) = 0;
            else
                [~,q] = min(abs(ALamb{p}(:,1)-Frequency));
                Frequency = ALamb{p}(q,1);
                ALambModes(p) = 1;
            end
        else
            ALambModes(p) = 1;
        end
        if  p == 10
            break
        end
    end
    if  any(ALambModes)
        String = append(String,newline,'A:   ',num2str(numel(find(ALambModes == 1))),'  [',num2str(ALambModes),']');
    end
end
if  ~isempty(SLamb{1})
    for p = 1:length(SLamb)
        q = find(SLamb{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(SLamb{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(SLamb{p}(end,1)) || (isempty(z) && Frequency < SLamb{p}(1,1)) || (~isempty(z) && Frequency < SLamb{p}(z,1))
                SLambModes(p) = 0;
            else
                [~,q] = min(abs(SLamb{p}(:,1)-Frequency));
                Frequency = SLamb{p}(q,1);
                SLambModes(p) = 1;
            end
        else
            SLambModes(p) = 1;
        end
        if  p == 10
            break
        end
    end
    if  any(SLambModes)
        String = append(String,newline,'S:   ',num2str(numel(find(SLambModes == 1))),'  [',num2str(SLambModes),']');
    end
end
if  ~isempty(BLamb{1})
    for p = 1:length(BLamb)
        q = find(BLamb{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(BLamb{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(BLamb{p}(end,1)) || (isempty(z) && Frequency < BLamb{p}(1,1)) || (~isempty(z) && Frequency < BLamb{p}(z,1))
                BLambModes(p) = 0;
            else
                [~,q] = min(abs(BLamb{p}(:,1)-Frequency));
                Frequency = BLamb{p}(q,1);
                BLambModes(p) = 1;
            end
        else
            BLambModes(p) = 1;
        end
        if  p == 10
            break
        end
    end
    if  any(BLambModes)
        String = append(String,newline,'B:   ',num2str(numel(find(BLambModes == 1))),'  [',num2str(BLambModes),']');
    end
end
if  ~isempty(AShear{1})
    for p = 1:length(AShear)
        q = find(AShear{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(AShear{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(AShear{p}(end,1)) || (isempty(z) && Frequency < AShear{p}(1,1)) || (~isempty(z) && Frequency < AShear{p}(z,1))
                AShearModes(p) = 0;
            else
                [~,q] = min(abs(AShear{p}(:,1)-Frequency));
                Frequency = AShear{p}(q,1);
                AShearModes(p) = 1;
            end
        else
            AShearModes(p) = 1;
        end
        if  p == 10
            break
        end
    end
    if  any(AShearModes)
        String = append(String,newline,'ASH: ',num2str(numel(find(AShearModes == 1))),'  [',num2str(AShearModes),']');
    end
end
if  ~isempty(SShear{1})
    for p = 1:length(SShear)
        q = find(SShear{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(SShear{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(SShear{p}(end,1)) || (isempty(z) && Frequency < SShear{p}(1,1)) || (~isempty(z) && Frequency < SShear{p}(z,1))
                SShearModes(p) = 0;
            else
                [~,q] = min(abs(SShear{p}(:,1)-Frequency));
                Frequency = SShear{p}(q,1);
                SShearModes(p) = 1;
            end
        else
            SShearModes(p) = 1; 
        end
        if  p == 10
            break
        end
    end
    if  any(SShearModes)
        String = append(String,newline,'SSH: ',num2str(numel(find(SShearModes == 1))),'  [',num2str(SShearModes),']');
    end
end
if  ~isempty(BShear{1})
    for p = 1:length(BShear)
        q = find(BShear{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(BShear{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(BShear{p}(end,1)) || (isempty(z) && Frequency < BShear{p}(1,1)) || (~isempty(z) && Frequency < BShear{p}(z,1))
                BShearModes(p) = 0;
            else
                [~,q] = min(abs(BShear{p}(:,1)-Frequency));
                Frequency = BShear{p}(q,1);
                BShearModes(p) = 1;
            end
        else
            BShearModes(p) = 1; 
        end
        if  p == 10
            break
        end
    end
    if  any(BShearModes)
        String = append(String,newline,'BSH: ',num2str(numel(find(BShearModes == 1))),'  [',num2str(BShearModes),']');
    end
end
String = append(String,newline);
% disp(String)