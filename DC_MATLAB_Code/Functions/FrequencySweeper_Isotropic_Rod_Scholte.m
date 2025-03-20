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
function [HLScholte,HFScholte] = FrequencySweeper_Isotropic_Rod_Scholte(Material,R,Fluid,SweepRange,FlexuralModeOrders)
Resolution = 1e-5; % (kHz)
PhaseVelocity = (1-1e-4)*Fluid.Velocity;

%#ok<*AGROW>
HLScholte = [];
HFScholte = [];
XRoughLScholte = [];
XRoughFScholte = zeros(FlexuralModeOrders,1);
if  Resolution < abs(SweepRange(1)-SweepRange(2))
    Bisections = ceil(log2(Resolution/(abs(SweepRange(1)-SweepRange(2))))/log2(2*.25));
else
    Bisections = 1;
end
R2 = R^2;
Density = Fluid.Density/Material.Density;
AngularFrequency = 2*pi*SweepRange*1e3;
k = AngularFrequency/PhaseVelocity;
k2 = k.^2;
kT2 = (AngularFrequency/Material.TransverseVelocity).^2;
y2 = kT2-k2;
xR = sqrt((AngularFrequency/Material.LongitudinalVelocity).^2-k2)*R;
yR = sqrt(kT2-k2)*R;
zR = sqrt((AngularFrequency/Fluid.Velocity).^2-k2)*R;
J0x = besselj(0,xR);
J0y = besselj(0,yR);
J1xR = xR.*besselj(1,xR);
J1yR = yR.*besselj(1,yR);
H0z = besselh(0,zR);
H1zR = -zR.*besselh(1,zR);
for i = 1:length(SweepRange)
    M(1,1) = R2*(.5*kT2(i)-k2(i))*J0x(i)-J1xR(i);
    M(1,2) = k(i)*(J1yR(i)-y2(i)*R2*J0y(i));
    M(1,3) = .5*kT2(i)*R2*Density*H0z(i);
    M(2,1) = -2*k(i)*J1xR(i);
    M(2,2) = (k2(i)-y2(i))*J1yR(i);
    M(3,1) = -J1xR(i);
    M(3,2) = k(i)*J1yR(i);
    M(3,3) = H1zR(i);
    Y(i) = abs(det(M));
    if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
        XRoughLScholte(end+1) = SweepRange(i-1);
    end
end
for n = 1:FlexuralModeOrders
    n2 = n^2;
    Jnx = besselj(n,xR);
    Jny = besselj(n,yR);
    dJnxR = n*Jnx-xR.*besselj(n+1,xR);
    dJnyR = n*Jny-yR.*besselj(n+1,yR);
    Hnz = besselh(n,zR);
    dHnzR = n*Hnz-zR.*besselh(n+1,zR);
    for i = 1:length(SweepRange)
        M(1,1) = dJnxR(i)+(.5*kT2(i)*R2-k2(i)*R2-n2)*Jnx(i);
        M(1,2) = -k(i)*(dJnyR(i)+(y2(i)*R2-n2)*Jny(i));
        M(1,3) = n*(dJnyR(i)-Jny(i));
        M(1,4) = .5*kT2(i)*R2*Density*Hnz(i);
        M(2,1) = 2*n*(dJnxR(i)-Jnx(i));
        M(2,2) = 2*k(i)*n*(Jny(i)-dJnyR(i));
        M(2,3) = 2*dJnyR(i)+(y2(i)*R2-2*n2)*Jny(i);
        M(3,1) = 2*k(i)*dJnxR(i);
        M(3,2) = (y2(i)-k2(i))*dJnyR(i);
        M(3,3) = -k(i)*n*Jny(i);
        M(4,1) = dJnxR(i);
        M(4,2) = -k(i)*dJnyR(i);
        M(4,3) = -n*Jny(i);
        M(4,4) = dHnzR(i);
        Y(i) = abs(det(M));
        if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
            XRoughFScholte(n,1+numel(find(XRoughFScholte(n,:) > 0))) = SweepRange(i-1);
        end
    end
    if  all(XRoughFScholte(n,:) == 0)
        XRoughFScholte(n:end,:) = [];
        break
    end
end
if  isempty(XRoughLScholte) && all(all(XRoughFScholte == 0))
    disp('No higher order Scholte modes found!')
    return
end
M = 0;
for j = 1:length(XRoughLScholte)
    Frequency = [XRoughLScholte(j)-(SweepRange(2)-SweepRange(1)) XRoughLScholte(j)+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
        AngularFrequency = 2*pi*Frequency*1e3;
        for i = 1:length(Frequency)
            k = AngularFrequency(i)/PhaseVelocity;
            k2 = k^2;
            kT2 = (AngularFrequency(i)/Material.TransverseVelocity)^2;
            y2 = kT2-k2;
            xR = sqrt((AngularFrequency(i)/Material.LongitudinalVelocity)^2-k2)*R;
            yR = sqrt(y2)*R;
            zR = sqrt((AngularFrequency(i)/Fluid.Velocity)^2-k2)*R;
            J0x = besselj(0,xR);
            J0y = besselj(0,yR);
            J1xR = xR*besselj(1,xR);
            J1yR = yR*besselj(1,yR);
            H0z = besselh(0,zR);
            H1zR = -zR*besselh(1,zR);
            M(1,1) = R2*(.5*kT2-k2)*J0x-J1xR;
            M(1,2) = k*(J1yR-y2*R2*J0y);
            M(1,3) = .5*kT2*R2*Density*H0z;
            M(2,1) = -2*k*J1xR;
            M(2,2) = (k2-y2)*J1yR;
            M(3,1) = -J1xR;
            M(3,2) = k*J1yR;
            M(3,3) = H1zR;
            Y(i) = abs(det(M));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    HLScholte(end+1) = Frequency(i-1);
                end
                Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                break
            end
        end
    end
end
HFScholte = zeros(size(XRoughFScholte));
for n = 1:height(XRoughFScholte)
    n2 = n^2;
    for j = 1:numel(find(XRoughFScholte(n,:) > 0))
        Frequency = [XRoughFScholte(n,j)-(SweepRange(2)-SweepRange(1)) XRoughFScholte(n,j)+(SweepRange(2)-SweepRange(1))];
        for o = 1:Bisections
            Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
            AngularFrequency = 2*pi*Frequency*1e3;
            for i = 1:length(Frequency)
                k = AngularFrequency(i)/PhaseVelocity;
                k2 = k^2;
                kT2 = (AngularFrequency(i)/Material.TransverseVelocity)^2;
                y2 = kT2-k2;
                xR = sqrt((AngularFrequency(i)/Material.LongitudinalVelocity)^2-k2)*R;
                yR = sqrt(y2)*R;
                zR = sqrt((AngularFrequency(i)/Fluid.Velocity)^2-k2)*R;
                Jnx = besselj(n,xR);
                Jny = besselj(n,yR);
                dJnxR = n*Jnx-xR*besselj(n+1,xR);
                dJnyR = n*Jny-yR*besselj(n+1,yR);
                Hnz = besselh(n,zR);
                dHnzR = n*Hnz-zR*besselh(n+1,zR);
                M(1,1) = dJnxR+(.5*kT2*R2-k2*R2-n2)*Jnx;
                M(1,2) = -k*(dJnyR+(y2*R2-n2)*Jny);
                M(1,3) = n*(dJnyR-Jny);
                M(1,4) = .5*kT2*R2*Density*Hnz;
                M(2,1) = 2*n*(dJnxR-Jnx);
                M(2,2) = 2*k*n*(Jny-dJnyR);
                M(2,3) = 2*dJnyR+(y2*R2-2*n2)*Jny;
                M(3,1) = 2*k*dJnxR;
                M(3,2) = (y2-k2)*dJnyR;
                M(3,3) = -k*n*Jny;
                M(4,1) = dJnxR;
                M(4,2) = -k*dJnyR;
                M(4,3) = -n*Jny;
                M(4,4) = dHnzR;
                Y(i) = abs(det(M));
                if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                    if  o == Bisections
                        HFScholte(n,j) = Frequency(i-1);
                    end
                    Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                    break
                end
            end
        end
    end
end
String = ['Frq. @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode              Frq.(kHz)'];
if  any(HFScholte)
    for i = 1:numel(find(HFScholte(1,:) > 0))
        if  i < 9
            if  HFScholte(1,i) < 1e2
                String = append(String,newline,'FScholte(1,',num2str(i+1),')       ',num2str(HFScholte(1,i),'%.3f'));
            elseif HFScholte(1,i) >= 1e2 && HFScholte(1,i) < 1e3
                String = append(String,newline,'FScholte(1,',num2str(i+1),')      ',num2str(HFScholte(1,i),'%.3f'));
            elseif HFScholte(1,i) >= 1e3 && HFScholte(1,i) < 1e4
                String = append(String,newline,'FScholte(1,',num2str(i+1),')     ',num2str(HFScholte(1,i),'%.3f'));
            elseif HFScholte(1,i) >= 1e4
                String = append(String,newline,'FScholte(1,',num2str(i+1),')    ',num2str(HFScholte(1,i),'%.3f'));
            end
        else
            if  HFScholte(1,i) < 1e2
                String = append(String,newline,'FScholte(1,',num2str(i+1),')      ',num2str(HFScholte(1,i),'%.3f'));
            elseif HFScholte(1,i) >= 1e2 && HFScholte(1,i) < 1e3
                String = append(String,newline,'FScholte(1,',num2str(i+1),')     ',num2str(HFScholte(1,i),'%.3f'));
            elseif HFScholte(1,i) >= 1e3 && HFScholte(1,i) < 1e4
                String = append(String,newline,'FScholte(1,',num2str(i+1),')    ',num2str(HFScholte(1,i),'%.3f'));
            elseif HFScholte(1,i) >= 1e4
                String = append(String,newline,'FScholte(1,',num2str(i+1),')   ',num2str(HFScholte(1,i),'%.3f'));
            end
        end
    end    
    for n = 2:height(HFScholte)
        for i = 1:numel(find(HFScholte(n,:) > 0))
            if  i < 10
                if  HFScholte(n,i) < 1e2
                    String = append(String,newline,'FScholte(',num2str(n),',',num2str(i),')       ',num2str(HFScholte(n,i),'%.3f'));
                elseif HFScholte(n,i) >= 1e2 && HFScholte(n,i) < 1e3
                    String = append(String,newline,'FScholte(',num2str(n),',',num2str(i),')      ',num2str(HFScholte(n,i),'%.3f'));
                elseif HFScholte(n,i) >= 1e3 && HFScholte(n,i) < 1e4
                    String = append(String,newline,'FScholte(',num2str(n),',',num2str(i),')     ',num2str(HFScholte(n,i),'%.3f'));
                elseif HFScholte(n,i) >= 1e4
                    String = append(String,newline,'FScholte(',num2str(n),',',num2str(i),')    ',num2str(HFScholte(n,i),'%.3f'));
                end
            else
                if  HFScholte(n,i) < 1e2
                    String = append(String,newline,'FScholte(',num2str(n),',',num2str(i),')      ',num2str(HFScholte(n,i),'%.3f'));
                elseif HFScholte(n,i) >= 1e2 && HFScholte(n,i) < 1e3
                    String = append(String,newline,'FScholte(',num2str(n),',',num2str(i),')     ',num2str(HFScholte(n,i),'%.3f'));
                elseif HFScholte(n,i) >= 1e3 && HFScholte(n,i) < 1e4
                    String = append(String,newline,'FScholte(',num2str(n),',',num2str(i),')    ',num2str(HFScholte(n,i),'%.3f'));
                elseif HFScholte(n,i) >= 1e4
                    String = append(String,newline,'FScholte(',num2str(n),',',num2str(i),')   ',num2str(HFScholte(n,i),'%.3f'));
                end
            end
        end
    end
end
String = append(String,newline);
if  any(HLScholte)
    for i = 1:length(HLScholte)
        if  i < 9
            if  HLScholte(i) < 1e2
                String = append(String,newline,'LScholte(0,',num2str(i),')       ',num2str(HLScholte(i),'%.3f'));
            elseif HLScholte(i) >= 1e2 && HLScholte(i) < 1e3
                String = append(String,newline,'LScholte(0,',num2str(i),')      ',num2str(HLScholte(i),'%.3f'));
            elseif HLScholte(i) >= 1e3 && HLScholte(i) < 1e4
                String = append(String,newline,'LScholte(0,',num2str(i),')     ',num2str(HLScholte(i),'%.3f'));
            elseif HLScholte(i) >= 1e4
                String = append(String,newline,'LScholte(0,',num2str(i),')    ',num2str(HLScholte(i),'%.3f'));
            end
        else
            if  HLScholte(i) < 1e2
                String = append(String,newline,'LScholte(0,',num2str(i),')      ',num2str(HLScholte(i),'%.3f'));
            elseif HLScholte(i) >= 1e2 && HLScholte(i) < 1e3
                String = append(String,newline,'LScholte(0,',num2str(i),')     ',num2str(HLScholte(i),'%.3f'));
            elseif HLScholte(i) >= 1e3 && HLScholte(i) < 1e4
                String = append(String,newline,'LScholte(0,',num2str(i),')    ',num2str(HLScholte(i),'%.3f'));
            elseif HLScholte(i) >= 1e4
                String = append(String,newline,'LScholte(0,',num2str(i),')   ',num2str(HLScholte(i),'%.3f'));
            end
        end
    end
end
String = append(String,newline,newline);
if  any(HFScholte)
    String = append(String,'FScholte: ',num2str(numel(find(HFScholte > 0))));
end
String = append(String,newline);
if  any(HLScholte)
    String = append(String,'LScholte: ',num2str(length(HLScholte)));
end
disp([String,newline,'----------------'])