%Editor : M. Taufiqul Huda (5022211007)
% Pembangkit diagram energy band sambungan PN pada kondisi equilibrium
% Dari buku Robert Pierret

% Konstanta
T = 300;        % suhu (Kelvin)
k = 8.617e-5;   % Konstanta Boltzmann (eV/K)
e0 = 8.85e-14;  % permitivitas free space (F/cm)
q = 1.602e-19;  % muatan elektron (coul, C)
KS = 11.8;      % konstanta dielektrik Si (silicon)
ni = 1.0e10;    % konsentrasi charge carrier intrinsic pada suhu 300K
EG = 1.12;      % Energy band gap Silicon (eV)

% Konstanta kontrol gambar
xleft = -3.5e-4;
xright = -xleft;

NA = input('Isikan nilai doping sisi-p (cm^-3), NA = ');
ND = input('Isikan nilai doping sisi-n (cm^-3), ND = ');
V = input('isikan nilai tegangan bias (volt), V = ');

% Komputasi
Vbi = k*T*log((NA*ND)/ni^2);

xN = sqrt(2*KS*e0/q*NA*(Vbi+V)/(ND*(NA+ND)));
xP = sqrt(2*KS*e0/q*ND*(Vbi-V)/(NA*(NA+ND)));
x = linspace(xleft, xright, 200);

Vx1 = ((Vbi+V)-q*ND.*(xN-x).^2/(2*KS*e0).*(x<=xN)).*(x>=0);
Vx2 = 0.5*q*NA.*(xP+x).^2/(KS*e0).*(x>=-xP & x<0);

Vx = Vx1 + Vx2;
VMAX = 5;
EFn = -Vx(1) + VMAX/2-k*T*log(NA/ni);
EFp = -Vx(200) + VMAX/2+k*T*log(ND/ni);

% Plot gambar
close

plot(x, -Vx+EG/2+VMAX/2);
axis([xleft xright 0 VMAX]);
axis('off'); hold on

plot(x, -Vx-EG/2+VMAX/2);
plot(x, -Vx+VMAX/2,'g:');
plot([xleft 0],[EFn EFn],'g');
plot([0 xright],[EFp EFp],'g');
plot([0 0],[0.15 VMAX-0.5], 'g--');

text(xleft*1.08,(-Vx(1)+EG/2+VMAX/2-.05),'Ec');
text(xright*1.02,(-Vx(200)+EG/2+VMAX/2-.05),'Ec');

text(xleft*1.08,(-Vx(1)-EG/2+VMAX/2-.05),'Ev');
text(xright*1.02,(-Vx(200)-EG/2+VMAX/2-.05),'Ev');

text(xleft*1.08,(-Vx(1)++VMAX/2-.05),'Ei');

text(xright*1.02, EFp-0.05,'EF');

text(.18, 0, 'p-side');
text(.47, 0, 'x=0');
text(.75, 0, 'n-side');
text(xleft*0.5,VMAX-0.3,'Tipe-N');
text(xright*0.5,VMAX-0.3,'Tipe-P');
hold off

