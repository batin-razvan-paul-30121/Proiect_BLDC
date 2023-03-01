%%
% Plotare date
t=batin.X.Data';
u=double(batin.Y(1,3).Data');
w=double(batin.Y(1,2).Data');
y=double(batin.Y(1,1).Data');
figure
subplot(3,1,1),plot(t, u, LineWidth=1.2), grid;
title('Intrarea u');
xlabel('Timp [s]', 'FontSize',14);
ylabel('u [%]', 'FontSize',14);
subplot(3,1,2),plot(t, w, 'r', LineWidth=1.2), grid;
title('Viteza unghiulara \omega');
xlabel('Timp [s]', 'FontSize',14);
ylabel('\omega [rad/s]', 'FontSize',14);
subplot(3,1,3),plot(t, y, 'g', LineWidth=1.2), grid;
title('Pozitia/Iesirea y');
xlabel('Timp [s]', 'FontSize',14);
ylabel('y [impulsuri]', 'FontSize',14);
figure
plot(t,[u*200 w y], LineWidth=1.2), grid;
legend('Intrarea u scalata', 'Viteza unghiulara \omega', 'Pozitia/Iesirea y');
title('Vizualizarea datelor experimentale');
xlabel('Timp [s]', 'FontSize',12);
ylabel('u [%], \omega [rad/s], y [impulsuri]', 'FontSize',12);


%%
%Indecsii regasiti de pe grafic folosind cursor_info.DataIndex
i1=3814;
i2=3959;
i3=4415;
i4=4520;

%Valorile medii pentru viteza si intrare
wst=mean(w(i3:i4));
ust=mean(u(i3:i4));

w0=mean(w(i1:i2));
u0=mean(u(i1:i2));

%Factorul de proportionalitate
K=(wst-w0)/(ust-u0); %239.72

%Timpul de urcare
i5=4112;
i6=4280;
tk=t(i5:i6);

%Algoritmul Matlab pentru regresie
xk=log(wst-w(i5:i6));
Areg=[sum(tk.^2) sum(tk);
    sum(tk) length(tk)];
breg=[sum(tk.*xk);
    sum(xk)];
sol=Areg\breg;
T=-1/sol(1);

%Calcul timp mort
i7=4066;
i8=4122;
Tm=t(i8)-t(i7);

H = tf(K, [T 1], 'iodelay', Tm);

N = round(Tm/(t(2)-t(1))); 
uN = [u(1)*ones(N,1); u(1:end-N)];
A = -1/T;
B = K/T;
C = 1;
D = 0;
figure
plot(t, u*200);
wsim = lsim(A,B,C,D,uN,t,w(1));
hold on
plot(t,[w wsim]), grid;
legend('Intrarea u scalata', 'Viteza unghiulara \omega', 'Viteza unghiulara simulata \omegasim');
xlabel('Timp [s]','FontSize',12);
ylabel('u [%], \omega [rad/s], \omegasim [rad/s]', 'FontSize',12);
%eroarea medie patratica
eMIN_w = norm(w-wsim)/norm(w-mean(w));

%%
% Autocorelatia pt modelul obtinut la treapta
e = w(i1:i4)-wsim(i1:i4);    
N = i4-i1; 
R0 = sum(e.^2)/N; 
N1 = 0;
for i = 1:N
    R(i) = sum([e(1:N-i)]'*e((i+1):N))/N;
    RN(i) = R(i)/R0;
    if abs(RN(i)) <= (2.17/sqrt(N))   % 2.17/sqrt(N) = 0.0817
        N1 = i;    % N1 = 125
        break
    end
end    
Nr_e = N-N1; % Nr_e = 581
R1 = sum(e(1:N-1).*e(2:N))/N;
Rn1 = R1/R0; % Rn1 = 0.9760
R2 = sum(e(1:N-2).*e(3:N))/N;
Rn2 = R2/R0 % Rn2 = 0.8517
R125 = sum(e(1:N-124).*e(125:N))/N;
Rn125 = R125/R0; % Rn125 = 0.0819
R126 = sum(e(1:N-125).*e(126:N))/N;
Rn126 = R126/R0; % Rn126 = 0.0810 < 2.17/sqrt(N) = 0.0817
%%
%Identificare neparametrica viteza-pozitie
Ki_yw = (y(i4)-y(i3))/wst/(t(i4)-t(i3))
H_yw = tf(Ki_yw, [1, 0]);
% Spatiul starilor:
Ayw = 0;
Byw = Ki_yw;
Cyw = 1;
Dyw = 0;
ysim_yw = lsim(Ayw, Byw, Cyw, Dyw, w, t, y(1));
figure
plot(t, [y, ysim_yw]), grid;
title('Pozitia si Pozitia simulata');
legend('y', 'y simulat')
xlabel('Timp [s]','FontSize',12);
ylabel('y [impulsuri], ysim [impulsuri]', 'FontSize', 12);
eMIN_yw = norm(y-ysim_yw)/norm(y-mean(y))% 10.74
%%
% Partea 2
% autocorelatia si intercorelatia
i9=1830;
i10=3348;
i11=4688;
i12=6199;

Te=t(2)-t(1);

data_id_w=iddata(w(i9:i10),u(i9:i10),Te);
data_vd_w=iddata(w(i11:i12),u(i11:i12),Te);
figure
plot(data_id_w), grid
title('Date identificare');
xlabel('Timp [s]','FontSize',12);
ylabel('Viteza-id si SPAB-id', 'FontSize', 12);
figure
plot(data_vd_w)
title('Date validare');
xlabel('Timp [s]','FontSize',12);
ylabel('Viteza-vd si SPAB-vd', 'FontSize', 12);

data_id_poz=iddata(y(i9:i10),w(i9:i10),Te);
data_vd_poz=iddata(y(i11:i12),w(i11:i12),Te);

%%
% armax fara decimare
mw_armax=armax(data_id_w,[1,1,1,1]);
mw_armax9=armax(data_id_w,[1,1,1,9])%mai bun fit

figure
resid(data_vd_w,mw_armax);
figure
compare(data_vd_w,mw_armax,mw_armax9);

mw_armax9_inceput=tf([0 3.986],[1 -0.9846],Te,'variable','z^-1');

Mw_armax9_inceput=d2c(mw_armax9_inceput,'zoh');

Hw_armax_c = tf(4782,[1 18.48],'iodelay',9*Te)

%am facut ss ca sa dau cond initiale nule
sys_w_armax = ss(Hw_armax_c);
figure
plot(t,200*u); 
wsim2 = lsim(sys_w_armax.A, sys_w_armax.B, sys_w_armax.C, sys_w_armax.D, u, t, w(1)/sys_w_armax.C);
hold on
plot(t,[w,wsim2]), grid;
legend('Intrarea u scalata', 'Viteza unghiulara \omega', 'Viteza unghiulara simulata cu armax \omegasim-armax');
xlabel('Timp [s]','FontSize',12);
ylabel('u [%], \omega [rad/s], \omegasim [rad/s]', 'FontSize',12);

eMIN_w2 = norm(w-wsim2)/norm(w-mean(w));

%%
%decimare date viteza si pozitie
%194.2334 apare de 8 ori
idx_id=[i9:i10]';
idx_vd=[i11:i12]';
%decimare
u_decimat_id=u(i9:8:i10);
w_decimat_id=w(i9:8:i10);
y_decimat_id=y(i9:8:i10);

u_decimat_vd=u(i11:8:i12);
w_decimat_vd=w(i11:8:i12);
y_decimat_vd=y(i11:8:i12);

%date id si validare la viteza cu decimare
dw_decimat_id=iddata(w_decimat_id,u_decimat_id,8*Te);
dw_decimat_vd=iddata(w_decimat_vd,u_decimat_vd,8*Te);

%date id si validare la pozitie cu decimare
dy_decimat_id=iddata(y_decimat_id,w_decimat_id,8*Te);
dy_decimat_vd=iddata(y_decimat_vd,w_decimat_vd,8*Te);

%%
%metoda oe viteza(cu decimare);
mwd_oe=oe(dw_decimat_id,[1, 1, 1]);
figure
resid(dw_decimat_vd,mwd_oe);
figure
compare(dw_decimat_vd,mwd_oe);

mwd_oe_inceput=tf(mwd_oe.B,mwd_oe.F,Te,'iodelay',1,'variable','z^-1');
Hw_oe=d2c(mwd_oe_inceput,'zoh');



sys_w_oe = ss(Hw_oe);
figure
plot(t,200*u); 
wsim3 = lsim(sys_w_oe.A, sys_w_oe.B, sys_w_oe.C, sys_w_oe.D, u, t, w(1)/sys_w_oe.C);
hold on
plot(t,[w,wsim3]), grid;
legend('Intrarea u scalata', 'Viteza unghiulara \omega', 'Viteza unghiulara simulata cu oe \omegasim-oe');
xlabel('Timp [s]','FontSize',12);
ylabel('u [%], \omega [rad/s], \omegasim [rad/s]', 'FontSize',12);

eMIN_w3 = norm(w-wsim3)/norm(w-mean(w))

%%
%armax pozitie cu decimare
myd_armax=armax(dy_decimat_id,[1,1,1,0]);
figure
resid(dy_decimat_vd,myd_armax);
figure
compare(dy_decimat_vd,myd_armax);%e validata

b1 = myd_armax.B;
Te_ARMAX = 0.00672; % Perioada de esnationare dupa decimare

% Constanta de integrare:
Ki_nou = b1/Te_ARMAX;

Hyw_c_2 = tf(Ki_nou, [1, 0]);
sys_y = ss(Hyw_c_2);
ysim2 = lsim(sys_y.A, sys_y.B, sys_y.C, sys_y.D, w, t, y(1)/sys_y.C);

figure
plot(t,[y, ysim2]), grid;
title('Pozitia si Pozitia simulata cu armax');
legend('y', 'y simulat-armax')
xlabel('Timp [s]','FontSize',12);
ylabel('y [impulsuri], ysim_armax [impulsuri]', 'FontSize', 12);

eMIN_y2 = norm(y-ysim2)/norm(y-mean(y));
%%
%metoda oe pozitie cu decimare
myd_oe=oe(dy_decimat_id,[1, 1, 0]);
figure
resid(dy_decimat_vd,myd_oe);
figure
compare(dy_decimat_vd,myd_oe);%e validata

Hyw_oe = tf(myd_oe.B, myd_oe.F, Te, 'variable', 'z^-1');
b1_oe = myd_oe.B;
Te_oe = 0.00672;
Ki_nou_oe = b1_oe/Te_oe;
Hyw_c_3 = tf(Ki_nou_oe, [1, 0]);
sys_y_oe = ss(Hyw_c_3);
ysim3 = lsim(sys_y_oe.A, sys_y_oe.B, sys_y_oe.C, sys_y_oe.D, w, t, y(1)/sys_y_oe.C);

figure
plot(t,[y, ysim3]), grid;
title('Pozitia si Pozitia simulata cu oe');
legend('y', 'y simulat-oe')
xlabel('Timp [s]','FontSize',12);
ylabel('y [impulsuri], ysim_oe [impulsuri]', 'FontSize', 12);

eMIN_y3 = norm(y-ysim3)/norm(y-mean(y));