%%ENGI9496 Project Simulation
%Jordan Mayo
%Feb/Mar 2021
%
%
%Use script to determine recommended collision avoidance boundaries for 
%encounters with uncompliant obstacles
%
%For this setup, both ships are the same, modelled after a test Laker 
%for which the NRC has completed testing to obtain the required
%coefficients. The ship can be changed by modifying the Ship parameter
%blocks

clear all

%Ship 1 Parameters
WLL1 = 182.8; %Water line length [m]
B1 = 26.7; %Beam width [m]
T1 = 5; %Draft [m]
XudotCoef1 = -0.1; %Xudot coefficient 
A_rudder1 = 25.4; %Rudder area [m^2]
Delta_max1 = 5*pi/180; %max rudder angle [rads]
u01 = 10; %Reference speed [knots]

SV1 = 14984; %Submerged volume [m^3]
WP1 = 4393; %Waterplane area [m^2]
MS1 = 120.2; %Midship section area [m^2]
LCB1 = 0.0319*WLL1; %Longitudinal centre of buoyancy relative to centre [m]
AT1 = B1*T1; %Transom area [m^2]
k2Coef1 = 1.5; %Factor based on physical geometry of hull;
WAA1 = 50; %Wetted area of appendages [m^2]
xG1 = 0; %Location of centre of gravity from ship centre, P [m]
hwind1 = 17.52; %height of ship above waterline [m]

WSA1_ax = WLL1*(B1+T1); %Wetted Surface area; may need to be checked [m^2]
Awind1_ax = hwind1*B1; %Area of ship exposed to the wind [m^2]
WSA1_tr = WLL1*T1; %Transverse area exposed to current [m^2]
Awind1_tr = WLL1*hwind1; %Transverse area exposed to wind [m^2]

%Drag coefficients
CWA1 = 1;
CWN1 = 1;
CCN1 = 1.2;

%Strip method parameters
n1 = 5; %number of segments
Lseg1 = WLL1/n1; %Segment lengtH

%Normal calcs
ncoord1 = zeros(1,n1);
ACseg1 = WSA1_tr/n1;
AWseg1 = Awind1_tr/n1;

%Ship 2 Parameters
WLL2 = 182.8; %Water line length [m]
B2 = 26.7; %Beam width [m]
T2 = 5; %Draft [m]
XudotCoef2 = -0.1; %Xudot coefficient 
A_rudder2 = 25.4; %Rudder area [m^2]
Delta_max2 = 5*pi/180; %max rudder angle [rads]
u02 = 10; %Reference speed [knots]

SV2 = 14984; %Submerged volume [m^3]
WP2 = 4393; %Waterplane area [m^2]
MS2 = 120.2; %Midship section area [m^2]
LCB2 = 0.0319*WLL1; %Longitudinal centre of buoyancy relative to centre [m]
AT2 = B1*T1; %Transom area [m^2]
k2Coef2 = 1.5; %Factor based on physical geometry of hull;
WAA2 = 50; %Wetted area of appendages [m^2]
xG2 = 0; %Location of centre of gravity from ship centre, P [m]
hwind2 = 17.52; %height of ship above waterline [m]

WSA2_ax = WLL2*(B2+T2); %Wetted Surface area; may need to be checked
Awind2_ax = hwind2*B2; %Area of ship exposed to the wind [m^2]
WSA2_tr = WLL2*T2; %Transverse area exposed to current [m^2]
Awind2_tr = WLL2*hwind2; %Transverse area exposed to wind [m^2]

%Drag coefficients
CWA2 = 1;
CWN2 = 1;
CCN2 = 1.2;

%Strip method parameters
n2 = 5; %number of segments
Lseg2 = WLL2/n2; %Segment length

%Normal calcs
ncoord2 = zeros(1,n2);
ACseg2 = WSA2_tr/n2;
AWseg2 = Awind2_tr/n2;

%Constants
rho = 1025; %Density salt water [kg/m^3]
nu = 0.00000118831; %Kinematic viscosity [m^2/s]
g = 9.81; %Gravity [m/s^2]

%Wind/Current Area calcs
DCseg1 = 0.5*rho*ACseg1*CCN1;
DWseg1 = 0.5*rho*AWseg1*CWN1;
for k=1:n1
    ncoord1(k)=-WLL1/2+k*Lseg1-Lseg1/2;    
end

DCseg2 = 0.5*rho*ACseg2*CCN2;
DWseg2 = 0.5*rho*AWseg2*CWN2;
for k=1:n2
    ncoord2(k)=-WLL2/2+k*Lseg2-Lseg2/2;    
end

%Ship 1 Calculations

CB1 = SV1/(WLL1*B1*T1); %Block Coefficient
m1 = CB1*WLL1*B1*T1*rho; %Mass of ship [kg]
Xudot1 = XudotCoef1*m1; %Xudot parameter
IPz1 = (1/12)*m1*(WLL1^2+B1^2); %moment of inertia about vertical axis through ship centre [kg*m^2]

%Dimensionless parameter conversions

u01 = u01*0.514444; %Convert reference speed to [m/s]
mp1 = m1/(0.5*rho*WLL1^3); %m1'
IPz1p = IPz1/(0.5*rho*WLL1^5); %IPz1'
xG1p = xG1/WLL1; %xG1'


Xudot1p = Xudot1/(0.5*rho*WLL1^3); %Xudot'

%Hydrodynamic Derivatives from Clarke (1982)
Yv1p = -pi*(T1/WLL1)^2*(1+0.3*CB1*B1/T1);
Yvdot1p = -pi*(T1/WLL1)^2*(1+0.16*CB1*B1/T1-5.1*(B1/T1)^2);

Yv1 = Yv1p*rho*WLL1^2*u01/2;
Yvdot1 = Yvdot1p*rho*WLL1^4/2;

Yr1p = -pi*(T1/WLL1)^2*(-1/2+2.2*B1/WLL1-0.08*B1/T1);
Yrdot1p = -pi*(T1/WLL1)^2*(0.67*B1/WLL1-0.0033*(B1/T1)^2);

Yr1 = Yr1p*rho*WLL1^3*u01/2;
Yrdot1 = Yrdot1p*rho*WLL1^4/2;

Nv1p = -pi*(T1/WLL1)^2*(1/2+2.4*T1/WLL1);
Nvdot1p = -pi*(T1/WLL1)^2*(1.1*B1/WLL1-0.041*B1/T1);

Nv1 = Nv1p*rho*WLL1^3*u01/2;
Nvdot1 = Nvdot1p*rho*WLL1^4/2;

Nr1p = -pi*(T1/WLL1)^2*(1/4+0.039*B1/T1-0.56*B1/WLL1);
Nrdot1p = -pi*(T1/WLL1)^2*(1/12+0.017*CB1*B1/T1-0.33*B1/WLL1);

Nr1 = Nr1p*rho*WLL1^4*u01/2;
Nrdot1 = Nrdot1p*rho*WLL1^5/2;

alph = 3; %Rudder factor from Clarke (1982)
Ydelta1p = A_rudder1/WLL1^2*alph;
Ndelta1p = -0.5*Ydelta1p;

Ydelta1 = Ydelta1p*rho*WLL1^2*u01^2/2;
Ndelta1 = Ndelta1p*rho*WLL1^3*u01^2/2;

%Total Water Resistance RT = RF(1+k1)+RAPP+RW+RTR+RA

CM1 = MS1/(B1*T1); %Midship section area coefficient
CWP1 = WP1/(B1*WLL1); %Waterplane area coefficient


CP1 = SV1/(WP1*WLL1); %Prismatic coefficient
CLCB1 = LCB1/WLL1*100; %Buoyancy coefficient

WA1 = WLL1*(2*T1+B1)*sqrt(CM1)*(0.453+0.4425*CB1-0.2862*CM1-0.003467*B1/T1+0.3696*CWP1); %Wetted area approximation [m^2]

%Find k1
LR1 = (1-CP1+0.06*CP1*CLCB1/(4*CP1-1))*WLL1; 

if T1/WLL1 > 0.05
    c121 = (T1/WLL1)^0.2228446;
elseif T1/WLL1 >0.02
    c121 = 48.2*(T1/WLL1-0.02)^2.078 +0.479948;
else
    c121 = 0.479948;
end

k1 = 0.93+c121*(B1/LR1)^0.92497*(0.95-CP1)^(-0.521448)*(1-CP1+0.0225*CLCB1)^0.6906 -1;

%RF and RAPP depend on speed

%RW, wave resistance
if B1/WLL1 < 0.11
    c71 = 0.229577*(B1/WLL1)^0.33333;
elseif B1/WLL1 < 0.25
    c71=B1/WLL1;
else
    c71 = 0.5-0.0625*WLL1/B1;
end

c51 = 1-0.8*AT1/(B1*T1*CM1);

if WLL1/B1<12
    lambda1 = 1.446*CP1 - 0.03*WLL1/B1;
else
    lambda1 = 1.446*CP1 - 0.36;
end

if CP1<0.8
    c161 = 8.07981*CP1 - 13.8673*CP1^2 + 6.984388*CP1^3;
else
    c161 = 1.73014 - 0.7067*CP1;
end

m11 = 0.0140407*WLL1/T1 - 1.75254*SV1^(1/3)/WLL1 - 4.79323*B1/WLL1 - c161;

if WLL1^3/SV1 < 512
    c151 = -1.69385;
elseif WLL1^3/SV1 < 1727
    c151 = -1.69385 + (WLL1/SV1^(1/3)-8)/2.36;
else
    c151 = 0;
end

d1=-0.9;

iE = 1+89*exp(-(WLL1/B1)^0.80856*(1-CWP1)^0.30484*(1-CP1-0.0225*CLCB1)^0.6367*(LR1/B1)^0.34574*(100*SV1/WLL1^3)^0.16302);
c11 = 2223105*c71^3.78613*(T1/B1)^1.07961*(90-iE)^(-1.37565);

%Wave Resistance depends on speed and is calculated in loop

%RTR, Resistance from transom stern calculated in loop

%RA, model ship correlation resistance

if T1/WLL1 < 0.04
    c41 = T1/WLL1;
else
    c41 = 0.04;
end

CA1 = 0.006*(WLL1+100)^(-0.16) - 0.00205 + 0.003*sqrt(WLL1/7.5)*CB1^4*(0.04-c41);

%Model Ship Correlation Resistance depends on speed

%Using First Order Nomoto model by assuming fore/aft symmetry for the
%vessel - valid assumption for this vessel

K = -(u01/WLL1)*Ndelta1p/Nr1p;
T = -(WLL1/u01)*(IPz1p-Nrdot1p)/Nr1p;

%Yr1p = 0;
%Yrdot1p = 0;
%Nv1p = 0;
%Nvdot1p = 0;

%Ship 2 Calculations

CB2 = SV2/(WLL2*B2*T2); %Block Coefficient
m2 = CB2*WLL2*B2*T2*rho; %Mass of ship [kg]
Xudot2 = XudotCoef2*m2; %Xudot parameter
IPz2 = (1/12)*m2*(WLL2^2+B2^2); %moment of inertia about vertical axis through ship centre [kg*m^2]

%Dimensionless parameter conversions

u02 = u02*0.514444; %Convert reference speed to [m/s]
mp2 = m2/(0.5*rho*WLL2^3); %m2'
IPz2p = IPz2/(0.5*rho*WLL2^5); %IPz2'
xG2p = xG2/WLL2; %xG2'

Xudot2p = Xudot2/(0.5*rho*WLL2^3); %Xudot'

%Hydrodynamic Derivatives from Clarke (1982)
Yv2p = -pi*(T2/WLL2)^2*(1+0.3*CB2*B2/T2);
Yvdot2p = -pi*(T2/WLL2)^2*(1+0.16*CB2*B2/T2-5.1*(B2/T2)^2);

Yv2 = Yv2p*rho*WLL2^2*u02/2;
Yvdot2 = Yvdot2p*rho*WLL2^4/2;

Yr2p = -pi*(T2/WLL2)^2*(-1/2+2.2*B2/WLL2-0.08*B2/T2);
Yrdot2p = -pi*(T2/WLL2)^2*(0.67*B2/WLL2-0.0033*(B2/T2)^2);

Yr2 = Yr2p*rho*WLL2^3*u02/2;
Yrdot2 = Yrdot2p*rho*WLL2^4/2;

Nv2p = -pi*(T2/WLL2)^2*(1/2+2.4*T2/WLL2);
Nvdot2p = -pi*(T2/WLL2)^2*(1.1*B2/WLL2-0.041*B2/T2);

Nv2 = Nv2p*rho*WLL2^3*u02/2;
Nvdot2 = Nvdot2p*rho*WLL2^4/2;

Nr2p = -pi*(T2/WLL2)^2*(1/4+0.039*B2/T2-0.56*B2/WLL2);
Nrdot2p = -pi*(T2/WLL2)^2*(1/12+0.017*CB2*B2/T2-0.33*B2/WLL2);

Nr2 = Nr2p*rho*WLL2^4*u02/2;
Nrdot2 = Nrdot2p*rho*WLL2^5/2;

alph2 = 3; %Rudder factor from Clarke (1982)
Ydelta2p = A_rudder2/WLL2^2*alph2;
Ndelta2p = -0.5*Ydelta2p;

Ydelta2 = Ydelta2p*rho*WLL2^2*u02^2/2;
Ndelta2 = Ndelta2p*rho*WLL2^3*u02^2/2;

%Total Water Resistance RT = RF(1+k1)+RAPP+RW+RTR+RA

CM2 = MS2/(B2*T2); %Midship section area coefficient
CWP2 = WP2/(B2*WLL2); %Waterplane area coefficient


CP2 = SV2/(WP2*WLL2); %Prismatic coefficient
CLCB2 = LCB2/WLL2*100; %Buoyancy coefficient

WA2 = WLL2*(2*T2+B2)*sqrt(CM2)*(0.453+0.4425*CB2-0.2862*CM2-0.003467*B2/T2+0.3696*CWP2); %Wetted area approximation [m^2]

%Find k2
LR2 = (1-CP2+0.06*CP2*CLCB2/(4*CP2-1))*WLL2; 

if T2/WLL2 > 0.05
    c122 = (T2/WLL2)^0.2228446;
elseif T2/WLL2 >0.02
    c122 = 48.2*(T2/WLL2-0.02)^2.078 +0.479948;
else
    c122 = 0.479948;
end

k2 = 0.93+c122*(B2/LR2)^0.92497*(0.95-CP2)^(-0.521448)*(1-CP2+0.0225*CLCB2)^0.6906 -1;

%RF and RAPP depend on speed

%RW, wave resistance
if B2/WLL2 < 0.11
    c72 = 0.229577*(B2/WLL2)^0.33333;
elseif B2/WLL2 < 0.25
    c72=B2/WLL2;
else
    c72 = 0.5-0.0625*WLL2/B2;
end

c52 = 1-0.8*AT2/(B2*T2*CM2);

if WLL2/B2<12
    lambda2 = 1.446*CP2 - 0.03*WLL2/B2;
else
    lambda2 = 1.446*CP2 - 0.36;
end

if CP2<0.8
    c162 = 8.07981*CP2 - 13.8673*CP2^2 + 6.984388*CP2^3;
else
    c162 = 1.73014 - 0.7067*CP2;
end

m12 = 0.0140407*WLL2/T2 - 1.75254*SV2^(1/3)/WLL2 - 4.79323*B2/WLL2 - c162;

if WLL2^3/SV2 < 512
    c152 = -1.69385;
elseif WLL2^3/SV2 < 1727
    c152 = -1.69385 + (WLL2/SV2^(1/3)-8)/2.36;
else
    c152 = 0;
end

d2=-0.9;

iE2 = 1+89*exp(-(WLL2/B2)^0.80856*(1-CWP2)^0.30484*(1-CP2-0.0225*CLCB2)^0.6367*(LR2/B2)^0.34574*(100*SV2/WLL2^3)^0.16302);
c12 = 2223105*c72^3.78613*(T2/B2)^1.07961*(90-iE2)^(-1.37565);

%Wave Resistance depends on speed and is calculated in loop

%RTR, Resistance from transom stern calculated in loop

%RA, model ship correlation resistance

if T2/WLL2 < 0.04
    c42 = T2/WLL2;
else
    c42 = 0.04;
end

CA2 = 0.006*(WLL2+100)^(-0.16) - 0.00205 + 0.003*sqrt(WLL2/7.5)*CB2^4*(0.04-c42);

%Model Ship Correlation Resistance depends on speed

%Using First Order Nomoto model by assuming fore/aft symmetry for the
%vessel - valid assumption for this vessel

%K2 = -(u02/WLL2)*Ndelta2p/Nr2p;
%T2 = -(WLL2/u02)*(IPz2p-Nrdot2p)/Nr2p;

%Yr2p = 0;
%Yrdot2p = 0;
%Nv2p = 0;
%Nvdot2p = 0;

%Initial Conditions
compliant = false; %Select compliant(true) or non-compliant(false) 

collPoint = 2500; %closest point of approach, select point on x axis

initheading1 = 0*pi/180; %Initial Heading of vessel 1 in rad
u1i = 1*u01; %Initial Speed of vessel 1 in m/s
delta1 = 0*pi/180; %Initial Rudder angle in rad
initRad1 = 2500; %arbitrary selection
x1init = -initRad1*cos(initheading1)+collPoint; %Initial Position controlled vessel
y1init = -initRad1*sin(initheading1);

%COLREG boundaries
overtakingmax = 22.5*pi/180;
overtakingmin = 337.5*pi/180;
headonmin = 157.5*pi/180; %157.5
headonmax = 202.5*pi/180; %202.5
crossing1max = 157.5*pi/180; %Intruder turns
crossing1min = 22.5*pi/180;
crossing2min = 202.5*pi/180; %Controlled vessel turns
crossing2max = 337.5*pi/180;

%Simulation Parameters
timestep = 0.1; %Sim timestep in s

%Environmental factors
VcurrMag = 0;%0.8267;
VcurrD = 0;%1.7676; %current direction in rad
Vcurrx = VcurrMag*cos(VcurrD);
Vcurry = VcurrMag*sin(VcurrD);
VwindMag = 0;%19.4542;
VwindD = 0;%4.2134; %wind direction in rad
Vwindx = VwindMag*cos(VwindD);
Vwindy = VwindMag*sin(VwindD);

emergzone1 = WLL1*2; %When CAP is within this boundary, turn at designated TCAP threshold

udot1 = 0;

rdot1 = 0;
r1(1) = 0;

vdot1 = 0;
v1 = 0;

x1 = zeros(1,1);
y1 = zeros(1,1);

udot2 = 0;

rdot2 = 0;
r2 = 0;

vdot2 = 0;
v2 = 0;

x2 = zeros(1,1);
y2 = zeros(1,1);

%Controller Params

%Surge Controller
kps1 = 1000;
kis1 = 50;
kds1 = 1;
es1 = 0;

kps2 = 1000;
kis2 = 50;
kds2 = 1;
es2 = 0;

%Rudder Controller
kpr1 = 0.0000001;
kir1 = 0.000001;
kdr1 = 0.01;
ep1 = 0;

kpr2 = 0.0000001;
kir2 = 0.000001;
kdr2 = 0.01;
ep2 = 0;

%For sweep analysis
startit = 28;
endit = 157;
stepit = 0.5;
iter = (endit-startit)/stepit;
TCAP = zeros(1,iter);
m = 1;

%for j = startit:stepit:endit
    
    initheading2 = 40*pi/180; %Heading of vessel 2 in rad
    u2i = 1*u02; %Initial speed of vessel 2 in m/s
    delta2 = 0*pi/180;
    initRad2 = (initRad1/u1i)*u2i;
    x2init = -initRad2*cos(initheading2)+collPoint; %Initial Position intruder
    y2init = -initRad2*sin(initheading2);
    
    %Path gen
    res = 50; %meters between pts

    path1 = zeros(2,collPoint/res);
    path2 = zeros(2,collPoint/res);

    for i = 1:1:2*collPoint/res
        path1(:,i) = [i*res;0];
        path2(:,i) = [x2init+res*i*cos(initheading2);y2init+res*i*sin(initheading2)];
    end
    
    TCAPturn1 = 172; %When TCAP goes below this number, turn sharply
    safeturn1 = false;
    tcap = 0;
    
    while safeturn1 == false
        
        phi_err1 = 0;
        phi_err2 = 0;
        phi_e1 = 0;
        phi_e2 = 0;
        sur_err1 = 0;
        sur_err2 = 0;
        sur_e1 = 0;
        sur_e2 = 0;

        flag = 0;
        z1 = 1; %Path index
        z2 = 1;
        distpath1 = 0;
        distpath2 = 0;
        distboat1 = 0;
        distboat2 = 0;
        relanglepath1 = 0;
        relanglepath2 = 0;

        checkclear1 = false;
        checkencount1 = false;
        i=2;
        x1 = zeros(1,1);
        y1 = zeros(1,1);
        x1(1) = x1init;
        y1(1) = y1init;
        delta1 = 0;
        rdot1 = 0;
        r1(1) = 0;
        udot1 = 0;
        u1 = u1i;
        v1 = 0;
        vdot1 = 0;
        heading1 = initheading1;

        x2 = zeros(1,1);
        y2 = zeros(1,1);
        x2(1) = x2init;
        y2(1) = y2init;
        delta2 = 0;
        rdot2 = 0;
        r2 = 0;
        udot2 = 0;
        vdot2 = 0;
        u2 = u2i;
        v2 = 0;
        heading2 = initheading2;

        relheading = heading1 - heading2;

        if relheading < 0
            relheading = relheading + 2*pi;
        end

        q=0; %index for trailing error measure
        sur_e1 = zeros(1,50);
        sur_e2 = zeros(1,50);

        while checkclear1 == false
            %rdot1 = (K*delta1-r1)/T; %Nomoto Model
            %r1 = r1 + rdot1*timestep;
            
            q = q+1;
            if q>50
                q=1;
            end
            
            if abs(r1(i-1))>0.0043633
                r1(i-1) = sign(r1(i-1))*0.0043633;
            end

            if abs(r2)>0.0043633
                r2 = sign(r2)*0.0043633;
            end

            RE1 = u1.*WLL1 / nu;
            CF1 = 0.075 ./ (log10(RE1)-2).^2;
            RF1 = 0.5 * rho * u1.^2*WA1.*CF1; %Friction Resistance force

            RAPP1 = 0.5 * rho * u1.^2 * k2Coef1 .* WAA1 * CF1; %Resistance due to appendages

            Fn1 = u1 ./ sqrt(g*WLL1);
            m21 = c151*CP1^2*exp(-0.1*Fn1.^(-2));
            RW1 = c11*c51*SV1*rho*g*exp(m11*Fn1.^d1 + m21.*cos(lambda1*Fn1.^(-2)) ); %Wave resistance

            FnT1 = u1./ sqrt(2*g*AT1/(B1+B1*CWP1));
            c61 = 0.2*(1-0.2*FnT1) .* (FnT1<5);
            RTR1 = 0.5 * rho .* u1.^2 * AT1 .* c61; %Transom stern resistance

            RA1 = 0.5*rho.*u1.^2*WAA1*CA1; %Correlation resistance

            RT = RF1+RAPP1+RW1+RTR1+RA1; %Total Surge resistance force
            if flag == 0
                Fp = RT;
            end

            CCA1 = RT/(0.5*rho*WSA1_ax*u1^2); %Drag Coefficient from current in axial direction

            DCA1 = 0.5*rho*WSA1_ax*CCA1;
            DWA1 = 0.5*1.2041*Awind1_ax*CWA1; %rho_air = 1.2041

            Vcurri1 = Vcurrx*cos(heading1)+Vcurry*sin(heading1);
            Vcurrj1 = Vcurrx*-sin(heading1)+Vcurry*cos(heading1);
            Vwindi1 = Vwindx*cos(heading1)+Vwindy*sin(heading1);
            Vwindj1 = Vwindx*-sin(heading1)+Vwindy*cos(heading1);

            Xwind1 = DWA1*(Vwindi1-u1)*abs(Vwindi1-u1);
            Xcurr1 = DCA1*(Vcurri1-u1)*abs(Vcurri1-u1);

            Ywind1 = zeros(1,n1);
            Ycurr1 = zeros(1,n1);

            for k=1:n1
                Ywind1(k) = DWseg1*(Vwindj1-v1-r1(i-1)*ncoord1(k))*abs(Vwindj1-v1-r1(i-1)*ncoord1(k));
                Ycurr1(k) = DCseg1*(Vcurrj1-v1-r1(i-1)*ncoord1(k))*abs(Vwindj1-v1-r1(i-1)*ncoord1(k));
            end

            Ywind1tot = sum(Ywind1);
            Ycurr1tot = sum(Ycurr1);

            Nwind1 = ncoord1*Ywind1';
            Ncurr1 = ncoord1*Ycurr1';

            Xnet1 = Fp+Xwind1+Xcurr1-RT;
            udot1 = (Fp+Xwind1+Xcurr1-RT)/(m1-Xudot1);
            u1 = u1+udot1*timestep;
            
            sur_err1 = u1i - u1;
            sur_e1(q) = sur_err1;
            sur_e1avg = sum(sur_e1)/50;
            input = kps1 * sur_err1 + kis1*sur_e1avg*timestep + kds1 * (0 - udot1);
            Fp = Fp + input;

            Fv1 = v1*Yv1+r1(i-1)*Yr1+rdot1*Yrdot1+delta1*Ydelta1+Ywind1tot+Ycurr1tot;
            vdot1 = Fv1/(m1-Yvdot1);
            v1 = v1 + vdot1*timestep;

            NH1 = v1*Nv1+vdot1*Nvdot1+r1(i-1)*Nr1;
            NR1 = -WLL1*delta1*Yr1/2;
            Ntot1 = NH1+NR1+Nwind1+Ncurr1;
            rdot1 = Ntot1/(IPz1-Nrdot1);
            r1(i+1) = r1(i-1) + rdot1*timestep;

            heading1 = heading1 + r1(i-1)*timestep;

            RE2 = u2.*WLL2 / nu;
            CF2 = 0.075 ./ (log10(RE2)-2).^2;
            RF2 = 0.5 * rho * u2.^2*WA2.*CF2; %Friction Resistance force

            RAPP2 = 0.5 * rho * u2.^2 * k2Coef2 .* WAA2 * CF2; %Resistance due to appendages

            Fn2 = u2 ./ sqrt(g*WLL2);
            m22 = c152*CP2^2*exp(-0.1*Fn2.^(-2));
            RW2 = c12*c52*SV2*rho*g*exp(m12*Fn2.^d2 + m22.*cos(lambda2*Fn2.^(-2)) ); %Wave resistance

            FnT2 = u2./ sqrt(2*g*AT2/(B2+B2*CWP2));
            c62 = 0.2*(1-0.2*FnT2) .* (FnT2<5);
            RTR2 = 0.5 * rho .* u2.^2 * AT2 .* c62; %Transom stern resistance

            RA2 = 0.5*rho.*u2.^2*WAA2*CA2; %Correlation resistance

            RT2 = RF2+RAPP2+RW2+RTR2+RA2; %Total Surge resistance force
            if flag == 0
                Fp2 = RT2;
            end
            CCA2 = RT2/(0.5*rho*WSA2_ax*u2^2);

            DCA2 = 0.5*rho*WSA2_ax*CCA2;
            DWA2 = 0.5*1.2041*Awind2_ax; %rho_air = 1.2041

            Vcurri2 = Vcurrx*cos(heading2)+Vcurry*sin(heading2);
            Vcurrj2 = Vcurrx*-sin(heading2)+Vcurry*cos(heading2);
            Vwindi2 = Vwindx*cos(heading2)+Vwindy*sin(heading2);
            Vwindj2 = Vwindx*-sin(heading2)+Vwindy*cos(heading2);

            Xwind2 = DWA2*(Vwindi2-u2)*abs(Vwindi2-u2);
            Xcurr2 = DCA2*(Vcurri2-u2)*abs(Vcurri2-u2);

            Ywind2 = zeros(1,n2);
            Ycurr2 = zeros(1,n2);

            for k=1:n2
                Ywind2(k) = DWseg2*(Vwindj2-v2-r2*ncoord2(k))*abs(Vwindj2-v2-r2*ncoord2(k));
                Ycurr2(k) = DCseg2*(Vcurrj2-v2-r2*ncoord2(k))*abs(Vwindj2-v2-r2*ncoord2(k));
            end

            Ywind2tot = sum(Ywind2);
            Ycurr2tot = sum(Ycurr2);

            Nwind2 = ncoord2*Ywind2';
            Ncurr2 = ncoord2*Ycurr2';      

            Xnet2 = Fp2+Xwind2+Xcurr2-RT2;    
            udot2 = (Fp2+Xwind2+Xcurr2-RT2)/(m2-Xudot2);
            u2 = u2+udot2*timestep;

            sur_err2 = u2i - u2;
            sur_e2(q) = sur_err2;
            sur_e2avg = sum(sur_e2)/50;
            input = kps2 * sur_err2 + kis2*sur_e2avg*timestep + kds2 * (0 - udot2);
            Fp2 = Fp2 + input;

            Fv2 = v2*Yv2+r2*Yr2+rdot2*Yrdot2+delta2*Ydelta2+Ywind2tot+Ycurr2tot;
            vdot2 = Fv2/(m2-Yvdot2);
            v2 = v2 + vdot2*timestep;

            NH2 = v2*Nv2+vdot2*Nvdot2+r2*Nr2;

            if checkencount1 == 1
                x=0;%do nothing
            end

            NR2 = -WLL2*delta2*Yr2/2;
            Ntot2 = NH2+NR2+Nwind2+Ncurr2;
            rdot2 = Ntot2/(IPz2-Nrdot2);
            r2 = r2 + rdot2*timestep;

            heading2 = heading2 + r2*timestep;

            x1(i) = x1(i-1)+u1*cos(heading1)*timestep+v1*sin(heading1)*timestep;
            y1(i) = y1(i-1)+u1*sin(heading1)*timestep+v1*cos(heading1)*timestep;

            x2(i) = x2(i-1)+u2*cos(heading2)*timestep+v2*sin(heading2)*timestep;
            y2(i) = y2(i-1)+u2*sin(heading2)*timestep+v2*cos(heading2)*timestep;

            vrelx = u1*cos(heading1)+v1*sin(heading1)-u2*cos(heading2)-v2*sin(heading2);
            vrely = u1*sin(heading1)+v1*cos(heading1)-u2*sin(heading2)-v2*cos(heading2);
            srelx = x1(i)-x2(i);
            srely = y1(i)-y2(i);

            vmag = sqrt(vrelx^2+vrely^2);
            drel = sqrt(srelx^2+srely^2);

            tcap = -(srelx*vrelx+srely*vrely)/vmag^2;
            dcap = sqrt((srelx+tcap*vrelx)^2+(srely+tcap*vrely)^2);
            cap = sqrt((srelx+vrelx*tcap)^2+(srely+vrely*tcap)^2);

            if cap<emergzone1 && tcap<TCAPturn1 && compliant==false
                checkencount1 = true;

                if relheading < overtakingmax || relheading >= 0
                    delta1 = Delta_max1;
                end

                if relheading > crossing1min && relheading < crossing1max
                    delta1 = -Delta_max1;
                end

                if relheading > headonmin && relheading < headonmax
                    delta1 = Delta_max1;
                end

                if relheading > crossing2min && relheading < crossing2max
                    delta1 = Delta_max1;
                end
                
                distpath2 = sqrt(((path2(1,z2)-collPoint)^2)+(path2(2,z2)^2));
                distboat2 = sqrt(((x2(i)-collPoint)^2)+(y2(i)^2));

                if distboat2<distpath2
                    z2 = z2+1;
                end

                relanglepath2 = atan((path2(2,z2)-y2(i))/(path2(1,z2)-x2(i)));
                
                if heading2>2*pi
                    heading2 = heading2 - 2*pi;
                end

                if heading2<-2*pi
                    heading2 = heading2 +2*pi;
                end
                
                phi_err2 = relanglepath2-heading2+initheading2;
                phi_e2(q) = phi_err2;
                phi_e2avg = sum(phi_e2)/50;
                
                input = kpr2 * phi_err2 + kir2*phi_e2avg*timestep + kdr2 * (0 - r2);
                delta2 = delta2 - input;
            end

            if cap<emergzone1 && tcap<TCAPturn1 && compliant==true
                checkencount1 = true;

                if relheading < overtakingmax || relheading > overtakingmin
                    delta1 = Delta_max1;
                end

                if relheading > crossing1min && relheading < crossing1max
                    delta2 = Delta_max2;
                end

                if relheading > headonmin && relheading < headonmax
                    delta1 = Delta_max1;
                    delta2 = Delta_max2;
                end

                if relheading > crossing2min && relheading < crossing2max
                    delta1 = Delta_max1;
                end        
            end

            if drel<emergzone1
                TCAPturn1 = TCAPturn1+1;
                safeturn1 = false;
                checkclear1 = true;
            end


            if checkencount1 == true && tcap<0
                checkclear1 = true;
                safeturn1 = true;
            end

            if checkencount1 == false
                distpath1 = sqrt(((path1(1,z1)-collPoint)^2)+(path1(2,z1)^2));
                distboat1 = sqrt(((x1(i)-collPoint)^2)+(y1(i)^2));
                distpath2 = sqrt(((path2(1,z2)-collPoint)^2)+(path2(2,z2)^2));
                distboat2 = sqrt(((x2(i)-collPoint)^2)+(y2(i)^2));

                if distboat1<distpath1
                    z1 = z1+1;
                end

                if distboat2<distpath2
                    z2 = z2+1;
                end

                relanglepath1 = atan((path1(2,z1)-y1(i))/(path1(1,z1)-x1(i)));
                relanglepath2 = atan((path2(2,z2)-y2(i))/(path2(1,z2)-x2(i)));

                if heading1>2*pi
                    heading1 = heading1 - 2*pi;
                end

                if heading1<-2*pi
                    heading1 = heading1 + 2*pi;
                end

                if heading2>2*pi
                    heading2 = heading2 - 2*pi;
                end

                if heading2<-2*pi
                    heading2 = heading2 +2*pi;
                end


                phi_err1 = relanglepath1-heading1;
                phi_e1(q) = phi_err1;
                phi_e1avg = sum(phi_e1)/50;

                %Need to figure out whats going on here - need this block for
                %head-on, need to remove it for crossing

%                 if heading2<-pi/2
%                     phi_err2 = relanglepath2-heading2-pi;
%                 elseif heading2>pi/2
%                     phi_err2 = relanglepath2-heading2+pi;
%                 else
%                     phi_err2 = relanglepath2-heading2;
%                 end             

                phi_err2 = relanglepath2-heading2;
                phi_e2(q) = phi_err2;
                phi_e2avg = sum(phi_e2)/50;

                input = kpr1 * phi_err1 + kir1*phi_e1avg*timestep + kdr1 * (0 - r1(i-1));
                delta1 = delta1 - input;
                input = kpr2 * phi_err2 + kir2*phi_e2avg*timestep + kdr2 * (0 - r2);
                delta2 = delta2 - input;
            end

            i = i+1;
            flag = 1;
            
%             if i==1600
%                 figure()
%                 plot(x1,y1,'b','LineWidth',5)
%                 hold on
%                 %plot(x1(i-1),y1(i-1),'x','LineWidth',5)
%                 %hold on
%                 plot(x2,y2,'r','LineWidth',5)
%                 hold on
%                 plot(x2(i-1),y2(i-1),'x','LineWidth',5)
%                 hold on
%                 th = 0:pi/16:2*pi;
%                 xcirc = emergzone1*cos(th)+x1(i-1);
%                 ycirc = emergzone1*sin(th)+y1(i-1);
%                 plot(xcirc,ycirc,'g','LineWidth',5)
%                 hold on
%                 rect1X = [-WLL1/2 WLL1/2 WLL1/2 -WLL1/2];
%                 rect1Y = [-B1/2 -B1/2 B1/2 B1/2];
%                 rect1Xrot = rect1X*cos(heading1) - rect1Y*sin(heading1);
%                 rect1Yrot = rect1X*sin(heading1) + rect1Y*cos(heading1);
%                 plot(rect1Xrot + x1(i-1), rect1Yrot+ y1(i-1),'g','LineWidth',3);
%                 xlabel("Distance [m]")
%                 ylabel("Distance [m]")
%                 %xlim([0 1000])
%                 %ylim([-1000 500])
%                 axis equal
%             end
%             
%             if i==4000
%                 figure()
%                 plot(x1,y1,'b','LineWidth',5)
%                 hold on
%                 %plot(x1(i-1),y1(i-1),'x','LineWidth',5)
%                 %hold on
%                 plot(x2,y2,'r','LineWidth',5)
%                 hold on
%                 plot(x2(i-1),y2(i-1),'x','LineWidth',5)
%                 hold on
%                 th = 0:pi/16:2*pi;
%                 xcirc = emergzone1*cos(th)+x1(i-1);
%                 ycirc = emergzone1*sin(th)+y1(i-1);
%                 plot(xcirc,ycirc,'g','LineWidth',5)
%                 hold on
%                 rect1X = [-WLL1/2 WLL1/2 WLL1/2 -WLL1/2];
%                 rect1Y = [-B1/2 -B1/2 B1/2 B1/2];
%                 rect1Xrot = rect1X*cos(heading1) - rect1Y*sin(heading1);
%                 rect1Yrot = rect1X*sin(heading1) + rect1Y*cos(heading1);
%                 plot(rect1Xrot + x1(i-1), rect1Yrot+ y1(i-1),'g','LineWidth',3);
%                 %xlim([0 1000])
%                 %ylim([-1000 500])
%                 xlabel("Distance [m]")
%                 ylabel("Distance [m]")
%                 axis equal
%             end
            
        end

    end
    TCAP(m) = TCAPturn1;
    m = m+1;
    
%end 

figure()
plot(x1,y1,'b','LineWidth',5)
hold on
%plot(x1(i-1),y1(i-1),'x','LineWidth',5)
%hold on
plot(x2,y2,'r','LineWidth',5)
hold on
plot(x2(i-1),y2(i-1),'x','LineWidth',5)
hold on
th = 0:pi/16:2*pi;
xcirc = emergzone1*cos(th)+x1(i-1);
ycirc = emergzone1*sin(th)+y1(i-1);
plot(xcirc,ycirc,'g','LineWidth',5)
hold on
rect1X = [-WLL1/2 WLL1/2 WLL1/2 -WLL1/2];
rect1Y = [-B1/2 -B1/2 B1/2 B1/2];
rect1Xrot = rect1X*cos(heading1) - rect1Y*sin(heading1);
rect1Yrot = rect1X*sin(heading1) + rect1Y*cos(heading1);
plot(rect1Xrot + x1(i-1), rect1Yrot+ y1(i-1),'g','LineWidth',3);
%xlim([0 1000])
%ylim([-1000 500])
axis equal


% TCAPPoints = zeros(2,iter+1);
% angleNorm1 = linspace((startit+180)*pi/180,(endit+180)*pi/180,iter+1);
% polarplot(angleNorm1,TCAP)
% %rlim([170 190])
% %thetalim([-25 25])





