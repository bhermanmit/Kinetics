% Solve simulataneous order equations
% Bryan Herman

% note this script already takes into accout
% that a_43 = 0, a_41 = a_31 and a_42 = a_32

% chose gamma for stability
% 0.231 for Kaps-Rentrop, 0.5 for Shampine
gam = 0.5;

% comPute used Polynomials
P1 = 1;
P2 = 1/2 - gam;
P3 = 1/3;
P4 = 1/6 - gam + gam^2;
P5 = 1/4;
P6 = 1/8 - 1/3*gam;
P7 = 1/12 - 1/3*gam;
P8 = 1/24 - 1/2*gam + 3/2*gam^2 - gam^3;
P9 = 1/5;
P14 = 1/20 - 1/4*gam;

% calculate free Parameters
a_2 = 2*gam;
a_3 = (P9-P5*a_2)/(P5-P3*a_2);
m_3 = 0;

% Eq. B.17
B17 = [1 1 1; 0 a_2^2 a_3^2; 0 a_2^3 a_3^3]\[P1; P3; P5];
m_1 = B17(1);
m_2 = B17(2);
m_4 = B17(3);

% Eq. B.18
B18 = [m_4*a_2^2 m_4*a_3^2; m_4*a_2^3 m_4*a_3^3]\[P7; P14];
beta_42 = B18(1);
beta_43 = B18(2);

% Eq. B.19
u = P8/(m_4*beta_43);

% Eq. B.20
B20 = [1 1 1; 0 a_2^2 a_3^2; 0 0 u]\[P1; P3; P4];
mh_1 = B20(1);
mh_2 = B20(2);
mh_3 = B20(3);

% Eq. B.21
B21 = [m_4*beta_42 m_4*beta_43; mh_2 mh_3]\[P4; P2];
beta_2 = B21(1);
beta_3 = B21(2);

% Eq. B.22
beta_32 = u/beta_2;

% Eq. B.23
a_32 = P6/(m_3*a_3*beta_2 + m_4*a_3*beta_2);

% Eq. B.24
beta_4 = (P2 - m_2*beta_2 - m_3*beta_3)/m_4;

% solve for rest of Parameters
a_21 = a_2;
a_31 = a_3 - a_32;
beta_21 = beta_2;
beta_31 = beta_3 - beta_32;
beta_41 = beta_4 - beta_43 - beta_42;
c_21 = beta_21 - a_21;
c_31 = beta_31 - a_31;
c_32 = beta_32 - a_32;
c_41 = beta_41 - a_31;
c_42 = beta_42 - a_32;
c_43 = beta_43;

% Derived coefficients
d_1 = 0;
d_2 = a_21;
d_3 = a_31 + a_32;
d_4 = a_31 + a_32;
b_1 = gam;
b_2 = c_21 + gam;
b_3 = c_31 + c_32 + gam;
b_4 = c_41 + c_42 + c_43 + gam;

% Note that
a_41 = a_31;
a_42 = a_32;
a_43 = 0;