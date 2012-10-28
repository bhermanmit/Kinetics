% Calculation of Coefficients

% Here we generate coefficient from Shampine --> Kaps form --> Num Rec form

gam = 0.5;
a_s21 = 1;
a_s31 = 24/25;
a_s32 = 3/25;
a_s41 = 24/25;
a_s42 = 3/25;
a_s43 = 0;

c_s21 = -4;
c_s31 = 186/25;
c_s32 = 6/5;
c_s41 = -56/125;
c_s42 = -27/125;
c_s43 = -1/5;

m_s1 = 19/18;
m_s2 = 1/4;
m_s3 = 25/216;
m_s4 = 125/216;

mp_s1 = 97/108;
mp_s2 = 11/72;
mp_s3 = 25/216;

b_s1 = gam;
b_s2 = gam + c_s21*b_s1;
b_s3 = gam + c_s31*b_s1 + c_s32*b_s2;
b_s4 = gam + c_s41*b_s1 + c_s42*b_s2 + c_s43*b_s3;

d_s1 = 0;
d_s2 = a_s21*b_s1/gam;
d_s3 = a_s31*b_s1/gam + a_s32*b_s2/gam;
d_s4 = a_s41*b_s1/gam + a_s42*b_s2/gam + a_s43*b_s3/gam;

e_s1 = m_s1 - mp_s1;
e_s2 = m_s2 - mp_s2;
e_s3 = m_s3 - mp_s3;
e_s4 = m_s4;

%% Compute all Kaps Rentrop Form
a_k21 = a_s21;
a_k31 = a_s31 + a_s32*c_s21;
a_k32 = a_s32;
a_k41 = a_s41 + a_s42*c_s21 + a_s43*c_s31 + a_s43*c_s32*c_s21;
a_k42 = a_s42 + a_s43*c_s32;
a_k43 = a_s43;

c_k21 = gam*c_s21;
c_k31 = gam*(c_s31 + c_s32*c_s21);
c_k32 = gam*c_s32;
c_k41 = gam*(c_s41 + c_s42*c_s21 + c_s43*c_s32*c_s21 + c_s43*c_s31);
c_k42 = gam*(c_s42 + c_s43*c_s32);
c_k43 = gam*c_s43;

m_k1 = m_s1 + m_s2*c_s21 + m_s3*c_s31 + m_s3*c_s32*c_s21 + m_s4*c_s41 + m_s4*c_s42*c_s21 + m_s4*c_s43*c_s31 + m_s4*c_s43*c_s32*c_s21;
m_k2 = m_s2 + m_s3*c_s32 + m_s4*c_s42 + m_s4*c_s43*c_s32;
m_k3 = m_s3 + m_s4*c_s43;
m_k4 = m_s4;

mp_k1 = mp_s1 + mp_s2*c_s21 + mp_s3*c_s31 + mp_s3*c_s32*c_s21;
mp_k2 = mp_s2 + mp_s3*c_s32;
mp_k3 = mp_s3;

e_k1 = m_k1 - mp_k1;
e_k2 = m_k2 - mp_k2;
e_k3 = m_k3 - mp_k3;
e_k4 = m_k4;

%% Compute all Numerical Recipes Form

a_n21 = a_k21/gam;
a_n31 = a_k31/gam - a_k32*c_k21/gam^2;
a_n32 = a_k32/gam;
a_n41 = a_k41/gam - a_k42*c_k21/gam^2 - a_k43*c_k31/gam^2 + a_k43*c_k32*c_k21/gam^3; 
a_n42 = a_k42/gam - a_k43*c_k32/gam^2;
a_n43 = a_k43/gam;

c_n21 = c_k21/gam^2;
c_n31 = c_k31/gam^2 - c_k32*c_k21/gam^3;
c_n32 = c_k32/gam^2;
c_n41 = c_k41/gam^2 - c_k42*c_k21/gam^3 - c_k43*c_k31/gam^3 + c_k43*c_k32*c_k21/gam^4;
c_n42 = c_k42/gam^2 - c_k43*c_k32/gam^3;
c_n43 = c_k43/gam^2;

m_n1 = m_k1/gam - m_k2*c_k21/gam^2 - m_k3*c_k31/gam^2 + m_k3*c_k32*c_k21/gam^3 - m_k4*c_k41/gam^2 + m_k4*c_k42*c_k21/gam^3 + m_k4*c_k43*c_k31/gam^3 - m_k4*c_k43*c_k32*c_k21/gam^4;
m_n2 = m_k2/gam - m_k3*c_k32/gam^2 - m_k4*c_k42/gam^2 + m_k4*c_k43*c_k32/gam^3;
m_n3 = m_k3/gam - m_k4*c_k43/gam^2;
m_n4 = m_k4/gam;

mp_n1 = mp_k1/gam - mp_k2*c_k21/gam^2 - mp_k3*c_k31/gam^2 + mp_k3*c_k32*c_k21/gam^3;
mp_n2 = mp_k2/gam - mp_k3*c_k32/gam^2;
mp_n3 = mp_k3/gam;

e_n1 = m_n1 - mp_n1;
e_n2 = m_n2 - mp_n2;
e_n3 = m_n3 - mp_n3;
e_n4 = m_n4;