% Conversion to numerically efficient form
% Bryan Herman

a_p21 = a_21/gam;
a_p31 = a_31/gam - a_32*c_21/gam^2;
a_p32 = a_32/gam;
a_p41 = a_41/gam - a_42*c_21/gam^2 - a_43*c_31/gam^2 + a_43*c_32*c_21/gam^3; 
a_p42 = a_42/gam - a_43*c_32/gam^2;
a_p43 = a_43/gam;

c_p21 = c_21/gam^2;
c_p31 = c_31/gam^2 - c_32*c_21/gam^3;
c_p32 = c_32/gam^2;
c_p41 = c_41/gam^2 - c_42*c_21/gam^3 - c_43*c_31/gam^3 + c_43*c_32*c_21/gam^4;
c_p42 = c_42/gam^2 - c_43*c_32/gam^3;
c_p43 = c_43/gam^2;

m_p1 = m_1/gam - m_2*c_21/gam^2 - m_3*c_31/gam^2 + m_3*c_32*c_21/gam^3 - m_4*c_41/gam^2 + m_4*c_42*c_21/gam^3 + m_4*c_43*c_31/gam^3 - m_4*c_43*c_32*c_21/gam^4;
m_p2 = m_2/gam - m_3*c_32/gam^2 - m_4*c_42/gam^2 + m_4*c_43*c_32/gam^3;
m_p3 = m_3/gam - m_4*c_43/gam^2;
m_p4 = m_4/gam;

mh_p1 = mh_1/gam - mh_2*c_21/gam^2 - mh_3*c_31/gam^2 + mh_3*c_32*c_21/gam^3;
mh_p2 = mh_2/gam - mh_3*c_32/gam^2;
mh_p3 = mh_3/gam;

e_p1 = m_p1 - mh_p1;
e_p2 = m_p2 - mh_p2;
e_p3 = m_p3 - mh_p3;
e_p4 = m_p4;