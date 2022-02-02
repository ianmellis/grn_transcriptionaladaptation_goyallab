function dydt = odefun1(t,y,rates)
r_onbasal_A1 = rates(13);
r_off_A1 = rates(22);
r_prodon_A1 = rates(5);
r_deg_A1 = rates(9);
r_bind_byA1_B1 = rates(26);
n_A1 = rates(32);
k_A1 = rates(28);
r_unbind_byA1_B1 = rates(36);
r_bound_byA1_B1 = rates(20);
r_off_B1 = rates(25);
r_prodon_B1 = rates(8);
r_deg_B1 = rates(12);

dydt = zeros(8,1);
dydt(1) = r_onbasal_A1*y(2) - r_off_A1*y(1);
dydt(2) = -r_onbasal_A1*y(2) + r_off_A1*y(1);
dydt(3) = r_prodon_A1*y(1) - r_deg_A1*y(3);
dydt(4) = (r_bind_byA1_B1*y(3)^n_A1/(k_A1^n_A1 + y(3)^n_A1))*y(5) - r_unbind_byA1_B1*y(4);
dydt(5) = -(r_bind_byA1_B1*y(3)^n_A1/(k_A1^n_A1 + y(3)^n_A1))*y(5) + r_unbind_byA1_B1*y(4);
dydt(6) = r_bound_byA1_B1*y(7) - r_off_B1*y(6);
dydt(7) = -r_bound_byA1_B1*y(7) + r_off_B1*y(6);
dydt(8) = r_prodon_B1*y(6) - r_deg_B1*y(8);
end