function dydt = odefun_wtwt(t,y,rates)
% diploid, wt/wt ODE, manually excludes paralog and nonsense

% vector: y: 
% y(1) = ona_1
% y(2) = offa_1
% y(3) = ona_2
% y(4) = offa_2
% y(5) = a 
% y(6) = bba_1
% y(7) = unb_1
% y(8) = bba_2
% y(9) = unb_2
% y(10) = onb_1
% y(11) = offb_1
% y(12) = onb_2
% y(13) = offb_2
% y(14) = b

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

dydt = zeros(14,1);
dydt(1) = r_onbasal_A1*y(2) - r_off_A1*y(1);
dydt(2) = -r_onbasal_A1*y(2) + r_off_A1*y(1);
dydt(3) = r_onbasal_A1*y(4) - r_off_A1*y(3);
dydt(4) = -r_onbasal_A1*y(4) + r_off_A1*y(3);
dydt(5) = r_prodon_A1*y(1) + r_prodon_A1*y(3) - r_deg_A1*y(5);
dydt(6) = (r_bind_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1))*y(7) - r_unbind_byA1_B1*y(6);
dydt(7) = -(r_bind_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1))*y(7) + r_unbind_byA1_B1*y(6);
dydt(8) = (r_bind_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1))*y(9) - r_unbind_byA1_B1*y(8);
dydt(9) = -(r_bind_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1))*y(9) + r_unbind_byA1_B1*y(8);
dydt(10) = r_bound_byA1_B1*y(11) - r_off_B1*y(10);
dydt(11) = -r_bound_byA1_B1*y(11) + r_off_B1*y(10);
dydt(12) = r_bound_byA1_B1*y(13) - r_off_B1*y(12);
dydt(13) = -r_bound_byA1_B1*y(13) + r_off_B1*y(12);
dydt(14) = r_prodon_B1*y(6) - r_deg_B1*y(8);
end