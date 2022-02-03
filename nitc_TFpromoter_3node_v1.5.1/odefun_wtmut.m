function dydt = odefun_wtmut(t,y,rates)
% diploid, wt/wt ODE, manually excludes paralog and nonsense

% vector: y: 
% y(1) = ona_1 (wt)
% y(2) = offa_1 (wt)
% y(3) = ona_2 (mut)
% y(4) = offa_2 (mut)
% y(5) = a 
% y(6) = anons
% y(7) = onap_1 
% y(8) = offap_1 
% y(9) = onap_2 
% y(10) = offap_2 
% y(11) = aprim 
% y(12) = bba_1
% y(13) = bbap_1
% y(14) = unb_1
% y(15) = bba_2
% y(16) = bbap_2
% y(17) = unb_2
% y(18) = onb_1
% y(19) = offb_1
% y(20) = onb_2
% y(21) = offb_2
% y(22) = b

r_onbasal_A1 = rates(13);
r_onbasal_Anonsense1 = rates(14);
r_nitc_byAnonsense1_A1 = rates(17);
r_nitc_byAnonsense1_Anonsense1 = rates(18);
r_nitc_byAnonsense1_Aprime1 = rates(19);
r_off_A1 = rates(22);
r_off_Anonsense1 = rates(23);
r_prodon_A1 = rates(5);
r_prodon_Anonsense1 = rates(6);
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
dydt(1) = y(2)*(r_onbasal_A1 + r_nitc_byAnonsense1_A1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) - r_off_A1*y(1);
dydt(2) = -y(2)*((r_onbasal_A1 + r_nitc_byAnonsense1_A1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1))) + r_off_A1*y(1);
dydt(3) = y(4)*(r_onbasal_Anonsense1* + r_nitc_byAnonsense1_Anonsense1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) - r_off_Anonsense1*y(3);
dydt(4) = -y(4)*((r_onbasal_Anonsense1 + r_nitc_byAnonsense1_Anonsense1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1))) + r_off_Anonsense1*y(3);
dydt(5) = r_prodon_A1*y(1) - r_deg_A1*y(5);
dydt(6) = r_prodon_Anonsense1*y(1) - r_deg_Anonsense1*y(6);
dydt(7) = y(8)*(r_onbasal_Aprime1 + r_nitc_byAnonsense1_A1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) - r_off_Aprime1*y(7);
dydt(8) = -y(8)*((r_onbasal_Aprime1 + r_nitc_byAnonsense1_A1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1))) + r_off_Aprime1*y(7);
dydt(9) = y(10)*(r_onbasal_Aprime1 + r_nitc_byAnonsense1_Aprime1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) - r_off_Aprime1*y(9);
dydt(10) = -y(10)*((r_onbasal_Aprime1 + r_nitc_byAnonsense1_Aprime1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1))) + r_off_Aprime1*y(9);
dydt(11) = r_prodon_Aprime1*y(7) + r_prodon_Aprime1*y(9) - r_deg_A1*y(11);
dydt(6) = (r_bind_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1))*y(7) - r_unbind_byA1_B1*y(6);
dydt(7) = -(r_bind_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1))*y(7) + r_unbind_byA1_B1*y(6);
dydt(8) = (r_bind_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1))*y(9) - r_unbind_byA1_B1*y(8);
dydt(9) = -(r_bind_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1))*y(9) + r_unbind_byA1_B1*y(8);
dydt(10) = r_bound_byA1_B1*y(11) - r_off_B1*y(10);
dydt(11) = -r_bound_byA1_B1*y(11) + r_off_B1*y(10);
dydt(12) = r_bound_byA1_B1*y(13) - r_off_B1*y(12);
dydt(13) = -r_bound_byA1_B1*y(13) + r_off_B1*y(12);
dydt(14) = r_prodon_B1*y(10) + r_prodon_B1*y(12) - r_deg_B1*y(14);
end