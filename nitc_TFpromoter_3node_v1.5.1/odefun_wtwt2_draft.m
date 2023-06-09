function dydt = odefun_wtwt2(t,y,rates)
% diploid, wt/wt ODE, manually excludes paralog and nonsense

r_prodbasal_A1 = rates(1);
r_prodbasal_Anonsense1 = rates(2);
r_prodbasal_Aprime1 = rates(3);
r_prodbasal_B1 = rates(4);
r_prodon_A1 = rates(5);
r_prodon_Anonsense1 = rates(6);
r_prodon_Aprime1 = rates(7);
r_prodon_B1 = rates(8);
r_deg_A1 = rates(9);
r_deg_Anonsense1 = rates(10);
r_deg_Aprime1 = rates(11);
r_deg_B1 = rates(12);
r_onbasal_A1 = rates(13);
r_onbasal_Anonsense1 = rates(14);
r_onbasal_Aprime1 = rates(15);
r_onbasal_B1 = rates(16);
r_nitc_byAnonsense1_A1 = rates(17);
r_nitc_byAnonsense1_Anonsense1 = rates(18);
r_nitc_byAnonsense1_Aprime1 = rates(19);
r_bound_byA1_B1 = rates(20);
r_bound_byAprime1_B1 = rates(21);
r_off_A1 = rates(22);
r_off_Anonsense1 = rates(23);
r_off_Aprime1 = rates(24);
r_off_B1 = rates(25);
r_bind_byA1_B1 = rates(26);
r_bind_byAprime1_B1 = rates(27);
k_A1 = rates(28);
k_Anonsense1 = rates(29);
k_Aprime1 = rates(30);
k_B1 = rates(31);
n_A1 = rates(32);
n_Anonsense1 = rates(33);
n_Aprime1 = rates(34);
n_B1 = rates(35);
r_unbind_byA1_B1 = rates(36);
r_unbind_byAprime1_B1 = rates(37);

% r_onbasal_A1 = rates(13);
% r_onbasal_Anonsense1 = rates(14);
% r_nitc_byAnonsense1_A1 = rates(17);
% r_nitc_byAnonsense1_Anonsense1 = rates(18);
% r_nitc_byAnonsense1_Aprime1 = rates(19);
% r_off_A1 = rates(22);
% r_off_Anonsense1 = rates(23);
% r_prodon_A1 = rates(5);
% r_prodon_Anonsense1 = rates(6);
% r_deg_A1 = rates(9);
% r_bind_byA1_B1 = rates(26);
% n_A1 = rates(32);
% k_A1 = rates(28);
% r_unbind_byA1_B1 = rates(36);
% r_bound_byA1_B1 = rates(20);
% r_off_B1 = rates(25);
% r_prodon_B1 = rates(8);
% r_deg_B1 = rates(12);

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

dydt = zeros(22,1);
% y(1) = ona_1 (wt)
dydt(1) = y(2)*(r_onbasal_A1 + r_nitc_byAnonsense1_A1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) ...
    - r_off_A1*y(1);
% y(2) = offa_1 (wt)
dydt(2) = -y(2)*((r_onbasal_A1 + r_nitc_byAnonsense1_A1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1))) ...
    + r_off_A1*y(1);
% y(3) = ona_2 (mut)
dydt(3) = y(4)*(r_onbasal_A1 + r_nitc_byAnonsense1_Anonsense1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) ...
    - r_off_A1*y(3);
% y(4) = offa_2 (mut)
dydt(4) = -y(4)*(r_onbasal_A1 + r_nitc_byAnonsense1_Anonsense1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) ...
    + r_off_A1*y(3);
% y(5) = a 
dydt(5) = r_prodon_A1*y(1) + r_prodon_A1*y(3) - r_deg_A1*y(5);
% y(6) = anons
dydt(6) = 0;
% y(7) = onap_1 
dydt(7) = y(8)*(r_onbasal_Aprime1 + r_nitc_byAnonsense1_A1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) ...
    - r_off_Aprime1*y(7);
% y(8) = offap_1 
dydt(8) = -y(8)*((r_onbasal_Aprime1 + r_nitc_byAnonsense1_A1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1))) ...
    + r_off_Aprime1*y(7);
% y(9) = onap_2 
dydt(9) = y(10)*(r_onbasal_Aprime1 + r_nitc_byAnonsense1_Aprime1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) ...
    - r_off_Aprime1*y(9);
% y(10) = offap_2 
dydt(10) = -y(10)*((r_onbasal_Aprime1 + r_nitc_byAnonsense1_Aprime1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1))) ...
    + r_off_Aprime1*y(9);
% y(11) = aprim 
dydt(11) = r_prodon_Aprime1*y(7) + r_prodon_Aprime1*y(9) - r_deg_Aprime1*y(11);
% y(12) = bba_1
dydt(12) = (r_bind_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1))*y(14) - r_unbind_byA1_B1*y(12);
% y(13) = bbap_1
dydt(13) = (r_bind_byAprime1_B1*y(11)^n_Aprime1/(k_Aprime1^n_Aprime1 + y(11)^n_Aprime1))*y(14) - r_unbind_byAprime1_B1*y(13);
% y(14) = unb_1
dydt(14) = -y(14)*((r_bind_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1)) + (r_bind_byAprime1_B1*y(11)^n_Aprime1/(k_Aprime1^n_Aprime1 + y(11)^n_Aprime1))) ...
    + r_unbind_byA1_B1*y(12) + r_unbind_byAprime1_B1*y(13);
% y(15) = bba_2
dydt(15) = (r_bind_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1))*y(17) - r_unbind_byA1_B1*y(15);
% y(16) = bbap_2
dydt(16) = (r_bind_byAprime1_B1*y(11)^n_Aprime1/(k_Aprime1^n_Aprime1 + y(11)^n_Aprime1))*y(17) - r_unbind_byAprime1_B1*y(16);
% y(17) = unb_2
dydt(17) = -y(17)*((r_bind_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1)) + (r_bind_byAprime1_B1*y(11)^n_Aprime1/(k_Aprime1^n_Aprime1 + y(11)^n_Aprime1))) ...
    + r_unbind_byA1_B1*y(15) + r_unbind_byAprime1_B1*y(16);
% y(18) = onb_1
dydt(18) = y(19)*(r_bound_byA1_B1*y(12) + r_bound_byAprime1_B1*y(13)) - r_off_B1*y(18);
% y(19) = offb_1
dydt(19) = -y(19)*(r_bound_byA1_B1*y(12) + r_bound_byAprime1_B1*y(13)) + r_off_B1*y(18);
% y(20) = onb_2
dydt(20) = y(21)*(r_bound_byA1_B1*y(15) + r_bound_byAprime1_B1*y(16)) - r_off_B1*y(20);
% y(21) = offb_2
dydt(21) = -y(21)*(r_bound_byA1_B1*y(15) + r_bound_byAprime1_B1*y(16)) + r_off_B1*y(20);
% y(22) = b
dydt(22) = r_prodon_B1*y(18) + r_prodon_B1*y(20) - r_deg_B1*y(22);
end