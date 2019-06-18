/** @file  DIRK_Jacobian_cBSSN_Reg.h
 *  @author Francesco Torsello
 *  @brief The Jacobian of the regularized cBSSN equations, needed by DIRK.
 *  @version 2019-06-18T14:56:08
 *  @image html DIRK_Jacobian_cBSSN_Reg.png
 */

Jacobian[1][1]=
	(-2 * exp(-2 * gconf(m,n)) * gconf_r(m,n) 
    /*1*/  * gAlp(m,n) * p(m,n)) / (gA(m,n) * Lt(m,n));

Jacobian[1][3]=
	-gAlp(m,n) / 6.;

Jacobian[1][5]=
	-((exp(-2 * gconf(m,n)) * gconf_r(m,n) * gAlp(m,n)
    /*2*/  * p(m,n)) / (pow2(gA(m,n)) * Lt(m,n)));

Jacobian[1][22]=
	gconf_r(m,n);

Jacobian[2][2]=
	(2 * exp(-2 * fconf(m,n)) * fconf_r(m,n) 
    /*1*/  * fAlp(m,n) * p(m,n)) / (fA(m,n) * Lt(m,n));

Jacobian[2][4]=
	-fAlp(m,n) / 6.;

Jacobian[2][7]=
	(exp(-2 * fconf(m,n)) * fconf_r(m,n) * fAlp(m,n) 
    /*1*/  * p(m,n)) / (pow2(fA(m,n)) * Lt(m,n));

Jacobian[2][22]=
	fconf_r(m,n);

Jacobian[3][1]=
	(exp(-4 * gconf(m,n)) * (-(k_g * exp(2 
    /*4*/  * (fconf(m,n) + gconf(m,n))) * fA(m,n) * gB(m,n) 
    /*3*/  * pow2(gA(m,n)) * r(m,n) * (gAlp(m,n) * ((2 + 4 
    /*6*/  * pow2(Lt(m,n))) * P_1_1(R(m,n)) + (3 - 6 
    /*6*/  * pow2(Lt(m,n))) * P_2_1(R(m,n)) + 2 * b_1) + 2 
    /*4*/  * fAlp(m,n) * (2 * P_1_2(R(m,n)) + b_2) * Lt(m,n)))
    /*2*/  - 2 * (k_g * exp(4 * gconf(m,n)) * fAlp(m,n) 
    /*3*/  * gB(m,n) * pow3(gA(m,n)) * r(m,n) * (P_2_1(R(m,n))
    /*4*/  + b_1) - 2 * (-(gAlp_r(m,n) * gA_r(m,n) * gB(m,n)
    /*5*/  * r(m,n)) + gA(m,n) * (2 * gAlp_r(m,n) * gB_r(m,n)
    /*5*/  * r(m,n) + gB(m,n) * (gAlp_rr(m,n) * r(m,n) + 2 
    /*6*/  * gAlp_r(m,n) * (1 + gconf_r(m,n) * r(m,n))))) 
    /*3*/  * Lt(m,n) + exp(2 * gconf(m,n)) * gAlp(m,n) 
    /*3*/  * gB(m,n) * pow2(gA(m,n)) * r(m,n) * (p(m,n) 
    /*4*/  * gtrK_r(m,n) + k_g * exp(2 * gconf(m,n)) * gA(m,n)
    /*4*/  * (P_1_0(R(m,n)) + b_0) * Lt(m,n))))) / (gB(m,n) 
    /*1*/  * pow3(gA(m,n)) * r(m,n) * Lt(m,n));

Jacobian[3][2]=
	(k_g * exp(-2 * gconf(m,n)) * (2 * exp(2 
    /*3*/  * gconf(m,n)) * gA(m,n) * (fAlp(m,n) 
    /*3*/  * (P_2_1(R(m,n)) + b_1) + gAlp(m,n) * (P_1_0(R(m,n))
    /*4*/  + b_0) * Lt(m,n)) + exp(2 * fconf(m,n)) * fA(m,n)
    /*2*/  * (gAlp(m,n) * ((2 + 4 * pow2(Lt(m,n))) 
    /*4*/  * P_1_1(R(m,n)) + (3 - 6 * pow2(Lt(m,n))) 
    /*4*/  * P_2_1(R(m,n)) + 2 * b_1) + 2 * fAlp(m,n) * (2 
    /*4*/  * P_1_2(R(m,n)) + b_2) * Lt(m,n)))) / (gA(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[3][3]=
	(2 * gAlp(m,n) * (gtrA(m,n) + gtrK(m,n))) / 3.;

Jacobian[3][5]=
	(exp(-4 * gconf(m,n)) * (-2 * exp(2 * gconf(m,n))
    /*2*/  * gAlp(m,n) * gB(m,n) * p(m,n) * pow2(gA(m,n)) 
    /*2*/  * r(m,n) * gtrK_r(m,n) - 6 * gAlp_r(m,n) * gA_r(m,n)
    /*2*/  * gB(m,n) * r(m,n) * Lt(m,n) + 4 * gA(m,n) * (2 
    /*3*/  * gAlp_r(m,n) * gB_r(m,n) * r(m,n) + gB(m,n) 
    /*3*/  * (gAlp_rr(m,n) * r(m,n) + 2 * gAlp_r(m,n) * (1 
    /*5*/  + gconf_r(m,n) * r(m,n)))) * Lt(m,n) - k_g * exp(2 
    /*3*/  * (fconf(m,n) + gconf(m,n))) * fA(m,n) * gB(m,n) 
    /*2*/  * pow2(gA(m,n)) * r(m,n) * (gAlp(m,n) * (2 
    /*4*/  * P_1_1(R(m,n)) + (1 - 2 * pow2(Lt(m,n))) 
    /*4*/  * P_2_1(R(m,n))) + 2 * fAlp(m,n) * P_1_2(R(m,n)) 
    /*3*/  * Lt(m,n)))) / (2. * Power(gA(m,n),4) * gB(m,n) 
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[3][6]=
	(exp(-4 * gconf(m,n)) * (-(k_g * exp(4 
    /*4*/  * gconf(m,n)) * fAlp(m,n) * gB(m,n) * pow2(gA(m,n))
    /*3*/  * R(m,n) * (P_1_2(R(m,n)) - b_2)) + (2 
    /*3*/  * gAlp_r(m,n) * gB_r(m,n) + k_g * exp(4 
    /*4*/  * gconf(m,n)) * gAlp(m,n) * gB(m,n) * pow2(gA(m,n))
    /*3*/  * R(m,n) * b_1) * Lt(m,n) + k_g * exp(2 
    /*3*/  * (fconf(m,n) + gconf(m,n))) * fA(m,n) * gA(m,n) 
    /*2*/  * gB(m,n) * R(m,n) * (gAlp(m,n) * ((-1 + 2 
    /*5*/  * pow2(Lt(m,n))) * P_1_2(R(m,n)) + b_2) + fAlp(m,n)
    /*3*/  * b_3 * Lt(m,n)))) / (pow2(gA(m,n)) 
    /*1*/  * pow2(gB(m,n)) * Lt(m,n));

Jacobian[3][7]=
	(k_g * exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*1*/  * (gAlp(m,n) * (2 * P_1_1(R(m,n)) + (1 - 2 
    /*4*/  * pow2(Lt(m,n))) * P_2_1(R(m,n))) + 2 * fAlp(m,n) 
    /*2*/  * P_1_2(R(m,n)) * Lt(m,n))) / (2. * gA(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[3][8]=
	-((k_g * exp(2 * fconf(m,n) - 4 * gconf(m,n)) 
    /*2*/  * (-(exp(2 * gconf(m,n)) * gA(m,n) * (fAlp(m,n) 
    /*5*/  * (P_1_2(R(m,n)) - b_2) - gAlp(m,n) * b_1 
    /*5*/  * Lt(m,n))) + exp(2 * fconf(m,n)) * fA(m,n) 
    /*3*/  * (gAlp(m,n) * ((-1 + 2 * pow2(Lt(m,n))) 
    /*5*/  * P_1_2(R(m,n)) + b_2) + fAlp(m,n) * b_3 
    /*4*/  * Lt(m,n)))) / (gA(m,n) * gB(m,n) * Lt(m,n)));

Jacobian[3][9]=
	2 * gA1(m,n) * gAlp(m,n);

Jacobian[3][10]=
	4 * gA2(m,n) * gAlp(m,n);

Jacobian[3][22]=
	gtrK_r(m,n);

Jacobian[4][1]=
	(k_f * exp(-2 * fconf(m,n)) * (-(exp(2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (fAlp(m,n) * ((-2 + 4 
    /*6*/  * pow2(Lt(m,n))) * P_1_1(R(m,n)) + (-3 + 2 
    /*6*/  * pow2(Lt(m,n))) * P_2_1(R(m,n)) - 2 * b_1) - 2 
    /*4*/  * gAlp(m,n) * (P_2_0(R(m,n)) + b_0) * Lt(m,n))) + 2
    /*2*/  * exp(2 * fconf(m,n)) * fA(m,n) * (gAlp(m,n) * (2
    /*4*/  * P_1_1(R(m,n)) + b_1) + fAlp(m,n) 
    /*3*/  * (P_1_2(R(m,n)) + b_2) * Lt(m,n)))) / (fA(m,n) 
    /*1*/  * pow2(R(m,n)) * Lt(m,n));

Jacobian[4][2]=
	(exp(-2 * (2 * fconf(m,n) + gconf(m,n))) * (-4 
    /*2*/  * exp(2 * gconf(m,n)) * fAlp_r(m,n) * fA_r(m,n) 
    /*2*/  * gB(m,n) * pow2(R(m,n)) * r(m,n) * Lt(m,n) + 4 
    /*2*/  * fA(m,n) * R(m,n) * (2 * exp(2 * fconf(m,n)) 
    /*3*/  * fAlp_r(m,n) * fB_r(m,n) * r(m,n) + exp(2 
    /*4*/  * gconf(m,n)) * gB(m,n) * (fAlp_rr(m,n) * r(m,n) + 2
    /*4*/  * fAlp_r(m,n) * (1 + fconf_r(m,n) * r(m,n))) 
    /*3*/  * R(m,n)) * Lt(m,n) + exp(2 * (fconf(m,n) 
    /*4*/  + gconf(m,n))) * gB(m,n) * pow2(fA(m,n)) * r(m,n) 
    /*2*/  * (fAlp(m,n) * (2 * p(m,n) * pow2(R(m,n)) 
    /*4*/  * ftrK_r(m,n) + k_f * exp(2 * gconf(m,n)) * gA(m,n)
    /*4*/  * ((-2 + 4 * pow2(Lt(m,n))) * P_1_1(R(m,n)) + (-3
    /*6*/  + 2 * pow2(Lt(m,n))) * P_2_1(R(m,n)) - 2 * b_1)) 
    /*3*/  - 2 * k_f * exp(2 * gconf(m,n)) * gAlp(m,n) 
    /*3*/  * gA(m,n) * (P_2_0(R(m,n)) + b_0) * Lt(m,n)) - 2 
    /*2*/  * k_f * exp(4 * fconf(m,n) + 2 * gconf(m,n)) 
    /*2*/  * gB(m,n) * pow3(fA(m,n)) * r(m,n) * (gAlp(m,n) * (2
    /*4*/  * P_1_1(R(m,n)) + b_1) + fAlp(m,n) 
    /*3*/  * (P_1_2(R(m,n)) + b_2) * Lt(m,n)))) / (gB(m,n) 
    /*1*/  * pow2(R(m,n)) * pow3(fA(m,n)) * r(m,n) * Lt(m,n));

Jacobian[4][4]=
	(2 * fAlp(m,n) * (ftrA(m,n) + ftrK(m,n))) / 3.;

Jacobian[4][5]=
	-(k_f * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*1*/  * (fAlp(m,n) * (2 * P_1_1(R(m,n)) + (-3 + 2 
    /*4*/  * pow2(Lt(m,n))) * P_2_1(R(m,n))) + 2 * gAlp(m,n) 
    /*2*/  * (P_1_0(R(m,n)) - P_2_0(R(m,n))) * Lt(m,n))) / (2.
    /*1*/  * fA(m,n) * pow2(R(m,n)) * Lt(m,n));

Jacobian[4][6]=
	-((k_f * exp(-2 * fconf(m,n)) * (gAlp(m,n) 
    /*3*/  * (exp(2 * fconf(m,n)) * fA(m,n) * (b_1 + 2 * R(m,n)
    /*5*/  * b_2) - exp(2 * gconf(m,n)) * gA(m,n) 
    /*4*/  * (P_1_0(R(m,n)) + b_0) * Lt(m,n)) + fAlp(m,n) 
    /*3*/  * (exp(2 * gconf(m,n)) * gA(m,n) * (2 * (-1 
    /*6*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) - b_1) + exp(2 
    /*5*/  * fconf(m,n)) * fA(m,n) * R(m,n) * b_3 * Lt(m,n))))
    /*1*/  / (fA(m,n) * gB(m,n) * pow2(R(m,n)) * Lt(m,n)));

Jacobian[4][7]=
	(exp(-2 * (2 * fconf(m,n) + gconf(m,n))) * (exp(2
    /*3*/  * (fconf(m,n) + gconf(m,n))) * fAlp(m,n) * gB(m,n)
    /*2*/  * pow2(fA(m,n)) * r(m,n) * (k_f * exp(2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (2 * P_1_1(R(m,n)) + (-3 
    /*5*/  + 2 * pow2(Lt(m,n))) * P_2_1(R(m,n))) + 2 * p(m,n) 
    /*3*/  * pow2(R(m,n)) * ftrK_r(m,n)) + 2 * R(m,n) * (-(k_f
    /*4*/  * exp(2 * fconf(m,n) + 4 * gconf(m,n)) * gAlp(m,n)
    /*4*/  * gA(m,n) * gB(m,n) * pow2(fA(m,n)) 
    /*4*/  * P_1_1(R(m,n)) * r(m,n)) - 3 * exp(2 * gconf(m,n))
    /*3*/  * fAlp_r(m,n) * fA_r(m,n) * gB(m,n) * r(m,n) 
    /*3*/  * R(m,n) + 2 * fA(m,n) * (2 * exp(2 * fconf(m,n)) 
    /*4*/  * fAlp_r(m,n) * fB_r(m,n) * r(m,n) + exp(2 
    /*5*/  * gconf(m,n)) * gB(m,n) * (fAlp_rr(m,n) * r(m,n) + 2
    /*5*/  * fAlp_r(m,n) * (1 + fconf_r(m,n) * r(m,n))) 
    /*4*/  * R(m,n))) * Lt(m,n))) / (2. * Power(fA(m,n),4) 
    /*1*/  * gB(m,n) * pow2(R(m,n)) * r(m,n) * Lt(m,n));

Jacobian[4][8]=
	(exp(-4 * gconf(m,n)) * (k_f * exp(2 * (fconf(m,n)
    /*4*/  + gconf(m,n))) * gAlp(m,n) * gB(m,n) 
    /*2*/  * pow2(fA(m,n)) * (b_1 + 2 * R(m,n) * b_2) + 2 
    /*2*/  * fAlp_r(m,n) * fB_r(m,n) * R(m,n) * Lt(m,n) - k_f 
    /*2*/  * exp(4 * gconf(m,n)) * fA(m,n) * gAlp(m,n) 
    /*2*/  * gA(m,n) * gB(m,n) * (P_1_0(R(m,n)) + b_0) 
    /*2*/  * Lt(m,n) + k_f * exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*2*/  * fA(m,n) * gB(m,n) * (exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * (2 * (-1 + pow2(Lt(m,n))) * P_1_1(R(m,n)) - b_1)
    /*3*/  + exp(2 * fconf(m,n)) * fA(m,n) * R(m,n) * b_3 
    /*3*/  * Lt(m,n)))) / (pow2(fA(m,n)) * pow2(gB(m,n)) 
    /*1*/  * pow3(R(m,n)) * Lt(m,n));

Jacobian[4][11]=
	2 * fA1(m,n) * fAlp(m,n);

Jacobian[4][12]=
	4 * fA2(m,n) * fAlp(m,n);

Jacobian[4][22]=
	ftrK_r(m,n);

Jacobian[5][1]=
	(-2 * exp(-2 * gconf(m,n)) * gAlp(m,n) * gA_r(m,n)
    /*1*/  * p(m,n)) / (gA(m,n) * Lt(m,n));

Jacobian[5][5]=
	gdet_pff(m,n) / (6. * gdet(m,n)) - gA1(m,n) 
    /*0*/  * gAlp(m,n) + gBet_r(m,n) + gAlp(m,n) * (gtrA(m,n) 
    /*1*/  / 3. - (exp(-2 * gconf(m,n)) * gA_r(m,n) * p(m,n)) 
    /*1*/  / (pow2(gA(m,n)) * Lt(m,n)));

Jacobian[5][9]=
	-(gAlp(m,n) * gA(m,n));

Jacobian[5][22]=
	gA_r(m,n);

Jacobian[6][1]=
	(-2 * exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) 
    /*1*/  * (gB(m,n) + gB_r(m,n) * r(m,n))) / (gA(m,n) 
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[6][5]=
	-((exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) 
    /*2*/  * (gB(m,n) + gB_r(m,n) * r(m,n))) / (pow2(gA(m,n)) 
    /*2*/  * r(m,n) * Lt(m,n)));

Jacobian[6][6]=
	gdet_pff(m,n) / (6. * gdet(m,n)) - gA2(m,n) 
    /*0*/  * gAlp(m,n) + gBetr(m,n) - (exp(-2 * gconf(m,n)) 
    /*1*/  * gAlp(m,n) * p(m,n)) / (gA(m,n) * Lt(m,n) * r(m,n))
    /*0*/  + (gAlp(m,n) * gtrA(m,n)) / 3. + (exp(-2 
    /*2*/  * gconf(m,n)) * gAlp(m,n) * p(m,n)) / (gA(m,n) 
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[6][10]=
	-(gAlp(m,n) * gB(m,n));

Jacobian[6][22]=
	gB_r(m,n) + gB(m,n) / r(m,n);

Jacobian[7][2]=
	(2 * exp(-2 * fconf(m,n)) * fAlp(m,n) * fA_r(m,n)
    /*1*/  * p(m,n)) / (fA(m,n) * Lt(m,n));

Jacobian[7][7]=
	fdet_pff(m,n) / (6. * fdet(m,n)) - fA1(m,n) 
    /*0*/  * fAlp(m,n) + fBet_r(m,n) + fAlp(m,n) * (ftrA(m,n) 
    /*1*/  / 3. + (exp(-2 * fconf(m,n)) * fA_r(m,n) * p(m,n)) 
    /*1*/  / (pow2(fA(m,n)) * Lt(m,n)));

Jacobian[7][11]=
	-(fAlp(m,n) * fA(m,n));

Jacobian[7][22]=
	fA_r(m,n);

Jacobian[8][2]=
	(2 * exp(-4 * fconf(m,n)) * fAlp(m,n) * p(m,n) 
    /*1*/  * (exp(2 * fconf(m,n)) * fB_r(m,n) * r(m,n) + exp(2
    /*3*/  * gconf(m,n)) * gB(m,n) * R(m,n))) / (fA(m,n) 
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[8][7]=
	(exp(-4 * fconf(m,n)) * fAlp(m,n) * p(m,n) 
    /*1*/  * (exp(2 * fconf(m,n)) * fB_r(m,n) * r(m,n) + exp(2
    /*3*/  * gconf(m,n)) * gB(m,n) * R(m,n))) 
    /*0*/  / (pow2(fA(m,n)) * r(m,n) * Lt(m,n));

Jacobian[8][8]=
	-(fA2(m,n) * fAlp(m,n)) + gBetr(m,n) 
    /*0*/  + (fdet_pff(m,n) / fdet(m,n) - (6 * exp(-2 
    /*3*/  * gconf(m,n)) * gAlp(m,n) * p(m,n)) / (gA(m,n) 
    /*2*/  * Lt(m,n) * r(m,n)) + fAlp(m,n) * (2 * ftrA(m,n) 
    /*2*/  - (6 * exp(-2 * fconf(m,n)) * p(m,n)) / (fA(m,n) 
    /*3*/  * r(m,n) * Lt(m,n)))) / 6.;

Jacobian[8][12]=
	-(exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*1*/  * fAlp(m,n) * gB(m,n) * R(m,n));

Jacobian[8][22]=
	fB_r(m,n) + (exp(-2 * fconf(m,n) + 2 * gconf(m,n))
    /*1*/  * gB(m,n) * R(m,n)) / r(m,n);

Jacobian[9][1]=
	(-2 * exp(-4 * gconf(m,n)) * (3 * exp(2 
    /*3*/  * gconf(m,n)) * gA1_r(m,n) * gAlp(m,n) * p(m,n) 
    /*2*/  * pow3(gA(m,n)) * pow3(gB(m,n)) * r(m,n) + 2 * k_g 
    /*2*/  * exp(4 * gconf(m,n)) * fAlp(m,n) * Power(gA(m,n),4)
    /*2*/  * pow3(gB(m,n)) * r(m,n) * R(m,n) * (b_2 + 2 
    /*3*/  * R(m,n) * b_3) + 8 * gAlp(m,n) * Power(gA(m,n),4) 
    /*2*/  * gB_r(m,n) * Lt(m,n) - 8 * gAlp(m,n) * gsig(m,n) 
    /*2*/  * gB_r(m,n) * pow2(gA(m,n)) * pow2(gB(m,n)) 
    /*2*/  * pow2(r(m,n)) * Lt(m,n) - 4 * gAlp(m,n) * gL(m,n) 
    /*2*/  * Power(gA(m,n),4) * pow3(gB(m,n)) * Lt(m,n) - 8 
    /*2*/  * gAlp(m,n) * gA(m,n) * gA_r(m,n) * pow3(gB(m,n)) 
    /*2*/  * Lt(m,n) + 8 * gconf_r(m,n) * gAlp(m,n) 
    /*2*/  * pow2(gA(m,n)) * pow3(gB(m,n)) * Lt(m,n) + 4 
    /*2*/  * gAlp_r(m,n) * pow2(gA(m,n)) * pow3(gB(m,n)) 
    /*2*/  * Lt(m,n) - 12 * gAlp(m,n) * gA(m,n) * gA_r(m,n) 
    /*2*/  * gB_r(m,n) * pow2(gB(m,n)) * r(m,n) * Lt(m,n) + 8 
    /*2*/  * gconf_r(m,n) * gAlp(m,n) * gB_r(m,n) 
    /*2*/  * pow2(gA(m,n)) * pow2(gB(m,n)) * r(m,n) * Lt(m,n) 
    /*2*/  + 4 * gAlp_r(m,n) * gB_r(m,n) * pow2(gA(m,n)) 
    /*2*/  * pow2(gB(m,n)) * r(m,n) * Lt(m,n) + 4 * gAlp(m,n) 
    /*2*/  * gB_rr(m,n) * pow2(gA(m,n)) * pow2(gB(m,n)) 
    /*2*/  * r(m,n) * Lt(m,n) - 4 * gAlp(m,n) * gB(m,n) 
    /*2*/  * pow2(gA(m,n)) * pow2(gB_r(m,n)) * r(m,n) * Lt(m,n)
    /*2*/  + 4 * gAlp(m,n) * gL_r(m,n) * Power(gA(m,n),4) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * Lt(m,n) + 8 
    /*2*/  * gconf_r(m,n) * gAlp(m,n) * gA(m,n) * gA_r(m,n) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * Lt(m,n) + 4 * gAlp_r(m,n)
    /*2*/  * gA(m,n) * gA_r(m,n) * pow3(gB(m,n)) * r(m,n) 
    /*2*/  * Lt(m,n) - 4 * gAlp(m,n) * gA(m,n) * gA_rr(m,n) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * Lt(m,n) - 8 
    /*2*/  * gconf_rr(m,n) * gAlp(m,n) * pow2(gA(m,n)) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * Lt(m,n) + 16 
    /*2*/  * gconf_r(m,n) * gAlp_r(m,n) * pow2(gA(m,n)) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * Lt(m,n) - 4 
    /*2*/  * gAlp_rr(m,n) * pow2(gA(m,n)) * pow3(gB(m,n)) 
    /*2*/  * r(m,n) * Lt(m,n) - 12 * gAlp(m,n) * gsig(m,n) 
    /*2*/  * pow2(gA(m,n)) * pow3(gB(m,n)) * r(m,n) * Lt(m,n) 
    /*2*/  + 16 * gAlp(m,n) * pow2(gconf_r(m,n)) 
    /*2*/  * pow2(gA(m,n)) * pow3(gB(m,n)) * r(m,n) * Lt(m,n) 
    /*2*/  + 12 * gAlp(m,n) * pow2(gA_r(m,n)) * pow3(gB(m,n)) 
    /*2*/  * r(m,n) * Lt(m,n) + 2 * k_g * exp(4 * gconf(m,n)) 
    /*2*/  * gAlp(m,n) * Power(gA(m,n),4) * pow3(gB(m,n)) 
    /*2*/  * r(m,n) * R(m,n) * b_1 * Lt(m,n) + 4 * k_g * exp(4
    /*3*/  * gconf(m,n)) * gAlp(m,n) * Power(gA(m,n),4) 
    /*2*/  * pow2(R(m,n)) * pow3(gB(m,n)) * r(m,n) * b_2 
    /*2*/  * Lt(m,n) - 2 * k_g * exp(2 * (fconf(m,n) 
    /*4*/  + gconf(m,n))) * fA(m,n) * pow3(gA(m,n)) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * (gAlp(m,n) * (2 * (-2 
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) - 3 * (-1 
    /*5*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)) - b_1) - fAlp(m,n)
    /*3*/  * (2 * P_1_2(R(m,n)) + b_2) * Lt(m,n)))) / (3. 
    /*1*/  * Power(gA(m,n),4) * pow3(gB(m,n)) * r(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[9][2]=
	(4 * k_g * exp(-2 * gconf(m,n)) * (-(exp(2 
    /*4*/  * fconf(m,n)) * fA(m,n) * (gAlp(m,n) * (2 * (-2 
    /*6*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) - 3 * (-1 
    /*6*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)) - b_1) - fAlp(m,n)
    /*4*/  * (2 * P_1_2(R(m,n)) + b_2) * Lt(m,n))) + exp(2 
    /*3*/  * gconf(m,n)) * gA(m,n) * R(m,n) * (fAlp(m,n) * (b_2
    /*4*/  + 2 * R(m,n) * b_3) + gAlp(m,n) * (b_1 + 2 
    /*4*/  * R(m,n) * b_2) * Lt(m,n)))) / (3. * gA(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[9][3]=
	(gAlp(m,n) * (3 * gA1(m,n) - gtrA(m,n))) / 3.;

Jacobian[9][5]=
	(exp(-4 * gconf(m,n)) * (-2 * gA(m,n) * gB(m,n) 
    /*2*/  * (3 * gAlp_r(m,n) * gA_r(m,n) * gB(m,n) * r(m,n) 
    /*3*/  + 2 * gA(m,n) * (gAlp_r(m,n) * gB_r(m,n) * r(m,n) 
    /*4*/  + gB(m,n) * (-(gAlp_rr(m,n) * r(m,n)) + gAlp_r(m,n)
    /*5*/  * (1 + 4 * gconf_r(m,n) * r(m,n))))) * Lt(m,n) - 2
    /*2*/  * k_g * exp(2 * (fconf(m,n) + gconf(m,n))) 
    /*2*/  * fA(m,n) * pow2(gB(m,n)) * pow3(gA(m,n)) * r(m,n) 
    /*2*/  * (gAlp(m,n) * (P_1_1(R(m,n)) + (-1 + pow2(Lt(m,n)))
    /*4*/  * P_2_1(R(m,n))) + fAlp(m,n) * P_1_2(R(m,n)) 
    /*3*/  * Lt(m,n)) + gAlp(m,n) * (-3 * exp(2 * gconf(m,n)) 
    /*3*/  * gA1_r(m,n) * p(m,n) * pow2(gB(m,n)) 
    /*3*/  * pow3(gA(m,n)) * r(m,n) - 24 * pow2(gA_r(m,n)) 
    /*3*/  * pow2(gB(m,n)) * r(m,n) * Lt(m,n) + 4 * gsig(m,n) 
    /*3*/  * gB(m,n) * pow2(gA(m,n)) * r(m,n) * (3 * gB(m,n) 
    /*4*/  + 2 * gB_r(m,n) * r(m,n)) * Lt(m,n) - 4 
    /*3*/  * pow2(gA(m,n)) * (gB(m,n) * (2 * gconf_r(m,n) 
    /*5*/  * gB_r(m,n) + gB_rr(m,n)) * r(m,n) - pow2(gB_r(m,n))
    /*4*/  * r(m,n) + 2 * pow2(gB(m,n)) * (gconf_r(m,n) 
    /*5*/  - gconf_rr(m,n) * r(m,n) + 2 * pow2(gconf_r(m,n)) 
    /*5*/  * r(m,n))) * Lt(m,n) + 6 * gA(m,n) * gB(m,n) * (3 
    /*4*/  * gA_r(m,n) * gB_r(m,n) * r(m,n) + gB(m,n) 
    /*4*/  * (gA_rr(m,n) * r(m,n) + gA_r(m,n) * (2 - 2 
    /*6*/  * gconf_r(m,n) * r(m,n)))) * Lt(m,n)))) / (3. 
    /*1*/  * Power(gA(m,n),5) * pow2(gB(m,n)) * r(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[9][6]=
	(-2 * exp(-4 * gconf(m,n)) * (k_g * exp(4 
    /*3*/  * gconf(m,n)) * fAlp(m,n) * pow3(gA(m,n)) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * R(m,n) * (b_2 + 2 
    /*3*/  * R(m,n) * b_3) + (gAlp_r(m,n) * gA(m,n) * gB_r(m,n)
    /*3*/  * pow2(gB(m,n)) * r(m,n) + gAlp(m,n) * (-3 
    /*4*/  * gA_r(m,n) * gB_r(m,n) * pow2(gB(m,n)) * r(m,n) 
    /*4*/  + gA(m,n) * gB(m,n) * r(m,n) * (gB(m,n) * (2 
    /*6*/  * gconf_r(m,n) * gB_r(m,n) + gB_rr(m,n)) - 2 
    /*5*/  * pow2(gB_r(m,n)) - 2 * gsig(m,n) * gB(m,n) 
    /*5*/  * gB_r(m,n) * r(m,n)) + pow3(gA(m,n)) * (6 
    /*5*/  * gB_r(m,n) + k_g * exp(4 * gconf(m,n)) 
    /*5*/  * pow3(gB(m,n)) * r(m,n) * R(m,n) * (b_1 + 2 
    /*6*/  * R(m,n) * b_2)))) * Lt(m,n) + k_g * exp(2 
    /*3*/  * (fconf(m,n) + gconf(m,n))) * fA(m,n) 
    /*2*/  * pow2(gA(m,n)) * pow3(gB(m,n)) * r(m,n) * R(m,n) 
    /*2*/  * (gAlp(m,n) * (2 * (-1 + pow2(Lt(m,n))) 
    /*4*/  * P_1_2(R(m,n)) - b_2) - fAlp(m,n) * b_3 
    /*3*/  * Lt(m,n)))) / (3. * Power(gB(m,n),4) 
    /*1*/  * pow3(gA(m,n)) * r(m,n) * Lt(m,n));

Jacobian[9][7]=
	(2 * k_g * exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*1*/  * (gAlp(m,n) * (P_1_1(R(m,n)) + (-1 + pow2(Lt(m,n)))
    /*3*/  * P_2_1(R(m,n))) + fAlp(m,n) * P_1_2(R(m,n)) 
    /*2*/  * Lt(m,n))) / (3. * gA(m,n) * Lt(m,n));

Jacobian[9][8]=
	(-2 * k_g * exp(2 * fconf(m,n) - 4 * gconf(m,n)) 
    /*1*/  * (exp(2 * gconf(m,n)) * gA(m,n) * (fAlp(m,n) * (2 
    /*4*/  * P_1_2(R(m,n)) + b_2) + gAlp(m,n) * (2 
    /*4*/  * P_1_1(R(m,n)) + b_1) * Lt(m,n)) - exp(2 
    /*3*/  * fconf(m,n)) * fA(m,n) * (gAlp(m,n) * (2 * (-1 
    /*5*/  + pow2(Lt(m,n))) * P_1_2(R(m,n)) - b_2) - fAlp(m,n)
    /*3*/  * b_3 * Lt(m,n)))) / (3. * gA(m,n) * gB(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[9][9]=
	gAlp(m,n) * (gtrA(m,n) + gtrK(m,n));

Jacobian[9][13]=
	(-2 * exp(-4 * gconf(m,n)) * gAlp(m,n)) / (3. 
    /*1*/  * r(m,n));

Jacobian[9][15]=
	(-2 * exp(-4 * gconf(m,n)) * gAlp(m,n) * (3 
    /*2*/  * gB(m,n) + 2 * gB_r(m,n) * r(m,n))) / (3. * gB(m,n)
    /*1*/  * pow2(gA(m,n)));

Jacobian[9][22]=
	gA1_r(m,n);

Jacobian[10][1]=
	(2 * exp(-4 * gconf(m,n)) * (-3 * exp(2 
    /*3*/  * gconf(m,n)) * gA2_r(m,n) * gAlp(m,n) * p(m,n) 
    /*2*/  * pow3(gA(m,n)) * pow3(gB(m,n)) * r(m,n) + k_g 
    /*2*/  * exp(4 * gconf(m,n)) * fAlp(m,n) * Power(gA(m,n),4)
    /*2*/  * pow3(gB(m,n)) * r(m,n) * R(m,n) * (b_2 + 2 
    /*3*/  * R(m,n) * b_3) + 4 * gAlp(m,n) * Power(gA(m,n),4) 
    /*2*/  * gB_r(m,n) * Lt(m,n) - 4 * gAlp(m,n) * gsig(m,n) 
    /*2*/  * gB_r(m,n) * pow2(gA(m,n)) * pow2(gB(m,n)) 
    /*2*/  * pow2(r(m,n)) * Lt(m,n) - 2 * gAlp(m,n) * gL(m,n) 
    /*2*/  * Power(gA(m,n),4) * pow3(gB(m,n)) * Lt(m,n) - 4 
    /*2*/  * gAlp(m,n) * gA(m,n) * gA_r(m,n) * pow3(gB(m,n)) 
    /*2*/  * Lt(m,n) + 4 * gconf_r(m,n) * gAlp(m,n) 
    /*2*/  * pow2(gA(m,n)) * pow3(gB(m,n)) * Lt(m,n) + 2 
    /*2*/  * gAlp_r(m,n) * pow2(gA(m,n)) * pow3(gB(m,n)) 
    /*2*/  * Lt(m,n) - 6 * gAlp(m,n) * gA(m,n) * gA_r(m,n) 
    /*2*/  * gB_r(m,n) * pow2(gB(m,n)) * r(m,n) * Lt(m,n) + 4 
    /*2*/  * gconf_r(m,n) * gAlp(m,n) * gB_r(m,n) 
    /*2*/  * pow2(gA(m,n)) * pow2(gB(m,n)) * r(m,n) * Lt(m,n) 
    /*2*/  + 2 * gAlp_r(m,n) * gB_r(m,n) * pow2(gA(m,n)) 
    /*2*/  * pow2(gB(m,n)) * r(m,n) * Lt(m,n) + 2 * gAlp(m,n) 
    /*2*/  * gB_rr(m,n) * pow2(gA(m,n)) * pow2(gB(m,n)) 
    /*2*/  * r(m,n) * Lt(m,n) - 2 * gAlp(m,n) * gB(m,n) 
    /*2*/  * pow2(gA(m,n)) * pow2(gB_r(m,n)) * r(m,n) * Lt(m,n)
    /*2*/  + 2 * gAlp(m,n) * gL_r(m,n) * Power(gA(m,n),4) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * Lt(m,n) + 4 
    /*2*/  * gconf_r(m,n) * gAlp(m,n) * gA(m,n) * gA_r(m,n) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * Lt(m,n) + 2 * gAlp_r(m,n)
    /*2*/  * gA(m,n) * gA_r(m,n) * pow3(gB(m,n)) * r(m,n) 
    /*2*/  * Lt(m,n) - 2 * gAlp(m,n) * gA(m,n) * gA_rr(m,n) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * Lt(m,n) - 4 
    /*2*/  * gconf_rr(m,n) * gAlp(m,n) * pow2(gA(m,n)) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * Lt(m,n) + 8 
    /*2*/  * gconf_r(m,n) * gAlp_r(m,n) * pow2(gA(m,n)) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * Lt(m,n) - 2 
    /*2*/  * gAlp_rr(m,n) * pow2(gA(m,n)) * pow3(gB(m,n)) 
    /*2*/  * r(m,n) * Lt(m,n) - 6 * gAlp(m,n) * gsig(m,n) 
    /*2*/  * pow2(gA(m,n)) * pow3(gB(m,n)) * r(m,n) * Lt(m,n) 
    /*2*/  + 8 * gAlp(m,n) * pow2(gconf_r(m,n)) * pow2(gA(m,n))
    /*2*/  * pow3(gB(m,n)) * r(m,n) * Lt(m,n) + 6 * gAlp(m,n)
    /*2*/  * pow2(gA_r(m,n)) * pow3(gB(m,n)) * r(m,n) 
    /*2*/  * Lt(m,n) + k_g * exp(4 * gconf(m,n)) * gAlp(m,n) 
    /*2*/  * Power(gA(m,n),4) * pow3(gB(m,n)) * r(m,n) * R(m,n)
    /*2*/  * b_1 * Lt(m,n) + 2 * k_g * exp(4 * gconf(m,n)) 
    /*2*/  * gAlp(m,n) * Power(gA(m,n),4) * pow2(R(m,n)) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * b_2 * Lt(m,n) - k_g 
    /*2*/  * exp(2 * (fconf(m,n) + gconf(m,n))) * fA(m,n) 
    /*2*/  * pow3(gA(m,n)) * pow3(gB(m,n)) * r(m,n) 
    /*2*/  * (gAlp(m,n) * (2 * (-2 + pow2(Lt(m,n))) 
    /*4*/  * P_1_1(R(m,n)) - 3 * (-1 + pow2(Lt(m,n))) 
    /*4*/  * P_2_1(R(m,n)) - b_1) - fAlp(m,n) * (2 
    /*4*/  * P_1_2(R(m,n)) + b_2) * Lt(m,n)))) / (3. 
    /*1*/  * Power(gA(m,n),4) * pow3(gB(m,n)) * r(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[10][2]=
	(-2 * k_g * exp(-2 * gconf(m,n)) * (-(exp(2 
    /*4*/  * fconf(m,n)) * fA(m,n) * (gAlp(m,n) * (2 * (-2 
    /*6*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) - 3 * (-1 
    /*6*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)) - b_1) - fAlp(m,n)
    /*4*/  * (2 * P_1_2(R(m,n)) + b_2) * Lt(m,n))) + exp(2 
    /*3*/  * gconf(m,n)) * gA(m,n) * R(m,n) * (fAlp(m,n) * (b_2
    /*4*/  + 2 * R(m,n) * b_3) + gAlp(m,n) * (b_1 + 2 
    /*4*/  * R(m,n) * b_2) * Lt(m,n)))) / (3. * gA(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[10][3]=
	(gAlp(m,n) * (3 * gA2(m,n) - gtrA(m,n))) / 3.;

Jacobian[10][5]=
	(exp(-4 * gconf(m,n)) * (gA(m,n) * gB(m,n) * (3 
    /*3*/  * gAlp_r(m,n) * gA_r(m,n) * gB(m,n) * r(m,n) + 2 
    /*3*/  * gA(m,n) * (gAlp_r(m,n) * gB_r(m,n) * r(m,n) 
    /*4*/  + gB(m,n) * (-(gAlp_rr(m,n) * r(m,n)) + gAlp_r(m,n)
    /*5*/  * (1 + 4 * gconf_r(m,n) * r(m,n))))) * Lt(m,n) 
    /*2*/  + k_g * exp(2 * (fconf(m,n) + gconf(m,n))) * fA(m,n)
    /*2*/  * pow2(gB(m,n)) * pow3(gA(m,n)) * r(m,n) 
    /*2*/  * (gAlp(m,n) * (P_1_1(R(m,n)) + (-1 + pow2(Lt(m,n)))
    /*4*/  * P_2_1(R(m,n))) + fAlp(m,n) * P_1_2(R(m,n)) 
    /*3*/  * Lt(m,n)) - gAlp(m,n) * (3 * exp(2 * gconf(m,n)) 
    /*3*/  * gA2_r(m,n) * p(m,n) * pow2(gB(m,n)) 
    /*3*/  * pow3(gA(m,n)) * r(m,n) - 12 * pow2(gA_r(m,n)) 
    /*3*/  * pow2(gB(m,n)) * r(m,n) * Lt(m,n) + 2 * gsig(m,n) 
    /*3*/  * gB(m,n) * pow2(gA(m,n)) * r(m,n) * (3 * gB(m,n) 
    /*4*/  + 2 * gB_r(m,n) * r(m,n)) * Lt(m,n) - 2 
    /*3*/  * pow2(gA(m,n)) * (gB(m,n) * (2 * gconf_r(m,n) 
    /*5*/  * gB_r(m,n) + gB_rr(m,n)) * r(m,n) - pow2(gB_r(m,n))
    /*4*/  * r(m,n) + 2 * pow2(gB(m,n)) * (gconf_r(m,n) 
    /*5*/  - gconf_rr(m,n) * r(m,n) + 2 * pow2(gconf_r(m,n)) 
    /*5*/  * r(m,n))) * Lt(m,n) + 3 * gA(m,n) * gB(m,n) * (3 
    /*4*/  * gA_r(m,n) * gB_r(m,n) * r(m,n) + gB(m,n) 
    /*4*/  * (gA_rr(m,n) * r(m,n) + gA_r(m,n) * (2 - 2 
    /*6*/  * gconf_r(m,n) * r(m,n)))) * Lt(m,n)))) / (3. 
    /*1*/  * Power(gA(m,n),5) * pow2(gB(m,n)) * r(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[10][6]=
	(exp(-4 * gconf(m,n)) * (k_g * exp(4 * gconf(m,n))
    /*2*/  * fAlp(m,n) * pow3(gA(m,n)) * pow3(gB(m,n)) 
    /*2*/  * r(m,n) * R(m,n) * (b_2 + 2 * R(m,n) * b_3) 
    /*2*/  + (gAlp_r(m,n) * gA(m,n) * gB_r(m,n) * pow2(gB(m,n))
    /*3*/  * r(m,n) + gAlp(m,n) * (-3 * gA_r(m,n) * gB_r(m,n)
    /*4*/  * pow2(gB(m,n)) * r(m,n) + gA(m,n) * gB(m,n) 
    /*4*/  * r(m,n) * (gB(m,n) * (2 * gconf_r(m,n) * gB_r(m,n)
    /*6*/  + gB_rr(m,n)) - 2 * pow2(gB_r(m,n)) - 2 
    /*5*/  * gsig(m,n) * gB(m,n) * gB_r(m,n) * r(m,n)) 
    /*4*/  + pow3(gA(m,n)) * (6 * gB_r(m,n) + k_g * exp(4 
    /*6*/  * gconf(m,n)) * pow3(gB(m,n)) * r(m,n) * R(m,n) 
    /*5*/  * (b_1 + 2 * R(m,n) * b_2)))) * Lt(m,n) + k_g 
    /*2*/  * exp(2 * (fconf(m,n) + gconf(m,n))) * fA(m,n) 
    /*2*/  * pow2(gA(m,n)) * pow3(gB(m,n)) * r(m,n) * R(m,n) 
    /*2*/  * (gAlp(m,n) * (2 * (-1 + pow2(Lt(m,n))) 
    /*4*/  * P_1_2(R(m,n)) - b_2) - fAlp(m,n) * b_3 
    /*3*/  * Lt(m,n)))) / (3. * Power(gB(m,n),4) 
    /*1*/  * pow3(gA(m,n)) * r(m,n) * Lt(m,n));

Jacobian[10][7]=
	-(k_g * exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*1*/  * (gAlp(m,n) * (P_1_1(R(m,n)) + (-1 + pow2(Lt(m,n)))
    /*3*/  * P_2_1(R(m,n))) + fAlp(m,n) * P_1_2(R(m,n)) 
    /*2*/  * Lt(m,n))) / (3. * gA(m,n) * Lt(m,n));

Jacobian[10][8]=
	(k_g * exp(2 * fconf(m,n) - 4 * gconf(m,n)) 
    /*1*/  * (exp(2 * gconf(m,n)) * gA(m,n) * (fAlp(m,n) * (2 
    /*4*/  * P_1_2(R(m,n)) + b_2) + gAlp(m,n) * (2 
    /*4*/  * P_1_1(R(m,n)) + b_1) * Lt(m,n)) - exp(2 
    /*3*/  * fconf(m,n)) * fA(m,n) * (gAlp(m,n) * (2 * (-1 
    /*5*/  + pow2(Lt(m,n))) * P_1_2(R(m,n)) - b_2) - fAlp(m,n)
    /*3*/  * b_3 * Lt(m,n)))) / (3. * gA(m,n) * gB(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[10][10]=
	gAlp(m,n) * (gtrA(m,n) + gtrK(m,n));

Jacobian[10][13]=
	(exp(-4 * gconf(m,n)) * gAlp(m,n)) / (3. * r(m,n));

Jacobian[10][15]=
	(exp(-4 * gconf(m,n)) * gAlp(m,n) * (3 * gB(m,n) 
    /*2*/  + 2 * gB_r(m,n) * r(m,n))) / (3. * gB(m,n) 
    /*1*/  * pow2(gA(m,n)));

Jacobian[10][22]=
	gA2_r(m,n);

Jacobian[11][1]=
	(-4 * k_f * exp(-2 * fconf(m,n)) * (-(exp(2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (fAlp(m,n) * (2 * (-1 
    /*6*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) + pow2(Lt(m,n)) 
    /*5*/  * P_2_1(R(m,n)) + b_1) + gAlp(m,n) * (P_2_0(R(m,n))
    /*5*/  + b_0) * Lt(m,n))) + exp(2 * fconf(m,n)) * fA(m,n)
    /*2*/  * (gAlp(m,n) * (P_1_1(R(m,n)) - b_1) + fAlp(m,n) 
    /*3*/  * (P_1_2(R(m,n)) - b_2) * Lt(m,n)))) / (3. * fA(m,n)
    /*1*/  * pow2(R(m,n)) * Lt(m,n));

Jacobian[11][2]=
	(-2 * exp(-4 * fconf(m,n) - 6 * gconf(m,n)) * (2 
    /*2*/  * exp(4 * gconf(m,n)) * fA(m,n) * pow2(gB(m,n)) 
    /*2*/  * (-(k_f * exp(4 * fconf(m,n) + 2 * gconf(m,n)) 
    /*4*/  * gAlp(m,n) * gB(m,n) * pow3(fA(m,n)) * r(m,n) 
    /*4*/  * (P_2_0(R(m,n)) + b_0)) + 2 * exp(2 * gconf(m,n)) 
    /*3*/  * fAlp_r(m,n) * fA_r(m,n) * gB(m,n) * pow3(R(m,n)) 
    /*3*/  * r(m,n) * Lt(m,n) + 2 * fA(m,n) * pow2(R(m,n)) 
    /*3*/  * (exp(2 * fconf(m,n)) * fAlp_r(m,n) * fB_r(m,n) 
    /*4*/  * r(m,n) + exp(2 * gconf(m,n)) * gB(m,n) 
    /*4*/  * (-(fAlp_rr(m,n) * r(m,n)) + fAlp_r(m,n) * (1 + 4 
    /*6*/  * fconf_r(m,n) * r(m,n))) * R(m,n)) * Lt(m,n) + k_f
    /*3*/  * exp(2 * fconf(m,n) + 4 * gconf(m,n)) * gAlp(m,n)
    /*3*/  * gA(m,n) * gB(m,n) * pow2(fA(m,n)) * r(m,n) 
    /*3*/  * R(m,n) * (P_2_0(R(m,n)) + b_0) * Lt(m,n)) 
    /*2*/  + fAlp(m,n) * (exp(2 * fconf(m,n) + 6 * gconf(m,n))
    /*3*/  * pow3(fA(m,n)) * pow3(gB(m,n)) * r(m,n) * (-3 
    /*4*/  * fA1_r(m,n) * p(m,n) * pow3(R(m,n)) + 2 * k_f 
    /*4*/  * exp(2 * gconf(m,n)) * gA(m,n) * (P_1_0(R(m,n)) - 2
    /*5*/  * P_2_0(R(m,n)) + 2 * pow2(Lt(m,n)) 
    /*5*/  * P_1_1(R(m,n)) * R(m,n) + pow2(Lt(m,n)) 
    /*5*/  * P_2_1(R(m,n)) * R(m,n) - b_0)) + 12 * exp(6 
    /*4*/  * gconf(m,n)) * pow2(fA_r(m,n)) * pow3(gB(m,n)) 
    /*3*/  * pow3(R(m,n)) * r(m,n) * Lt(m,n) + 4 * exp(4 
    /*4*/  * gconf(m,n)) * fA(m,n) * pow2(gB(m,n)) 
    /*3*/  * pow2(R(m,n)) * (-3 * exp(2 * fconf(m,n)) 
    /*4*/  * fA_r(m,n) * fB_r(m,n) * r(m,n) + exp(2 
    /*5*/  * gconf(m,n)) * gB(m,n) * (-(fA_rr(m,n) * r(m,n)) 
    /*5*/  + 2 * fA_r(m,n) * (-1 + fconf_r(m,n) * r(m,n))) 
    /*4*/  * R(m,n)) * Lt(m,n) - 4 * exp(2 * gconf(m,n)) 
    /*3*/  * gB(m,n) * pow2(fA(m,n)) * R(m,n) * (exp(4 
    /*5*/  * fconf(m,n)) * pow2(fB_r(m,n)) * r(m,n) - 2 * exp(4
    /*5*/  * gconf(m,n)) * pow2(gB(m,n)) * pow2(R(m,n)) 
    /*4*/  * (fconf_r(m,n) - fconf_rr(m,n) * r(m,n) + 2 
    /*5*/  * pow2(fconf_r(m,n)) * r(m,n)) - exp(2 * (fconf(m,n)
    /*6*/  + gconf(m,n))) * (2 * fconf_r(m,n) * fB_r(m,n) 
    /*5*/  + fB_rr(m,n)) * gB(m,n) * r(m,n) * R(m,n) + exp(2 
    /*5*/  * gconf(m,n)) * fsig(m,n) * gB(m,n) * r(m,n) 
    /*4*/  * R(m,n) * (2 * exp(2 * fconf(m,n)) * fB_r(m,n) 
    /*5*/  * r(m,n) + 3 * exp(2 * gconf(m,n)) * gB(m,n) 
    /*5*/  * R(m,n))) * Lt(m,n) - 2 * Power(fA(m,n),4) * (-4 
    /*4*/  * exp(6 * fconf(m,n)) * fB_r(m,n) + exp(6 
    /*5*/  * gconf(m,n)) * pow3(gB(m,n)) * (2 * fL(m,n) 
    /*5*/  * pow3(R(m,n)) - 2 * fL_r(m,n) * pow3(R(m,n)) 
    /*5*/  * r(m,n) + k_f * exp(4 * fconf(m,n)) * P_2_1(R(m,n))
    /*5*/  * r(m,n) + k_f * exp(4 * fconf(m,n)) * r(m,n) 
    /*5*/  * b_1)) * Lt(m,n)))) / (3. * Power(fA(m,n),4) 
    /*1*/  * pow3(gB(m,n)) * pow3(R(m,n)) * r(m,n) * Lt(m,n));

Jacobian[11][4]=
	(fAlp(m,n) * (3 * fA1(m,n) - ftrA(m,n))) / 3.;

Jacobian[11][5]=
	(2 * k_f * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*1*/  * (fAlp(m,n) * (-P_1_1(R(m,n)) + pow2(Lt(m,n)) 
    /*3*/  * P_2_1(R(m,n))) + gAlp(m,n) * (-P_1_0(R(m,n)) 
    /*3*/  + P_2_0(R(m,n))) * Lt(m,n))) / (3. * fA(m,n) 
    /*1*/  * pow2(R(m,n)) * Lt(m,n));

Jacobian[11][6]=
	(-2 * k_f * exp(-2 * fconf(m,n)) * (-(exp(2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (fAlp(m,n) * ((-1 + 2 
    /*6*/  * pow2(Lt(m,n))) * P_1_1(R(m,n)) + b_1) + gAlp(m,n)
    /*4*/  * (P_1_0(R(m,n)) + b_0) * Lt(m,n))) + exp(2 
    /*3*/  * fconf(m,n)) * fA(m,n) * (gAlp(m,n) 
    /*3*/  * (P_1_1(R(m,n)) - b_1) + fAlp(m,n) * (P_1_2(R(m,n))
    /*4*/  - b_2) * Lt(m,n)))) / (3. * fA(m,n) * gB(m,n) 
    /*1*/  * pow2(R(m,n)) * Lt(m,n));

Jacobian[11][7]=
	(exp(-4 * (fconf(m,n) + gconf(m,n))) * (-2 * exp(2
    /*3*/  * gconf(m,n)) * fA(m,n) * gB(m,n) * R(m,n) * (k_f
    /*3*/  * exp(2 * fconf(m,n) + 4 * gconf(m,n)) * gAlp(m,n)
    /*3*/  * gA(m,n) * gB(m,n) * pow2(fA(m,n)) 
    /*3*/  * P_1_1(R(m,n)) * r(m,n) + 3 * exp(2 * gconf(m,n)) 
    /*3*/  * fAlp_r(m,n) * fA_r(m,n) * gB(m,n) * r(m,n) 
    /*3*/  * R(m,n) + 2 * fA(m,n) * (exp(2 * fconf(m,n)) 
    /*4*/  * fAlp_r(m,n) * fB_r(m,n) * r(m,n) + exp(2 
    /*5*/  * gconf(m,n)) * gB(m,n) * (-(fAlp_rr(m,n) * r(m,n))
    /*5*/  + fAlp_r(m,n) * (1 + 4 * fconf_r(m,n) * r(m,n))) 
    /*4*/  * R(m,n))) * Lt(m,n) + fAlp(m,n) * (exp(2 
    /*4*/  * fconf(m,n) + 4 * gconf(m,n)) * pow2(gB(m,n)) 
    /*3*/  * pow3(fA(m,n)) * (3 * fA1_r(m,n) * p(m,n) 
    /*4*/  * pow2(R(m,n)) + 2 * k_f * exp(2 * gconf(m,n)) 
    /*4*/  * gA(m,n) * (P_1_1(R(m,n)) - pow2(Lt(m,n)) 
    /*5*/  * P_2_1(R(m,n)))) * r(m,n) - 24 * exp(4 
    /*4*/  * gconf(m,n)) * pow2(fA_r(m,n)) * pow2(gB(m,n)) 
    /*3*/  * pow2(R(m,n)) * r(m,n) * Lt(m,n) - 6 * exp(2 
    /*4*/  * gconf(m,n)) * fA(m,n) * gB(m,n) * R(m,n) * (-3 
    /*4*/  * exp(2 * fconf(m,n)) * fA_r(m,n) * fB_r(m,n) 
    /*4*/  * r(m,n) + exp(2 * gconf(m,n)) * gB(m,n) 
    /*4*/  * (-(fA_rr(m,n) * r(m,n)) + 2 * fA_r(m,n) * (-1 
    /*6*/  + fconf_r(m,n) * r(m,n))) * R(m,n)) * Lt(m,n) + 4 
    /*3*/  * pow2(fA(m,n)) * (exp(4 * fconf(m,n)) 
    /*4*/  * pow2(fB_r(m,n)) * r(m,n) - 2 * exp(4 * gconf(m,n))
    /*4*/  * pow2(gB(m,n)) * pow2(R(m,n)) * (fconf_r(m,n) 
    /*5*/  - fconf_rr(m,n) * r(m,n) + 2 * pow2(fconf_r(m,n)) 
    /*5*/  * r(m,n)) - exp(2 * (fconf(m,n) + gconf(m,n))) * (2
    /*5*/  * fconf_r(m,n) * fB_r(m,n) + fB_rr(m,n)) * gB(m,n)
    /*4*/  * r(m,n) * R(m,n) + exp(2 * gconf(m,n)) 
    /*4*/  * fsig(m,n) * gB(m,n) * r(m,n) * R(m,n) * (2 * exp(2
    /*6*/  * fconf(m,n)) * fB_r(m,n) * r(m,n) + 3 * exp(2 
    /*6*/  * gconf(m,n)) * gB(m,n) * R(m,n))) * Lt(m,n)))) 
    /*0*/  / (3. * Power(fA(m,n),5) * pow2(gB(m,n)) 
    /*1*/  * pow2(R(m,n)) * r(m,n) * Lt(m,n));

Jacobian[11][8]=
	(2 * exp(-8 * gconf(m,n)) * (3 * exp(4 
    /*3*/  * gconf(m,n)) * fAlp(m,n) * fA_r(m,n) * fB_r(m,n) 
    /*2*/  * pow2(gB(m,n)) * pow2(R(m,n)) * r(m,n) * Lt(m,n) 
    /*2*/  + exp(2 * gconf(m,n)) * fA(m,n) * gB(m,n) * r(m,n) 
    /*2*/  * R(m,n) * (-(exp(2 * gconf(m,n)) * fAlp_r(m,n) 
    /*4*/  * fB_r(m,n) * gB(m,n) * R(m,n)) + fAlp(m,n) * (2 
    /*4*/  * exp(2 * fconf(m,n)) * pow2(fB_r(m,n)) - exp(2 
    /*5*/  * gconf(m,n)) * (2 * fconf_r(m,n) * fB_r(m,n) 
    /*5*/  + fB_rr(m,n)) * gB(m,n) * R(m,n) + 2 * exp(2 
    /*5*/  * gconf(m,n)) * fsig(m,n) * fB_r(m,n) * gB(m,n) 
    /*4*/  * r(m,n) * R(m,n))) * Lt(m,n) - k_f * exp(8 
    /*3*/  * gconf(m,n)) * gA(m,n) * pow2(fA(m,n)) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * (fAlp(m,n) * (-2 
    /*4*/  * pow2(Lt(m,n)) * P_1_0(R(m,n)) + (-1 + 2 
    /*5*/  * pow2(Lt(m,n))) * P_2_0(R(m,n)) - b_0) + gAlp(m,n)
    /*3*/  * R(m,n) * (P_1_0(R(m,n)) + b_0) * Lt(m,n)) 
    /*2*/  + exp(2 * fconf(m,n)) * pow3(fA(m,n)) * (k_f * exp(6
    /*4*/  * gconf(m,n)) * gAlp(m,n) * pow3(gB(m,n)) * r(m,n)
    /*3*/  * (P_2_0(R(m,n)) + b_0) + fAlp(m,n) * (-6 * exp(2
    /*5*/  * fconf(m,n)) * fB_r(m,n) + k_f * exp(6 
    /*5*/  * gconf(m,n)) * pow3(gB(m,n)) * r(m,n) 
    /*4*/  * (P_2_1(R(m,n)) + b_1)) * Lt(m,n)))) / (3. 
    /*1*/  * Power(gB(m,n),4) * pow3(fA(m,n)) * r(m,n) 
    /*1*/  * Power(R(m,n),4) * Lt(m,n));

Jacobian[11][11]=
	fAlp(m,n) * (ftrA(m,n) + ftrK(m,n));

Jacobian[11][14]=
	(-2 * exp(-4 * fconf(m,n)) * fAlp(m,n)) / (3. 
    /*1*/  * r(m,n));

Jacobian[11][16]=
	(2 * exp(-4 * fconf(m,n)) * fAlp(m,n) * (-3 - (2 
    /*3*/  * exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fB_r(m,n) 
    /*3*/  * r(m,n)) / (gB(m,n) * R(m,n)))) / (3. 
    /*1*/  * pow2(fA(m,n)));

Jacobian[11][22]=
	fA1_r(m,n);

Jacobian[12][1]=
	(2 * k_f * exp(-2 * fconf(m,n)) * (-(exp(2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (fAlp(m,n) * (2 * (-1 
    /*6*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) + pow2(Lt(m,n)) 
    /*5*/  * P_2_1(R(m,n)) + b_1) + gAlp(m,n) * (P_2_0(R(m,n))
    /*5*/  + b_0) * Lt(m,n))) + exp(2 * fconf(m,n)) * fA(m,n)
    /*2*/  * (gAlp(m,n) * (P_1_1(R(m,n)) - b_1) + fAlp(m,n) 
    /*3*/  * (P_1_2(R(m,n)) - b_2) * Lt(m,n)))) / (3. * fA(m,n)
    /*1*/  * pow2(R(m,n)) * Lt(m,n));

Jacobian[12][2]=
	(2 * exp(-4 * fconf(m,n) - 6 * gconf(m,n)) 
    /*1*/  * (exp(4 * gconf(m,n)) * fA(m,n) * pow2(gB(m,n)) 
    /*2*/  * (-(k_f * exp(4 * fconf(m,n) + 2 * gconf(m,n)) 
    /*4*/  * gAlp(m,n) * gB(m,n) * pow3(fA(m,n)) * r(m,n) 
    /*4*/  * (P_2_0(R(m,n)) + b_0)) + 2 * exp(2 * gconf(m,n)) 
    /*3*/  * fAlp_r(m,n) * fA_r(m,n) * gB(m,n) * pow3(R(m,n)) 
    /*3*/  * r(m,n) * Lt(m,n) + 2 * fA(m,n) * pow2(R(m,n)) 
    /*3*/  * (exp(2 * fconf(m,n)) * fAlp_r(m,n) * fB_r(m,n) 
    /*4*/  * r(m,n) + exp(2 * gconf(m,n)) * gB(m,n) 
    /*4*/  * (-(fAlp_rr(m,n) * r(m,n)) + fAlp_r(m,n) * (1 + 4 
    /*6*/  * fconf_r(m,n) * r(m,n))) * R(m,n)) * Lt(m,n) + k_f
    /*3*/  * exp(2 * fconf(m,n) + 4 * gconf(m,n)) * gAlp(m,n)
    /*3*/  * gA(m,n) * gB(m,n) * pow2(fA(m,n)) * r(m,n) 
    /*3*/  * R(m,n) * (P_2_0(R(m,n)) + b_0) * Lt(m,n)) 
    /*2*/  + fAlp(m,n) * (exp(2 * fconf(m,n) + 6 * gconf(m,n))
    /*3*/  * pow3(fA(m,n)) * pow3(gB(m,n)) * r(m,n) * (3 
    /*4*/  * fA2_r(m,n) * p(m,n) * pow3(R(m,n)) + k_f * exp(2 
    /*5*/  * gconf(m,n)) * gA(m,n) * (P_1_0(R(m,n)) - 2 
    /*5*/  * P_2_0(R(m,n)) + 2 * pow2(Lt(m,n)) * P_1_1(R(m,n))
    /*5*/  * R(m,n) + pow2(Lt(m,n)) * P_2_1(R(m,n)) * R(m,n)
    /*5*/  - b_0)) + 6 * exp(6 * gconf(m,n)) 
    /*3*/  * pow2(fA_r(m,n)) * pow3(gB(m,n)) * pow3(R(m,n)) 
    /*3*/  * r(m,n) * Lt(m,n) + 2 * exp(4 * gconf(m,n)) 
    /*3*/  * fA(m,n) * pow2(gB(m,n)) * pow2(R(m,n)) * (-3 
    /*4*/  * exp(2 * fconf(m,n)) * fA_r(m,n) * fB_r(m,n) 
    /*4*/  * r(m,n) + exp(2 * gconf(m,n)) * gB(m,n) 
    /*4*/  * (-(fA_rr(m,n) * r(m,n)) + 2 * fA_r(m,n) * (-1 
    /*6*/  + fconf_r(m,n) * r(m,n))) * R(m,n)) * Lt(m,n) - 2 
    /*3*/  * exp(2 * gconf(m,n)) * gB(m,n) * pow2(fA(m,n)) 
    /*3*/  * R(m,n) * (exp(4 * fconf(m,n)) * pow2(fB_r(m,n)) 
    /*4*/  * r(m,n) - 2 * exp(4 * gconf(m,n)) * pow2(gB(m,n)) 
    /*4*/  * pow2(R(m,n)) * (fconf_r(m,n) - fconf_rr(m,n) 
    /*5*/  * r(m,n) + 2 * pow2(fconf_r(m,n)) * r(m,n)) - exp(2
    /*5*/  * (fconf(m,n) + gconf(m,n))) * (2 * fconf_r(m,n) 
    /*5*/  * fB_r(m,n) + fB_rr(m,n)) * gB(m,n) * r(m,n) 
    /*4*/  * R(m,n) + exp(2 * gconf(m,n)) * fsig(m,n) * gB(m,n)
    /*4*/  * r(m,n) * R(m,n) * (2 * exp(2 * fconf(m,n)) 
    /*5*/  * fB_r(m,n) * r(m,n) + 3 * exp(2 * gconf(m,n)) 
    /*5*/  * gB(m,n) * R(m,n))) * Lt(m,n) - Power(fA(m,n),4) 
    /*3*/  * (-4 * exp(6 * fconf(m,n)) * fB_r(m,n) + exp(6 
    /*5*/  * gconf(m,n)) * pow3(gB(m,n)) * (2 * fL(m,n) 
    /*5*/  * pow3(R(m,n)) - 2 * fL_r(m,n) * pow3(R(m,n)) 
    /*5*/  * r(m,n) + k_f * exp(4 * fconf(m,n)) * P_2_1(R(m,n))
    /*5*/  * r(m,n) + k_f * exp(4 * fconf(m,n)) * r(m,n) 
    /*5*/  * b_1)) * Lt(m,n)))) / (3. * Power(fA(m,n),4) 
    /*1*/  * pow3(gB(m,n)) * pow3(R(m,n)) * r(m,n) * Lt(m,n));

Jacobian[12][4]=
	(fAlp(m,n) * (3 * fA2(m,n) - ftrA(m,n))) / 3.;

Jacobian[12][5]=
	(k_f * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*1*/  * (fAlp(m,n) * (P_1_1(R(m,n)) - pow2(Lt(m,n)) 
    /*3*/  * P_2_1(R(m,n))) + gAlp(m,n) * (P_1_0(R(m,n)) 
    /*3*/  - P_2_0(R(m,n))) * Lt(m,n))) / (3. * fA(m,n) 
    /*1*/  * pow2(R(m,n)) * Lt(m,n));

Jacobian[12][6]=
	(k_f * exp(-2 * fconf(m,n)) * (-(exp(2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (fAlp(m,n) * ((-1 + 2 
    /*6*/  * pow2(Lt(m,n))) * P_1_1(R(m,n)) + b_1) + gAlp(m,n)
    /*4*/  * (P_1_0(R(m,n)) + b_0) * Lt(m,n))) + exp(2 
    /*3*/  * fconf(m,n)) * fA(m,n) * (gAlp(m,n) 
    /*3*/  * (P_1_1(R(m,n)) - b_1) + fAlp(m,n) * (P_1_2(R(m,n))
    /*4*/  - b_2) * Lt(m,n)))) / (3. * fA(m,n) * gB(m,n) 
    /*1*/  * pow2(R(m,n)) * Lt(m,n));

Jacobian[12][7]=
	(exp(-4 * (fconf(m,n) + gconf(m,n))) * (exp(2 
    /*3*/  * gconf(m,n)) * fA(m,n) * gB(m,n) * R(m,n) * (k_f 
    /*3*/  * exp(2 * fconf(m,n) + 4 * gconf(m,n)) * gAlp(m,n) 
    /*3*/  * gA(m,n) * gB(m,n) * pow2(fA(m,n)) * P_1_1(R(m,n))
    /*3*/  * r(m,n) + 3 * exp(2 * gconf(m,n)) * fAlp_r(m,n) 
    /*3*/  * fA_r(m,n) * gB(m,n) * r(m,n) * R(m,n) + 2 
    /*3*/  * fA(m,n) * (exp(2 * fconf(m,n)) * fAlp_r(m,n) 
    /*4*/  * fB_r(m,n) * r(m,n) + exp(2 * gconf(m,n)) * gB(m,n)
    /*4*/  * (-(fAlp_rr(m,n) * r(m,n)) + fAlp_r(m,n) * (1 + 4
    /*6*/  * fconf_r(m,n) * r(m,n))) * R(m,n))) * Lt(m,n) 
    /*2*/  + fAlp(m,n) * (-(exp(2 * fconf(m,n) + 4 
    /*5*/  * gconf(m,n)) * pow2(gB(m,n)) * pow3(fA(m,n)) * (-3
    /*5*/  * fA2_r(m,n) * p(m,n) * pow2(R(m,n)) + k_f * exp(2
    /*6*/  * gconf(m,n)) * gA(m,n) * (P_1_1(R(m,n)) 
    /*6*/  - pow2(Lt(m,n)) * P_2_1(R(m,n)))) * r(m,n)) + 12 
    /*3*/  * exp(4 * gconf(m,n)) * pow2(fA_r(m,n)) 
    /*3*/  * pow2(gB(m,n)) * pow2(R(m,n)) * r(m,n) * Lt(m,n) 
    /*3*/  + 3 * exp(2 * gconf(m,n)) * fA(m,n) * gB(m,n) 
    /*3*/  * R(m,n) * (-3 * exp(2 * fconf(m,n)) * fA_r(m,n) 
    /*4*/  * fB_r(m,n) * r(m,n) + exp(2 * gconf(m,n)) * gB(m,n)
    /*4*/  * (-(fA_rr(m,n) * r(m,n)) + 2 * fA_r(m,n) * (-1 
    /*6*/  + fconf_r(m,n) * r(m,n))) * R(m,n)) * Lt(m,n) - 2 
    /*3*/  * pow2(fA(m,n)) * (exp(4 * fconf(m,n)) 
    /*4*/  * pow2(fB_r(m,n)) * r(m,n) - 2 * exp(4 * gconf(m,n))
    /*4*/  * pow2(gB(m,n)) * pow2(R(m,n)) * (fconf_r(m,n) 
    /*5*/  - fconf_rr(m,n) * r(m,n) + 2 * pow2(fconf_r(m,n)) 
    /*5*/  * r(m,n)) - exp(2 * (fconf(m,n) + gconf(m,n))) * (2
    /*5*/  * fconf_r(m,n) * fB_r(m,n) + fB_rr(m,n)) * gB(m,n)
    /*4*/  * r(m,n) * R(m,n) + exp(2 * gconf(m,n)) 
    /*4*/  * fsig(m,n) * gB(m,n) * r(m,n) * R(m,n) * (2 * exp(2
    /*6*/  * fconf(m,n)) * fB_r(m,n) * r(m,n) + 3 * exp(2 
    /*6*/  * gconf(m,n)) * gB(m,n) * R(m,n))) * Lt(m,n)))) 
    /*0*/  / (3. * Power(fA(m,n),5) * pow2(gB(m,n)) 
    /*1*/  * pow2(R(m,n)) * r(m,n) * Lt(m,n));

Jacobian[12][8]=
	(exp(-8 * gconf(m,n)) * (-3 * exp(4 * gconf(m,n))
    /*2*/  * fAlp(m,n) * fA_r(m,n) * fB_r(m,n) 
    /*2*/  * pow2(gB(m,n)) * pow2(R(m,n)) * r(m,n) * Lt(m,n) 
    /*2*/  - exp(2 * gconf(m,n)) * fA(m,n) * gB(m,n) * r(m,n) 
    /*2*/  * R(m,n) * (-(exp(2 * gconf(m,n)) * fAlp_r(m,n) 
    /*4*/  * fB_r(m,n) * gB(m,n) * R(m,n)) + fAlp(m,n) * (2 
    /*4*/  * exp(2 * fconf(m,n)) * pow2(fB_r(m,n)) - exp(2 
    /*5*/  * gconf(m,n)) * (2 * fconf_r(m,n) * fB_r(m,n) 
    /*5*/  + fB_rr(m,n)) * gB(m,n) * R(m,n) + 2 * exp(2 
    /*5*/  * gconf(m,n)) * fsig(m,n) * fB_r(m,n) * gB(m,n) 
    /*4*/  * r(m,n) * R(m,n))) * Lt(m,n) + k_f * exp(8 
    /*3*/  * gconf(m,n)) * gA(m,n) * pow2(fA(m,n)) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * (fAlp(m,n) * (-2 
    /*4*/  * pow2(Lt(m,n)) * P_1_0(R(m,n)) + (-1 + 2 
    /*5*/  * pow2(Lt(m,n))) * P_2_0(R(m,n)) - b_0) + gAlp(m,n)
    /*3*/  * R(m,n) * (P_1_0(R(m,n)) + b_0) * Lt(m,n)) 
    /*2*/  - exp(2 * fconf(m,n)) * pow3(fA(m,n)) * (k_f * exp(6
    /*4*/  * gconf(m,n)) * gAlp(m,n) * pow3(gB(m,n)) * r(m,n)
    /*3*/  * (P_2_0(R(m,n)) + b_0) + fAlp(m,n) * (-6 * exp(2
    /*5*/  * fconf(m,n)) * fB_r(m,n) + k_f * exp(6 
    /*5*/  * gconf(m,n)) * pow3(gB(m,n)) * r(m,n) 
    /*4*/  * (P_2_1(R(m,n)) + b_1)) * Lt(m,n)))) / (3. 
    /*1*/  * Power(gB(m,n),4) * pow3(fA(m,n)) * r(m,n) 
    /*1*/  * Power(R(m,n),4) * Lt(m,n));

Jacobian[12][12]=
	fAlp(m,n) * (ftrA(m,n) + ftrK(m,n));

Jacobian[12][14]=
	(exp(-4 * fconf(m,n)) * fAlp(m,n)) / (3. * r(m,n));

Jacobian[12][16]=
	(exp(-4 * fconf(m,n)) * fAlp(m,n) * (3 + (2 
    /*3*/  * exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fB_r(m,n) 
    /*3*/  * r(m,n)) / (gB(m,n) * R(m,n)))) / (3. 
    /*1*/  * pow2(fA(m,n)));

Jacobian[12][22]=
	fA2_r(m,n);

Jacobian[13][1]=
	2 * gAlp(m,n) * (-4 * k_g * exp(4 * gconf(m,n)) 
    /*1*/  * gj(m,n) - (exp(-2 * gconf(m,n)) * p(m,n) 
    /*2*/  * (gA(m,n) * (-2 + gL_r(m,n) * pow2(gB(m,n)) 
    /*4*/  * pow2(r(m,n))) + 4 * k_g * exp(2 * (fconf(m,n) 
    /*5*/  + gconf(m,n))) * fA(m,n) * pow2(gB(m,n)) 
    /*3*/  * pow2(r(m,n)) * P_1_2(R(m,n)) * R(m,n) * Lt(m,n)))
    /*1*/  / (pow2(gA(m,n)) * pow2(gB(m,n)) * pow2(r(m,n)) 
    /*2*/  * Lt(m,n)));

Jacobian[13][2]=
	(-4 * k_g * exp(2 * fconf(m,n)) * fA(m,n) 
    /*1*/  * gAlp(m,n) * p(m,n) * (2 * P_1_1(R(m,n)) - 3 
    /*2*/  * P_2_1(R(m,n)))) / pow2(gA(m,n));

Jacobian[13][5]=
	(exp(-2 * gconf(m,n)) * (gdet(m,n) 
    /*2*/  * gdet_pff_r(m,n) * exp(2 * gconf(m,n)) * gA(m,n) 
    /*2*/  * pow2(gB(m,n)) * pow2(r(m,n)) * Lt(m,n) 
    /*2*/  - gdet_pff(m,n) * gdet_r(m,n) * exp(2 * gconf(m,n))
    /*2*/  * gA(m,n) * pow2(gB(m,n)) * pow2(r(m,n)) * Lt(m,n)
    /*2*/  - pow2(gdet(m,n)) * (-2 * exp(2 * gconf(m,n)) 
    /*3*/  * gA(m,n) * pow2(gB(m,n)) * pow2(r(m,n)) * (6 
    /*4*/  * gA1(m,n) * gAlp_r(m,n) - 3 * gBet_rr(m,n) - 2 
    /*4*/  * gAlp_r(m,n) * gtrA(m,n)) * Lt(m,n) + gAlp(m,n) 
    /*3*/  * (3 * p(m,n) * pow2(gA(m,n)) * (-2 + gL_r(m,n) 
    /*5*/  * pow2(gB(m,n)) * pow2(r(m,n))) + 12 * exp(2 
    /*5*/  * gconf(m,n)) * gAsig(m,n) * gA_r(m,n) 
    /*4*/  * pow2(gB(m,n)) * Power(r(m,n),4) * Lt(m,n) + 4 
    /*4*/  * exp(2 * gconf(m,n)) * gA(m,n) * gB(m,n) 
    /*4*/  * pow2(r(m,n)) * (2 * gAsig(m,n) * pow2(r(m,n)) * (6
    /*6*/  * gconf_r(m,n) * gB(m,n) + gB_r(m,n) + gsig(m,n) 
    /*6*/  * gB(m,n) * r(m,n)) + gB(m,n) * (3 * k_g * exp(2 
    /*7*/  * fconf(m,n)) * fA(m,n) * p(m,n) * P_2_1(R(m,n)) - 2
    /*6*/  * (gtrA_r(m,n) + gtrK_r(m,n)))) * Lt(m,n))))) 
    /*0*/  / (3. * Power(gA(m,n),4) * pow2(gdet(m,n)) 
    /*1*/  * pow2(gB(m,n)) * pow2(r(m,n)) * Lt(m,n));

Jacobian[13][6]=
	(-4 * exp(-2 * gconf(m,n)) * (-3 * exp(2 
    /*3*/  * gconf(m,n)) * Lt(m,n) * pow2(gA(m,n)) * (gBet(m,n)
    /*3*/  - gBet_r(m,n) * r(m,n)) * Lt(m,n) + gAlp(m,n) 
    /*2*/  * (exp(2 * gconf(m,n)) * gB(m,n) * Lt(m,n) 
    /*3*/  * pow2(r(m,n)) * (gAsig(m,n) * gB_r(m,n) 
    /*4*/  * pow2(r(m,n)) + 3 * k_g * exp(2 * fconf(m,n)) 
    /*4*/  * fA(m,n) * gB(m,n) * p(m,n) * P_1_2(R(m,n)) 
    /*4*/  * R(m,n)) * Lt(m,n) + 3 * gA(m,n) * p(m,n) 
    /*3*/  * (-Lt(m,n) + Lt(m,n))))) / (3. * Lt(m,n) 
    /*1*/  * pow2(gA(m,n)) * pow2(r(m,n)) * pow3(gB(m,n)) 
    /*1*/  * Lt(m,n));

Jacobian[13][7]=
	(2 * k_g * exp(2 * fconf(m,n)) * gAlp(m,n) 
    /*1*/  * p(m,n) * P_2_1(R(m,n))) / pow2(gA(m,n));

Jacobian[13][8]=
	(4 * k_g * exp(4 * fconf(m,n) - 2 * gconf(m,n)) 
    /*1*/  * fA(m,n) * gAlp(m,n) * p(m,n) * P_1_2(R(m,n))) 
    /*0*/  / (gB(m,n) * pow2(gA(m,n)));

Jacobian[13][9]=
	(-2 * gAlp_r(m,n)) / pow2(gA(m,n));

Jacobian[13][13]=
	-gdet_pff(m,n) / (3. * gdet(m,n)) - gBet_r(m,n);

Jacobian[13][15]=
	(4 * gAsig(m,n) * gAlp(m,n) * pow3(r(m,n))) / (3.
    /*1*/  * pow2(gA(m,n)));

Jacobian[13][17]=
	(4 * gAlp(m,n) * pow2(r(m,n)) * (gA_r(m,n) 
    /*2*/  * gB(m,n) + gA(m,n) * (6 * gconf_r(m,n) * gB(m,n) 
    /*3*/  + gB_r(m,n)) + gsig(m,n) * gA(m,n) * gB(m,n) 
    /*2*/  * r(m,n))) / (3. * gB(m,n) * pow3(gA(m,n)));

Jacobian[13][22]=
	gL_r(m,n) - 2 / (pow2(gB(m,n)) * pow2(r(m,n)));

Jacobian[14][1]=
	(-4 * k_f * exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*1*/  * gA(m,n) * p(m,n) * (2 * P_1_1(R(m,n)) 
    /*2*/  + P_2_1(R(m,n)))) / (pow2(fA(m,n)) * pow2(R(m,n)));

Jacobian[14][2]=
	2 * fAlp(m,n) * (-4 * k_f * exp(4 * fconf(m,n)) 
    /*1*/  * fj(m,n) + (p(m,n) * ((4 * k_f * exp(2 
    /*5*/  * gconf(m,n)) * gA(m,n) * P_1_1(R(m,n))) 
    /*3*/  / pow2(R(m,n)) + (exp(-2 * fconf(m,n)) * fA(m,n) 
    /*4*/  * (fL_r(m,n) - (2 * exp(4 * fconf(m,n) - 4 
    /*7*/  * gconf(m,n))) / (pow2(gB(m,n)) * pow2(r(m,n)) 
    /*6*/  * pow2(R(m,n))))) / Lt(m,n))) / pow2(fA(m,n)));

Jacobian[14][5]=
	(-2 * k_f * exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*1*/  * p(m,n) * P_2_1(R(m,n))) / (pow2(fA(m,n)) 
    /*1*/  * pow2(R(m,n)));

Jacobian[14][6]=
	(-4 * k_f * exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*1*/  * gA(m,n) * p(m,n) * P_1_1(R(m,n))) / (gB(m,n) 
    /*1*/  * pow2(fA(m,n)) * pow2(R(m,n)));

Jacobian[14][7]=
	(4 * fAsig(m,n) * fAlp(m,n) * pow2(r(m,n)) * (-3 
    /*2*/  * fA_r(m,n) - 2 * fsig(m,n) * fA(m,n) * r(m,n) 
    /*2*/  + fA(m,n) * (-12 * fconf_r(m,n) - (2 * exp(2 
    /*5*/  * fconf(m,n) - 2 * gconf(m,n)) * fB_r(m,n)) 
    /*3*/  / (gB(m,n) * R(m,n)))) + fA(m,n) * (fdet_pff_r(m,n)
    /*2*/  / fdet(m,n) + 12 * fA1(m,n) * fAlp_r(m,n) - 6 
    /*2*/  * fBet_rr(m,n) - (fdet_pff(m,n) * fdet_r(m,n)) 
    /*2*/  / pow2(fdet(m,n)) - 4 * fAlp_r(m,n) * ftrA(m,n) 
    /*2*/  + fAlp(m,n) * ((12 * k_f * exp(2 * gconf(m,n)) 
    /*4*/  * gA(m,n) * p(m,n) * P_2_1(R(m,n))) / pow2(R(m,n)) 
    /*3*/  + 8 * (ftrA_r(m,n) + ftrK_r(m,n)) + (3 * exp(-2 
    /*5*/  * fconf(m,n)) * fA(m,n) * p(m,n) * (fL_r(m,n) - (2 
    /*6*/  * exp(4 * fconf(m,n) - 4 * gconf(m,n))) 
    /*5*/  / (pow2(gB(m,n)) * pow2(r(m,n)) * pow2(R(m,n))))) 
    /*3*/  / Lt(m,n)))) / (3. * Power(fA(m,n),4));

Jacobian[14][8]=
	(4 * exp(2 * fconf(m,n) - 8 * gconf(m,n)) * (3 
    /*2*/  * exp(4 * fconf(m,n)) * pow2(fA(m,n)) * (exp(2 
    /*4*/  * gconf(m,n)) * gBet(m,n) * gA(m,n) * Lt(m,n) 
    /*3*/  - gAlp(m,n) * p(m,n) - exp(2 * gconf(m,n)) 
    /*3*/  * fBet_r(m,n) * gA(m,n) * Lt(m,n) * r(m,n)) 
    /*2*/  * Lt(m,n) + exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*2*/  * gA(m,n) * Lt(m,n) * (-3 * exp(2 * fconf(m,n)) 
    /*3*/  * fA(m,n) * p(m,n) + exp(2 * gconf(m,n)) * gB(m,n) 
    /*3*/  * pow2(r(m,n)) * (3 * k_f * exp(4 * gconf(m,n)) 
    /*4*/  * gA(m,n) * gB(m,n) * p(m,n) * P_1_1(R(m,n)) - exp(2
    /*5*/  * fconf(m,n)) * fAsig(m,n) * fB_r(m,n) 
    /*4*/  * pow2(r(m,n)) * R(m,n)) * Lt(m,n)))) / (3. 
    /*1*/  * gA(m,n) * Lt(m,n) * pow2(fA(m,n)) * pow2(r(m,n)) 
    /*1*/  * pow3(gB(m,n)) * pow3(R(m,n)) * Lt(m,n));

Jacobian[14][11]=
	(-2 * fAlp_r(m,n)) / pow2(fA(m,n));

Jacobian[14][14]=
	-fdet_pff(m,n) / (3. * fdet(m,n)) - fBet_r(m,n);

Jacobian[14][16]=
	(4 * fAsig(m,n) * fAlp(m,n) * pow3(r(m,n))) / (3.
    /*1*/  * pow2(fA(m,n)));

Jacobian[14][18]=
	(4 * fAlp(m,n) * pow2(r(m,n)) * (fA_r(m,n) 
    /*2*/  + fsig(m,n) * fA(m,n) * r(m,n) + fA(m,n) * (6 
    /*3*/  * fconf_r(m,n) + (exp(2 * fconf(m,n) - 2 
    /*5*/  * gconf(m,n)) * fB_r(m,n)) / (gB(m,n) * R(m,n))))) 
    /*0*/  / (3. * pow3(fA(m,n)));

Jacobian[14][22]=
	fL_r(m,n) - (2 * exp(4 * fconf(m,n) - 4 
    /*2*/  * gconf(m,n))) / (pow2(gB(m,n)) * pow2(r(m,n)) 
    /*1*/  * pow2(R(m,n)));

Jacobian[15][1]=
	(-2 * exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) 
    /*1*/  * (2 * pow2(gA(m,n)) + pow2(gB(m,n)) * pow2(r(m,n))
    /*2*/  * (2 * gsig(m,n) + gsig_r(m,n) * r(m,n)))) 
    /*0*/  / (gA(m,n) * pow2(gB(m,n)) * pow3(r(m,n)) * Lt(m,n));

Jacobian[15][5]=
	(exp(-2 * gconf(m,n)) * (4 * exp(2 * gconf(m,n)) 
    /*2*/  * Lt(m,n) * pow3(gA(m,n)) * (gBet(m,n) - gBet_r(m,n)
    /*3*/  * r(m,n)) * Lt(m,n) + gAlp(m,n) * (-(Lt(m,n) 
    /*4*/  * p(m,n) * pow2(gB(m,n)) * pow2(r(m,n)) * (2 
    /*5*/  * gsig(m,n) + gsig_r(m,n) * r(m,n))) + 2 * p(m,n) 
    /*3*/  * pow2(gA(m,n)) * (Lt(m,n) - 2 * Lt(m,n)) + 4 
    /*3*/  * exp(2 * gconf(m,n)) * gAsig(m,n) * Lt(m,n) 
    /*3*/  * pow3(gA(m,n)) * pow3(r(m,n)) * Lt(m,n)))) 
    /*0*/  / (Lt(m,n) * pow2(gA(m,n)) * pow2(gB(m,n)) 
    /*1*/  * pow3(r(m,n)) * Lt(m,n));

Jacobian[15][6]=
	(-4 * exp(-2 * gconf(m,n)) * gA(m,n) * (exp(2 
    /*3*/  * gconf(m,n)) * gA(m,n) * Lt(m,n) * (gBet(m,n) 
    /*3*/  - gBet_r(m,n) * r(m,n)) * Lt(m,n) + gAlp(m,n) 
    /*2*/  * (-(p(m,n) * Lt(m,n)) + Lt(m,n) * (p(m,n) + exp(2 
    /*5*/  * gconf(m,n)) * gAsig(m,n) * gA(m,n) * pow3(r(m,n))
    /*4*/  * Lt(m,n))))) / (Lt(m,n) * pow3(gB(m,n)) 
    /*1*/  * pow3(r(m,n)) * Lt(m,n));

Jacobian[15][15]=
	2 * (gBetr(m,n) + (exp(-2 * gconf(m,n)) 
    /*2*/  * gAlp(m,n) * p(m,n) * (Lt(m,n) - Lt(m,n))) 
    /*1*/  / (gA(m,n) * Lt(m,n) * r(m,n) * Lt(m,n)));

Jacobian[15][17]=
	(2 * gAlp(m,n) * pow2(gA(m,n))) / pow2(gB(m,n));

Jacobian[15][22]=
	gsig_r(m,n) + (2 * pow2(gA(m,n))) / (pow2(gB(m,n))
    /*1*/  * pow3(r(m,n))) + (2 * gsig(m,n)) / r(m,n);

Jacobian[16][2]=
	(2 * exp(-2 * (fconf(m,n) + 2 * gconf(m,n))) 
    /*1*/  * fAlp(m,n) * p(m,n) * (2 * exp(4 * fconf(m,n)) 
    /*2*/  * pow2(fA(m,n)) + exp(4 * gconf(m,n)) 
    /*2*/  * pow2(gB(m,n)) * pow2(r(m,n)) * pow2(R(m,n)) * (2 
    /*3*/  * fsig(m,n) + fsig_r(m,n) * r(m,n)))) / (fA(m,n) 
    /*1*/  * pow2(gB(m,n)) * pow2(R(m,n)) * pow3(r(m,n)) 
    /*1*/  * Lt(m,n));

Jacobian[16][7]=
	(exp(-2 * (fconf(m,n) + 3 * gconf(m,n))) * (4 
    /*2*/  * exp(6 * fconf(m,n)) * pow3(fA(m,n)) * (exp(2 
    /*4*/  * gconf(m,n)) * gBet(m,n) * gA(m,n) * Lt(m,n) 
    /*3*/  - gAlp(m,n) * p(m,n) - exp(2 * gconf(m,n)) 
    /*3*/  * fBet_r(m,n) * gA(m,n) * Lt(m,n) * r(m,n)) 
    /*2*/  * Lt(m,n) + exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*2*/  * gA(m,n) * Lt(m,n) * (-2 * exp(4 * fconf(m,n)) 
    /*3*/  * p(m,n) * pow2(fA(m,n)) + exp(4 * gconf(m,n)) 
    /*3*/  * p(m,n) * pow2(gB(m,n)) * pow2(r(m,n)) 
    /*3*/  * pow2(R(m,n)) * (2 * fsig(m,n) + fsig_r(m,n) 
    /*4*/  * r(m,n)) + 4 * exp(6 * fconf(m,n)) * fAsig(m,n) 
    /*3*/  * pow3(fA(m,n)) * pow3(r(m,n)) * Lt(m,n)))) 
    /*0*/  / (gA(m,n) * Lt(m,n) * pow2(fA(m,n)) * pow2(gB(m,n))
    /*1*/  * pow2(R(m,n)) * pow3(r(m,n)) * Lt(m,n));

Jacobian[16][8]=
	(-4 * exp(4 * fconf(m,n) - 8 * gconf(m,n)) 
    /*1*/  * fA(m,n) * (exp(2 * fconf(m,n)) * fA(m,n) * (exp(2
    /*4*/  * gconf(m,n)) * gBet(m,n) * gA(m,n) * Lt(m,n) 
    /*3*/  - gAlp(m,n) * p(m,n) - exp(2 * gconf(m,n)) 
    /*3*/  * fBet_r(m,n) * gA(m,n) * Lt(m,n) * r(m,n)) 
    /*2*/  * Lt(m,n) - exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*2*/  * gA(m,n) * Lt(m,n) * (p(m,n) - exp(2 * fconf(m,n))
    /*3*/  * fAsig(m,n) * fA(m,n) * pow3(r(m,n)) * Lt(m,n))))
    /*0*/  / (gA(m,n) * Lt(m,n) * pow3(gB(m,n)) 
    /*1*/  * pow3(r(m,n)) * pow3(R(m,n)) * Lt(m,n));

Jacobian[16][16]=
	2 * (gBetr(m,n) + (p(m,n) * (-((exp(-2 
    /*6*/  * gconf(m,n)) * gAlp(m,n)) / (gA(m,n) * Lt(m,n))) 
    /*3*/  - (exp(-2 * fconf(m,n)) * fAlp(m,n)) / (fA(m,n) 
    /*4*/  * Lt(m,n)))) / r(m,n));

Jacobian[16][18]=
	(2 * exp(4 * fconf(m,n) - 4 * gconf(m,n)) 
    /*1*/  * fAlp(m,n) * pow2(fA(m,n))) / (pow2(gB(m,n)) 
    /*1*/  * pow2(R(m,n)));

Jacobian[16][22]=
	fsig_r(m,n) + (2 * exp(4 * fconf(m,n) - 4 
    /*2*/  * gconf(m,n)) * pow2(fA(m,n))) / (pow2(gB(m,n)) 
    /*1*/  * pow2(R(m,n)) * pow3(r(m,n))) + (2 * fsig(m,n)) 
    /*0*/  / r(m,n);

Jacobian[17][1]=
	(-2 * exp(-4 * gconf(m,n)) * (2 * exp(2 
    /*3*/  * gconf(m,n)) * gAsig(m,n) * gAlp(m,n) * gB(m,n) 
    /*2*/  * p(m,n) * pow2(r(m,n)) * pow3(gA(m,n)) + exp(2 
    /*3*/  * gconf(m,n)) * gAsig_r(m,n) * gAlp(m,n) * gB(m,n) 
    /*2*/  * p(m,n) * pow3(gA(m,n)) * pow3(r(m,n)) + k_g 
    /*2*/  * exp(4 * gconf(m,n)) * fAlp(m,n) * Power(gA(m,n),4)
    /*2*/  * gB(m,n) * r(m,n) * R(m,n) * b_2 + 2 * k_g 
    /*2*/  * exp(4 * gconf(m,n)) * fAlp(m,n) * Power(gA(m,n),4)
    /*2*/  * gB(m,n) * pow2(R(m,n)) * r(m,n) * b_3 - 2 
    /*2*/  * gAlp(m,n) * gL(m,n) * Power(gA(m,n),4) * gB(m,n) 
    /*2*/  * Lt(m,n) + 4 * gconf_r(m,n) * gAlp(m,n) * gB(m,n) 
    /*2*/  * pow2(gA(m,n)) * Lt(m,n) + 2 * gAlp_r(m,n) 
    /*2*/  * gB(m,n) * pow2(gA(m,n)) * Lt(m,n) + 2 * gAlp(m,n)
    /*2*/  * gsig_r(m,n) * gB(m,n) * pow2(gA(m,n)) 
    /*2*/  * pow2(r(m,n)) * Lt(m,n) - 8 * gAlp(m,n) * gsig(m,n)
    /*2*/  * gB_r(m,n) * pow2(gA(m,n)) * pow2(r(m,n)) 
    /*2*/  * Lt(m,n) + 4 * gAlp(m,n) * gsig_r(m,n) 
    /*2*/  * pow2(r(m,n)) * pow3(gB(m,n)) * Lt(m,n) - 2 
    /*2*/  * gAlp(m,n) * gL(m,n) * gsig(m,n) * pow2(gA(m,n)) 
    /*2*/  * pow2(r(m,n)) * pow3(gB(m,n)) * Lt(m,n) + gAlp(m,n)
    /*2*/  * gsig_rr(m,n) * pow3(gB(m,n)) * pow3(r(m,n)) 
    /*2*/  * Lt(m,n) + 2 * gAlp(m,n) * pow2(gsig(m,n)) 
    /*2*/  * pow3(gB(m,n)) * pow3(r(m,n)) * Lt(m,n) - gAlp(m,n)
    /*2*/  * gL(m,n) * gsig_r(m,n) * pow2(gA(m,n)) 
    /*2*/  * pow3(gB(m,n)) * pow3(r(m,n)) * Lt(m,n) + 2 
    /*2*/  * gAlp(m,n) * gL_r(m,n) * Power(gA(m,n),4) * gB(m,n)
    /*2*/  * r(m,n) * Lt(m,n) + 4 * gconf_r(m,n) * gAlp(m,n)
    /*2*/  * gA(m,n) * gA_r(m,n) * gB(m,n) * r(m,n) * Lt(m,n)
    /*2*/  + 2 * gAlp_r(m,n) * gA(m,n) * gA_r(m,n) * gB(m,n)
    /*2*/  * r(m,n) * Lt(m,n) - 8 * gAlp(m,n) * gA(m,n) 
    /*2*/  * gA_r(m,n) * gB_r(m,n) * r(m,n) * Lt(m,n) - 4 
    /*2*/  * gconf_rr(m,n) * gAlp(m,n) * gB(m,n) 
    /*2*/  * pow2(gA(m,n)) * r(m,n) * Lt(m,n) + 8 
    /*2*/  * gconf_r(m,n) * gAlp_r(m,n) * gB(m,n) 
    /*2*/  * pow2(gA(m,n)) * r(m,n) * Lt(m,n) - 2 
    /*2*/  * gAlp_rr(m,n) * gB(m,n) * pow2(gA(m,n)) * r(m,n) 
    /*2*/  * Lt(m,n) + 4 * gconf_r(m,n) * gAlp(m,n) * gB_r(m,n)
    /*2*/  * pow2(gA(m,n)) * r(m,n) * Lt(m,n) + 2 
    /*2*/  * gAlp_r(m,n) * gB_r(m,n) * pow2(gA(m,n)) * r(m,n) 
    /*2*/  * Lt(m,n) + 8 * gAlp(m,n) * gB(m,n) 
    /*2*/  * pow2(gconf_r(m,n)) * pow2(gA(m,n)) * r(m,n) 
    /*2*/  * Lt(m,n) + 6 * gAlp(m,n) * gB(m,n) 
    /*2*/  * pow2(gA_r(m,n)) * r(m,n) * Lt(m,n) + k_g * exp(4 
    /*3*/  * gconf(m,n)) * gAlp(m,n) * Power(gA(m,n),4) 
    /*2*/  * gB(m,n) * r(m,n) * R(m,n) * b_1 * Lt(m,n) + 2 
    /*2*/  * k_g * exp(4 * gconf(m,n)) * gAlp(m,n) 
    /*2*/  * Power(gA(m,n),4) * gB(m,n) * pow2(R(m,n)) * r(m,n)
    /*2*/  * b_2 * Lt(m,n) - k_g * exp(2 * (fconf(m,n) 
    /*4*/  + gconf(m,n))) * fA(m,n) * gB(m,n) * pow3(gA(m,n)) 
    /*2*/  * r(m,n) * (gAlp(m,n) * (2 * (-2 + pow2(Lt(m,n))) 
    /*4*/  * P_1_1(R(m,n)) - 3 * (-1 + pow2(Lt(m,n))) 
    /*4*/  * P_2_1(R(m,n)) - b_1) - fAlp(m,n) * (2 
    /*4*/  * P_1_2(R(m,n)) + b_2) * Lt(m,n)))) 
    /*0*/  / (Power(gA(m,n),4) * gB(m,n) * pow3(r(m,n)) 
    /*1*/  * Lt(m,n));

Jacobian[17][2]=
	(2 * k_g * exp(-2 * gconf(m,n)) * (-(exp(2 
    /*4*/  * fconf(m,n)) * fA(m,n) * (gAlp(m,n) * (2 * (-2 
    /*6*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) - 3 * (-1 
    /*6*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)) - b_1) - fAlp(m,n)
    /*4*/  * (2 * P_1_2(R(m,n)) + b_2) * Lt(m,n))) + exp(2 
    /*3*/  * gconf(m,n)) * gA(m,n) * R(m,n) * (fAlp(m,n) * (b_2
    /*4*/  + 2 * R(m,n) * b_3) + gAlp(m,n) * (b_1 + 2 
    /*4*/  * R(m,n) * b_2) * Lt(m,n)))) / (gA(m,n) 
    /*1*/  * pow2(r(m,n)) * Lt(m,n));

Jacobian[17][3]=
	gAsig(m,n) * gAlp(m,n);

Jacobian[17][5]=
	(exp(-4 * gconf(m,n)) * (-2 * exp(2 * gconf(m,n))
    /*2*/  * gAsig(m,n) * gAlp(m,n) * gB(m,n) * p(m,n) 
    /*2*/  * pow2(r(m,n)) * pow3(gA(m,n)) - exp(2 * gconf(m,n))
    /*2*/  * gAsig_r(m,n) * gAlp(m,n) * gB(m,n) * p(m,n) 
    /*2*/  * pow3(gA(m,n)) * pow3(r(m,n)) - 4 * gconf_r(m,n) 
    /*2*/  * gAlp(m,n) * gB(m,n) * pow2(gA(m,n)) * Lt(m,n) - 2
    /*2*/  * gAlp_r(m,n) * gB(m,n) * pow2(gA(m,n)) * Lt(m,n)
    /*2*/  - 2 * gAlp(m,n) * gsig_r(m,n) * gB(m,n) 
    /*2*/  * pow2(gA(m,n)) * pow2(r(m,n)) * Lt(m,n) + 8 
    /*2*/  * gAlp(m,n) * gsig(m,n) * gB_r(m,n) * pow2(gA(m,n))
    /*2*/  * pow2(r(m,n)) * Lt(m,n) - 8 * gAlp(m,n) 
    /*2*/  * gsig_r(m,n) * pow2(r(m,n)) * pow3(gB(m,n)) 
    /*2*/  * Lt(m,n) + 2 * gAlp(m,n) * gL(m,n) * gsig(m,n) 
    /*2*/  * pow2(gA(m,n)) * pow2(r(m,n)) * pow3(gB(m,n)) 
    /*2*/  * Lt(m,n) - 2 * gAlp(m,n) * gsig_rr(m,n) 
    /*2*/  * pow3(gB(m,n)) * pow3(r(m,n)) * Lt(m,n) - 4 
    /*2*/  * gAlp(m,n) * pow2(gsig(m,n)) * pow3(gB(m,n)) 
    /*2*/  * pow3(r(m,n)) * Lt(m,n) + gAlp(m,n) * gL(m,n) 
    /*2*/  * gsig_r(m,n) * pow2(gA(m,n)) * pow3(gB(m,n)) 
    /*2*/  * pow3(r(m,n)) * Lt(m,n) - 6 * gconf_r(m,n) 
    /*2*/  * gAlp(m,n) * gA(m,n) * gA_r(m,n) * gB(m,n) * r(m,n)
    /*2*/  * Lt(m,n) - 3 * gAlp_r(m,n) * gA(m,n) * gA_r(m,n)
    /*2*/  * gB(m,n) * r(m,n) * Lt(m,n) + 12 * gAlp(m,n) 
    /*2*/  * gA(m,n) * gA_r(m,n) * gB_r(m,n) * r(m,n) * Lt(m,n)
    /*2*/  + 4 * gconf_rr(m,n) * gAlp(m,n) * gB(m,n) 
    /*2*/  * pow2(gA(m,n)) * r(m,n) * Lt(m,n) - 8 
    /*2*/  * gconf_r(m,n) * gAlp_r(m,n) * gB(m,n) 
    /*2*/  * pow2(gA(m,n)) * r(m,n) * Lt(m,n) + 2 
    /*2*/  * gAlp_rr(m,n) * gB(m,n) * pow2(gA(m,n)) * r(m,n) 
    /*2*/  * Lt(m,n) - 4 * gconf_r(m,n) * gAlp(m,n) * gB_r(m,n)
    /*2*/  * pow2(gA(m,n)) * r(m,n) * Lt(m,n) - 2 
    /*2*/  * gAlp_r(m,n) * gB_r(m,n) * pow2(gA(m,n)) * r(m,n) 
    /*2*/  * Lt(m,n) - 8 * gAlp(m,n) * gB(m,n) 
    /*2*/  * pow2(gconf_r(m,n)) * pow2(gA(m,n)) * r(m,n) 
    /*2*/  * Lt(m,n) - 12 * gAlp(m,n) * gB(m,n) 
    /*2*/  * pow2(gA_r(m,n)) * r(m,n) * Lt(m,n) - k_g * exp(2 
    /*3*/  * (fconf(m,n) + gconf(m,n))) * fA(m,n) * gB(m,n) 
    /*2*/  * pow3(gA(m,n)) * r(m,n) * (gAlp(m,n) 
    /*3*/  * (P_1_1(R(m,n)) + (-1 + pow2(Lt(m,n))) 
    /*4*/  * P_2_1(R(m,n))) + fAlp(m,n) * P_1_2(R(m,n)) 
    /*3*/  * Lt(m,n)))) / (Power(gA(m,n),5) * gB(m,n) 
    /*1*/  * pow3(r(m,n)) * Lt(m,n));

Jacobian[17][6]=
	(exp(-4 * gconf(m,n)) * (-(k_g * exp(4 
    /*4*/  * gconf(m,n)) * fAlp(m,n) * Power(gA(m,n),4) 
    /*3*/  * gB(m,n) * R(m,n) * (b_2 + 2 * R(m,n) * b_3)) 
    /*2*/  + (-(gAlp_r(m,n) * gB_r(m,n) * pow2(gA(m,n))) 
    /*3*/  + gAlp(m,n) * (4 * gA(m,n) * gA_r(m,n) * gB_r(m,n) 
    /*4*/  + 2 * pow2(gsig(m,n)) * pow2(r(m,n)) * pow3(gB(m,n))
    /*4*/  - pow2(gA(m,n)) * (2 * gconf_r(m,n) * gB_r(m,n) 
    /*5*/  + gL(m,n) * gsig_r(m,n) * pow2(r(m,n)) 
    /*5*/  * pow3(gB(m,n))) - 2 * gsig(m,n) * pow2(gA(m,n)) 
    /*4*/  * (-2 * gB_r(m,n) + gL(m,n) * pow3(gB(m,n))) 
    /*4*/  * r(m,n) + pow3(gB(m,n)) * r(m,n) * (4 * gsig_r(m,n)
    /*5*/  + gsig_rr(m,n) * r(m,n)) - k_g * exp(4 
    /*5*/  * gconf(m,n)) * Power(gA(m,n),4) * gB(m,n) * R(m,n)
    /*4*/  * (b_1 + 2 * R(m,n) * b_2))) * Lt(m,n) - k_g 
    /*2*/  * exp(2 * (fconf(m,n) + gconf(m,n))) * fA(m,n) 
    /*2*/  * gB(m,n) * pow3(gA(m,n)) * R(m,n) * (gAlp(m,n) * (2
    /*4*/  * (-1 + pow2(Lt(m,n))) * P_1_2(R(m,n)) - b_2) 
    /*3*/  - fAlp(m,n) * b_3 * Lt(m,n)))) / (Power(gA(m,n),4) 
    /*1*/  * pow2(gB(m,n)) * pow2(r(m,n)) * Lt(m,n));

Jacobian[17][7]=
	(k_g * exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*1*/  * (gAlp(m,n) * (P_1_1(R(m,n)) + (-1 + pow2(Lt(m,n)))
    /*3*/  * P_2_1(R(m,n))) + fAlp(m,n) * P_1_2(R(m,n)) 
    /*2*/  * Lt(m,n))) / (gA(m,n) * pow2(r(m,n)) * Lt(m,n));

Jacobian[17][8]=
	-((k_g * exp(2 * fconf(m,n) - 4 * gconf(m,n)) 
    /*2*/  * (exp(2 * gconf(m,n)) * gA(m,n) * (fAlp(m,n) * (2 
    /*5*/  * P_1_2(R(m,n)) + b_2) + gAlp(m,n) * (2 
    /*5*/  * P_1_1(R(m,n)) + b_1) * Lt(m,n)) - exp(2 
    /*4*/  * fconf(m,n)) * fA(m,n) * (gAlp(m,n) * (2 * (-1 
    /*6*/  + pow2(Lt(m,n))) * P_1_2(R(m,n)) - b_2) - fAlp(m,n)
    /*4*/  * b_3 * Lt(m,n)))) / (gA(m,n) * gB(m,n) 
    /*2*/  * pow2(r(m,n)) * Lt(m,n)));

Jacobian[17][13]=
	-(exp(-4 * gconf(m,n)) * gAlp(m,n) * (2 
    /*2*/  * pow2(gA(m,n)) + pow2(gB(m,n)) * pow2(r(m,n)) * (2
    /*3*/  * gsig(m,n) + gsig_r(m,n) * r(m,n)))) / (2. 
    /*1*/  * pow2(gA(m,n)) * pow3(r(m,n)));

Jacobian[17][15]=
	(exp(-4 * gconf(m,n)) * gAlp(m,n) 
    /*1*/  * (-(pow2(gA(m,n)) * (4 * gB_r(m,n) + gL(m,n) 
    /*4*/  * pow3(gB(m,n)))) + 2 * gsig(m,n) * pow3(gB(m,n)) 
    /*2*/  * r(m,n))) / (Power(gA(m,n),4) * gB(m,n) * r(m,n));

Jacobian[17][17]=
	2 * gBetr(m,n) + (exp(-2 * gconf(m,n)) * gAlp(m,n)
    /*1*/  * (-2 * p(m,n) * Lt(m,n) + Lt(m,n) * (2 * p(m,n) 
    /*3*/  + exp(2 * gconf(m,n)) * gA(m,n) * r(m,n) 
    /*3*/  * (gtrA(m,n) + gtrK(m,n)) * Lt(m,n)))) / (gA(m,n) 
    /*1*/  * Lt(m,n) * r(m,n) * Lt(m,n));

Jacobian[17][22]=
	gAsig_r(m,n) + (2 * gAsig(m,n)) / r(m,n);

Jacobian[18][1]=
	(-2 * k_f * exp(-2 * fconf(m,n)) * (-(exp(2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (fAlp(m,n) * (2 * (-1 
    /*6*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) + pow2(Lt(m,n)) 
    /*5*/  * P_2_1(R(m,n)) + b_1) + gAlp(m,n) * (P_2_0(R(m,n))
    /*5*/  + b_0) * Lt(m,n))) + exp(2 * fconf(m,n)) * fA(m,n)
    /*2*/  * (gAlp(m,n) * (P_1_1(R(m,n)) - b_1) + fAlp(m,n) 
    /*3*/  * (P_1_2(R(m,n)) - b_2) * Lt(m,n)))) / (fA(m,n) 
    /*1*/  * pow2(r(m,n)) * pow2(R(m,n)) * Lt(m,n));

Jacobian[18][2]=
	(2 * exp(-2 * (4 * fconf(m,n) + gconf(m,n))) 
    /*1*/  * (-(exp(2 * gconf(m,n)) * fAlp(m,n) * gB(m,n) 
    /*3*/  * pow2(R(m,n)) * r(m,n) * (6 * exp(4 * fconf(m,n)) 
    /*4*/  * pow2(fA_r(m,n)) + 2 * exp(4 * gconf(m,n)) 
    /*4*/  * pow2(fsig(m,n)) * pow2(gB(m,n)) * pow2(r(m,n)) 
    /*4*/  * pow2(R(m,n)) + exp(4 * gconf(m,n)) * pow2(gB(m,n))
    /*4*/  * pow2(R(m,n)) * r(m,n) * (4 * fsig_r(m,n) 
    /*5*/  + fsig_rr(m,n) * r(m,n))) * Lt(m,n)) - 2 * exp(4 
    /*3*/  * fconf(m,n)) * fA(m,n) * fA_r(m,n) * r(m,n) 
    /*2*/  * R(m,n) * (exp(2 * gconf(m,n)) * fAlp_r(m,n) 
    /*3*/  * gB(m,n) * R(m,n) + 2 * fAlp(m,n) * (-2 * exp(2 
    /*5*/  * fconf(m,n)) * fB_r(m,n) + exp(2 * gconf(m,n)) 
    /*4*/  * fconf_r(m,n) * gB(m,n) * R(m,n))) * Lt(m,n) 
    /*2*/  + pow2(fA(m,n)) * R(m,n) * (-2 * exp(4 * fconf(m,n))
    /*3*/  * (exp(2 * fconf(m,n)) * fAlp_r(m,n) * fB_r(m,n) 
    /*4*/  * r(m,n) + exp(2 * gconf(m,n)) * gB(m,n) 
    /*4*/  * (-(fAlp_rr(m,n) * r(m,n)) + fAlp_r(m,n) * (1 + 4 
    /*6*/  * fconf_r(m,n) * r(m,n))) * R(m,n)) + fAlp(m,n) 
    /*3*/  * (exp(6 * gconf(m,n)) * fL(m,n) * pow2(r(m,n)) 
    /*4*/  * pow3(gB(m,n)) * pow3(R(m,n)) * (2 * fsig(m,n) 
    /*5*/  + fsig_r(m,n) * r(m,n)) - 2 * exp(4 * fconf(m,n)) 
    /*4*/  * (-2 * exp(2 * fconf(m,n)) * fB_r(m,n) * r(m,n) 
    /*5*/  * (-fconf_r(m,n) + 2 * fsig(m,n) * r(m,n)) + exp(2 
    /*6*/  * gconf(m,n)) * gB(m,n) * (2 * fconf_r(m,n) 
    /*6*/  + fsig_r(m,n) * pow2(r(m,n)) - 2 * fconf_rr(m,n) 
    /*6*/  * r(m,n) + 4 * pow2(fconf_r(m,n)) * r(m,n)) 
    /*5*/  * R(m,n)))) * Lt(m,n) - exp(6 * fconf(m,n) + 2 
    /*3*/  * gconf(m,n)) * gB(m,n) * pow3(fA(m,n)) * r(m,n) 
    /*2*/  * (fAlp(m,n) * (-(p(m,n) * pow2(R(m,n)) * r(m,n) 
    /*5*/  * (2 * fAsig(m,n) + fAsig_r(m,n) * r(m,n))) + k_f 
    /*4*/  * exp(2 * gconf(m,n)) * gA(m,n) * (2 * (-1 
    /*6*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) + pow2(Lt(m,n)) 
    /*5*/  * P_2_1(R(m,n)) + b_1)) + k_f * exp(2 * gconf(m,n))
    /*3*/  * gAlp(m,n) * gA(m,n) * (P_2_0(R(m,n)) + b_0) 
    /*3*/  * Lt(m,n)) + exp(4 * fconf(m,n) + 2 * gconf(m,n)) 
    /*2*/  * Power(fA(m,n),4) * gB(m,n) * (k_f * exp(4 
    /*4*/  * fconf(m,n)) * gAlp(m,n) * r(m,n) * (P_1_1(R(m,n))
    /*4*/  - b_1) + fAlp(m,n) * (2 * fL(m,n) * pow2(R(m,n)) 
    /*4*/  - 2 * fL_r(m,n) * pow2(R(m,n)) * r(m,n) + k_f 
    /*4*/  * exp(4 * fconf(m,n)) * P_1_2(R(m,n)) * r(m,n) - k_f
    /*4*/  * exp(4 * fconf(m,n)) * r(m,n) * b_2) * Lt(m,n))))
    /*0*/  / (Power(fA(m,n),4) * gB(m,n) * pow2(R(m,n)) 
    /*1*/  * pow3(r(m,n)) * Lt(m,n));

Jacobian[18][4]=
	fAsig(m,n) * fAlp(m,n);

Jacobian[18][5]=
	(k_f * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*1*/  * (fAlp(m,n) * (-P_1_1(R(m,n)) + pow2(Lt(m,n)) 
    /*3*/  * P_2_1(R(m,n))) + gAlp(m,n) * (-P_1_0(R(m,n)) 
    /*3*/  + P_2_0(R(m,n))) * Lt(m,n))) / (fA(m,n) 
    /*1*/  * pow2(r(m,n)) * pow2(R(m,n)) * Lt(m,n));

Jacobian[18][6]=
	-((k_f * exp(-2 * fconf(m,n)) * (-(exp(2 
    /*5*/  * gconf(m,n)) * gA(m,n) * (fAlp(m,n) * ((-1 + 2 
    /*7*/  * pow2(Lt(m,n))) * P_1_1(R(m,n)) + b_1) + gAlp(m,n)
    /*5*/  * (P_1_0(R(m,n)) + b_0) * Lt(m,n))) + exp(2 
    /*4*/  * fconf(m,n)) * fA(m,n) * (gAlp(m,n) 
    /*4*/  * (P_1_1(R(m,n)) - b_1) + fAlp(m,n) * (P_1_2(R(m,n))
    /*5*/  - b_2) * Lt(m,n)))) / (fA(m,n) * gB(m,n) 
    /*2*/  * pow2(r(m,n)) * pow2(R(m,n)) * Lt(m,n)));

Jacobian[18][7]=
	(exp(-2 * (4 * fconf(m,n) + gconf(m,n))) 
    /*1*/  * (-(exp(4 * fconf(m,n)) * fA(m,n) * R(m,n) * (k_f 
    /*4*/  * exp(2 * fconf(m,n) + 4 * gconf(m,n)) * gAlp(m,n) 
    /*4*/  * gA(m,n) * gB(m,n) * pow2(fA(m,n)) * P_1_1(R(m,n))
    /*4*/  * r(m,n) + 3 * exp(2 * gconf(m,n)) * fAlp_r(m,n) 
    /*4*/  * fA_r(m,n) * gB(m,n) * r(m,n) * R(m,n) + 2 
    /*4*/  * fA(m,n) * (exp(2 * fconf(m,n)) * fAlp_r(m,n) 
    /*5*/  * fB_r(m,n) * r(m,n) + exp(2 * gconf(m,n)) * gB(m,n)
    /*5*/  * (-(fAlp_rr(m,n) * r(m,n)) + fAlp_r(m,n) * (1 + 4
    /*7*/  * fconf_r(m,n) * r(m,n))) * R(m,n))) * Lt(m,n)) 
    /*2*/  + fAlp(m,n) * (exp(6 * fconf(m,n) + 2 * gconf(m,n))
    /*3*/  * gB(m,n) * pow3(fA(m,n)) * r(m,n) * (k_f * exp(2
    /*5*/  * gconf(m,n)) * gA(m,n) * (P_1_1(R(m,n)) 
    /*5*/  - pow2(Lt(m,n)) * P_2_1(R(m,n))) + p(m,n) 
    /*4*/  * pow2(R(m,n)) * r(m,n) * (2 * fAsig(m,n) 
    /*5*/  + fAsig_r(m,n) * r(m,n))) - 2 * exp(2 * gconf(m,n))
    /*3*/  * gB(m,n) * pow2(R(m,n)) * r(m,n) * (6 * exp(4 
    /*5*/  * fconf(m,n)) * pow2(fA_r(m,n)) + 2 * exp(4 
    /*5*/  * gconf(m,n)) * pow2(fsig(m,n)) * pow2(gB(m,n)) 
    /*4*/  * pow2(r(m,n)) * pow2(R(m,n)) + exp(4 * gconf(m,n))
    /*4*/  * pow2(gB(m,n)) * pow2(R(m,n)) * r(m,n) * (4 
    /*5*/  * fsig_r(m,n) + fsig_rr(m,n) * r(m,n))) * Lt(m,n) 
    /*3*/  + 6 * exp(4 * fconf(m,n)) * fA(m,n) * fA_r(m,n) 
    /*3*/  * r(m,n) * R(m,n) * (2 * exp(2 * fconf(m,n)) 
    /*4*/  * fB_r(m,n) - exp(2 * gconf(m,n)) * fconf_r(m,n) 
    /*4*/  * gB(m,n) * R(m,n)) * Lt(m,n) + pow2(fA(m,n)) 
    /*3*/  * R(m,n) * (exp(6 * gconf(m,n)) * fL(m,n) 
    /*4*/  * pow2(r(m,n)) * pow3(gB(m,n)) * pow3(R(m,n)) * (2 
    /*5*/  * fsig(m,n) + fsig_r(m,n) * r(m,n)) - 2 * exp(4 
    /*5*/  * fconf(m,n)) * (-2 * exp(2 * fconf(m,n)) 
    /*5*/  * fB_r(m,n) * r(m,n) * (-fconf_r(m,n) + 2 
    /*6*/  * fsig(m,n) * r(m,n)) + exp(2 * gconf(m,n)) 
    /*5*/  * gB(m,n) * (2 * fconf_r(m,n) + fsig_r(m,n) 
    /*6*/  * pow2(r(m,n)) - 2 * fconf_rr(m,n) * r(m,n) + 4 
    /*6*/  * pow2(fconf_r(m,n)) * r(m,n)) * R(m,n))) 
    /*3*/  * Lt(m,n)))) / (Power(fA(m,n),5) * gB(m,n) 
    /*1*/  * pow2(R(m,n)) * pow3(r(m,n)) * Lt(m,n));

Jacobian[18][8]=
	(exp(-6 * fconf(m,n) - 4 * gconf(m,n)) * (4 
    /*2*/  * exp(6 * fconf(m,n)) * fAlp(m,n) * fA(m,n) 
    /*2*/  * fA_r(m,n) * fB_r(m,n) * R(m,n) * Lt(m,n) 
    /*2*/  - pow2(fA(m,n)) * (exp(6 * fconf(m,n)) * fAlp_r(m,n)
    /*3*/  * fB_r(m,n) + fAlp(m,n) * (-2 * exp(6 
    /*5*/  * fconf(m,n)) * fB_r(m,n) * (-fconf_r(m,n) + 2 
    /*5*/  * fsig(m,n) * r(m,n)) + exp(6 * gconf(m,n)) 
    /*4*/  * fL(m,n) * pow3(gB(m,n)) * pow3(R(m,n)) * r(m,n) 
    /*4*/  * (2 * fsig(m,n) + fsig_r(m,n) * r(m,n)))) * R(m,n)
    /*2*/  * Lt(m,n) + exp(6 * gconf(m,n)) * fAlp(m,n) 
    /*2*/  * pow3(gB(m,n)) * r(m,n) * (4 * fsig_r(m,n) 
    /*3*/  + fsig_rr(m,n) * r(m,n) + 2 * pow2(fsig(m,n)) 
    /*3*/  * r(m,n)) * Power(R(m,n),4) * Lt(m,n) - k_f * exp(6
    /*3*/  * fconf(m,n) + 4 * gconf(m,n)) * gA(m,n) * gB(m,n)
    /*2*/  * pow3(fA(m,n)) * (fAlp(m,n) * ((-1 + 2 
    /*5*/  * pow2(Lt(m,n))) * P_1_1(R(m,n)) + b_1) + gAlp(m,n)
    /*3*/  * (P_1_0(R(m,n)) + b_0) * Lt(m,n)) + k_f * exp(8 
    /*3*/  * fconf(m,n) + 2 * gconf(m,n)) * Power(fA(m,n),4) 
    /*2*/  * gB(m,n) * (gAlp(m,n) * (P_1_1(R(m,n)) - b_1) 
    /*3*/  + fAlp(m,n) * (P_1_2(R(m,n)) - b_2) * Lt(m,n)))) 
    /*0*/  / (Power(fA(m,n),4) * pow2(gB(m,n)) * pow2(r(m,n)) 
    /*1*/  * pow3(R(m,n)) * Lt(m,n));

Jacobian[18][14]=
	-(exp(-8 * fconf(m,n)) * fAlp(m,n) * (2 * exp(4 
    /*3*/  * fconf(m,n)) * pow2(fA(m,n)) + exp(4 * gconf(m,n))
    /*2*/  * pow2(gB(m,n)) * pow2(r(m,n)) * pow2(R(m,n)) * (2
    /*3*/  * fsig(m,n) + fsig_r(m,n) * r(m,n)))) / (2. 
    /*1*/  * pow2(fA(m,n)) * pow3(r(m,n)));

Jacobian[18][16]=
	(exp(-2 * (4 * fconf(m,n) + gconf(m,n))) 
    /*1*/  * fAlp(m,n) * (-(pow2(fA(m,n)) * (4 * exp(6 
    /*5*/  * fconf(m,n)) * fB_r(m,n) + exp(6 * gconf(m,n)) 
    /*4*/  * fL(m,n) * pow3(gB(m,n)) * pow3(R(m,n)))) + 2 
    /*2*/  * exp(6 * gconf(m,n)) * fsig(m,n) * pow3(gB(m,n)) 
    /*2*/  * pow3(R(m,n)) * r(m,n))) / (Power(fA(m,n),4) 
    /*1*/  * gB(m,n) * r(m,n) * R(m,n));

Jacobian[18][18]=
	2 * gBetr(m,n) - (2 * exp(-2 * gconf(m,n)) 
    /*1*/  * gAlp(m,n) * p(m,n)) / (gA(m,n) * Lt(m,n) * r(m,n))
    /*0*/  + fAlp(m,n) * (ftrA(m,n) + ftrK(m,n) - (2 * exp(-2
    /*3*/  * fconf(m,n)) * p(m,n)) / (fA(m,n) * r(m,n) 
    /*2*/  * Lt(m,n)));

Jacobian[18][22]=
	fAsig_r(m,n) + (2 * fAsig(m,n)) / r(m,n);

Jacobian[19][1]=
	(-2 * exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) 
    /*1*/  * (2 * pfD(m,n) + pfD_r(m,n) * r(m,n))) / (gA(m,n) 
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[19][5]=
	-((exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) * (2
    /*3*/  * pfD(m,n) + pfD_r(m,n) * r(m,n))) 
    /*1*/  / (pow2(gA(m,n)) * r(m,n) * Lt(m,n)));

Jacobian[19][19]=
	gBet_r(m,n) + 2 * gBetr(m,n) - gAlp_r(m,n) 
    /*0*/  * pfv(m,n) + gAlp(m,n) * (-pfv_r(m,n) - (2 
    /*2*/  * pfv(m,n)) / r(m,n) + (2 * exp(-2 * gconf(m,n)) 
    /*2*/  * p(m,n) * (Lt(m,n) - Lt(m,n))) / (gA(m,n) * Lt(m,n)
    /*2*/  * r(m,n) * Lt(m,n)));

Jacobian[19][22]=
	pfD_r(m,n) + (2 * pfD(m,n)) / r(m,n);

Jacobian[20][1]=
	(-2 * exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) 
    /*1*/  * (2 * pfS(m,n) + pfS_r(m,n) * r(m,n))) / (gA(m,n) 
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[20][5]=
	-((exp(-2 * gconf(m,n)) * gAlp(m,n) * (pfS_r(m,n)
    /*3*/  * p(m,n) * r(m,n) + pfS(m,n) * (2 * p(m,n) + exp(2
    /*5*/  * gconf(m,n)) * pfv(m,n) * gA_r(m,n) * r(m,n) 
    /*4*/  * Lt(m,n)))) / (pow2(gA(m,n)) * r(m,n) * Lt(m,n)));

Jacobian[20][19]=
	-gAlp_r(m,n);

Jacobian[20][20]=
	2 * gBet_r(m,n) + 2 * gBetr(m,n) - gAlp_r(m,n) 
    /*0*/  * pfv(m,n) + (exp(-2 * gconf(m,n)) * gAlp(m,n) * (-2
    /*2*/  * p(m,n) * Lt(m,n) + Lt(m,n) * (2 * p(m,n) + exp(2
    /*4*/  * gconf(m,n)) * (-(pfv_r(m,n) * gA(m,n) * r(m,n))
    /*4*/  + pfv(m,n) * (gA_r(m,n) * r(m,n) + 2 * gA(m,n) 
    /*5*/  * (-1 + gconf_r(m,n) * r(m,n)))) * Lt(m,n)))) 
    /*0*/  / (gA(m,n) * Lt(m,n) * r(m,n) * Lt(m,n));

Jacobian[20][21]=
	-gAlp_r(m,n);

Jacobian[20][22]=
	pfS_r(m,n) + (2 * pfS(m,n)) / r(m,n);

Jacobian[21][1]=
	(exp(-4 * gconf(m,n)) * (-2 * exp(2 * gconf(m,n))
    /*2*/  * gAlp(m,n) * gA(m,n) * p(m,n) * (2 * pftau(m,n) 
    /*3*/  + pftau_r(m,n) * r(m,n)) + 4 * gAlp_r(m,n) 
    /*2*/  * pfS(m,n) * r(m,n) * Lt(m,n))) / (pow2(gA(m,n)) 
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[21][5]=
	(exp(-4 * gconf(m,n)) * (-(exp(2 * gconf(m,n)) 
    /*3*/  * gAlp(m,n) * gA(m,n) * p(m,n) * (2 * pftau(m,n) 
    /*4*/  + pftau_r(m,n) * r(m,n))) + 2 * gAlp_r(m,n) 
    /*2*/  * pfS(m,n) * r(m,n) * Lt(m,n))) / (pow3(gA(m,n)) 
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[21][20]=
	gK1(m,n) * gAlp(m,n) * pfv(m,n) - (exp(-4 
    /*2*/  * gconf(m,n)) * gAlp_r(m,n)) / pow2(gA(m,n));

Jacobian[21][21]=
	gBet_r(m,n) + 2 * gBetr(m,n) - gAlp_r(m,n) 
    /*0*/  * pfv(m,n) + gAlp(m,n) * (-pfv_r(m,n) - (2 
    /*2*/  * pfv(m,n)) / r(m,n) + (2 * exp(-2 * gconf(m,n)) 
    /*2*/  * p(m,n) * (Lt(m,n) - Lt(m,n))) / (gA(m,n) * Lt(m,n)
    /*2*/  * r(m,n) * Lt(m,n)));

Jacobian[21][22]=
	pftau_r(m,n) + (2 * pftau(m,n)) / r(m,n);

Jacobian[22][22]=
	q_r(m,n);

Jacobian[22][23]=
	0.75;

Jacobian[23][1]=
	(-2 * exp(-2 * gconf(m,n)) * p(m,n) 
    /*1*/  * pow3(gAlp(m,n)) * (gA(m,n) * gA_rr(m,n) 
    /*2*/  * pow2(gB(m,n)) * pow2(r(m,n)) - pow2(gA_r(m,n)) 
    /*2*/  * pow2(gB(m,n)) * pow2(r(m,n)) - 2 * pow2(gA(m,n)) 
    /*2*/  * (pow2(gB(m,n)) - gB(m,n) * gB_rr(m,n) 
    /*3*/  * pow2(r(m,n)) + pow2(gB_r(m,n)) * pow2(r(m,n))) + 2
    /*2*/  * gL(m,n) * gA_r(m,n) * pow2(gB(m,n)) 
    /*2*/  * pow2(r(m,n)) * pow3(gA(m,n)) + Power(gA(m,n),4) 
    /*2*/  * (-6 + 3 * gL_r(m,n) * pow2(gB(m,n)) * pow2(r(m,n))
    /*3*/  + 4 * gL(m,n) * gB(m,n) * r(m,n) * (gB(m,n) 
    /*4*/  + gB_r(m,n) * r(m,n))))) / (3. * Power(gA(m,n),5) 
    /*1*/  * pow2(gB(m,n)) * pow2(r(m,n)) * Lt(m,n));

Jacobian[23][5]=
	-(exp(-2 * gconf(m,n)) * pow2(gAlp(m,n)) * (exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * Lt(m,n) * (gBet(m,n) 
    /*3*/  * (3 * gA(m,n) * gA_rr(m,n) * pow2(gB(m,n)) 
    /*4*/  * pow2(r(m,n)) - 4 * pow2(gA_r(m,n)) * pow2(gB(m,n))
    /*4*/  * pow2(r(m,n)) - 4 * pow2(gA(m,n)) 
    /*4*/  * (pow2(gB(m,n)) - gB(m,n) * gB_rr(m,n) 
    /*5*/  * pow2(r(m,n)) + pow2(gB_r(m,n)) * pow2(r(m,n)))) 
    /*3*/  + gA(m,n) * gB(m,n) * r(m,n) * (4 * gBet_r(m,n) 
    /*4*/  * gA(m,n) * gB(m,n) - 8 * gA1(m,n) * gAlp_r(m,n) 
    /*4*/  * gA(m,n) * gB(m,n) * r(m,n) + 8 * gA2(m,n) 
    /*4*/  * gAlp_r(m,n) * gA(m,n) * gB(m,n) * r(m,n) + 8 
    /*4*/  * gBet_rr(m,n) * gA(m,n) * gB(m,n) * r(m,n) + 3 
    /*4*/  * gBet_r(m,n) * gA_r(m,n) * gB(m,n) * r(m,n) + 4 
    /*4*/  * gBet_r(m,n) * gA(m,n) * gB_r(m,n) * r(m,n) + 2 
    /*4*/  * gL(m,n) * gB(m,n) * hShiD(gA(m,n)) * pow2(gA(m,n))
    /*4*/  * r(m,n))) * Lt(m,n) + gAlp(m,n) 
    /*2*/  * (Power(gA(m,n),4) * Lt(m,n) * p(m,n) * (-6 + 3 
    /*4*/  * gL_r(m,n) * pow2(gB(m,n)) * pow2(r(m,n)) + 4 
    /*4*/  * gL(m,n) * gB(m,n) * r(m,n) * (gB(m,n) + gB_r(m,n)
    /*5*/  * r(m,n))) + gA(m,n) * gA_rr(m,n) * p(m,n) 
    /*3*/  * pow2(gB(m,n)) * pow2(r(m,n)) * (4 * Lt(m,n) - 3 
    /*4*/  * Lt(m,n)) + p(m,n) * pow2(gA_r(m,n)) 
    /*3*/  * pow2(gB(m,n)) * pow2(r(m,n)) * (-5 * Lt(m,n) + 4 
    /*4*/  * Lt(m,n)) + 4 * gB(m,n) * Lt(m,n) * pow2(r(m,n)) 
    /*3*/  * pow3(gA(m,n)) * (2 * exp(2 * gconf(m,n)) 
    /*4*/  * gAsig(m,n) * pow2(r(m,n)) * (6 * gconf_r(m,n) 
    /*5*/  * gB(m,n) + gB_r(m,n) + gsig(m,n) * gB(m,n) 
    /*5*/  * r(m,n)) * Lt(m,n) + gB(m,n) * (gL(m,n) * gA_r(m,n)
    /*5*/  * p(m,n) - 2 * exp(2 * gconf(m,n)) * (gtrA_r(m,n)
    /*6*/  + gtrK_r(m,n)) * Lt(m,n))) + pow2(gA(m,n)) * (2 
    /*4*/  * gB(m,n) * gB_rr(m,n) * p(m,n) * pow2(r(m,n)) * (3
    /*5*/  * Lt(m,n) - 2 * Lt(m,n)) + 2 * p(m,n) 
    /*4*/  * pow2(gB_r(m,n)) * pow2(r(m,n)) * (-3 * Lt(m,n) + 2
    /*5*/  * Lt(m,n)) + pow2(gB(m,n)) * (4 * p(m,n) * Lt(m,n)
    /*5*/  - 6 * Lt(m,n) * (p(m,n) - 2 * exp(2 * gconf(m,n))
    /*6*/  * gAsig(m,n) * gA_r(m,n) * Power(r(m,n),4) 
    /*6*/  * Lt(m,n))))))) / (3. * Power(gA(m,n),6) * Lt(m,n) 
    /*1*/  * pow2(gB(m,n)) * pow2(r(m,n)) * Lt(m,n));

Jacobian[23][6]=
	(-2 * exp(-2 * gconf(m,n)) * pow2(gAlp(m,n)) 
    /*1*/  * (exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n) 
    /*2*/  * (gBet(m,n) * ((gB(m,n) * gB_rr(m,n) - 2 
    /*5*/  * pow2(gB_r(m,n))) * pow2(r(m,n)) + 2 
    /*4*/  * pow2(gA(m,n)) * (-3 + gL(m,n) * gB(m,n) 
    /*5*/  * gB_r(m,n) * pow2(r(m,n)))) + gBet_r(m,n) * r(m,n)
    /*3*/  * (6 * pow2(gA(m,n)) + gB(m,n) * gB_r(m,n) 
    /*4*/  * r(m,n))) * Lt(m,n) + gAlp(m,n) * (-(p(m,n) 
    /*4*/  * (-(gB(m,n) * gB_rr(m,n)) + 2 * pow2(gB_r(m,n))) 
    /*4*/  * pow2(r(m,n)) * (Lt(m,n) - Lt(m,n))) + 2 * p(m,n) 
    /*3*/  * pow2(gA(m,n)) * (-3 + gL(m,n) * gB(m,n) 
    /*4*/  * gB_r(m,n) * pow2(r(m,n))) * (Lt(m,n) - Lt(m,n)) 
    /*3*/  + 2 * exp(2 * gconf(m,n)) * gAsig(m,n) * gA(m,n) 
    /*3*/  * gB(m,n) * gB_r(m,n) * Lt(m,n) * Power(r(m,n),4) 
    /*3*/  * Lt(m,n)))) / (3. * Lt(m,n) * pow2(r(m,n)) 
    /*1*/  * pow3(gA(m,n)) * pow3(gB(m,n)) * Lt(m,n));

Jacobian[23][9]=
	(-4 * gAlp_r(m,n) * pow2(gAlp(m,n))) / (3. 
    /*1*/  * pow2(gA(m,n)));

Jacobian[23][10]=
	(4 * gAlp_r(m,n) * pow2(gAlp(m,n))) / (3. 
    /*1*/  * pow2(gA(m,n)));

Jacobian[23][13]=
	(exp(-2 * gconf(m,n)) * pow2(gAlp(m,n)) * (2 
    /*2*/  * gAlp(m,n) * p(m,n) * (gA_r(m,n) * gB(m,n) 
    /*3*/  * Lt(m,n) * r(m,n) + 2 * gA(m,n) * (gB_r(m,n) 
    /*4*/  * Lt(m,n) * r(m,n) + gB(m,n) * (Lt(m,n) - Lt(m,n))))
    /*2*/  + exp(2 * gconf(m,n)) * gA(m,n) * (4 * gBetr(m,n)
    /*3*/  * gA(m,n) * gB(m,n) + gB(m,n) * (-(gBet_r(m,n) 
    /*5*/  * gA(m,n)) + 2 * hShiD(gA(m,n))) + 4 * gA(m,n) 
    /*3*/  * hShiD(gB(m,n))) * Lt(m,n) * r(m,n) * Lt(m,n))) 
    /*0*/  / (3. * gB(m,n) * Lt(m,n) * pow2(gA(m,n)) * r(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[23][15]=
	(4 * gAsig(m,n) * pow3(gAlp(m,n)) * pow3(r(m,n)))
    /*0*/  / (3. * pow2(gA(m,n)));

Jacobian[23][17]=
	(4 * pow2(r(m,n)) * pow3(gAlp(m,n)) * (gA_r(m,n) 
    /*2*/  * gB(m,n) + gA(m,n) * (6 * gconf_r(m,n) * gB(m,n) 
    /*3*/  + gB_r(m,n)) + gsig(m,n) * gA(m,n) * gB(m,n) 
    /*2*/  * r(m,n))) / (3. * gB(m,n) * pow3(gA(m,n)));

Jacobian[23][22]=
	(3 * Bq_r(m,n) * Power(gA(m,n),4) * pow2(gB(m,n))
    /*1*/  * pow2(r(m,n)) + pow2(gAlp(m,n)) * (gA(m,n) 
    /*2*/  * gA_rr(m,n) * pow2(gB(m,n)) * pow2(r(m,n)) 
    /*2*/  - pow2(gA_r(m,n)) * pow2(gB(m,n)) * pow2(r(m,n)) - 2
    /*2*/  * pow2(gA(m,n)) * (pow2(gB(m,n)) - gB(m,n) 
    /*3*/  * gB_rr(m,n) * pow2(r(m,n)) + pow2(gB_r(m,n)) 
    /*3*/  * pow2(r(m,n))) + 2 * gL(m,n) * gA_r(m,n) 
    /*2*/  * pow2(gB(m,n)) * pow2(r(m,n)) * pow3(gA(m,n)) 
    /*2*/  + Power(gA(m,n),4) * (-6 + 4 * gL(m,n) * gB(m,n) 
    /*3*/  * r(m,n) * (gB(m,n) + gB_r(m,n) * r(m,n))))) / (3. 
    /*1*/  * Power(gA(m,n),4) * pow2(gB(m,n)) * pow2(r(m,n)));

Jacobian[23][23]=
	-eta;

