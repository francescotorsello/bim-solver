/** @file  DIRK_Jacobian_cBSSN_RegRed.h
 *  @author Francesco Torsello
 *  @brief The Jacobian of the regularized reduced cBSSN equations, needed by DIRK.
 *  @version 2019-06-18T14:56:13
 *  @image html DIRK_Jacobian_cBSSN_RegRed.png
 */

Jacobian[1][1]=
	(-2 * exp(-2 * gconf(m,n)) * gDconf(m,n) 
    /*1*/  * gAlp(m,n) * p(m,n)) / (gA(m,n) * Lt(m,n));

Jacobian[1][3]=
	gBet(m,n) + gAlp(m,n) * p(m,n) * (-(exp(-2 
    /*3*/  * gconf(m,n)) / (gA(m,n) * Lt(m,n))) + exp(-2 
    /*2*/  * gconf(m,n)) / (gA(m,n) * Lt(m,n)));

Jacobian[1][5]=
	-gAlp(m,n) / 6.;

Jacobian[1][7]=
	-((exp(-2 * gconf(m,n)) * gDconf(m,n) * gAlp(m,n)
    /*2*/  * p(m,n)) / (pow2(gA(m,n)) * Lt(m,n)));

Jacobian[2][2]=
	(2 * exp(-2 * fconf(m,n)) * fDconf(m,n) 
    /*1*/  * fAlp(m,n) * p(m,n)) / (fA(m,n) * Lt(m,n));

Jacobian[2][4]=
	gBet(m,n) - (exp(-2 * gconf(m,n)) * gAlp(m,n) 
    /*1*/  * p(m,n)) / (gA(m,n) * Lt(m,n)) - (exp(-2 
    /*2*/  * fconf(m,n)) * fAlp(m,n) * p(m,n)) / (fA(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[2][6]=
	-fAlp(m,n) / 6.;

Jacobian[2][11]=
	(exp(-2 * fconf(m,n)) * fDconf(m,n) * fAlp(m,n) 
    /*1*/  * p(m,n)) / (pow2(fA(m,n)) * Lt(m,n));

Jacobian[3][1]=
	(-2 * exp(-2 * gconf(m,n)) * gDconf(m,n) 
    /*1*/  * gDAlp(m,n) * p(m,n)) / (gA(m,n) * Lt(m,n)) 
    /*0*/  + gAlp(m,n) * ((-2 * exp(-2 * gconf(m,n)) 
    /*2*/  * gDconf(m,n) * p_r(m,n)) / (gA(m,n) * Lt(m,n)) + (2
    /*2*/  * exp(-2 * gconf(m,n)) * p(m,n) * (gDA(m,n) 
    /*3*/  * gDconf(m,n) * Lt(m,n) - gDconf_r(m,n) * Lt(m,n) 
    /*3*/  + 2 * pow2(gDconf(m,n)) * Lt(m,n) + gDconf(m,n) 
    /*3*/  * Lt_r(m,n))) / (gA(m,n) * pow2(Lt(m,n))));

Jacobian[3][3]=
	gBet_r(m,n) + p(m,n) * (-((exp(-2 * gconf(m,n)) 
    /*3*/  * gAlp_r(m,n)) / (gA(m,n) * Lt(m,n))) + (exp(-2 
    /*3*/  * gconf(m,n)) * gDAlp(m,n)) / (gA(m,n) * Lt(m,n))) 
    /*0*/  + gAlp(m,n) * (-((exp(-2 * gconf(m,n)) * p_r(m,n)) 
    /*2*/  / (gA(m,n) * Lt(m,n))) + (exp(-2 * gconf(m,n)) 
    /*2*/  * p_r(m,n)) / (gA(m,n) * Lt(m,n)) + p(m,n) * (((2 
    /*4*/  * exp(-2 * gconf(m,n)) * gconf_r(m,n)) / gA(m,n) 
    /*3*/  + (exp(-2 * gconf(m,n)) * gA_r(m,n)) 
    /*3*/  / pow2(gA(m,n))) / Lt(m,n) + (exp(-2 * gconf(m,n)) 
    /*3*/  * Lt_r(m,n)) / (gA(m,n) * pow2(Lt(m,n))) - (exp(-2 
    /*4*/  * gconf(m,n)) * (gDA(m,n) * Lt(m,n) + 4 
    /*4*/  * gDconf(m,n) * Lt(m,n) + Lt_r(m,n))) / (gA(m,n) 
    /*3*/  * pow2(Lt(m,n)))));

Jacobian[3][5]=
	-gDAlp(m,n) / 6.;

Jacobian[3][7]=
	-((exp(-2 * gconf(m,n)) * gDconf(m,n) * gDAlp(m,n)
    /*2*/  * p(m,n)) / (pow2(gA(m,n)) * Lt(m,n))) + gAlp(m,n)
    /*0*/  * (-((exp(-2 * gconf(m,n)) * gDconf(m,n) 
    /*3*/  * p_r(m,n)) / (pow2(gA(m,n)) * Lt(m,n))) + (exp(-2 
    /*3*/  * gconf(m,n)) * p(m,n) * (gDA(m,n) * gDconf(m,n) 
    /*3*/  * Lt(m,n) - gDconf_r(m,n) * Lt(m,n) + 2 
    /*3*/  * pow2(gDconf(m,n)) * Lt(m,n) + gDconf(m,n) 
    /*3*/  * Lt_r(m,n))) / (pow2(gA(m,n)) * pow2(Lt(m,n))));

Jacobian[3][9]=
	-((exp(-2 * gconf(m,n)) * gDconf(m,n) * gAlp(m,n)
    /*2*/  * p(m,n)) / (gA(m,n) * Lt(m,n)));

Jacobian[4][2]=
	(2 * exp(-2 * fconf(m,n)) * fDconf(m,n) 
    /*1*/  * fDAlp(m,n) * p(m,n)) / (fA(m,n) * Lt(m,n)) 
    /*0*/  + fAlp(m,n) * ((2 * exp(-2 * fconf(m,n)) 
    /*2*/  * fDconf(m,n) * p_r(m,n)) / (fA(m,n) * Lt(m,n)) - (2
    /*2*/  * exp(-2 * fconf(m,n)) * p(m,n) * (fDA(m,n) 
    /*3*/  * fDconf(m,n) * Lt(m,n) - fDconf_r(m,n) * Lt(m,n) 
    /*3*/  + 2 * pow2(fDconf(m,n)) * Lt(m,n) + fDconf(m,n) 
    /*3*/  * Lt_r(m,n))) / (fA(m,n) * pow2(Lt(m,n))));

Jacobian[4][4]=
	gBet_r(m,n) + gAlp(m,n) * (-((exp(-2 * gconf(m,n))
    /*3*/  * p_r(m,n)) / (gA(m,n) * Lt(m,n))) + p(m,n) * (((2
    /*4*/  * exp(-2 * gconf(m,n)) * gconf_r(m,n)) / gA(m,n) 
    /*3*/  + (exp(-2 * gconf(m,n)) * gA_r(m,n)) 
    /*3*/  / pow2(gA(m,n))) / Lt(m,n) + (exp(-2 * gconf(m,n)) 
    /*3*/  * Lt_r(m,n)) / (gA(m,n) * pow2(Lt(m,n))))) + p(m,n)
    /*0*/  * (-((exp(-2 * gconf(m,n)) * gAlp_r(m,n)) 
    /*2*/  / (gA(m,n) * Lt(m,n))) - (exp(-2 * fconf(m,n)) 
    /*2*/  * fDAlp(m,n)) / (fA(m,n) * Lt(m,n))) + fAlp(m,n) 
    /*0*/  * (-((exp(-2 * fconf(m,n)) * p_r(m,n)) / (fA(m,n) 
    /*3*/  * Lt(m,n))) + (exp(-2 * fconf(m,n)) * p(m,n) 
    /*2*/  * (fDA(m,n) * Lt(m,n) + 4 * fDconf(m,n) * Lt(m,n) 
    /*3*/  + Lt_r(m,n))) / (fA(m,n) * pow2(Lt(m,n))));

Jacobian[4][6]=
	-fDAlp(m,n) / 6.;

Jacobian[4][11]=
	(exp(-2 * fconf(m,n)) * fDconf(m,n) * fDAlp(m,n) 
    /*1*/  * p(m,n)) / (pow2(fA(m,n)) * Lt(m,n)) + fAlp(m,n) 
    /*0*/  * ((exp(-2 * fconf(m,n)) * fDconf(m,n) * p_r(m,n)) 
    /*1*/  / (pow2(fA(m,n)) * Lt(m,n)) - (exp(-2 * fconf(m,n))
    /*2*/  * p(m,n) * (fDA(m,n) * fDconf(m,n) * Lt(m,n) 
    /*3*/  - fDconf_r(m,n) * Lt(m,n) + 2 * pow2(fDconf(m,n)) 
    /*3*/  * Lt(m,n) + fDconf(m,n) * Lt_r(m,n))) 
    /*1*/  / (pow2(fA(m,n)) * pow2(Lt(m,n))));

Jacobian[4][13]=
	(exp(-2 * fconf(m,n)) * fDconf(m,n) * fAlp(m,n) 
    /*1*/  * p(m,n)) / (fA(m,n) * Lt(m,n));

Jacobian[5][1]=
	(exp(-4 * gconf(m,n)) * (4 * gDAlp_r(m,n) * r(m,n)
    /*2*/  + gDAlp(m,n) * (8 - 4 * gDA(m,n) * r(m,n) + 8 
    /*3*/  * gDB(m,n) * r(m,n) + 8 * gDconf(m,n) * r(m,n)))) 
    /*0*/  / (pow2(gA(m,n)) * r(m,n)) + k_g * (fAlp(m,n) * ((-2
    /*3*/  * exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) 
    /*3*/  * (2 * P_1_2(R(m,n)) + b_2)) / gA(m,n) - (2 
    /*3*/  * (P_2_1(R(m,n)) + b_1)) / Lt(m,n)) + gAlp(m,n) 
    /*1*/  * (-2 * (P_1_0(R(m,n)) + b_0) - (exp(2 * fconf(m,n)
    /*4*/  - 2 * gconf(m,n)) * fA(m,n) * ((2 + 4 
    /*5*/  * pow2(Lt(m,n))) * P_1_1(R(m,n)) + (3 - 6 
    /*5*/  * pow2(Lt(m,n))) * P_2_1(R(m,n)) + 2 * b_1)) 
    /*2*/  / (gA(m,n) * Lt(m,n)))) - (2 * exp(-2 * gconf(m,n))
    /*1*/  * gAlp(m,n) * p(m,n) * gtrK_r(m,n)) / (gA(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[5][2]=
	k_g * (fAlp(m,n) * ((2 * exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * fA(m,n) * (2 * P_1_2(R(m,n)) + b_2))
    /*2*/  / gA(m,n) + (2 * (P_2_1(R(m,n)) + b_1)) / Lt(m,n))
    /*1*/  + gAlp(m,n) * (2 * (P_1_0(R(m,n)) + b_0) + (exp(2
    /*4*/  * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * ((2 + 4
    /*5*/  * pow2(Lt(m,n))) * P_1_1(R(m,n)) + (3 - 6 
    /*5*/  * pow2(Lt(m,n))) * P_2_1(R(m,n)) + 2 * b_1)) 
    /*2*/  / (gA(m,n) * Lt(m,n))));

Jacobian[5][3]=
	(-2 * exp(-4 * gconf(m,n)) * gDAlp(m,n)) 
    /*0*/  / pow2(gA(m,n));

Jacobian[5][5]=
	(2 * gAlp(m,n) * (gtrA(m,n) + gtrK(m,n))) / 3.;

Jacobian[5][7]=
	(exp(-4 * gconf(m,n)) * (2 * gDAlp_r(m,n) * r(m,n)
    /*2*/  + gDAlp(m,n) * (4 - 2 * gDA(m,n) * r(m,n) + 4 
    /*3*/  * gDB(m,n) * r(m,n) + 4 * gDconf(m,n) * r(m,n)))) 
    /*0*/  / (pow3(gA(m,n)) * r(m,n)) + k_g * (-((exp(2 
    /*4*/  * fconf(m,n) - 2 * gconf(m,n)) * fAlp(m,n) * fA(m,n)
    /*3*/  * P_1_2(R(m,n))) / pow2(gA(m,n))) - (exp(2 
    /*3*/  * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * gAlp(m,n)
    /*2*/  * (2 * P_1_1(R(m,n)) + (1 - 2 * pow2(Lt(m,n))) 
    /*3*/  * P_2_1(R(m,n)))) / (2. * pow2(gA(m,n)) * Lt(m,n)))
    /*0*/  - (exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) 
    /*1*/  * gtrK_r(m,n)) / (pow2(gA(m,n)) * Lt(m,n));

Jacobian[5][8]=
	k_g * (fAlp(m,n) * ((exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * fA(m,n) * R(m,n) * b_3) / (gA(m,n) 
    /*3*/  * gB(m,n)) + (R(m,n) * (-P_1_2(R(m,n)) + b_2)) 
    /*2*/  / (gB(m,n) * Lt(m,n))) + gAlp(m,n) * ((R(m,n) * b_1)
    /*2*/  / gB(m,n) + (exp(2 * fconf(m,n) - 2 * gconf(m,n))
    /*3*/  * fA(m,n) * R(m,n) * ((-1 + 2 * pow2(Lt(m,n))) 
    /*4*/  * P_1_2(R(m,n)) + b_2)) / (gA(m,n) * gB(m,n) 
    /*3*/  * Lt(m,n))));

Jacobian[5][9]=
	(exp(-4 * gconf(m,n)) * gDAlp(m,n)) 
    /*0*/  / pow2(gA(m,n));

Jacobian[5][10]=
	(-2 * exp(-4 * gconf(m,n)) * gDAlp(m,n)) 
    /*0*/  / pow2(gA(m,n));

Jacobian[5][11]=
	k_g * ((exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*2*/  * fAlp(m,n) * P_1_2(R(m,n))) / gA(m,n) + (exp(2 
    /*3*/  * fconf(m,n) - 2 * gconf(m,n)) * gAlp(m,n) * (2 
    /*3*/  * P_1_1(R(m,n)) + (1 - 2 * pow2(Lt(m,n))) 
    /*3*/  * P_2_1(R(m,n)))) / (2. * gA(m,n) * Lt(m,n)));

Jacobian[5][12]=
	k_g * (fAlp(m,n) * (-((exp(4 * fconf(m,n) - 4 
    /*5*/  * gconf(m,n)) * fA(m,n) * b_3) / (gA(m,n) 
    /*4*/  * gB(m,n))) + (exp(2 * fconf(m,n) - 2 * gconf(m,n))
    /*3*/  * (P_1_2(R(m,n)) - b_2)) / (gB(m,n) * Lt(m,n))) 
    /*1*/  + gAlp(m,n) * (-((exp(2 * fconf(m,n) - 2 
    /*5*/  * gconf(m,n)) * b_1) / gB(m,n)) - (exp(4 
    /*4*/  * fconf(m,n) - 4 * gconf(m,n)) * fA(m,n) * ((-1 + 2
    /*5*/  * pow2(Lt(m,n))) * P_1_2(R(m,n)) + b_2)) 
    /*2*/  / (gA(m,n) * gB(m,n) * Lt(m,n))));

Jacobian[5][15]=
	2 * gA1(m,n) * gAlp(m,n);

Jacobian[5][16]=
	4 * gA2(m,n) * gAlp(m,n);

Jacobian[6][1]=
	k_f * (gAlp(m,n) * ((2 * exp(-2 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (P_2_0(R(m,n)) + b_0)) 
    /*2*/  / (fA(m,n) * pow2(R(m,n))) + (2 * (2 * P_1_1(R(m,n))
    /*4*/  + b_1)) / (pow2(R(m,n)) * Lt(m,n))) + fAlp(m,n) 
    /*1*/  * ((2 * (P_1_2(R(m,n)) + b_2)) / pow2(R(m,n)) 
    /*2*/  + (exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) 
    /*3*/  * ((2 - 4 * pow2(Lt(m,n))) * P_1_1(R(m,n)) + (3 - 2
    /*5*/  * pow2(Lt(m,n))) * P_2_1(R(m,n)) + 2 * b_1)) 
    /*2*/  / (fA(m,n) * pow2(R(m,n)) * Lt(m,n))));

Jacobian[6][2]=
	(exp(-4 * fconf(m,n)) * (4 * fDAlp_r(m,n) * r(m,n)
    /*2*/  + fDAlp(m,n) * (8 - 4 * fDA(m,n) * r(m,n) + 8 
    /*3*/  * fDB(m,n) * r(m,n) + 8 * fDconf(m,n) * r(m,n)))) 
    /*0*/  / (pow2(fA(m,n)) * r(m,n)) + k_f * (fAlp(m,n) * ((-2
    /*3*/  * (P_1_2(R(m,n)) + b_2)) / pow2(R(m,n)) + (exp(-2
    /*4*/  * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) * ((-2 
    /*5*/  + 4 * pow2(Lt(m,n))) * P_1_1(R(m,n)) + (-3 + 2 
    /*5*/  * pow2(Lt(m,n))) * P_2_1(R(m,n)) - 2 * b_1)) 
    /*2*/  / (fA(m,n) * pow2(R(m,n)) * Lt(m,n))) + gAlp(m,n) 
    /*1*/  * ((-2 * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*3*/  * gA(m,n) * (P_2_0(R(m,n)) + b_0)) / (fA(m,n) 
    /*3*/  * pow2(R(m,n))) - (2 * (2 * P_1_1(R(m,n)) + b_1)) 
    /*2*/  / (pow2(R(m,n)) * Lt(m,n)))) + (2 * exp(-2 
    /*2*/  * fconf(m,n)) * fAlp(m,n) * p(m,n) * ftrK_r(m,n)) 
    /*0*/  / (fA(m,n) * Lt(m,n));

Jacobian[6][4]=
	(-2 * exp(-4 * fconf(m,n)) * fDAlp(m,n)) 
    /*0*/  / pow2(fA(m,n));

Jacobian[6][6]=
	(2 * fAlp(m,n) * (ftrA(m,n) + ftrK(m,n))) / 3.;

Jacobian[6][7]=
	k_f * ((exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*2*/  * gAlp(m,n) * (-P_1_0(R(m,n)) + P_2_0(R(m,n)))) 
    /*1*/  / (fA(m,n) * pow2(R(m,n))) + (exp(-2 * fconf(m,n) 
    /*3*/  + 2 * gconf(m,n)) * fAlp(m,n) * (-2 * P_1_1(R(m,n))
    /*3*/  + (3 - 2 * pow2(Lt(m,n))) * P_2_1(R(m,n)))) / (2.
    /*2*/  * fA(m,n) * pow2(R(m,n)) * Lt(m,n)));

Jacobian[6][8]=
	k_f * (fAlp(m,n) * (-(b_3 / (gB(m,n) * R(m,n))) 
    /*2*/  + (exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) 
    /*3*/  * (-2 * (-1 + pow2(Lt(m,n))) * P_1_1(R(m,n)) + b_1))
    /*2*/  / (fA(m,n) * gB(m,n) * pow2(R(m,n)) * Lt(m,n))) 
    /*1*/  + gAlp(m,n) * ((exp(-2 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (P_1_0(R(m,n)) + b_0)) 
    /*2*/  / (fA(m,n) * gB(m,n) * pow2(R(m,n))) - (b_1 + 2 
    /*3*/  * R(m,n) * b_2) / (gB(m,n) * pow2(R(m,n)) 
    /*3*/  * Lt(m,n))));

Jacobian[6][11]=
	(exp(-4 * fconf(m,n)) * (2 * fDAlp_r(m,n) * r(m,n)
    /*2*/  + fDAlp(m,n) * (4 - 2 * fDA(m,n) * r(m,n) + 4 
    /*3*/  * fDB(m,n) * r(m,n) + 4 * fDconf(m,n) * r(m,n)))) 
    /*0*/  / (pow3(fA(m,n)) * r(m,n)) + k_f * (-((exp(-2 
    /*4*/  * fconf(m,n) + 2 * gconf(m,n)) * gAlp(m,n) * gA(m,n)
    /*3*/  * P_1_1(R(m,n))) / (pow2(fA(m,n)) * R(m,n))) 
    /*1*/  + (exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * fAlp(m,n)
    /*2*/  * gA(m,n) * (2 * P_1_1(R(m,n)) + (-3 + 2 
    /*4*/  * pow2(Lt(m,n))) * P_2_1(R(m,n)))) / (2. 
    /*2*/  * pow2(fA(m,n)) * pow2(R(m,n)) * Lt(m,n))) + (exp(-2
    /*2*/  * fconf(m,n)) * fAlp(m,n) * p(m,n) * ftrK_r(m,n))
    /*0*/  / (pow2(fA(m,n)) * Lt(m,n));

Jacobian[6][12]=
	k_f * (fAlp(m,n) * ((exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * b_3) / (gB(m,n) * pow2(R(m,n))) 
    /*2*/  + (gA(m,n) * (2 * (-1 + pow2(Lt(m,n))) 
    /*4*/  * P_1_1(R(m,n)) - b_1)) / (fA(m,n) * gB(m,n) 
    /*3*/  * pow3(R(m,n)) * Lt(m,n))) + gAlp(m,n) * (-((gA(m,n)
    /*4*/  * (P_1_0(R(m,n)) + b_0)) / (fA(m,n) * gB(m,n) 
    /*4*/  * pow3(R(m,n)))) + (exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * (b_1 + 2 * R(m,n) * b_2)) / (gB(m,n)
    /*3*/  * pow3(R(m,n)) * Lt(m,n))));

Jacobian[6][13]=
	(exp(-4 * fconf(m,n)) * fDAlp(m,n)) 
    /*0*/  / pow2(fA(m,n));

Jacobian[6][14]=
	(-2 * exp(-4 * fconf(m,n)) * fDAlp(m,n)) 
    /*0*/  / pow2(fA(m,n));

Jacobian[6][17]=
	2 * fA1(m,n) * fAlp(m,n);

Jacobian[6][18]=
	4 * fA2(m,n) * fAlp(m,n);

Jacobian[7][1]=
	(-2 * exp(-2 * gconf(m,n)) * gDA(m,n) * gAlp(m,n)
    /*1*/  * p(m,n)) / Lt(m,n);

Jacobian[7][7]=
	gdet_pff(m,n) / (6. * gdet(m,n)) + gDA(m,n) 
    /*0*/  * gBet(m,n) + gBet_r(m,n) + gAlp(m,n) * (-((exp(-2 
    /*4*/  * gconf(m,n)) * gDA(m,n) * p(m,n)) / (gA(m,n) 
    /*3*/  * Lt(m,n))) + (-3 * gA1(m,n) + gtrA(m,n)) / 3.);

Jacobian[7][9]=
	gBet(m,n) * gA(m,n) + gAlp(m,n) * p(m,n) 
    /*0*/  * (-(exp(-2 * gconf(m,n)) / Lt(m,n)) + exp(-2 
    /*2*/  * gconf(m,n)) / Lt(m,n));

Jacobian[7][15]=
	-(gAlp(m,n) * gA(m,n));

Jacobian[8][1]=
	(-2 * exp(-2 * gconf(m,n)) * gAlp(m,n) * gB(m,n) 
    /*1*/  * p(m,n) * (1 + gDB(m,n) * r(m,n))) / (gA(m,n) 
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[8][7]=
	-((exp(-2 * gconf(m,n)) * gAlp(m,n) * gB(m,n) 
    /*2*/  * p(m,n) * (1 + gDB(m,n) * r(m,n))) / (pow2(gA(m,n))
    /*2*/  * r(m,n) * Lt(m,n)));

Jacobian[8][8]=
	gdet_pff(m,n) / (6. * gdet(m,n)) + gDB(m,n) 
    /*0*/  * gBet(m,n) + gBetr(m,n) + gAlp(m,n) * ((-3 
    /*2*/  * gA2(m,n) + gtrA(m,n)) / 3. + p(m,n) * (-((exp(-2 
    /*5*/  * gconf(m,n)) * (1 + gDB(m,n) * r(m,n))) / (gA(m,n)
    /*4*/  * Lt(m,n) * r(m,n))) + (exp(-2 * gconf(m,n)) * (1
    /*4*/  + gDB(m,n) * r(m,n))) / (gA(m,n) * r(m,n) 
    /*3*/  * Lt(m,n))));

Jacobian[8][10]=
	gBet(m,n) * gB(m,n) + gAlp(m,n) * p(m,n) 
    /*0*/  * (-((exp(-2 * gconf(m,n)) * gB(m,n)) / (gA(m,n) 
    /*3*/  * Lt(m,n))) + (exp(-2 * gconf(m,n)) * gB(m,n)) 
    /*1*/  / (gA(m,n) * Lt(m,n)));

Jacobian[8][16]=
	-(gAlp(m,n) * gB(m,n));

Jacobian[9][1]=
	(-2 * exp(-2 * gconf(m,n)) * gDA(m,n) * gDAlp(m,n)
    /*1*/  * p(m,n)) / (gA(m,n) * Lt(m,n)) + gAlp(m,n) * ((-2
    /*2*/  * exp(-2 * gconf(m,n)) * gDA(m,n) * p_r(m,n)) 
    /*1*/  / (gA(m,n) * Lt(m,n)) + (2 * exp(-2 * gconf(m,n)) 
    /*2*/  * p(m,n) * (-(gDA_r(m,n) * Lt(m,n)) + pow2(gDA(m,n))
    /*3*/  * Lt(m,n) + gDA(m,n) * (2 * gDconf(m,n) * Lt(m,n)
    /*4*/  + Lt_r(m,n)))) / (gA(m,n) * pow2(Lt(m,n))));

Jacobian[9][3]=
	(-2 * exp(-2 * gconf(m,n)) * gDA(m,n) * gAlp(m,n)
    /*1*/  * p(m,n)) / (gA(m,n) * Lt(m,n));

Jacobian[9][7]=
	-((exp(-2 * gconf(m,n)) * gDA(m,n) * gDAlp(m,n) 
    /*2*/  * p(m,n)) / (pow2(gA(m,n)) * Lt(m,n))) + gAlp(m,n) 
    /*0*/  * (-((exp(-2 * gconf(m,n)) * gDA(m,n) * p_r(m,n)) 
    /*2*/  / (pow2(gA(m,n)) * Lt(m,n))) + (exp(-2 * gconf(m,n))
    /*2*/  * p(m,n) * (-(gDA_r(m,n) * Lt(m,n)) 
    /*3*/  + pow2(gDA(m,n)) * Lt(m,n) + gDA(m,n) * (2 
    /*4*/  * gDconf(m,n) * Lt(m,n) + Lt_r(m,n)))) 
    /*1*/  / (pow2(gA(m,n)) * pow2(Lt(m,n))));

Jacobian[9][9]=
	gBet_r(m,n) + p(m,n) * (-((exp(-2 * gconf(m,n)) 
    /*3*/  * gAlp_r(m,n)) / (gA(m,n) * Lt(m,n))) + (exp(-2 
    /*3*/  * gconf(m,n)) * gDAlp(m,n)) / (gA(m,n) * Lt(m,n))) 
    /*0*/  + gAlp(m,n) * (-((exp(-2 * gconf(m,n)) * p_r(m,n)) 
    /*2*/  / (gA(m,n) * Lt(m,n))) + (exp(-2 * gconf(m,n)) 
    /*2*/  * p_r(m,n)) / (gA(m,n) * Lt(m,n)) + p(m,n) * (((2 
    /*4*/  * exp(-2 * gconf(m,n)) * gconf_r(m,n)) / gA(m,n) 
    /*3*/  + (exp(-2 * gconf(m,n)) * gA_r(m,n)) 
    /*3*/  / pow2(gA(m,n))) / Lt(m,n) + (exp(-2 * gconf(m,n)) 
    /*3*/  * Lt_r(m,n)) / (gA(m,n) * pow2(Lt(m,n))) - (exp(-2 
    /*4*/  * gconf(m,n)) * (2 * gDA(m,n) * Lt(m,n) + 2 
    /*4*/  * gDconf(m,n) * Lt(m,n) + Lt_r(m,n))) / (gA(m,n) 
    /*3*/  * pow2(Lt(m,n)))));

Jacobian[9][15]=
	-gDAlp(m,n);

Jacobian[10][1]=
	(-2 * exp(-2 * gconf(m,n)) * gDAlp(m,n) * p(m,n) 
    /*1*/  * (1 + gDB(m,n) * r(m,n))) / (gA(m,n) * r(m,n) 
    /*1*/  * Lt(m,n)) + gAlp(m,n) * ((-2 * exp(-2 * gconf(m,n))
    /*2*/  * p_r(m,n) * (1 + gDB(m,n) * r(m,n))) / (gA(m,n) 
    /*2*/  * r(m,n) * Lt(m,n)) + (2 * exp(-2 * gconf(m,n)) 
    /*2*/  * p(m,n) * ((1 - gDB_r(m,n) * pow2(r(m,n)) 
    /*4*/  + gDA(m,n) * r(m,n) * (1 + gDB(m,n) * r(m,n)) + 2 
    /*4*/  * gDconf(m,n) * r(m,n) * (1 + gDB(m,n) * r(m,n))) 
    /*3*/  * Lt(m,n) + r(m,n) * (1 + gDB(m,n) * r(m,n)) 
    /*3*/  * Lt_r(m,n))) / (gA(m,n) * pow2(r(m,n)) 
    /*2*/  * pow2(Lt(m,n))));

Jacobian[10][3]=
	(-2 * exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) 
    /*1*/  * (1 + gDB(m,n) * r(m,n))) / (gA(m,n) * r(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[10][7]=
	-((exp(-2 * gconf(m,n)) * gDAlp(m,n) * p(m,n) * (1
    /*3*/  + gDB(m,n) * r(m,n))) / (pow2(gA(m,n)) * r(m,n) 
    /*2*/  * Lt(m,n))) + gAlp(m,n) * (-((exp(-2 * gconf(m,n)) 
    /*3*/  * p_r(m,n) * (1 + gDB(m,n) * r(m,n))) 
    /*2*/  / (pow2(gA(m,n)) * r(m,n) * Lt(m,n))) + (exp(-2 
    /*3*/  * gconf(m,n)) * p(m,n) * ((1 - gDB_r(m,n) 
    /*4*/  * pow2(r(m,n)) + gDA(m,n) * r(m,n) * (1 + gDB(m,n) 
    /*5*/  * r(m,n)) + 2 * gDconf(m,n) * r(m,n) * (1 + gDB(m,n)
    /*5*/  * r(m,n))) * Lt(m,n) + r(m,n) * (1 + gDB(m,n) 
    /*4*/  * r(m,n)) * Lt_r(m,n))) / (pow2(gA(m,n)) 
    /*2*/  * pow2(r(m,n)) * pow2(Lt(m,n))));

Jacobian[10][9]=
	-((exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) * (1
    /*3*/  + gDB(m,n) * r(m,n))) / (gA(m,n) * r(m,n) 
    /*2*/  * Lt(m,n)));

Jacobian[10][10]=
	gBet_r(m,n) + p(m,n) * (-((exp(-2 * gconf(m,n)) 
    /*3*/  * gAlp_r(m,n)) / (gA(m,n) * Lt(m,n))) + (exp(-2 
    /*3*/  * gconf(m,n)) * gDAlp(m,n)) / (gA(m,n) * Lt(m,n))) 
    /*0*/  + gAlp(m,n) * (-((exp(-2 * gconf(m,n)) * p_r(m,n)) 
    /*2*/  / (gA(m,n) * Lt(m,n))) + (exp(-2 * gconf(m,n)) 
    /*2*/  * p_r(m,n)) / (gA(m,n) * Lt(m,n)) + p(m,n) * (((2 
    /*4*/  * exp(-2 * gconf(m,n)) * gconf_r(m,n)) / gA(m,n) 
    /*3*/  + (exp(-2 * gconf(m,n)) * gA_r(m,n)) 
    /*3*/  / pow2(gA(m,n))) / Lt(m,n) + (exp(-2 * gconf(m,n)) 
    /*3*/  * Lt_r(m,n)) / (gA(m,n) * pow2(Lt(m,n))) - (exp(-2 
    /*4*/  * gconf(m,n)) * (gDA(m,n) * Lt(m,n) + 2 
    /*4*/  * gDconf(m,n) * Lt(m,n) + Lt_r(m,n))) / (gA(m,n) 
    /*3*/  * pow2(Lt(m,n)))));

Jacobian[10][16]=
	-gDAlp(m,n);

Jacobian[11][2]=
	(2 * exp(-2 * fconf(m,n)) * fDA(m,n) * fAlp(m,n) 
    /*1*/  * p(m,n)) / Lt(m,n);

Jacobian[11][11]=
	fdet_pff(m,n) / (6. * fdet(m,n)) + fBet_r(m,n) 
    /*0*/  + fDA(m,n) * gBet(m,n) - (exp(-2 * gconf(m,n)) 
    /*1*/  * fDA(m,n) * gAlp(m,n) * p(m,n)) / (gA(m,n) 
    /*1*/  * Lt(m,n)) + (fAlp(m,n) * (-3 * fA1(m,n) 
    /*2*/  + ftrA(m,n))) / 3.;

Jacobian[11][13]=
	fA(m,n) * gBet(m,n) - (exp(-2 * gconf(m,n)) 
    /*1*/  * fA(m,n) * gAlp(m,n) * p(m,n)) / (gA(m,n) 
    /*1*/  * Lt(m,n)) - (exp(-2 * fconf(m,n)) * fAlp(m,n) 
    /*1*/  * p(m,n)) / Lt(m,n);

Jacobian[11][17]=
	-(fAlp(m,n) * fA(m,n));

Jacobian[12][2]=
	(2 * exp(-4 * fconf(m,n) + 2 * gconf(m,n)) 
    /*1*/  * fAlp(m,n) * gB(m,n) * p(m,n) * (1 + fDB(m,n) 
    /*2*/  * r(m,n)) * R(m,n)) / (fA(m,n) * r(m,n) * Lt(m,n));

Jacobian[12][11]=
	(exp(-4 * fconf(m,n) + 2 * gconf(m,n)) * fAlp(m,n)
    /*1*/  * gB(m,n) * p(m,n) * (1 + fDB(m,n) * r(m,n)) 
    /*1*/  * R(m,n)) / (pow2(fA(m,n)) * r(m,n) * Lt(m,n));

Jacobian[12][12]=
	fdet_pff(m,n) / (6. * fdet(m,n)) + fDB(m,n) 
    /*0*/  * gBet(m,n) + gBetr(m,n) - (exp(-2 * gconf(m,n)) 
    /*1*/  * gAlp(m,n) * p(m,n) * (1 + fDB(m,n) * r(m,n))) 
    /*0*/  / (gA(m,n) * Lt(m,n) * r(m,n)) + fAlp(m,n) * ((-3 
    /*2*/  * fA2(m,n) + ftrA(m,n)) / 3. - (exp(-2 * fconf(m,n))
    /*2*/  * p(m,n) * (1 + fDB(m,n) * r(m,n))) / (fA(m,n) 
    /*2*/  * r(m,n) * Lt(m,n)));

Jacobian[12][14]=
	exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gBet(m,n)
    /*0*/  * gB(m,n) * R(m,n) - (exp(-2 * fconf(m,n)) 
    /*1*/  * gAlp(m,n) * gB(m,n) * p(m,n) * R(m,n)) / (gA(m,n)
    /*1*/  * Lt(m,n)) - (exp(-4 * fconf(m,n) + 2 
    /*2*/  * gconf(m,n)) * fAlp(m,n) * gB(m,n) * p(m,n) 
    /*1*/  * R(m,n)) / (fA(m,n) * Lt(m,n));

Jacobian[12][18]=
	-(exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*1*/  * fAlp(m,n) * gB(m,n) * R(m,n));

Jacobian[13][2]=
	(2 * exp(-2 * fconf(m,n)) * fDA(m,n) * fDAlp(m,n)
    /*1*/  * p(m,n)) / (fA(m,n) * Lt(m,n)) + fAlp(m,n) * ((2
    /*2*/  * exp(-2 * fconf(m,n)) * fDA(m,n) * p_r(m,n)) 
    /*1*/  / (fA(m,n) * Lt(m,n)) - (2 * exp(-2 * fconf(m,n)) 
    /*2*/  * p(m,n) * (-(fDA_r(m,n) * Lt(m,n)) + pow2(fDA(m,n))
    /*3*/  * Lt(m,n) + fDA(m,n) * (2 * fDconf(m,n) * Lt(m,n)
    /*4*/  + Lt_r(m,n)))) / (fA(m,n) * pow2(Lt(m,n))));

Jacobian[13][4]=
	(2 * exp(-2 * fconf(m,n)) * fDA(m,n) * fAlp(m,n) 
    /*1*/  * p(m,n)) / (fA(m,n) * Lt(m,n));

Jacobian[13][11]=
	(exp(-2 * fconf(m,n)) * fDA(m,n) * fDAlp(m,n) 
    /*1*/  * p(m,n)) / (pow2(fA(m,n)) * Lt(m,n)) + fAlp(m,n) 
    /*0*/  * ((exp(-2 * fconf(m,n)) * fDA(m,n) * p_r(m,n)) 
    /*1*/  / (pow2(fA(m,n)) * Lt(m,n)) + (exp(-2 * fconf(m,n))
    /*2*/  * p(m,n) * (fDA_r(m,n) * Lt(m,n) - pow2(fDA(m,n))
    /*3*/  * Lt(m,n) - fDA(m,n) * (2 * fDconf(m,n) * Lt(m,n)
    /*4*/  + Lt_r(m,n)))) / (pow2(fA(m,n)) * pow2(Lt(m,n))));

Jacobian[13][13]=
	gBet_r(m,n) + gAlp(m,n) * (-((exp(-2 * gconf(m,n))
    /*3*/  * p_r(m,n)) / (gA(m,n) * Lt(m,n))) + p(m,n) * (((2
    /*4*/  * exp(-2 * gconf(m,n)) * gconf_r(m,n)) / gA(m,n) 
    /*3*/  + (exp(-2 * gconf(m,n)) * gA_r(m,n)) 
    /*3*/  / pow2(gA(m,n))) / Lt(m,n) + (exp(-2 * gconf(m,n)) 
    /*3*/  * Lt_r(m,n)) / (gA(m,n) * pow2(Lt(m,n))))) + p(m,n)
    /*0*/  * (-((exp(-2 * gconf(m,n)) * gAlp_r(m,n)) 
    /*2*/  / (gA(m,n) * Lt(m,n))) - (exp(-2 * fconf(m,n)) 
    /*2*/  * fDAlp(m,n)) / (fA(m,n) * Lt(m,n))) + fAlp(m,n) 
    /*0*/  * (-((exp(-2 * fconf(m,n)) * p_r(m,n)) / (fA(m,n) 
    /*3*/  * Lt(m,n))) + (exp(-2 * fconf(m,n)) * p(m,n) * (2 
    /*3*/  * fDA(m,n) * Lt(m,n) + 2 * fDconf(m,n) * Lt(m,n) 
    /*3*/  + Lt_r(m,n))) / (fA(m,n) * pow2(Lt(m,n))));

Jacobian[13][17]=
	-fDAlp(m,n);

Jacobian[14][2]=
	(2 * exp(-2 * fconf(m,n)) * fDAlp(m,n) * p(m,n) 
    /*1*/  * (1 + fDB(m,n) * r(m,n))) / (fA(m,n) * r(m,n) 
    /*1*/  * Lt(m,n)) + fAlp(m,n) * ((2 * exp(-2 * fconf(m,n))
    /*2*/  * p_r(m,n) * (1 + fDB(m,n) * r(m,n))) / (fA(m,n) 
    /*2*/  * r(m,n) * Lt(m,n)) - (2 * exp(-2 * fconf(m,n)) 
    /*2*/  * p(m,n) * ((1 - fDB_r(m,n) * pow2(r(m,n)) 
    /*4*/  + fDA(m,n) * r(m,n) * (1 + fDB(m,n) * r(m,n)) + 2 
    /*4*/  * fDconf(m,n) * r(m,n) * (1 + fDB(m,n) * r(m,n))) 
    /*3*/  * Lt(m,n) + r(m,n) * (1 + fDB(m,n) * r(m,n)) 
    /*3*/  * Lt_r(m,n))) / (fA(m,n) * pow2(r(m,n)) 
    /*2*/  * pow2(Lt(m,n))));

Jacobian[14][4]=
	(2 * exp(-2 * fconf(m,n)) * fAlp(m,n) * p(m,n) 
    /*1*/  * (1 + fDB(m,n) * r(m,n))) / (fA(m,n) * r(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[14][11]=
	(exp(-2 * fconf(m,n)) * fDAlp(m,n) * p(m,n) * (1 
    /*2*/  + fDB(m,n) * r(m,n))) / (pow2(fA(m,n)) * r(m,n) 
    /*1*/  * Lt(m,n)) + fAlp(m,n) * ((exp(-2 * fconf(m,n)) 
    /*2*/  * p_r(m,n) * (1 + fDB(m,n) * r(m,n))) 
    /*1*/  / (pow2(fA(m,n)) * r(m,n) * Lt(m,n)) + (exp(-2 
    /*3*/  * fconf(m,n)) * p(m,n) * (-((1 - fDB_r(m,n) 
    /*5*/  * pow2(r(m,n)) + fDA(m,n) * r(m,n) * (1 + fDB(m,n) 
    /*6*/  * r(m,n)) + 2 * fDconf(m,n) * r(m,n) * (1 + fDB(m,n)
    /*6*/  * r(m,n))) * Lt(m,n)) - r(m,n) * (1 + fDB(m,n) 
    /*4*/  * r(m,n)) * Lt_r(m,n))) / (pow2(fA(m,n)) 
    /*2*/  * pow2(r(m,n)) * pow2(Lt(m,n))));

Jacobian[14][13]=
	(exp(-2 * fconf(m,n)) * fAlp(m,n) * p(m,n) * (1 
    /*2*/  + fDB(m,n) * r(m,n))) / (fA(m,n) * r(m,n) * Lt(m,n));

Jacobian[14][14]=
	gBet_r(m,n) + gAlp(m,n) * (-((exp(-2 * gconf(m,n))
    /*3*/  * p_r(m,n)) / (gA(m,n) * Lt(m,n))) + p(m,n) * (((2
    /*4*/  * exp(-2 * gconf(m,n)) * gconf_r(m,n)) / gA(m,n) 
    /*3*/  + (exp(-2 * gconf(m,n)) * gA_r(m,n)) 
    /*3*/  / pow2(gA(m,n))) / Lt(m,n) + (exp(-2 * gconf(m,n)) 
    /*3*/  * Lt_r(m,n)) / (gA(m,n) * pow2(Lt(m,n))))) + p(m,n)
    /*0*/  * (-((exp(-2 * gconf(m,n)) * gAlp_r(m,n)) 
    /*2*/  / (gA(m,n) * Lt(m,n))) - (exp(-2 * fconf(m,n)) 
    /*2*/  * fDAlp(m,n)) / (fA(m,n) * Lt(m,n))) + fAlp(m,n) 
    /*0*/  * (-((exp(-2 * fconf(m,n)) * p_r(m,n)) / (fA(m,n) 
    /*3*/  * Lt(m,n))) + (exp(-2 * fconf(m,n)) * p(m,n) 
    /*2*/  * (fDA(m,n) * Lt(m,n) + 2 * fDconf(m,n) * Lt(m,n) 
    /*3*/  + Lt_r(m,n))) / (fA(m,n) * pow2(Lt(m,n))));

Jacobian[14][18]=
	-fDAlp(m,n);

Jacobian[15][1]=
	(-8 * exp(-4 * gconf(m,n)) * (-(gDAlp_r(m,n) 
    /*3*/  * r(m,n)) + gDAlp(m,n) * (1 + gDA(m,n) * r(m,n) 
    /*3*/  + gDB(m,n) * r(m,n) + 4 * gDconf(m,n) * r(m,n)))) 
    /*0*/  / (3. * pow2(gA(m,n)) * r(m,n)) + k_g * (gAlp(m,n) 
    /*1*/  * ((-4 * (3 * P_1_0(R(m,n)) - 2 * P_2_0(R(m,n)) 
    /*4*/  + b_0)) / 3. + (4 * exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * fA(m,n) * (2 * (-2 + pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) - 3 * (-1 + pow2(Lt(m,n))) 
    /*4*/  * P_2_1(R(m,n)) - b_1)) / (3. * gA(m,n) * Lt(m,n)))
    /*1*/  + fAlp(m,n) * ((-4 * exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * fA(m,n) * (2 * P_1_2(R(m,n)) + b_2))
    /*2*/  / (3. * gA(m,n)) - (4 * (3 * P_1_1(R(m,n)) - 2 
    /*4*/  * P_2_1(R(m,n)) + b_1)) / (3. * Lt(m,n)))) 
    /*0*/  + gAlp(m,n) * ((-8 * exp(-4 * gconf(m,n)) * (2 
    /*3*/  * gDB(m,n) + pow2(gB(m,n)) * (-gL(m,n) + gL_r(m,n) 
    /*4*/  * r(m,n)))) / (3. * pow2(gB(m,n)) * r(m,n)) + (8 
    /*2*/  * exp(-4 * gconf(m,n)) * (-2 * pow2(gDA(m,n)) 
    /*3*/  * r(m,n) - 4 * pow2(gDconf(m,n)) * r(m,n) - 2 
    /*3*/  * gDconf(m,n) * (1 + gDB(m,n) * r(m,n)) + gDA(m,n) 
    /*3*/  * (2 + 3 * gDB(m,n) * r(m,n) - 2 * gDconf(m,n) 
    /*4*/  * r(m,n)) + r(m,n) * (gDA_r(m,n) - gDB_r(m,n) + 2 
    /*4*/  * gDconf_r(m,n) + gsig(m,n) * (3 + 2 * gDB(m,n) 
    /*5*/  * r(m,n))))) / (3. * pow2(gA(m,n)) * r(m,n)) - (2 
    /*2*/  * exp(-2 * gconf(m,n)) * gA1_r(m,n) * p(m,n)) 
    /*1*/  / (gA(m,n) * Lt(m,n)));

Jacobian[15][2]=
	k_g * (gAlp(m,n) * ((4 * R(m,n) * (b_1 + 2 
    /*4*/  * R(m,n) * b_2)) / 3. - (4 * exp(2 * fconf(m,n) - 2
    /*4*/  * gconf(m,n)) * fA(m,n) * (2 * (-2 
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) - 3 * (-1 
    /*5*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)) - b_1)) / (3. 
    /*3*/  * gA(m,n) * Lt(m,n))) + fAlp(m,n) * ((4 * exp(2 
    /*4*/  * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * (2 
    /*4*/  * P_1_2(R(m,n)) + b_2)) / (3. * gA(m,n)) + (4 
    /*3*/  * R(m,n) * (b_2 + 2 * R(m,n) * b_3)) / (3. 
    /*3*/  * Lt(m,n))));

Jacobian[15][3]=
	(8 * exp(-4 * gconf(m,n)) * gDAlp(m,n)) / (3. 
    /*1*/  * pow2(gA(m,n))) + (4 * exp(-4 * gconf(m,n)) 
    /*1*/  * gAlp(m,n) * (1 + gDA(m,n) * r(m,n) + gDB(m,n) 
    /*2*/  * r(m,n) + 4 * gDconf(m,n) * r(m,n))) / (3. 
    /*1*/  * pow2(gA(m,n)) * r(m,n));

Jacobian[15][5]=
	gAlp(m,n) * (gA1(m,n) - gtrA(m,n) / 3.);

Jacobian[15][7]=
	(-4 * exp(-4 * gconf(m,n)) * (-(gDAlp_r(m,n) 
    /*3*/  * r(m,n)) + gDAlp(m,n) * (1 + gDA(m,n) * r(m,n) 
    /*3*/  + gDB(m,n) * r(m,n) + 4 * gDconf(m,n) * r(m,n)))) 
    /*0*/  / (3. * pow3(gA(m,n)) * r(m,n)) + gAlp(m,n) * ((4 
    /*2*/  * exp(-4 * gconf(m,n)) * (-2 * pow2(gDA(m,n)) 
    /*3*/  * r(m,n) - 4 * pow2(gDconf(m,n)) * r(m,n) - 2 
    /*3*/  * gDconf(m,n) * (1 + gDB(m,n) * r(m,n)) + gDA(m,n) 
    /*3*/  * (2 + 3 * gDB(m,n) * r(m,n) - 2 * gDconf(m,n) 
    /*4*/  * r(m,n)) + r(m,n) * (gDA_r(m,n) - gDB_r(m,n) + 2 
    /*4*/  * gDconf_r(m,n) + gsig(m,n) * (3 + 2 * gDB(m,n) 
    /*5*/  * r(m,n))))) / (3. * pow3(gA(m,n)) * r(m,n)) 
    /*1*/  - (exp(-2 * gconf(m,n)) * gA1_r(m,n) * p(m,n)) 
    /*1*/  / (pow2(gA(m,n)) * Lt(m,n))) + k_g * ((-2 * exp(2 
    /*3*/  * fconf(m,n) - 2 * gconf(m,n)) * fAlp(m,n) * fA(m,n)
    /*2*/  * P_1_2(R(m,n))) / (3. * pow2(gA(m,n))) - (2 
    /*2*/  * exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) 
    /*2*/  * gAlp(m,n) * (P_1_1(R(m,n)) + (-1 + pow2(Lt(m,n)))
    /*3*/  * P_2_1(R(m,n)))) / (3. * pow2(gA(m,n)) 
    /*2*/  * Lt(m,n)));

Jacobian[15][8]=
	(-8 * exp(-4 * gconf(m,n)) * gDB(m,n) * gAlp(m,n))
    /*0*/  / (3. * pow3(gB(m,n)) * r(m,n)) + k_g * (gAlp(m,n)
    /*1*/  * ((-2 * R(m,n) * (b_1 + 2 * R(m,n) * b_2)) / (3.
    /*3*/  * gB(m,n)) - (2 * exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * fA(m,n) * R(m,n) * (2 * (-1 
    /*5*/  + pow2(Lt(m,n))) * P_1_2(R(m,n)) - b_2)) / (3. 
    /*3*/  * gA(m,n) * gB(m,n) * Lt(m,n))) + fAlp(m,n) * ((2 
    /*3*/  * exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) 
    /*3*/  * R(m,n) * b_3) / (3. * gA(m,n) * gB(m,n)) - (2 
    /*3*/  * R(m,n) * (b_2 + 2 * R(m,n) * b_3)) / (3. * gB(m,n)
    /*3*/  * Lt(m,n))));

Jacobian[15][9]=
	(2 * exp(-4 * gconf(m,n)) * gDAlp(m,n)) / (3. 
    /*1*/  * pow2(gA(m,n))) + (2 * exp(-4 * gconf(m,n)) 
    /*1*/  * gAlp(m,n) * (-2 + 4 * gDA(m,n) * r(m,n) - 3 
    /*2*/  * gDB(m,n) * r(m,n) + 2 * gDconf(m,n) * r(m,n))) 
    /*0*/  / (3. * pow2(gA(m,n)) * r(m,n));

Jacobian[15][10]=
	(2 * exp(-4 * gconf(m,n)) * gDAlp(m,n)) / (3. 
    /*1*/  * pow2(gA(m,n))) + gAlp(m,n) * ((4 * exp(-4 
    /*3*/  * gconf(m,n))) / (3. * pow2(gB(m,n)) * r(m,n)) - (2
    /*2*/  * exp(-4 * gconf(m,n)) * (3 * gDA(m,n) - 2 
    /*3*/  * gDconf(m,n) + 2 * gsig(m,n) * r(m,n))) / (3. 
    /*2*/  * pow2(gA(m,n))));

Jacobian[15][11]=
	k_g * ((2 * exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*2*/  * fAlp(m,n) * P_1_2(R(m,n))) / (3. * gA(m,n)) + (2 
    /*2*/  * exp(2 * fconf(m,n) - 2 * gconf(m,n)) * gAlp(m,n) 
    /*2*/  * (P_1_1(R(m,n)) + (-1 + pow2(Lt(m,n))) 
    /*3*/  * P_2_1(R(m,n)))) / (3. * gA(m,n) * Lt(m,n)));

Jacobian[15][12]=
	k_g * (gAlp(m,n) * ((-2 * exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * (2 * P_1_1(R(m,n)) + b_1)) / (3. 
    /*3*/  * gB(m,n)) + (2 * exp(4 * fconf(m,n) - 4 
    /*4*/  * gconf(m,n)) * fA(m,n) * (2 * (-1 + pow2(Lt(m,n)))
    /*4*/  * P_1_2(R(m,n)) - b_2)) / (3. * gA(m,n) * gB(m,n)
    /*3*/  * Lt(m,n))) + fAlp(m,n) * ((-2 * exp(4 
    /*4*/  * fconf(m,n) - 4 * gconf(m,n)) * fA(m,n) * b_3) 
    /*2*/  / (3. * gA(m,n) * gB(m,n)) - (2 * exp(2 * fconf(m,n)
    /*4*/  - 2 * gconf(m,n)) * (2 * P_1_2(R(m,n)) + b_2)) 
    /*2*/  / (3. * gB(m,n) * Lt(m,n))));

Jacobian[15][15]=
	gAlp(m,n) * (gtrA(m,n) + gtrK(m,n));

Jacobian[15][19]=
	(-2 * exp(-4 * gconf(m,n)) * gAlp(m,n)) / (3. 
    /*1*/  * r(m,n));

Jacobian[15][21]=
	(-2 * exp(-4 * gconf(m,n)) * gAlp(m,n) * (3 + 2 
    /*2*/  * gDB(m,n) * r(m,n))) / (3. * pow2(gA(m,n)));

Jacobian[16][1]=
	(4 * exp(-4 * gconf(m,n)) * (-(gDAlp_r(m,n) 
    /*3*/  * r(m,n)) + gDAlp(m,n) * (1 + gDA(m,n) * r(m,n) 
    /*3*/  + gDB(m,n) * r(m,n) + 4 * gDconf(m,n) * r(m,n)))) 
    /*0*/  / (3. * pow2(gA(m,n)) * r(m,n)) + k_g * (gAlp(m,n) 
    /*1*/  * ((2 * (3 * P_1_0(R(m,n)) - 2 * P_2_0(R(m,n)) 
    /*4*/  + b_0)) / 3. - (2 * exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * fA(m,n) * (2 * (-2 + pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) - 3 * (-1 + pow2(Lt(m,n))) 
    /*4*/  * P_2_1(R(m,n)) - b_1)) / (3. * gA(m,n) * Lt(m,n)))
    /*1*/  + fAlp(m,n) * ((2 * exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * fA(m,n) * (2 * P_1_2(R(m,n)) + b_2))
    /*2*/  / (3. * gA(m,n)) + (2 * (3 * P_1_1(R(m,n)) - 2 
    /*4*/  * P_2_1(R(m,n)) + b_1)) / (3. * Lt(m,n)))) 
    /*0*/  + gAlp(m,n) * ((exp(-4 * gconf(m,n)) * (8 * gDB(m,n)
    /*3*/  - 4 * pow2(gB(m,n)) * (gL(m,n) - gL_r(m,n) 
    /*4*/  * r(m,n)))) / (3. * pow2(gB(m,n)) * r(m,n)) + (4 
    /*2*/  * exp(-4 * gconf(m,n)) * (2 * pow2(gDA(m,n)) 
    /*3*/  * r(m,n) + 4 * pow2(gDconf(m,n)) * r(m,n) + 2 
    /*3*/  * gDconf(m,n) * (1 + gDB(m,n) * r(m,n)) + gDA(m,n) 
    /*3*/  * (-2 - 3 * gDB(m,n) * r(m,n) + 2 * gDconf(m,n) 
    /*4*/  * r(m,n)) - r(m,n) * (gDA_r(m,n) - gDB_r(m,n) + 2 
    /*4*/  * gDconf_r(m,n) + gsig(m,n) * (3 + 2 * gDB(m,n) 
    /*5*/  * r(m,n))))) / (3. * pow2(gA(m,n)) * r(m,n)) - (2 
    /*2*/  * exp(-2 * gconf(m,n)) * gA2_r(m,n) * p(m,n)) 
    /*1*/  / (gA(m,n) * Lt(m,n)));

Jacobian[16][2]=
	k_g * (gAlp(m,n) * ((-2 * R(m,n) * (b_1 + 2 
    /*4*/  * R(m,n) * b_2)) / 3. + (2 * exp(2 * fconf(m,n) - 2
    /*4*/  * gconf(m,n)) * fA(m,n) * (2 * (-2 
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) - 3 * (-1 
    /*5*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)) - b_1)) / (3. 
    /*3*/  * gA(m,n) * Lt(m,n))) + fAlp(m,n) * ((-2 * exp(2 
    /*4*/  * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * (2 
    /*4*/  * P_1_2(R(m,n)) + b_2)) / (3. * gA(m,n)) - (2 
    /*3*/  * R(m,n) * (b_2 + 2 * R(m,n) * b_3)) / (3. 
    /*3*/  * Lt(m,n))));

Jacobian[16][3]=
	(-4 * exp(-4 * gconf(m,n)) * gDAlp(m,n)) / (3. 
    /*1*/  * pow2(gA(m,n))) - (2 * exp(-4 * gconf(m,n)) 
    /*1*/  * gAlp(m,n) * (1 + gDA(m,n) * r(m,n) + gDB(m,n) 
    /*2*/  * r(m,n) + 4 * gDconf(m,n) * r(m,n))) / (3. 
    /*1*/  * pow2(gA(m,n)) * r(m,n));

Jacobian[16][5]=
	gAlp(m,n) * (gA2(m,n) - gtrA(m,n) / 3.);

Jacobian[16][7]=
	(2 * exp(-4 * gconf(m,n)) * (-(gDAlp_r(m,n) 
    /*3*/  * r(m,n)) + gDAlp(m,n) * (1 + gDA(m,n) * r(m,n) 
    /*3*/  + gDB(m,n) * r(m,n) + 4 * gDconf(m,n) * r(m,n)))) 
    /*0*/  / (3. * pow3(gA(m,n)) * r(m,n)) + gAlp(m,n) 
    /*0*/  * ((exp(-4 * gconf(m,n)) * (4 * pow2(gDA(m,n)) 
    /*3*/  * r(m,n) + 8 * pow2(gDconf(m,n)) * r(m,n) + 4 
    /*3*/  * gDconf(m,n) * (1 + gDB(m,n) * r(m,n)) + gDA(m,n) 
    /*3*/  * (-4 - 6 * gDB(m,n) * r(m,n) + 4 * gDconf(m,n) 
    /*4*/  * r(m,n)) - 2 * r(m,n) * (gDA_r(m,n) - gDB_r(m,n) 
    /*4*/  + 2 * gDconf_r(m,n) + gsig(m,n) * (3 + 2 * gDB(m,n)
    /*5*/  * r(m,n))))) / (3. * pow3(gA(m,n)) * r(m,n)) 
    /*1*/  - (exp(-2 * gconf(m,n)) * gA2_r(m,n) * p(m,n)) 
    /*1*/  / (pow2(gA(m,n)) * Lt(m,n))) + k_g * ((exp(2 
    /*3*/  * fconf(m,n) - 2 * gconf(m,n)) * fAlp(m,n) * fA(m,n)
    /*2*/  * P_1_2(R(m,n))) / (3. * pow2(gA(m,n))) + (exp(2 
    /*3*/  * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * gAlp(m,n)
    /*2*/  * (P_1_1(R(m,n)) + (-1 + pow2(Lt(m,n))) 
    /*3*/  * P_2_1(R(m,n)))) / (3. * pow2(gA(m,n)) * Lt(m,n)));

Jacobian[16][8]=
	(4 * exp(-4 * gconf(m,n)) * gDB(m,n) * gAlp(m,n))
    /*0*/  / (3. * pow3(gB(m,n)) * r(m,n)) + k_g * (gAlp(m,n)
    /*1*/  * ((R(m,n) * (b_1 + 2 * R(m,n) * b_2)) / (3. 
    /*3*/  * gB(m,n)) + (exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*3*/  * fA(m,n) * R(m,n) * (2 * (-1 + pow2(Lt(m,n))) 
    /*4*/  * P_1_2(R(m,n)) - b_2)) / (3. * gA(m,n) * gB(m,n) 
    /*3*/  * Lt(m,n))) + fAlp(m,n) * (-(exp(2 * fconf(m,n) - 2
    /*4*/  * gconf(m,n)) * fA(m,n) * R(m,n) * b_3) / (3. 
    /*3*/  * gA(m,n) * gB(m,n)) + (R(m,n) * (b_2 + 2 * R(m,n) 
    /*4*/  * b_3)) / (3. * gB(m,n) * Lt(m,n))));

Jacobian[16][9]=
	-(exp(-4 * gconf(m,n)) * gDAlp(m,n)) / (3. 
    /*1*/  * pow2(gA(m,n))) + (exp(-4 * gconf(m,n)) * gAlp(m,n)
    /*1*/  * (2 - 4 * gDA(m,n) * r(m,n) + 3 * gDB(m,n) 
    /*2*/  * r(m,n) - 2 * gDconf(m,n) * r(m,n))) / (3. 
    /*1*/  * pow2(gA(m,n)) * r(m,n));

Jacobian[16][10]=
	-(exp(-4 * gconf(m,n)) * gDAlp(m,n)) / (3. 
    /*1*/  * pow2(gA(m,n))) + gAlp(m,n) * ((-2 * exp(-4 
    /*3*/  * gconf(m,n))) / (3. * pow2(gB(m,n)) * r(m,n)) 
    /*1*/  + (exp(-4 * gconf(m,n)) * (3 * gDA(m,n) - 2 
    /*3*/  * gDconf(m,n) + 2 * gsig(m,n) * r(m,n))) / (3. 
    /*2*/  * pow2(gA(m,n))));

Jacobian[16][11]=
	k_g * (-(exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*2*/  * fAlp(m,n) * P_1_2(R(m,n))) / (3. * gA(m,n)) 
    /*1*/  - (exp(2 * fconf(m,n) - 2 * gconf(m,n)) * gAlp(m,n)
    /*2*/  * (P_1_1(R(m,n)) + (-1 + pow2(Lt(m,n))) 
    /*3*/  * P_2_1(R(m,n)))) / (3. * gA(m,n) * Lt(m,n)));

Jacobian[16][12]=
	k_g * (fAlp(m,n) * ((exp(4 * fconf(m,n) - 4 
    /*4*/  * gconf(m,n)) * fA(m,n) * b_3) / (3. * gA(m,n) 
    /*3*/  * gB(m,n)) + (exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*3*/  * (2 * P_1_2(R(m,n)) + b_2)) / (3. * gB(m,n) 
    /*3*/  * Lt(m,n))) + gAlp(m,n) * ((exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * (2 * P_1_1(R(m,n)) + b_1)) / (3. 
    /*3*/  * gB(m,n)) + (exp(4 * fconf(m,n) - 4 * gconf(m,n)) 
    /*3*/  * fA(m,n) * (-2 * (-1 + pow2(Lt(m,n))) 
    /*4*/  * P_1_2(R(m,n)) + b_2)) / (3. * gA(m,n) * gB(m,n) 
    /*3*/  * Lt(m,n))));

Jacobian[16][16]=
	gAlp(m,n) * (gtrA(m,n) + gtrK(m,n));

Jacobian[16][19]=
	(exp(-4 * gconf(m,n)) * gAlp(m,n)) / (3. * r(m,n));

Jacobian[16][21]=
	(exp(-4 * gconf(m,n)) * gAlp(m,n) * (3 + 2 
    /*2*/  * gDB(m,n) * r(m,n))) / (3. * pow2(gA(m,n)));

Jacobian[17][1]=
	k_f * (gAlp(m,n) * ((4 * exp(-2 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (P_2_0(R(m,n)) + b_0)) 
    /*2*/  / (3. * fA(m,n) * pow2(R(m,n))) - (4 
    /*3*/  * (P_1_1(R(m,n)) - b_1)) / (3. * pow2(R(m,n)) 
    /*3*/  * Lt(m,n))) + fAlp(m,n) * ((-4 * (P_1_2(R(m,n)) 
    /*4*/  - b_2)) / (3. * pow2(R(m,n))) + (4 * exp(-2 
    /*4*/  * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) * (2 * (-1
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) + pow2(Lt(m,n)) 
    /*4*/  * P_2_1(R(m,n)) + b_1)) / (3. * fA(m,n) 
    /*3*/  * pow2(R(m,n)) * Lt(m,n))));

Jacobian[17][2]=
	(-8 * exp(-4 * fconf(m,n)) * (-(fDAlp_r(m,n) 
    /*3*/  * r(m,n)) + fDAlp(m,n) * (1 + fDA(m,n) * r(m,n) 
    /*3*/  + fDB(m,n) * r(m,n) + 4 * fDconf(m,n) * r(m,n)))) 
    /*0*/  / (3. * pow2(fA(m,n)) * r(m,n)) + k_f * (gAlp(m,n) 
    /*1*/  * ((-4 * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*3*/  * gA(m,n) * (P_2_0(R(m,n)) + b_0)) / (3. * fA(m,n) 
    /*3*/  * pow2(R(m,n))) + (4 * (P_1_1(R(m,n)) - b_1)) / (3.
    /*3*/  * pow2(R(m,n)) * Lt(m,n))) + fAlp(m,n) * ((4 
    /*3*/  * (P_1_2(R(m,n)) - b_2)) / (3. * pow2(R(m,n))) - (4
    /*3*/  * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n)
    /*3*/  * (2 * (-1 + pow2(Lt(m,n))) * P_1_1(R(m,n)) 
    /*4*/  + pow2(Lt(m,n)) * P_2_1(R(m,n)) + b_1)) / (3. 
    /*3*/  * fA(m,n) * pow2(R(m,n)) * Lt(m,n)))) + fAlp(m,n) 
    /*0*/  * ((8 * (-(exp(-4 * fconf(m,n)) * fL_r(m,n)) 
    /*3*/  + (exp(-4 * fconf(m,n)) * fL(m,n)) / r(m,n) - (2 
    /*4*/  * exp(-4 * gconf(m,n)) * fDB(m,n)) / (pow2(gB(m,n))
    /*4*/  * pow2(R(m,n)) * r(m,n)))) / 3. + (8 * exp(-4 
    /*3*/  * fconf(m,n)) * (-2 * pow2(fDA(m,n)) * r(m,n) - 4 
    /*3*/  * pow2(fDconf(m,n)) * r(m,n) - 2 * fDconf(m,n) * (1
    /*4*/  + fDB(m,n) * r(m,n)) + fDA(m,n) * (2 + 3 
    /*4*/  * fDB(m,n) * r(m,n) - 2 * fDconf(m,n) * r(m,n)) 
    /*3*/  + r(m,n) * (fDA_r(m,n) - fDB_r(m,n) + 2 
    /*4*/  * fDconf_r(m,n) + fsig(m,n) * (3 + 2 * fDB(m,n) 
    /*5*/  * r(m,n))))) / (3. * pow2(fA(m,n)) * r(m,n)) + (2 
    /*2*/  * exp(-2 * fconf(m,n)) * fA1_r(m,n) * p(m,n)) 
    /*1*/  / (fA(m,n) * Lt(m,n)));

Jacobian[17][4]=
	(8 * exp(-4 * fconf(m,n)) * fDAlp(m,n)) / (3. 
    /*1*/  * pow2(fA(m,n))) + (4 * exp(-4 * fconf(m,n)) 
    /*1*/  * fAlp(m,n) * (1 + fDA(m,n) * r(m,n) + fDB(m,n) 
    /*2*/  * r(m,n) + 4 * fDconf(m,n) * r(m,n))) / (3. 
    /*1*/  * pow2(fA(m,n)) * r(m,n));

Jacobian[17][6]=
	fAlp(m,n) * (fA1(m,n) - ftrA(m,n) / 3.);

Jacobian[17][7]=
	k_f * ((-2 * exp(-2 * fconf(m,n) + 2 * gconf(m,n))
    /*2*/  * gAlp(m,n) * (P_1_0(R(m,n)) - P_2_0(R(m,n)))) 
    /*1*/  / (3. * fA(m,n) * pow2(R(m,n))) + (2 * exp(-2 
    /*3*/  * fconf(m,n) + 2 * gconf(m,n)) * fAlp(m,n) 
    /*2*/  * (-P_1_1(R(m,n)) + pow2(Lt(m,n)) * P_2_1(R(m,n))))
    /*1*/  / (3. * fA(m,n) * pow2(R(m,n)) * Lt(m,n)));

Jacobian[17][8]=
	k_f * (gAlp(m,n) * ((2 * exp(-2 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (P_1_0(R(m,n)) + b_0)) 
    /*2*/  / (3. * fA(m,n) * gB(m,n) * pow2(R(m,n))) - (2 
    /*3*/  * (P_1_1(R(m,n)) - b_1)) / (3. * gB(m,n) 
    /*3*/  * pow2(R(m,n)) * Lt(m,n))) + fAlp(m,n) * ((-2 
    /*3*/  * (P_1_2(R(m,n)) - b_2)) / (3. * gB(m,n) 
    /*3*/  * pow2(R(m,n))) + (2 * exp(-2 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * gA(m,n) * ((-1 + 2 * pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) + b_1)) / (3. * fA(m,n) * gB(m,n)
    /*3*/  * pow2(R(m,n)) * Lt(m,n))));

Jacobian[17][11]=
	(-4 * exp(-4 * fconf(m,n)) * (-(fDAlp_r(m,n) 
    /*3*/  * r(m,n)) + fDAlp(m,n) * (1 + fDA(m,n) * r(m,n) 
    /*3*/  + fDB(m,n) * r(m,n) + 4 * fDconf(m,n) * r(m,n)))) 
    /*0*/  / (3. * pow3(fA(m,n)) * r(m,n)) + fAlp(m,n) * ((4 
    /*2*/  * exp(-4 * fconf(m,n)) * (-2 * pow2(fDA(m,n)) 
    /*3*/  * r(m,n) - 4 * pow2(fDconf(m,n)) * r(m,n) - 2 
    /*3*/  * fDconf(m,n) * (1 + fDB(m,n) * r(m,n)) + fDA(m,n) 
    /*3*/  * (2 + 3 * fDB(m,n) * r(m,n) - 2 * fDconf(m,n) 
    /*4*/  * r(m,n)) + r(m,n) * (fDA_r(m,n) - fDB_r(m,n) + 2 
    /*4*/  * fDconf_r(m,n) + fsig(m,n) * (3 + 2 * fDB(m,n) 
    /*5*/  * r(m,n))))) / (3. * pow3(fA(m,n)) * r(m,n)) 
    /*1*/  + (exp(-2 * fconf(m,n)) * fA1_r(m,n) * p(m,n)) 
    /*1*/  / (pow2(fA(m,n)) * Lt(m,n))) + k_f * ((-2 * exp(-2 
    /*3*/  * fconf(m,n) + 2 * gconf(m,n)) * gAlp(m,n) * gA(m,n)
    /*2*/  * P_1_1(R(m,n))) / (3. * pow2(fA(m,n)) * R(m,n)) 
    /*1*/  + (2 * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*2*/  * fAlp(m,n) * gA(m,n) * (P_1_1(R(m,n)) 
    /*3*/  - pow2(Lt(m,n)) * P_2_1(R(m,n)))) / (3. 
    /*2*/  * pow2(fA(m,n)) * pow2(R(m,n)) * Lt(m,n)));

Jacobian[17][12]=
	(-8 * exp(2 * fconf(m,n) - 6 * gconf(m,n)) 
    /*1*/  * fDB(m,n) * fAlp(m,n)) / (3. * pow3(gB(m,n)) 
    /*1*/  * pow3(R(m,n)) * r(m,n)) + k_f * (gAlp(m,n) * ((-2 
    /*3*/  * gA(m,n) * (P_1_0(R(m,n)) + b_0)) / (3. * fA(m,n) 
    /*3*/  * gB(m,n) * pow3(R(m,n))) + (2 * exp(2 * fconf(m,n)
    /*4*/  - 2 * gconf(m,n)) * (P_1_1(R(m,n)) - b_1)) / (3. 
    /*3*/  * gB(m,n) * pow3(R(m,n)) * Lt(m,n))) + fAlp(m,n) 
    /*1*/  * ((2 * exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*3*/  * (P_1_2(R(m,n)) - b_2)) / (3. * gB(m,n) 
    /*3*/  * pow3(R(m,n))) - (2 * gA(m,n) * ((-1 + 2 
    /*5*/  * pow2(Lt(m,n))) * P_1_1(R(m,n)) + b_1)) / (3. 
    /*3*/  * fA(m,n) * gB(m,n) * pow3(R(m,n)) * Lt(m,n))));

Jacobian[17][13]=
	(2 * exp(-4 * fconf(m,n)) * fDAlp(m,n)) / (3. 
    /*1*/  * pow2(fA(m,n))) + (2 * exp(-4 * fconf(m,n)) 
    /*1*/  * fAlp(m,n) * (-2 + 4 * fDA(m,n) * r(m,n) - 3 
    /*2*/  * fDB(m,n) * r(m,n) + 2 * fDconf(m,n) * r(m,n))) 
    /*0*/  / (3. * pow2(fA(m,n)) * r(m,n));

Jacobian[17][14]=
	(2 * exp(-4 * fconf(m,n)) * fDAlp(m,n)) / (3. 
    /*1*/  * pow2(fA(m,n))) + fAlp(m,n) * ((4 * exp(-4 
    /*3*/  * gconf(m,n))) / (3. * pow2(gB(m,n)) * pow2(R(m,n))
    /*2*/  * r(m,n)) - (2 * exp(-4 * fconf(m,n)) * (3 
    /*3*/  * fDA(m,n) - 2 * fDconf(m,n) + 2 * fsig(m,n) 
    /*3*/  * r(m,n))) / (3. * pow2(fA(m,n))));

Jacobian[17][17]=
	fAlp(m,n) * (ftrA(m,n) + ftrK(m,n));

Jacobian[17][20]=
	(-2 * exp(-4 * fconf(m,n)) * fAlp(m,n)) / (3. 
    /*1*/  * r(m,n));

Jacobian[17][22]=
	(-2 * exp(-4 * fconf(m,n)) * fAlp(m,n) * (3 + 2 
    /*2*/  * fDB(m,n) * r(m,n))) / (3. * pow2(fA(m,n)));

Jacobian[18][1]=
	k_f * (gAlp(m,n) * ((-2 * exp(-2 * fconf(m,n) + 2
    /*4*/  * gconf(m,n)) * gA(m,n) * (P_2_0(R(m,n)) + b_0)) 
    /*2*/  / (3. * fA(m,n) * pow2(R(m,n))) + (2 
    /*3*/  * (P_1_1(R(m,n)) - b_1)) / (3. * pow2(R(m,n)) 
    /*3*/  * Lt(m,n))) + fAlp(m,n) * ((2 * (P_1_2(R(m,n)) 
    /*4*/  - b_2)) / (3. * pow2(R(m,n))) - (2 * exp(-2 
    /*4*/  * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) * (2 * (-1
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) + pow2(Lt(m,n)) 
    /*4*/  * P_2_1(R(m,n)) + b_1)) / (3. * fA(m,n) 
    /*3*/  * pow2(R(m,n)) * Lt(m,n))));

Jacobian[18][2]=
	(4 * exp(-4 * fconf(m,n)) * (-(fDAlp_r(m,n) 
    /*3*/  * r(m,n)) + fDAlp(m,n) * (1 + fDA(m,n) * r(m,n) 
    /*3*/  + fDB(m,n) * r(m,n) + 4 * fDconf(m,n) * r(m,n)))) 
    /*0*/  / (3. * pow2(fA(m,n)) * r(m,n)) + k_f * (gAlp(m,n) 
    /*1*/  * ((2 * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*3*/  * gA(m,n) * (P_2_0(R(m,n)) + b_0)) / (3. * fA(m,n) 
    /*3*/  * pow2(R(m,n))) - (2 * (P_1_1(R(m,n)) - b_1)) / (3.
    /*3*/  * pow2(R(m,n)) * Lt(m,n))) + fAlp(m,n) * ((-2 
    /*3*/  * (P_1_2(R(m,n)) - b_2)) / (3. * pow2(R(m,n))) + (2
    /*3*/  * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n)
    /*3*/  * (2 * (-1 + pow2(Lt(m,n))) * P_1_1(R(m,n)) 
    /*4*/  + pow2(Lt(m,n)) * P_2_1(R(m,n)) + b_1)) / (3. 
    /*3*/  * fA(m,n) * pow2(R(m,n)) * Lt(m,n)))) + fAlp(m,n) 
    /*0*/  * ((4 * (exp(-4 * fconf(m,n)) * fL_r(m,n) - (exp(-4
    /*5*/  * fconf(m,n)) * fL(m,n)) / r(m,n) + (2 * exp(-4 
    /*5*/  * gconf(m,n)) * fDB(m,n)) / (pow2(gB(m,n)) 
    /*4*/  * pow2(R(m,n)) * r(m,n)))) / 3. + (4 * exp(-4 
    /*3*/  * fconf(m,n)) * (2 * pow2(fDA(m,n)) * r(m,n) + 4 
    /*3*/  * pow2(fDconf(m,n)) * r(m,n) + 2 * fDconf(m,n) * (1
    /*4*/  + fDB(m,n) * r(m,n)) + fDA(m,n) * (-2 - 3 
    /*4*/  * fDB(m,n) * r(m,n) + 2 * fDconf(m,n) * r(m,n)) 
    /*3*/  - r(m,n) * (fDA_r(m,n) - fDB_r(m,n) + 2 
    /*4*/  * fDconf_r(m,n) + fsig(m,n) * (3 + 2 * fDB(m,n) 
    /*5*/  * r(m,n))))) / (3. * pow2(fA(m,n)) * r(m,n)) + (2 
    /*2*/  * exp(-2 * fconf(m,n)) * fA2_r(m,n) * p(m,n)) 
    /*1*/  / (fA(m,n) * Lt(m,n)));

Jacobian[18][4]=
	(-4 * exp(-4 * fconf(m,n)) * fDAlp(m,n)) / (3. 
    /*1*/  * pow2(fA(m,n))) - (2 * exp(-4 * fconf(m,n)) 
    /*1*/  * fAlp(m,n) * (1 + fDA(m,n) * r(m,n) + fDB(m,n) 
    /*2*/  * r(m,n) + 4 * fDconf(m,n) * r(m,n))) / (3. 
    /*1*/  * pow2(fA(m,n)) * r(m,n));

Jacobian[18][6]=
	fAlp(m,n) * (fA2(m,n) - ftrA(m,n) / 3.);

Jacobian[18][7]=
	k_f * ((exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*2*/  * gAlp(m,n) * (P_1_0(R(m,n)) - P_2_0(R(m,n)))) / (3.
    /*2*/  * fA(m,n) * pow2(R(m,n))) + (exp(-2 * fconf(m,n) 
    /*3*/  + 2 * gconf(m,n)) * fAlp(m,n) * (P_1_1(R(m,n)) 
    /*3*/  - pow2(Lt(m,n)) * P_2_1(R(m,n)))) / (3. * fA(m,n) 
    /*2*/  * pow2(R(m,n)) * Lt(m,n)));

Jacobian[18][8]=
	k_f * (gAlp(m,n) * (-(exp(-2 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (P_1_0(R(m,n)) + b_0)) 
    /*2*/  / (3. * fA(m,n) * gB(m,n) * pow2(R(m,n))) 
    /*2*/  + (P_1_1(R(m,n)) - b_1) / (3. * gB(m,n) 
    /*3*/  * pow2(R(m,n)) * Lt(m,n))) + fAlp(m,n) 
    /*1*/  * ((P_1_2(R(m,n)) - b_2) / (3. * gB(m,n) 
    /*3*/  * pow2(R(m,n))) + (exp(-2 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * gA(m,n) * ((1 - 2 * pow2(Lt(m,n))) 
    /*4*/  * P_1_1(R(m,n)) - b_1)) / (3. * fA(m,n) * gB(m,n) 
    /*3*/  * pow2(R(m,n)) * Lt(m,n))));

Jacobian[18][11]=
	(2 * exp(-4 * fconf(m,n)) * (-(fDAlp_r(m,n) 
    /*3*/  * r(m,n)) + fDAlp(m,n) * (1 + fDA(m,n) * r(m,n) 
    /*3*/  + fDB(m,n) * r(m,n) + 4 * fDconf(m,n) * r(m,n)))) 
    /*0*/  / (3. * pow3(fA(m,n)) * r(m,n)) + fAlp(m,n) 
    /*0*/  * ((exp(-4 * fconf(m,n)) * (4 * pow2(fDA(m,n)) 
    /*3*/  * r(m,n) + 8 * pow2(fDconf(m,n)) * r(m,n) + 4 
    /*3*/  * fDconf(m,n) * (1 + fDB(m,n) * r(m,n)) + fDA(m,n) 
    /*3*/  * (-4 - 6 * fDB(m,n) * r(m,n) + 4 * fDconf(m,n) 
    /*4*/  * r(m,n)) - 2 * r(m,n) * (fDA_r(m,n) - fDB_r(m,n) 
    /*4*/  + 2 * fDconf_r(m,n) + fsig(m,n) * (3 + 2 * fDB(m,n)
    /*5*/  * r(m,n))))) / (3. * pow3(fA(m,n)) * r(m,n)) 
    /*1*/  + (exp(-2 * fconf(m,n)) * fA2_r(m,n) * p(m,n)) 
    /*1*/  / (pow2(fA(m,n)) * Lt(m,n))) + k_f * ((exp(-2 
    /*3*/  * fconf(m,n) + 2 * gconf(m,n)) * gAlp(m,n) * gA(m,n)
    /*2*/  * P_1_1(R(m,n))) / (3. * pow2(fA(m,n)) * R(m,n)) 
    /*1*/  + (exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * fAlp(m,n)
    /*2*/  * gA(m,n) * (-P_1_1(R(m,n)) + pow2(Lt(m,n)) 
    /*3*/  * P_2_1(R(m,n)))) / (3. * pow2(fA(m,n)) 
    /*2*/  * pow2(R(m,n)) * Lt(m,n)));

Jacobian[18][12]=
	(4 * exp(2 * fconf(m,n) - 6 * gconf(m,n)) 
    /*1*/  * fDB(m,n) * fAlp(m,n)) / (3. * pow3(gB(m,n)) 
    /*1*/  * pow3(R(m,n)) * r(m,n)) + k_f * (gAlp(m,n) 
    /*1*/  * ((gA(m,n) * (P_1_0(R(m,n)) + b_0)) / (3. * fA(m,n)
    /*3*/  * gB(m,n) * pow3(R(m,n))) + (exp(2 * fconf(m,n) 
    /*4*/  - 2 * gconf(m,n)) * (-P_1_1(R(m,n)) + b_1)) / (3. 
    /*3*/  * gB(m,n) * pow3(R(m,n)) * Lt(m,n))) + fAlp(m,n) 
    /*1*/  * ((exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*3*/  * (-P_1_2(R(m,n)) + b_2)) / (3. * gB(m,n) 
    /*3*/  * pow3(R(m,n))) + (gA(m,n) * ((-1 + 2 
    /*5*/  * pow2(Lt(m,n))) * P_1_1(R(m,n)) + b_1)) / (3. 
    /*3*/  * fA(m,n) * gB(m,n) * pow3(R(m,n)) * Lt(m,n))));

Jacobian[18][13]=
	-(exp(-4 * fconf(m,n)) * fDAlp(m,n)) / (3. 
    /*1*/  * pow2(fA(m,n))) + (exp(-4 * fconf(m,n)) * fAlp(m,n)
    /*1*/  * (2 - 4 * fDA(m,n) * r(m,n) + 3 * fDB(m,n) 
    /*2*/  * r(m,n) - 2 * fDconf(m,n) * r(m,n))) / (3. 
    /*1*/  * pow2(fA(m,n)) * r(m,n));

Jacobian[18][14]=
	-(exp(-4 * fconf(m,n)) * fDAlp(m,n)) / (3. 
    /*1*/  * pow2(fA(m,n))) + fAlp(m,n) * ((-2 * exp(-4 
    /*3*/  * gconf(m,n))) / (3. * pow2(gB(m,n)) * pow2(R(m,n))
    /*2*/  * r(m,n)) + (exp(-4 * fconf(m,n)) * (3 * fDA(m,n)
    /*3*/  - 2 * fDconf(m,n) + 2 * fsig(m,n) * r(m,n))) / (3.
    /*2*/  * pow2(fA(m,n))));

Jacobian[18][18]=
	fAlp(m,n) * (ftrA(m,n) + ftrK(m,n));

Jacobian[18][20]=
	(exp(-4 * fconf(m,n)) * fAlp(m,n)) / (3. * r(m,n));

Jacobian[18][22]=
	(exp(-4 * fconf(m,n)) * fAlp(m,n) * (3 + 2 
    /*2*/  * fDB(m,n) * r(m,n))) / (3. * pow2(fA(m,n)));

Jacobian[19][1]=
	k_g * gAlp(m,n) * (-8 * exp(4 * gconf(m,n)) 
    /*1*/  * gj(m,n) - (8 * exp(2 * fconf(m,n)) * fA(m,n) 
    /*2*/  * p(m,n) * P_1_2(R(m,n)) * R(m,n)) / pow2(gA(m,n)))
    /*0*/  + (exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) * (4
    /*2*/  - 2 * gL_r(m,n) * pow2(gB(m,n)) * pow2(r(m,n)))) 
    /*0*/  / (gA(m,n) * pow2(gB(m,n)) * pow2(r(m,n)) * Lt(m,n));

Jacobian[19][2]=
	(-4 * k_g * exp(2 * fconf(m,n)) * fA(m,n) 
    /*1*/  * gAlp(m,n) * p(m,n) * (2 * P_1_1(R(m,n)) - 3 
    /*2*/  * P_2_1(R(m,n)))) / pow2(gA(m,n));

Jacobian[19][3]=
	(8 * gAsig(m,n) * gAlp(m,n) * pow2(r(m,n))) 
    /*0*/  / pow2(gA(m,n));

Jacobian[19][7]=
	(-4 * k_g * exp(2 * fconf(m,n)) * fA(m,n) 
    /*1*/  * gAlp(m,n) * p(m,n) * P_2_1(R(m,n))) 
    /*0*/  / pow3(gA(m,n)) + (gdet_pff_r(m,n) / gdet(m,n) + 12
    /*1*/  * gA1(m,n) * gDAlp(m,n) - 6 * gBet_rr(m,n) 
    /*1*/  - (gdet_pff(m,n) * gdet_r(m,n)) / pow2(gdet(m,n)) 
    /*1*/  - 4 * gDAlp(m,n) * gtrA(m,n)) / (3. * pow3(gA(m,n)))
    /*0*/  + gAlp(m,n) * ((-8 * (gAsig(m,n) * pow2(r(m,n)) 
    /*3*/  * (gDA(m,n) + gDB(m,n) + 6 * gDconf(m,n) + gsig(m,n)
    /*4*/  * r(m,n)) - gtrA_r(m,n) - gtrK_r(m,n))) / (3. 
    /*2*/  * pow3(gA(m,n))) + (exp(-2 * gconf(m,n)) * p(m,n) 
    /*2*/  * (2 - gL_r(m,n) * pow2(gB(m,n)) * pow2(r(m,n)))) 
    /*1*/  / (pow2(gA(m,n)) * pow2(gB(m,n)) * pow2(r(m,n)) 
    /*2*/  * Lt(m,n)));

Jacobian[19][8]=
	(4 * (gBet(m,n) - gBet_r(m,n) * r(m,n))) 
    /*0*/  / (pow2(r(m,n)) * pow3(gB(m,n))) - (4 * k_g * exp(2
    /*2*/  * fconf(m,n)) * fA(m,n) * gAlp(m,n) * p(m,n) 
    /*1*/  * P_1_2(R(m,n)) * R(m,n)) / (gB(m,n) 
    /*1*/  * pow2(gA(m,n))) + gAlp(m,n) * p(m,n) * ((-4 
    /*2*/  * exp(-2 * gconf(m,n))) / (gA(m,n) * Lt(m,n) 
    /*2*/  * pow2(r(m,n)) * pow3(gB(m,n))) + (4 * exp(-2 
    /*3*/  * gconf(m,n))) / (gA(m,n) * pow2(r(m,n)) 
    /*2*/  * pow3(gB(m,n)) * Lt(m,n)));

Jacobian[19][9]=
	(4 * gAsig(m,n) * gAlp(m,n) * pow2(r(m,n))) / (3.
    /*1*/  * pow2(gA(m,n)));

Jacobian[19][10]=
	(4 * gAsig(m,n) * gAlp(m,n) * pow2(r(m,n))) / (3.
    /*1*/  * pow2(gA(m,n)));

Jacobian[19][11]=
	(2 * k_g * exp(2 * fconf(m,n)) * gAlp(m,n) 
    /*1*/  * p(m,n) * P_2_1(R(m,n))) / pow2(gA(m,n));

Jacobian[19][12]=
	(4 * k_g * exp(4 * fconf(m,n) - 2 * gconf(m,n)) 
    /*1*/  * fA(m,n) * gAlp(m,n) * p(m,n) * P_1_2(R(m,n))) 
    /*0*/  / (gB(m,n) * pow2(gA(m,n)));

Jacobian[19][15]=
	(-2 * gDAlp(m,n)) / pow2(gA(m,n));

Jacobian[19][19]=
	-gdet_pff(m,n) / (3. * gdet(m,n)) - gBet_r(m,n);

Jacobian[19][21]=
	(4 * gAsig(m,n) * gAlp(m,n) * pow3(r(m,n))) / (3.
    /*1*/  * pow2(gA(m,n)));

Jacobian[19][23]=
	(4 * gAlp(m,n) * pow2(r(m,n)) * (gDA(m,n) 
    /*2*/  + gDB(m,n) + 6 * gDconf(m,n) + gsig(m,n) * r(m,n)))
    /*0*/  / (3. * pow2(gA(m,n)));

Jacobian[20][1]=
	(-4 * k_f * exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*1*/  * gA(m,n) * p(m,n) * (2 * P_1_1(R(m,n)) 
    /*2*/  + P_2_1(R(m,n)))) / (pow2(fA(m,n)) * pow2(R(m,n)));

Jacobian[20][2]=
	k_f * fAlp(m,n) * (-8 * exp(4 * fconf(m,n)) 
    /*1*/  * fj(m,n) + (8 * exp(2 * gconf(m,n)) * gA(m,n) 
    /*2*/  * p(m,n) * P_1_1(R(m,n))) / (pow2(fA(m,n)) 
    /*2*/  * pow2(R(m,n)))) + (2 * exp(-2 * fconf(m,n)) 
    /*1*/  * fAlp(m,n) * p(m,n) * (fL_r(m,n) - (2 * exp(4 
    /*4*/  * fconf(m,n) - 4 * gconf(m,n))) / (pow2(gB(m,n)) 
    /*3*/  * pow2(r(m,n)) * pow2(R(m,n))))) / (fA(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[20][4]=
	(8 * fAsig(m,n) * fAlp(m,n) * pow2(r(m,n))) 
    /*0*/  / pow2(fA(m,n));

Jacobian[20][7]=
	(-2 * k_f * exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*1*/  * p(m,n) * P_2_1(R(m,n))) / (pow2(fA(m,n)) 
    /*1*/  * pow2(R(m,n)));

Jacobian[20][8]=
	(-4 * k_f * exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*1*/  * gA(m,n) * p(m,n) * P_1_1(R(m,n))) / (gB(m,n) 
    /*1*/  * pow2(fA(m,n)) * pow2(R(m,n)));

Jacobian[20][11]=
	(4 * k_f * exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*1*/  * gA(m,n) * p(m,n) * P_2_1(R(m,n))) / (pow2(R(m,n))
    /*1*/  * pow3(fA(m,n))) + (fdet_pff_r(m,n) / fdet(m,n) 
    /*1*/  + 12 * fA1(m,n) * fDAlp(m,n) - 6 * fBet_rr(m,n) 
    /*1*/  - (fdet_pff(m,n) * fdet_r(m,n)) / pow2(fdet(m,n)) 
    /*1*/  - 4 * fDAlp(m,n) * ftrA(m,n)) / (3. * pow3(fA(m,n)))
    /*0*/  + fAlp(m,n) * ((-8 * (fAsig(m,n) * pow2(r(m,n)) 
    /*3*/  * (fDA(m,n) + fDB(m,n) + 6 * fDconf(m,n) + fsig(m,n)
    /*4*/  * r(m,n)) - ftrA_r(m,n) - ftrK_r(m,n))) / (3. 
    /*2*/  * pow3(fA(m,n))) + (exp(-2 * fconf(m,n)) * p(m,n) 
    /*2*/  * (fL_r(m,n) - (2 * exp(4 * fconf(m,n) - 4 
    /*5*/  * gconf(m,n))) / (pow2(gB(m,n)) * pow2(r(m,n)) 
    /*4*/  * pow2(R(m,n))))) / (pow2(fA(m,n)) * Lt(m,n)));

Jacobian[20][12]=
	(-4 * exp(6 * fconf(m,n) - 8 * gconf(m,n)) 
    /*1*/  * gAlp(m,n) * p(m,n)) / (gA(m,n) * Lt(m,n) 
    /*1*/  * pow2(r(m,n)) * pow3(gB(m,n)) * pow3(R(m,n))) + (4
    /*1*/  * k_f * exp(2 * fconf(m,n)) * fAlp(m,n) * gA(m,n)
    /*1*/  * p(m,n) * P_1_1(R(m,n))) / (gB(m,n) 
    /*1*/  * pow2(fA(m,n)) * pow3(R(m,n))) + (4 * exp(6 
    /*2*/  * fconf(m,n) - 6 * gconf(m,n)) * (gBet(m,n) 
    /*2*/  - fBet_r(m,n) * r(m,n))) / (pow2(r(m,n)) 
    /*1*/  * pow3(gB(m,n)) * pow3(R(m,n))) - (4 * exp(4 
    /*2*/  * fconf(m,n) - 6 * gconf(m,n)) * fAlp(m,n) * p(m,n))
    /*0*/  / (fA(m,n) * pow2(r(m,n)) * pow3(gB(m,n)) 
    /*1*/  * pow3(R(m,n)) * Lt(m,n));

Jacobian[20][13]=
	(4 * fAsig(m,n) * fAlp(m,n) * pow2(r(m,n))) / (3.
    /*1*/  * pow2(fA(m,n)));

Jacobian[20][14]=
	(4 * fAsig(m,n) * fAlp(m,n) * pow2(r(m,n))) / (3.
    /*1*/  * pow2(fA(m,n)));

Jacobian[20][17]=
	(-2 * fDAlp(m,n)) / pow2(fA(m,n));

Jacobian[20][20]=
	-fdet_pff(m,n) / (3. * fdet(m,n)) - fBet_r(m,n);

Jacobian[20][22]=
	(4 * fAsig(m,n) * fAlp(m,n) * pow3(r(m,n))) / (3.
    /*1*/  * pow2(fA(m,n)));

Jacobian[20][24]=
	(4 * fAlp(m,n) * pow2(r(m,n)) * (fDA(m,n) 
    /*2*/  + fDB(m,n) + 6 * fDconf(m,n) + fsig(m,n) * r(m,n)))
    /*0*/  / (3. * pow2(fA(m,n)));

Jacobian[21][1]=
	gAlp(m,n) * p(m,n) * ((-4 * exp(-2 * gconf(m,n)) 
    /*2*/  * gA(m,n)) / (pow2(gB(m,n)) * pow3(r(m,n)) 
    /*2*/  * Lt(m,n)) - (2 * exp(-2 * gconf(m,n)) * (2 
    /*3*/  * gsig(m,n) + gsig_r(m,n) * r(m,n))) / (gA(m,n) 
    /*2*/  * r(m,n) * Lt(m,n)));

Jacobian[21][7]=
	(4 * gA(m,n) * (gBet(m,n) - gBet_r(m,n) * r(m,n)))
    /*0*/  / (pow2(gB(m,n)) * pow3(r(m,n))) + gAlp(m,n) * ((4
    /*2*/  * gAsig(m,n) * gA(m,n)) / pow2(gB(m,n)) + p(m,n) 
    /*1*/  * ((-4 * exp(-2 * gconf(m,n))) / (Lt(m,n) 
    /*3*/  * pow2(gB(m,n)) * pow3(r(m,n))) + (2 * exp(-2 
    /*4*/  * gconf(m,n))) / (pow2(gB(m,n)) * pow3(r(m,n)) 
    /*3*/  * Lt(m,n)) - (exp(-2 * gconf(m,n)) * (2 * gsig(m,n)
    /*4*/  + gsig_r(m,n) * r(m,n))) / (pow2(gA(m,n)) * r(m,n)
    /*3*/  * Lt(m,n))));

Jacobian[21][8]=
	(4 * pow2(gA(m,n)) * (-gBet(m,n) + gBet_r(m,n) 
    /*2*/  * r(m,n))) / (pow3(gB(m,n)) * pow3(r(m,n))) 
    /*0*/  + gAlp(m,n) * ((-4 * gAsig(m,n) * pow2(gA(m,n))) 
    /*1*/  / pow3(gB(m,n)) + p(m,n) * ((4 * exp(-2 
    /*4*/  * gconf(m,n)) * gA(m,n)) / (Lt(m,n) * pow3(gB(m,n))
    /*3*/  * pow3(r(m,n))) - (4 * exp(-2 * gconf(m,n)) 
    /*3*/  * gA(m,n)) / (pow3(gB(m,n)) * pow3(r(m,n)) 
    /*3*/  * Lt(m,n))));

Jacobian[21][21]=
	2 * gBetr(m,n) + gAlp(m,n) * p(m,n) * ((-2 
    /*2*/  * exp(-2 * gconf(m,n))) / (gA(m,n) * Lt(m,n) 
    /*2*/  * r(m,n)) + (2 * exp(-2 * gconf(m,n))) / (gA(m,n) 
    /*2*/  * r(m,n) * Lt(m,n)));

Jacobian[21][23]=
	(2 * gAlp(m,n) * pow2(gA(m,n))) / pow2(gB(m,n));

Jacobian[22][2]=
	fAlp(m,n) * p(m,n) * ((4 * exp(2 * fconf(m,n) - 4
    /*3*/  * gconf(m,n)) * fA(m,n)) / (pow2(gB(m,n)) 
    /*2*/  * pow2(R(m,n)) * pow3(r(m,n)) * Lt(m,n)) + (2 
    /*2*/  * exp(-2 * fconf(m,n)) * (2 * fsig(m,n) 
    /*3*/  + fsig_r(m,n) * r(m,n))) / (fA(m,n) * r(m,n) 
    /*2*/  * Lt(m,n)));

Jacobian[22][11]=
	(-4 * exp(4 * fconf(m,n) - 6 * gconf(m,n)) 
    /*1*/  * fA(m,n) * gAlp(m,n) * p(m,n)) / (gA(m,n) * Lt(m,n)
    /*1*/  * pow2(gB(m,n)) * pow2(R(m,n)) * pow3(r(m,n))) 
    /*0*/  + (4 * exp(4 * fconf(m,n) - 4 * gconf(m,n)) 
    /*1*/  * fA(m,n) * (gBet(m,n) - fBet_r(m,n) * r(m,n))) 
    /*0*/  / (pow2(gB(m,n)) * pow2(R(m,n)) * pow3(r(m,n))) 
    /*0*/  + fAlp(m,n) * ((4 * exp(4 * fconf(m,n) - 4 
    /*3*/  * gconf(m,n)) * fAsig(m,n) * fA(m,n)) 
    /*1*/  / (pow2(gB(m,n)) * pow2(R(m,n))) + p(m,n) * ((-2 
    /*3*/  * exp(2 * fconf(m,n) - 4 * gconf(m,n))) 
    /*2*/  / (pow2(gB(m,n)) * pow2(R(m,n)) * pow3(r(m,n)) 
    /*3*/  * Lt(m,n)) + (exp(-2 * fconf(m,n)) * (2 * fsig(m,n)
    /*4*/  + fsig_r(m,n) * r(m,n))) / (pow2(fA(m,n)) * r(m,n)
    /*3*/  * Lt(m,n))));

Jacobian[22][12]=
	(4 * exp(6 * fconf(m,n) - 8 * gconf(m,n)) 
    /*1*/  * gAlp(m,n) * p(m,n) * pow2(fA(m,n))) / (gA(m,n) 
    /*1*/  * Lt(m,n) * pow3(gB(m,n)) * pow3(r(m,n)) 
    /*1*/  * pow3(R(m,n))) + (4 * exp(6 * fconf(m,n) - 6 
    /*2*/  * gconf(m,n)) * pow2(fA(m,n)) * (-gBet(m,n) 
    /*2*/  + fBet_r(m,n) * r(m,n))) / (pow3(gB(m,n)) 
    /*1*/  * pow3(r(m,n)) * pow3(R(m,n))) + fAlp(m,n) * ((-4 
    /*2*/  * exp(6 * fconf(m,n) - 6 * gconf(m,n)) * fAsig(m,n)
    /*2*/  * pow2(fA(m,n))) / (pow3(gB(m,n)) * pow3(R(m,n)))
    /*1*/  + (4 * exp(4 * fconf(m,n) - 6 * gconf(m,n)) 
    /*2*/  * fA(m,n) * p(m,n)) / (pow3(gB(m,n)) * pow3(r(m,n))
    /*2*/  * pow3(R(m,n)) * Lt(m,n)));

Jacobian[22][22]=
	2 * gBetr(m,n) - (2 * exp(-2 * gconf(m,n)) 
    /*1*/  * gAlp(m,n) * p(m,n)) / (gA(m,n) * Lt(m,n) * r(m,n))
    /*0*/  - (2 * exp(-2 * fconf(m,n)) * fAlp(m,n) * p(m,n))
    /*0*/  / (fA(m,n) * r(m,n) * Lt(m,n));

Jacobian[22][24]=
	(2 * exp(4 * fconf(m,n) - 4 * gconf(m,n)) 
    /*1*/  * fAlp(m,n) * pow2(fA(m,n))) / (pow2(gB(m,n)) 
    /*1*/  * pow2(R(m,n)));

Jacobian[23][1]=
	(-4 * exp(-4 * gconf(m,n)) * (-(gDAlp_r(m,n) 
    /*3*/  * r(m,n)) + gDAlp(m,n) * (1 + gDA(m,n) * r(m,n) 
    /*3*/  + gDB(m,n) * r(m,n) + 4 * gDconf(m,n) * r(m,n)))) 
    /*0*/  / (pow2(gA(m,n)) * pow3(r(m,n))) + k_g * (gAlp(m,n)
    /*1*/  * ((-2 * (3 * P_1_0(R(m,n)) - 2 * P_2_0(R(m,n)) 
    /*4*/  + b_0)) / pow2(r(m,n)) + (2 * exp(2 * fconf(m,n) - 2
    /*4*/  * gconf(m,n)) * fA(m,n) * (2 * (-2 
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) - 3 * (-1 
    /*5*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)) - b_1)) / (gA(m,n)
    /*3*/  * pow2(r(m,n)) * Lt(m,n))) + fAlp(m,n) * ((-2 
    /*3*/  * exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) 
    /*3*/  * (2 * P_1_2(R(m,n)) + b_2)) / (gA(m,n) 
    /*3*/  * pow2(r(m,n))) - (2 * (3 * P_1_1(R(m,n)) - 2 
    /*4*/  * P_2_1(R(m,n)) + b_1)) / (pow2(r(m,n)) * Lt(m,n))))
    /*0*/  + gAlp(m,n) * ((4 * exp(-4 * gconf(m,n)) 
    /*2*/  * (gL(m,n) - gL_r(m,n) * r(m,n))) / pow3(r(m,n)) 
    /*1*/  - (2 * exp(-4 * gconf(m,n)) * pow2(gB(m,n)) * (4 
    /*3*/  * gsig_r(m,n) + gsig_rr(m,n) * r(m,n) + 2 
    /*3*/  * pow2(gsig(m,n)) * r(m,n))) / (Power(gA(m,n),4) 
    /*2*/  * r(m,n)) + (2 * exp(-4 * gconf(m,n)) * (-6 
    /*3*/  * pow2(gDA(m,n)) * r(m,n) - 8 * pow2(gDconf(m,n)) 
    /*3*/  * r(m,n) - 4 * gDconf(m,n) * (1 + gDB(m,n) * r(m,n))
    /*3*/  + gDA(m,n) * (8 * gDB(m,n) * r(m,n) - 4 
    /*4*/  * gDconf(m,n) * r(m,n)) + r(m,n) * (4 
    /*4*/  * gDconf_r(m,n) + gLD(gsig(m,n)) * pow2(gB(m,n)) 
    /*4*/  * pow2(r(m,n)) + 8 * gDB(m,n) * gsig(m,n) * r(m,n) 
    /*4*/  - 2 * gsig_r(m,n) * r(m,n) + 2 * gL(m,n) * gsig(m,n)
    /*4*/  * pow2(gB(m,n)) * r(m,n)))) / (pow2(gA(m,n)) 
    /*2*/  * pow3(r(m,n))) - (2 * exp(-2 * gconf(m,n)) * p(m,n)
    /*2*/  * (2 * gAsig(m,n) + gAsig_r(m,n) * r(m,n))) 
    /*1*/  / (gA(m,n) * r(m,n) * Lt(m,n)));

Jacobian[23][2]=
	k_g * (gAlp(m,n) * ((2 * R(m,n) * (b_1 + 2 
    /*4*/  * R(m,n) * b_2)) / pow2(r(m,n)) - (2 * exp(2 
    /*4*/  * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * (2 * (-2
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) - 3 * (-1 
    /*5*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)) - b_1)) / (gA(m,n)
    /*3*/  * pow2(r(m,n)) * Lt(m,n))) + fAlp(m,n) * ((2 
    /*3*/  * exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) 
    /*3*/  * (2 * P_1_2(R(m,n)) + b_2)) / (gA(m,n) 
    /*3*/  * pow2(r(m,n))) + (2 * R(m,n) * (b_2 + 2 * R(m,n) 
    /*4*/  * b_3)) / (pow2(r(m,n)) * Lt(m,n))));

Jacobian[23][3]=
	(4 * exp(-4 * gconf(m,n)) * gDAlp(m,n)) 
    /*0*/  / (pow2(gA(m,n)) * pow2(r(m,n))) + (2 * exp(-4 
    /*2*/  * gconf(m,n)) * gAlp(m,n) * (1 + gDA(m,n) * r(m,n) 
    /*2*/  + gDB(m,n) * r(m,n) + 4 * gDconf(m,n) * r(m,n))) 
    /*0*/  / (pow2(gA(m,n)) * pow3(r(m,n)));

Jacobian[23][5]=
	gAsig(m,n) * gAlp(m,n);

Jacobian[23][7]=
	(-2 * exp(-4 * gconf(m,n)) * (-(gDAlp_r(m,n) 
    /*3*/  * r(m,n)) + gDAlp(m,n) * (1 + gDA(m,n) * r(m,n) 
    /*3*/  + gDB(m,n) * r(m,n) + 4 * gDconf(m,n) * r(m,n)))) 
    /*0*/  / (pow3(gA(m,n)) * pow3(r(m,n))) + k_g * (-((exp(2 
    /*4*/  * fconf(m,n) - 2 * gconf(m,n)) * fAlp(m,n) * fA(m,n)
    /*3*/  * P_1_2(R(m,n))) / (pow2(gA(m,n)) * pow2(r(m,n))))
    /*1*/  - (exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n)
    /*2*/  * gAlp(m,n) * (P_1_1(R(m,n)) + (-1 
    /*4*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)))) / (pow2(gA(m,n))
    /*2*/  * pow2(r(m,n)) * Lt(m,n))) + gAlp(m,n) * ((-2 
    /*2*/  * exp(-4 * gconf(m,n)) * pow2(gB(m,n)) * (4 
    /*3*/  * gsig_r(m,n) + gsig_rr(m,n) * r(m,n) + 2 
    /*3*/  * pow2(gsig(m,n)) * r(m,n))) / (Power(gA(m,n),5) 
    /*2*/  * r(m,n)) + (exp(-4 * gconf(m,n)) * (-6 
    /*3*/  * pow2(gDA(m,n)) * r(m,n) - 8 * pow2(gDconf(m,n)) 
    /*3*/  * r(m,n) - 4 * gDconf(m,n) * (1 + gDB(m,n) * r(m,n))
    /*3*/  + gDA(m,n) * (8 * gDB(m,n) * r(m,n) - 4 
    /*4*/  * gDconf(m,n) * r(m,n)) + r(m,n) * (4 
    /*4*/  * gDconf_r(m,n) + gLD(gsig(m,n)) * pow2(gB(m,n)) 
    /*4*/  * pow2(r(m,n)) + 8 * gDB(m,n) * gsig(m,n) * r(m,n) 
    /*4*/  - 2 * gsig_r(m,n) * r(m,n) + 2 * gL(m,n) * gsig(m,n)
    /*4*/  * pow2(gB(m,n)) * r(m,n)))) / (pow3(gA(m,n)) 
    /*2*/  * pow3(r(m,n))) - (exp(-2 * gconf(m,n)) * p(m,n) 
    /*2*/  * (2 * gAsig(m,n) + gAsig_r(m,n) * r(m,n))) 
    /*1*/  / (pow2(gA(m,n)) * r(m,n) * Lt(m,n)));

Jacobian[23][8]=
	gAlp(m,n) * (-((exp(-4 * gconf(m,n)) * gB(m,n) 
    /*3*/  * (2 * gL(m,n) * gsig(m,n) + gLD(gsig(m,n)) 
    /*4*/  * r(m,n))) / (pow2(gA(m,n)) * r(m,n))) + (exp(-4 
    /*3*/  * gconf(m,n)) * gB(m,n) * (4 * gsig_r(m,n) 
    /*3*/  + gsig_rr(m,n) * r(m,n) + 2 * pow2(gsig(m,n)) 
    /*3*/  * r(m,n))) / (Power(gA(m,n),4) * r(m,n))) + k_g 
    /*0*/  * (gAlp(m,n) * (-((R(m,n) * (b_1 + 2 * R(m,n) 
    /*5*/  * b_2)) / (gB(m,n) * pow2(r(m,n)))) - (exp(2 
    /*4*/  * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * R(m,n) 
    /*3*/  * (2 * (-1 + pow2(Lt(m,n))) * P_1_2(R(m,n)) - b_2))
    /*2*/  / (gA(m,n) * gB(m,n) * pow2(r(m,n)) * Lt(m,n))) 
    /*1*/  + fAlp(m,n) * ((exp(2 * fconf(m,n) - 2 * gconf(m,n))
    /*3*/  * fA(m,n) * R(m,n) * b_3) / (gA(m,n) * gB(m,n) 
    /*3*/  * pow2(r(m,n))) - (R(m,n) * (b_2 + 2 * R(m,n) 
    /*4*/  * b_3)) / (gB(m,n) * pow2(r(m,n)) * Lt(m,n))));

Jacobian[23][9]=
	(exp(-4 * gconf(m,n)) * gDAlp(m,n)) 
    /*0*/  / (pow2(gA(m,n)) * pow2(r(m,n))) + (2 * exp(-4 
    /*2*/  * gconf(m,n)) * (3 * gDA(m,n) - 2 * gDB(m,n) 
    /*2*/  + gDconf(m,n)) * gAlp(m,n)) / (pow2(gA(m,n)) 
    /*1*/  * pow2(r(m,n)));

Jacobian[23][10]=
	(exp(-4 * gconf(m,n)) * gDAlp(m,n)) 
    /*0*/  / (pow2(gA(m,n)) * pow2(r(m,n))) + (exp(-4 
    /*2*/  * gconf(m,n)) * gAlp(m,n) * (-4 * gDA(m,n) + 2 
    /*2*/  * gDconf(m,n) - 4 * gsig(m,n) * r(m,n))) 
    /*0*/  / (pow2(gA(m,n)) * pow2(r(m,n)));

Jacobian[23][11]=
	k_g * ((exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*2*/  * fAlp(m,n) * P_1_2(R(m,n))) / (gA(m,n) 
    /*2*/  * pow2(r(m,n))) + (exp(2 * fconf(m,n) - 2 
    /*3*/  * gconf(m,n)) * gAlp(m,n) * (P_1_1(R(m,n)) + (-1 
    /*4*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)))) / (gA(m,n) 
    /*2*/  * pow2(r(m,n)) * Lt(m,n)));

Jacobian[23][12]=
	k_g * (gAlp(m,n) * (-((exp(2 * fconf(m,n) - 2 
    /*5*/  * gconf(m,n)) * (2 * P_1_1(R(m,n)) + b_1)) 
    /*3*/  / (gB(m,n) * pow2(r(m,n)))) + (exp(4 * fconf(m,n) 
    /*4*/  - 4 * gconf(m,n)) * fA(m,n) * (2 * (-1 
    /*5*/  + pow2(Lt(m,n))) * P_1_2(R(m,n)) - b_2)) / (gA(m,n)
    /*3*/  * gB(m,n) * pow2(r(m,n)) * Lt(m,n))) + fAlp(m,n) 
    /*1*/  * (-((exp(4 * fconf(m,n) - 4 * gconf(m,n)) * fA(m,n)
    /*4*/  * b_3) / (gA(m,n) * gB(m,n) * pow2(r(m,n)))) 
    /*2*/  - (exp(2 * fconf(m,n) - 2 * gconf(m,n)) * (2 
    /*4*/  * P_1_2(R(m,n)) + b_2)) / (gB(m,n) * pow2(r(m,n)) 
    /*3*/  * Lt(m,n))));

Jacobian[23][19]=
	gAlp(m,n) * (-(exp(-4 * gconf(m,n)) 
    /*2*/  / pow3(r(m,n))) - (exp(-4 * gconf(m,n)) 
    /*2*/  * pow2(gB(m,n)) * (2 * gsig(m,n) + gsig_r(m,n) 
    /*3*/  * r(m,n))) / (2. * pow2(gA(m,n)) * r(m,n)));

Jacobian[23][21]=
	gAlp(m,n) * ((2 * exp(-4 * gconf(m,n)) * gsig(m,n)
    /*2*/  * pow2(gB(m,n))) / Power(gA(m,n),4) - (exp(-4 
    /*3*/  * gconf(m,n)) * (4 * gDB(m,n) + gL(m,n) 
    /*3*/  * pow2(gB(m,n)))) / (pow2(gA(m,n)) * r(m,n)));

Jacobian[23][23]=
	2 * gBetr(m,n) + gAlp(m,n) * (gtrA(m,n) 
    /*1*/  + gtrK(m,n) + p(m,n) * ((-2 * exp(-2 * gconf(m,n)))
    /*2*/  / (gA(m,n) * Lt(m,n) * r(m,n)) + (2 * exp(-2 
    /*4*/  * gconf(m,n))) / (gA(m,n) * r(m,n) * Lt(m,n))));

Jacobian[24][1]=
	k_f * (gAlp(m,n) * ((2 * exp(-2 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (P_2_0(R(m,n)) + b_0)) 
    /*2*/  / (fA(m,n) * pow2(r(m,n)) * pow2(R(m,n))) - (2 
    /*3*/  * (P_1_1(R(m,n)) - b_1)) / (pow2(r(m,n)) 
    /*3*/  * pow2(R(m,n)) * Lt(m,n))) + fAlp(m,n) * ((-2 
    /*3*/  * (P_1_2(R(m,n)) - b_2)) / (pow2(r(m,n)) 
    /*3*/  * pow2(R(m,n))) + (2 * exp(-2 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (2 * (-1 + pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) + pow2(Lt(m,n)) * P_2_1(R(m,n)) 
    /*4*/  + b_1)) / (fA(m,n) * pow2(r(m,n)) * pow2(R(m,n)) 
    /*3*/  * Lt(m,n))));

Jacobian[24][2]=
	(-4 * exp(-4 * fconf(m,n)) * (-(fDAlp_r(m,n) 
    /*3*/  * r(m,n)) + fDAlp(m,n) * (1 + fDA(m,n) * r(m,n) 
    /*3*/  + fDB(m,n) * r(m,n) + 4 * fDconf(m,n) * r(m,n)))) 
    /*0*/  / (pow2(fA(m,n)) * pow3(r(m,n))) + k_f * (gAlp(m,n)
    /*1*/  * ((-2 * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*3*/  * gA(m,n) * (P_2_0(R(m,n)) + b_0)) / (fA(m,n) 
    /*3*/  * pow2(r(m,n)) * pow2(R(m,n))) + (2 * (P_1_1(R(m,n))
    /*4*/  - b_1)) / (pow2(r(m,n)) * pow2(R(m,n)) * Lt(m,n)))
    /*1*/  + fAlp(m,n) * ((2 * (P_1_2(R(m,n)) - b_2)) 
    /*2*/  / (pow2(r(m,n)) * pow2(R(m,n))) - (2 * exp(-2 
    /*4*/  * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) * (2 * (-1
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) + pow2(Lt(m,n)) 
    /*4*/  * P_2_1(R(m,n)) + b_1)) / (fA(m,n) * pow2(r(m,n)) 
    /*3*/  * pow2(R(m,n)) * Lt(m,n)))) + fAlp(m,n) * ((4 
    /*2*/  * exp(-4 * fconf(m,n)) * (fL(m,n) - fL_r(m,n) 
    /*3*/  * r(m,n))) / pow3(r(m,n)) - (2 * exp(-8 * fconf(m,n)
    /*3*/  + 4 * gconf(m,n)) * pow2(gB(m,n)) * pow2(R(m,n)) 
    /*2*/  * (4 * fsig_r(m,n) + fsig_rr(m,n) * r(m,n) + 2 
    /*3*/  * pow2(fsig(m,n)) * r(m,n))) / (Power(fA(m,n),4) 
    /*2*/  * r(m,n)) + (2 * exp(-8 * fconf(m,n)) * (4 * exp(4 
    /*4*/  * fconf(m,n)) * fDA(m,n) * (2 * fDB(m,n) 
    /*4*/  - fDconf(m,n)) * r(m,n) - 6 * exp(4 * fconf(m,n)) 
    /*3*/  * pow2(fDA(m,n)) * r(m,n) - 8 * exp(4 * fconf(m,n))
    /*3*/  * pow2(fDconf(m,n)) * r(m,n) - 4 * exp(4 
    /*4*/  * fconf(m,n)) * fDconf(m,n) * (1 + fDB(m,n) 
    /*4*/  * r(m,n)) + r(m,n) * (4 * exp(4 * fconf(m,n)) 
    /*4*/  * fDconf_r(m,n) + exp(4 * gconf(m,n)) 
    /*4*/  * fLD(fsig(m,n)) * pow2(gB(m,n)) * pow2(r(m,n)) 
    /*4*/  * pow2(R(m,n)) + 8 * exp(4 * fconf(m,n)) * fDB(m,n)
    /*4*/  * fsig(m,n) * r(m,n) - 2 * exp(4 * fconf(m,n)) 
    /*4*/  * fsig_r(m,n) * r(m,n) + 2 * exp(4 * gconf(m,n)) 
    /*4*/  * fL(m,n) * fsig(m,n) * pow2(gB(m,n)) * pow2(R(m,n))
    /*4*/  * r(m,n)))) / (pow2(fA(m,n)) * pow3(r(m,n))) + (2
    /*2*/  * exp(-2 * fconf(m,n)) * p(m,n) * (2 * fAsig(m,n)
    /*3*/  + fAsig_r(m,n) * r(m,n))) / (fA(m,n) * r(m,n) 
    /*2*/  * Lt(m,n)));

Jacobian[24][4]=
	(4 * exp(-4 * fconf(m,n)) * fDAlp(m,n)) 
    /*0*/  / (pow2(fA(m,n)) * pow2(r(m,n))) + (2 * exp(-4 
    /*2*/  * fconf(m,n)) * fAlp(m,n) * (1 + fDA(m,n) * r(m,n) 
    /*2*/  + fDB(m,n) * r(m,n) + 4 * fDconf(m,n) * r(m,n))) 
    /*0*/  / (pow2(fA(m,n)) * pow3(r(m,n)));

Jacobian[24][6]=
	fAsig(m,n) * fAlp(m,n);

Jacobian[24][7]=
	k_f * ((exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*2*/  * gAlp(m,n) * (-P_1_0(R(m,n)) + P_2_0(R(m,n)))) 
    /*1*/  / (fA(m,n) * pow2(r(m,n)) * pow2(R(m,n))) + (exp(-2
    /*3*/  * fconf(m,n) + 2 * gconf(m,n)) * fAlp(m,n) 
    /*2*/  * (-P_1_1(R(m,n)) + pow2(Lt(m,n)) * P_2_1(R(m,n))))
    /*1*/  / (fA(m,n) * pow2(r(m,n)) * pow2(R(m,n)) 
    /*2*/  * Lt(m,n)));

Jacobian[24][8]=
	k_f * (gAlp(m,n) * ((exp(-2 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * gA(m,n) * (P_1_0(R(m,n)) + b_0)) 
    /*2*/  / (fA(m,n) * gB(m,n) * pow2(r(m,n)) * pow2(R(m,n)))
    /*2*/  + (-P_1_1(R(m,n)) + b_1) / (gB(m,n) * pow2(r(m,n))
    /*3*/  * pow2(R(m,n)) * Lt(m,n))) + fAlp(m,n) 
    /*1*/  * ((-P_1_2(R(m,n)) + b_2) / (gB(m,n) * pow2(r(m,n))
    /*3*/  * pow2(R(m,n))) + (exp(-2 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * gA(m,n) * ((-1 + 2 * pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) + b_1)) / (fA(m,n) * gB(m,n) 
    /*3*/  * pow2(r(m,n)) * pow2(R(m,n)) * Lt(m,n))));

Jacobian[24][11]=
	(-2 * exp(-4 * fconf(m,n)) * (-(fDAlp_r(m,n) 
    /*3*/  * r(m,n)) + fDAlp(m,n) * (1 + fDA(m,n) * r(m,n) 
    /*3*/  + fDB(m,n) * r(m,n) + 4 * fDconf(m,n) * r(m,n)))) 
    /*0*/  / (pow3(fA(m,n)) * pow3(r(m,n))) + k_f * (-((exp(-2
    /*4*/  * fconf(m,n) + 2 * gconf(m,n)) * gAlp(m,n) 
    /*3*/  * gA(m,n) * P_1_1(R(m,n))) / (pow2(fA(m,n)) 
    /*3*/  * pow2(r(m,n)) * R(m,n))) + (exp(-2 * fconf(m,n) + 2
    /*3*/  * gconf(m,n)) * fAlp(m,n) * gA(m,n) 
    /*2*/  * (P_1_1(R(m,n)) - pow2(Lt(m,n)) * P_2_1(R(m,n)))) 
    /*1*/  / (pow2(fA(m,n)) * pow2(r(m,n)) * pow2(R(m,n)) 
    /*2*/  * Lt(m,n))) + fAlp(m,n) * ((-2 * exp(-8 * fconf(m,n)
    /*3*/  + 4 * gconf(m,n)) * pow2(gB(m,n)) * pow2(R(m,n)) 
    /*2*/  * (4 * fsig_r(m,n) + fsig_rr(m,n) * r(m,n) + 2 
    /*3*/  * pow2(fsig(m,n)) * r(m,n))) / (Power(fA(m,n),5) 
    /*2*/  * r(m,n)) + (exp(-8 * fconf(m,n)) * (4 * exp(4 
    /*4*/  * fconf(m,n)) * fDA(m,n) * (2 * fDB(m,n) 
    /*4*/  - fDconf(m,n)) * r(m,n) - 6 * exp(4 * fconf(m,n)) 
    /*3*/  * pow2(fDA(m,n)) * r(m,n) - 8 * exp(4 * fconf(m,n))
    /*3*/  * pow2(fDconf(m,n)) * r(m,n) - 4 * exp(4 
    /*4*/  * fconf(m,n)) * fDconf(m,n) * (1 + fDB(m,n) 
    /*4*/  * r(m,n)) + r(m,n) * (4 * exp(4 * fconf(m,n)) 
    /*4*/  * fDconf_r(m,n) + exp(4 * gconf(m,n)) 
    /*4*/  * fLD(fsig(m,n)) * pow2(gB(m,n)) * pow2(r(m,n)) 
    /*4*/  * pow2(R(m,n)) + 8 * exp(4 * fconf(m,n)) * fDB(m,n)
    /*4*/  * fsig(m,n) * r(m,n) - 2 * exp(4 * fconf(m,n)) 
    /*4*/  * fsig_r(m,n) * r(m,n) + 2 * exp(4 * gconf(m,n)) 
    /*4*/  * fL(m,n) * fsig(m,n) * pow2(gB(m,n)) * pow2(R(m,n))
    /*4*/  * r(m,n)))) / (pow3(fA(m,n)) * pow3(r(m,n))) 
    /*1*/  + (exp(-2 * fconf(m,n)) * p(m,n) * (2 * fAsig(m,n) 
    /*3*/  + fAsig_r(m,n) * r(m,n))) / (pow2(fA(m,n)) * r(m,n)
    /*2*/  * Lt(m,n)));

Jacobian[24][12]=
	fAlp(m,n) * (-((exp(-6 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * gB(m,n) * (2 * fL(m,n) * fsig(m,n) 
    /*4*/  + fLD(fsig(m,n)) * r(m,n)) * R(m,n)) 
    /*2*/  / (pow2(fA(m,n)) * r(m,n))) + (exp(-6 * fconf(m,n) 
    /*3*/  + 2 * gconf(m,n)) * gB(m,n) * (4 * fsig_r(m,n) 
    /*3*/  + fsig_rr(m,n) * r(m,n) + 2 * pow2(fsig(m,n)) 
    /*3*/  * r(m,n)) * R(m,n)) / (Power(fA(m,n),4) * r(m,n))) 
    /*0*/  + k_f * (gAlp(m,n) * (-((gA(m,n) * (P_1_0(R(m,n)) 
    /*5*/  + b_0)) / (fA(m,n) * gB(m,n) * pow2(r(m,n)) 
    /*4*/  * pow3(R(m,n)))) + (exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * (P_1_1(R(m,n)) - b_1)) / (gB(m,n) 
    /*3*/  * pow2(r(m,n)) * pow3(R(m,n)) * Lt(m,n))) 
    /*1*/  + fAlp(m,n) * ((exp(2 * fconf(m,n) - 2 * gconf(m,n))
    /*3*/  * (P_1_2(R(m,n)) - b_2)) / (gB(m,n) * pow2(r(m,n))
    /*3*/  * pow3(R(m,n))) + (gA(m,n) * ((1 - 2 
    /*5*/  * pow2(Lt(m,n))) * P_1_1(R(m,n)) - b_1)) / (fA(m,n)
    /*3*/  * gB(m,n) * pow2(r(m,n)) * pow3(R(m,n)) 
    /*3*/  * Lt(m,n))));

Jacobian[24][13]=
	(exp(-4 * fconf(m,n)) * fDAlp(m,n)) 
    /*0*/  / (pow2(fA(m,n)) * pow2(r(m,n))) + (2 * exp(-4 
    /*2*/  * fconf(m,n)) * (3 * fDA(m,n) - 2 * fDB(m,n) 
    /*2*/  + fDconf(m,n)) * fAlp(m,n)) / (pow2(fA(m,n)) 
    /*1*/  * pow2(r(m,n)));

Jacobian[24][14]=
	(exp(-4 * fconf(m,n)) * fDAlp(m,n)) 
    /*0*/  / (pow2(fA(m,n)) * pow2(r(m,n))) + (exp(-4 
    /*2*/  * fconf(m,n)) * fAlp(m,n) * (-4 * fDA(m,n) + 2 
    /*2*/  * fDconf(m,n) - 4 * fsig(m,n) * r(m,n))) 
    /*0*/  / (pow2(fA(m,n)) * pow2(r(m,n)));

Jacobian[24][20]=
	fAlp(m,n) * (-(exp(-4 * fconf(m,n)) 
    /*2*/  / pow3(r(m,n))) - (exp(-8 * fconf(m,n) + 4 
    /*3*/  * gconf(m,n)) * pow2(gB(m,n)) * pow2(R(m,n)) * (2 
    /*3*/  * fsig(m,n) + fsig_r(m,n) * r(m,n))) / (2. 
    /*2*/  * pow2(fA(m,n)) * r(m,n)));

Jacobian[24][22]=
	fAlp(m,n) * ((2 * exp(-8 * fconf(m,n) + 4 
    /*3*/  * gconf(m,n)) * fsig(m,n) * pow2(gB(m,n)) 
    /*2*/  * pow2(R(m,n))) / Power(fA(m,n),4) - (4 * exp(-4 
    /*3*/  * fconf(m,n)) * fDB(m,n) + exp(-8 * fconf(m,n) + 4 
    /*3*/  * gconf(m,n)) * fL(m,n) * pow2(gB(m,n)) 
    /*2*/  * pow2(R(m,n))) / (pow2(fA(m,n)) * r(m,n)));

Jacobian[24][24]=
	2 * gBetr(m,n) - (2 * exp(-2 * gconf(m,n)) 
    /*1*/  * gAlp(m,n) * p(m,n)) / (gA(m,n) * Lt(m,n) * r(m,n))
    /*0*/  + fAlp(m,n) * (ftrA(m,n) + ftrK(m,n) - (2 * exp(-2
    /*3*/  * fconf(m,n)) * p(m,n)) / (fA(m,n) * r(m,n) 
    /*2*/  * Lt(m,n)));

Jacobian[25][1]=
	(-2 * exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) 
    /*1*/  * (2 * pfD(m,n) + pfD_r(m,n) * r(m,n))) / (gA(m,n) 
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[25][7]=
	-((exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) * (2
    /*3*/  * pfD(m,n) + pfD_r(m,n) * r(m,n))) 
    /*1*/  / (pow2(gA(m,n)) * r(m,n) * Lt(m,n)));

Jacobian[25][27]=
	gBet_r(m,n) + 2 * gBetr(m,n) - gDAlp(m,n) 
    /*0*/  * pfv(m,n) + gAlp(m,n) * (-pfv_r(m,n) - (2 
    /*2*/  * pfv(m,n)) / r(m,n) + p(m,n) * ((-2 * exp(-2 
    /*4*/  * gconf(m,n))) / (gA(m,n) * Lt(m,n) * r(m,n)) + (2 
    /*3*/  * exp(-2 * gconf(m,n))) / (gA(m,n) * r(m,n) 
    /*3*/  * Lt(m,n))));

Jacobian[26][1]=
	(-2 * exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) 
    /*1*/  * (2 * pfS(m,n) + pfS_r(m,n) * r(m,n))) / (gA(m,n) 
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[26][3]=
	2 * gAlp(m,n) * pfS(m,n) * pfv(m,n);

Jacobian[26][7]=
	-((exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) * (2
    /*3*/  * pfS(m,n) + pfS_r(m,n) * r(m,n))) 
    /*1*/  / (pow2(gA(m,n)) * r(m,n) * Lt(m,n)));

Jacobian[26][9]=
	gAlp(m,n) * pfS(m,n) * pfv(m,n);

Jacobian[26][27]=
	-gDAlp(m,n);

Jacobian[26][28]=
	2 * gBet_r(m,n) + 2 * gBetr(m,n) - gDAlp(m,n) 
    /*0*/  * pfv(m,n) + gAlp(m,n) * ((-(pfv_r(m,n) * r(m,n)) 
    /*2*/  + pfv(m,n) * (-2 + gDA(m,n) * r(m,n) + 2 
    /*3*/  * gDconf(m,n) * r(m,n))) / r(m,n) + p(m,n) * ((-2 
    /*3*/  * exp(-2 * gconf(m,n))) / (gA(m,n) * Lt(m,n) 
    /*3*/  * r(m,n)) + (2 * exp(-2 * gconf(m,n))) / (gA(m,n) 
    /*3*/  * r(m,n) * Lt(m,n))));

Jacobian[26][29]=
	-gDAlp(m,n);

Jacobian[27][1]=
	(4 * exp(-4 * gconf(m,n)) * gDAlp(m,n) * pfS(m,n))
    /*0*/  / pow2(gA(m,n)) - (2 * exp(-2 * gconf(m,n)) 
    /*1*/  * gAlp(m,n) * p(m,n) * (2 * pftau(m,n) 
    /*2*/  + pftau_r(m,n) * r(m,n))) / (gA(m,n) * r(m,n) 
    /*1*/  * Lt(m,n));

Jacobian[27][7]=
	(2 * exp(-4 * gconf(m,n)) * gDAlp(m,n) * pfS(m,n))
    /*0*/  / pow3(gA(m,n)) - (exp(-2 * gconf(m,n)) 
    /*1*/  * gAlp(m,n) * p(m,n) * (2 * pftau(m,n) 
    /*2*/  + pftau_r(m,n) * r(m,n))) / (pow2(gA(m,n)) * r(m,n)
    /*1*/  * Lt(m,n));

Jacobian[27][28]=
	gK1(m,n) * gAlp(m,n) * pfv(m,n) - (exp(-4 
    /*2*/  * gconf(m,n)) * gDAlp(m,n)) / pow2(gA(m,n));

Jacobian[27][29]=
	gBet_r(m,n) + 2 * gBetr(m,n) - gDAlp(m,n) 
    /*0*/  * pfv(m,n) + gAlp(m,n) * (-pfv_r(m,n) - (2 
    /*2*/  * pfv(m,n)) / r(m,n) + p(m,n) * ((-2 * exp(-2 
    /*4*/  * gconf(m,n))) / (gA(m,n) * Lt(m,n) * r(m,n)) + (2 
    /*3*/  * exp(-2 * gconf(m,n))) / (gA(m,n) * r(m,n) 
    /*3*/  * Lt(m,n))));

Jacobian[29][1]=
	k_g * pow3(gAlp(m,n)) * (-8 * exp(4 * gconf(m,n))
    /*1*/  * gj(m,n) - (8 * exp(2 * fconf(m,n)) * fA(m,n) 
    /*2*/  * p(m,n) * P_1_2(R(m,n)) * R(m,n)) / pow2(gA(m,n)))
    /*0*/  + p(m,n) * pow3(gAlp(m,n)) * ((-2 * exp(-2 
    /*3*/  * gconf(m,n)) * (-2 + gDA_r(m,n) * pow2(r(m,n)) + 2
    /*3*/  * gDB_r(m,n) * pow2(r(m,n)))) / (3. * pow2(r(m,n))
    /*2*/  * pow3(gA(m,n)) * Lt(m,n)) - (2 * exp(-2 
    /*3*/  * gconf(m,n)) * (-6 + 3 * gL_r(m,n) * pow2(gB(m,n))
    /*3*/  * pow2(r(m,n)) + 2 * gL(m,n) * pow2(gB(m,n)) 
    /*3*/  * r(m,n) * (2 + gDA(m,n) * r(m,n) + 2 * gDB(m,n) 
    /*4*/  * r(m,n)))) / (3. * gA(m,n) * pow2(gB(m,n)) 
    /*2*/  * pow2(r(m,n)) * Lt(m,n)));

Jacobian[29][2]=
	(-4 * k_g * exp(2 * fconf(m,n)) * fA(m,n) * p(m,n)
    /*1*/  * pow3(gAlp(m,n)) * (2 * P_1_1(R(m,n)) - 3 
    /*2*/  * P_2_1(R(m,n)))) / pow2(gA(m,n));

Jacobian[29][3]=
	(8 * gAsig(m,n) * pow2(r(m,n)) * pow3(gAlp(m,n)))
    /*0*/  / pow2(gA(m,n));

Jacobian[29][7]=
	(-4 * k_g * exp(2 * fconf(m,n)) * fA(m,n) * p(m,n)
    /*1*/  * pow3(gAlp(m,n)) * P_2_1(R(m,n))) / pow3(gA(m,n))
    /*0*/  - (2 * pow2(gAlp(m,n)) * (-2 * gBet(m,n) 
    /*2*/  + gDAD(gBet(m,n)) * pow2(r(m,n)) + 2 
    /*2*/  * gDBD(gBet(m,n)) * pow2(r(m,n)) - 4 * gA1(m,n) 
    /*2*/  * gDAlp(m,n) * pow2(r(m,n)) + 4 * gA2(m,n) 
    /*2*/  * gDAlp(m,n) * pow2(r(m,n)) + 4 * gBet_rr(m,n) 
    /*2*/  * pow2(r(m,n)) + gShiD(gDA(m,n)) * pow2(r(m,n)) + 2
    /*2*/  * gShiD(gDB(m,n)) * pow2(r(m,n)) + 2 * gBet_r(m,n)
    /*2*/  * r(m,n))) / (3. * pow2(r(m,n)) * pow3(gA(m,n))) 
    /*0*/  + pow3(gAlp(m,n)) * ((-8 * (gAsig(m,n) 
    /*3*/  * pow2(r(m,n)) * (gDA(m,n) + gDB(m,n) + 6 
    /*4*/  * gDconf(m,n) + gsig(m,n) * r(m,n)) - gtrA_r(m,n) 
    /*3*/  - gtrK_r(m,n))) / (3. * pow3(gA(m,n))) + p(m,n) 
    /*1*/  * ((2 * exp(-2 * gconf(m,n)) * (-2 + gDA_r(m,n) 
    /*4*/  * pow2(r(m,n)) + 2 * gDB_r(m,n) * pow2(r(m,n)))) 
    /*2*/  / (3. * Power(gA(m,n),4) * Lt(m,n) * pow2(r(m,n))) 
    /*2*/  - (exp(-2 * gconf(m,n)) * (-2 + gDA_r(m,n) 
    /*4*/  * pow2(r(m,n)) + 2 * gDB_r(m,n) * pow2(r(m,n)))) 
    /*2*/  / (Power(gA(m,n),4) * pow2(r(m,n)) * Lt(m,n)) 
    /*2*/  + (exp(-2 * gconf(m,n)) * (6 - 3 * gL_r(m,n) 
    /*4*/  * pow2(gB(m,n)) * pow2(r(m,n)) - 2 * gL(m,n) 
    /*4*/  * pow2(gB(m,n)) * r(m,n) * (2 + gDA(m,n) * r(m,n) 
    /*5*/  + 2 * gDB(m,n) * r(m,n)))) / (3. * pow2(gA(m,n)) 
    /*3*/  * pow2(gB(m,n)) * pow2(r(m,n)) * Lt(m,n))));

Jacobian[29][8]=
	(4 * pow2(gAlp(m,n)) * (gBet(m,n) - gBet_r(m,n) 
    /*2*/  * r(m,n))) / (pow2(r(m,n)) * pow3(gB(m,n))) - (4 
    /*1*/  * k_g * exp(2 * fconf(m,n)) * fA(m,n) * p(m,n) 
    /*1*/  * pow3(gAlp(m,n)) * P_1_2(R(m,n)) * R(m,n)) 
    /*0*/  / (gB(m,n) * pow2(gA(m,n))) + p(m,n) 
    /*0*/  * pow3(gAlp(m,n)) * ((-4 * exp(-2 * gconf(m,n))) 
    /*1*/  / (gA(m,n) * Lt(m,n) * pow2(r(m,n)) * pow3(gB(m,n)))
    /*1*/  + (4 * exp(-2 * gconf(m,n))) / (gA(m,n) 
    /*2*/  * pow2(r(m,n)) * pow3(gB(m,n)) * Lt(m,n)));

Jacobian[29][9]=
	pow2(gAlp(m,n)) * ((2 * gBet(m,n) * gL(m,n)) / 3.
    /*1*/  + gBet_r(m,n) / (3. * pow2(gA(m,n)))) 
    /*0*/  + pow3(gAlp(m,n)) * ((4 * gAsig(m,n) * pow2(r(m,n)))
    /*1*/  / (3. * pow2(gA(m,n))) + p(m,n) * ((-2 * exp(-2 
    /*4*/  * gconf(m,n)) * gL(m,n)) / (3. * gA(m,n) * Lt(m,n))
    /*2*/  + (2 * exp(-2 * gconf(m,n)) * gL(m,n)) / (3. 
    /*3*/  * gA(m,n) * Lt(m,n))));

Jacobian[29][10]=
	pow2(gAlp(m,n)) * ((4 * gBet(m,n) * gL(m,n)) / 3.
    /*1*/  + (2 * gBet_r(m,n)) / (3. * pow2(gA(m,n)))) 
    /*0*/  + pow3(gAlp(m,n)) * ((4 * gAsig(m,n) * pow2(r(m,n)))
    /*1*/  / (3. * pow2(gA(m,n))) + p(m,n) * ((-4 * exp(-2 
    /*4*/  * gconf(m,n)) * gL(m,n)) / (3. * gA(m,n) * Lt(m,n))
    /*2*/  + (4 * exp(-2 * gconf(m,n)) * gL(m,n)) / (3. 
    /*3*/  * gA(m,n) * Lt(m,n))));

Jacobian[29][11]=
	(2 * k_g * exp(2 * fconf(m,n)) * p(m,n) 
    /*1*/  * pow3(gAlp(m,n)) * P_2_1(R(m,n))) / pow2(gA(m,n));

Jacobian[29][12]=
	(4 * k_g * exp(4 * fconf(m,n) - 2 * gconf(m,n)) 
    /*1*/  * fA(m,n) * p(m,n) * pow3(gAlp(m,n)) 
    /*1*/  * P_1_2(R(m,n))) / (gB(m,n) * pow2(gA(m,n)));

Jacobian[29][15]=
	(-4 * gDAlp(m,n) * pow2(gAlp(m,n))) / (3. 
    /*1*/  * pow2(gA(m,n)));

Jacobian[29][16]=
	(4 * gDAlp(m,n) * pow2(gAlp(m,n))) / (3. 
    /*1*/  * pow2(gA(m,n)));

Jacobian[29][19]=
	((2 * gDA(m,n) * gBet(m,n) + 4 * gDB(m,n) 
    /*2*/  * gBet(m,n) - gBet_r(m,n) + 4 * gBetr(m,n)) 
    /*1*/  * pow2(gAlp(m,n))) / 3. + p(m,n) * pow3(gAlp(m,n)) 
    /*0*/  * ((-2 * exp(-2 * gconf(m,n)) * (2 + gDA(m,n) 
    /*3*/  * r(m,n) + 2 * gDB(m,n) * r(m,n))) / (3. * gA(m,n) 
    /*2*/  * Lt(m,n) * r(m,n)) + (2 * exp(-2 * gconf(m,n)) * (2
    /*3*/  + gDA(m,n) * r(m,n) + 2 * gDB(m,n) * r(m,n))) 
    /*1*/  / (3. * gA(m,n) * r(m,n) * Lt(m,n)));

Jacobian[29][21]=
	(4 * gAsig(m,n) * pow3(gAlp(m,n)) * pow3(r(m,n)))
    /*0*/  / (3. * pow2(gA(m,n)));

Jacobian[29][23]=
	(4 * pow2(r(m,n)) * pow3(gAlp(m,n)) * (gDA(m,n) 
    /*2*/  + gDB(m,n) + 6 * gDconf(m,n) + gsig(m,n) * r(m,n)))
    /*0*/  / (3. * pow2(gA(m,n)));

