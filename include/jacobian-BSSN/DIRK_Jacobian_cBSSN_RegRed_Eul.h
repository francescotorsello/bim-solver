/** @file  DIRK_Jacobian_cBSSN_RegRed_Eul.h
 *  @author Francesco Torsello
 *  @brief The Jacobian of the Eulerian regularized reduced cBSSN equations, needed by DIRK.
 *  @version 2019-06-20T14:45:06
 *  @image html DIRK_Jacobian_cBSSN_RegRed_Eul.png
 */

Jacobian[0][0]=
	(-2 * gDconf(m,n) * gAlp(m,n) * p(m,n))
    /*0*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*1*/  * Lt(m,n));

Jacobian[0][2]=
	gBet(m,n) + gAlp(m,n) * p(m,n) * (-(1 / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n))) + 1
    /*1*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*2*/  * Lt(m,n)));

Jacobian[0][4]=
	-gAlp(m,n) / 6.;

Jacobian[0][6]=
	-((gDconf(m,n) * gAlp(m,n) * p(m,n)) / (TINY_Real
    /*2*/  + exp(2 * gconf(m,n)) * Power(gA(m,n),2)
    /*2*/  * Lt(m,n)));

Jacobian[0][25]=
	gDconf(m,n);

Jacobian[1][1]=
	(2 * fDconf(m,n) * fAlp(m,n) * p(m,n))
    /*0*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*1*/  * Lt(m,n));

Jacobian[1][3]=
	gBet(m,n) - (gAlp(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n))
    /*0*/  - (fAlp(m,n) * p(m,n)) / (TINY_Real + exp(2
    /*2*/  * fconf(m,n)) * fA(m,n) * Lt(m,n));

Jacobian[1][5]=
	-fAlp(m,n) / 6.;

Jacobian[1][10]=
	(fDconf(m,n) * fAlp(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(2 * fconf(m,n)) * Power(fA(m,n),2) * Lt(m,n));

Jacobian[1][25]=
	fDconf(m,n);

Jacobian[2][0]=
	(-2 * gDconf(m,n) * gDAlp(m,n) * p(m,n))
    /*0*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*1*/  * Lt(m,n)) + gAlp(m,n) * ((-2 * gDconf(m,n)
    /*2*/  * p_r(m,n)) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * gA(m,n) * Lt(m,n)) + (2 * p(m,n) * (gDA(m,n)
    /*3*/  * gDconf(m,n) * Lt(m,n) - gDconf_r(m,n) * Lt(m,n)
    /*3*/  + 2 * pow2(gDconf(m,n)) * Lt(m,n) + gDconf(m,n)
    /*3*/  * Lt_r(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * gA(m,n) * Power(Lt(m,n),2)));

Jacobian[2][2]=
	gBet_r(m,n) + p(m,n) * (-(gAlp_r(m,n) / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n)))
    /*1*/  + gDAlp(m,n) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * gA(m,n) * Lt(m,n))) + gAlp(m,n) * (-(p_r(m,n)
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * Lt(m,n))) + p_r(m,n) / (TINY_Real + exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * Lt(m,n)) + p(m,n) * (((2
    /*4*/  * gconf_r(m,n)) / (TINY_Real + exp(2 * gconf(m,n))
    /*4*/  * gA(m,n)) + gA_r(m,n) / (TINY_Real + exp(2
    /*5*/  * gconf(m,n)) * Power(gA(m,n),2))) / (TINY_Real
    /*3*/  + Lt(m,n)) + Lt_r(m,n) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * Power(Lt(m,n),2))
    /*2*/  + (-(gDA(m,n) * Lt(m,n)) - 4 * gDconf(m,n) * Lt(m,n)
    /*3*/  - Lt_r(m,n)) / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * Power(Lt(m,n),2))));

Jacobian[2][4]=
	-(gDAlp(m,n) / (6 + TINY_Real));

Jacobian[2][6]=
	-((gDconf(m,n) * gDAlp(m,n) * p(m,n)) / (TINY_Real
    /*2*/  + exp(2 * gconf(m,n)) * Power(gA(m,n),2)
    /*2*/  * Lt(m,n))) + gAlp(m,n) * (-((gDconf(m,n)
    /*3*/  * p_r(m,n)) / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * Power(gA(m,n),2) * Lt(m,n))) + (p(m,n) * (gDA(m,n)
    /*3*/  * gDconf(m,n) * Lt(m,n) - gDconf_r(m,n) * Lt(m,n)
    /*3*/  + 2 * pow2(gDconf(m,n)) * Lt(m,n) + gDconf(m,n)
    /*3*/  * Lt_r(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2) * Power(Lt(m,n),2)));

Jacobian[2][8]=
	-((gDconf(m,n) * gAlp(m,n) * p(m,n)) / (TINY_Real
    /*2*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n)));

Jacobian[2][25]=
	gDconf_r(m,n);

Jacobian[3][1]=
	(2 * fDconf(m,n) * fDAlp(m,n) * p(m,n))
    /*0*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*1*/  * Lt(m,n)) + fAlp(m,n) * ((2 * fDconf(m,n)
    /*2*/  * p_r(m,n)) / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * Lt(m,n)) - (2 * p(m,n) * (fDA(m,n)
    /*3*/  * fDconf(m,n) * Lt(m,n) - fDconf_r(m,n) * Lt(m,n)
    /*3*/  + 2 * pow2(fDconf(m,n)) * Lt(m,n) + fDconf(m,n)
    /*3*/  * Lt_r(m,n))) / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * Power(Lt(m,n),2)));

Jacobian[3][3]=
	gBet_r(m,n) + gAlp(m,n) * ((((2 * gconf_r(m,n))
    /*3*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n))
    /*3*/  + gA_r(m,n) / (TINY_Real + exp(2 * gconf(m,n))
    /*4*/  * Power(gA(m,n),2))) / (TINY_Real + Lt(m,n))
    /*2*/  + Lt_r(m,n) / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * Power(Lt(m,n),2))) * p(m,n) - p_r(m,n)
    /*1*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*2*/  * Lt(m,n))) + p(m,n) * (-(gAlp_r(m,n) / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n)))
    /*1*/  - fDAlp(m,n) / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * Lt(m,n))) + fAlp(m,n) * (-(p_r(m,n)
    /*2*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*3*/  * Lt(m,n))) + (p(m,n) * (fDA(m,n) * Lt(m,n) + 4
    /*3*/  * fDconf(m,n) * Lt(m,n) + Lt_r(m,n))) / (TINY_Real
    /*2*/  + exp(2 * fconf(m,n)) * fA(m,n) * Power(Lt(m,n),2)));

Jacobian[3][5]=
	-(fDAlp(m,n) / (6 + TINY_Real));

Jacobian[3][10]=
	(fDconf(m,n) * fDAlp(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(2 * fconf(m,n)) * Power(fA(m,n),2) * Lt(m,n))
    /*0*/  + fAlp(m,n) * ((fDconf(m,n) * p_r(m,n))
    /*1*/  / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * Power(fA(m,n),2) * Lt(m,n)) + (p(m,n)
    /*2*/  * (-(fDA(m,n) * fDconf(m,n) * Lt(m,n))
    /*3*/  + fDconf_r(m,n) * Lt(m,n) - 2 * pow2(fDconf(m,n))
    /*3*/  * Lt(m,n) - fDconf(m,n) * Lt_r(m,n))) / (TINY_Real
    /*2*/  + exp(2 * fconf(m,n)) * Power(fA(m,n),2)
    /*2*/  * Power(Lt(m,n),2)));

Jacobian[3][12]=
	(fDconf(m,n) * fAlp(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(2 * fconf(m,n)) * fA(m,n) * Lt(m,n));

Jacobian[3][25]=
	fDconf_r(m,n);

Jacobian[4][0]=
	(4 * gDAlp_r(m,n) * r(m,n) + gDAlp(m,n) * (8 - 4
    /*2*/  * gDA(m,n) * r(m,n) + 8 * gDB(m,n) * r(m,n) + 8
    /*2*/  * gDconf(m,n) * r(m,n))) / (TINY_Real + exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2) * r(m,n)) - (2
    /*1*/  * gAlp(m,n) * p(m,n) * gtrK_r(m,n)) / (TINY_Real
    /*1*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n)) + k_g
    /*0*/  * (fAlp(m,n) * ((-2 * exp(2 * fconf(m,n)) * fA(m,n)
    /*3*/  * (2 * P_1_2(R(m,n)) + b_2)) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n)) - (2 * (P_2_1(R(m,n))
    /*4*/  + b_1)) / (TINY_Real + Lt(m,n))) + gAlp(m,n) * (-2
    /*2*/  * (P_1_0(R(m,n)) + b_0) - (exp(2 * fconf(m,n))
    /*3*/  * fA(m,n) * ((2 + 4 * pow2(Lt(m,n))) * P_1_1(R(m,n))
    /*4*/  + (3 - 6 * pow2(Lt(m,n))) * P_2_1(R(m,n)) + 2
    /*4*/  * b_1)) / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * Lt(m,n))));

Jacobian[4][1]=
	k_g * (fAlp(m,n) * ((2 * exp(2 * fconf(m,n))
    /*3*/  * fA(m,n) * (2 * P_1_2(R(m,n)) + b_2)) / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * gA(m,n)) + (2
    /*3*/  * (P_2_1(R(m,n)) + b_1)) / (TINY_Real + Lt(m,n)))
    /*1*/  + gAlp(m,n) * (2 * (P_1_0(R(m,n)) + b_0) + (exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * ((2 + 4 * pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) + (3 - 6 * pow2(Lt(m,n)))
    /*4*/  * P_2_1(R(m,n)) + 2 * b_1)) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * Lt(m,n))));

Jacobian[4][2]=
	(-2 * gDAlp(m,n)) / (TINY_Real + exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2));

Jacobian[4][4]=
	(2 * gAlp(m,n) * gtrK(m,n)) / (3 + TINY_Real);

Jacobian[4][6]=
	(2 * gDAlp_r(m,n) * r(m,n) + gDAlp(m,n) * (4 - 2
    /*2*/  * gDA(m,n) * r(m,n) + 4 * gDB(m,n) * r(m,n) + 4
    /*2*/  * gDconf(m,n) * r(m,n))) / (TINY_Real + exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),3) * r(m,n))
    /*0*/  - (gAlp(m,n) * p(m,n) * gtrK_r(m,n)) / (TINY_Real
    /*1*/  + exp(2 * gconf(m,n)) * Power(gA(m,n),2) * Lt(m,n))
    /*0*/  + k_g * (-((exp(2 * fconf(m,n)) * fAlp(m,n)
    /*3*/  * fA(m,n) * P_1_2(R(m,n))) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * Power(gA(m,n),2))) - (exp(2
    /*3*/  * fconf(m,n)) * fA(m,n) * gAlp(m,n) * (2
    /*3*/  * P_1_1(R(m,n)) + (1 - 2 * pow2(Lt(m,n)))
    /*3*/  * P_2_1(R(m,n)))) / (TINY_Real + 2 * exp(2
    /*3*/  * gconf(m,n)) * Power(gA(m,n),2) * Lt(m,n)));

Jacobian[4][7]=
	k_g * (fAlp(m,n) * ((exp(2 * fconf(m,n)) * fA(m,n)
    /*3*/  * R(m,n) * b_3) / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * gB(m,n)) + (R(m,n) * (-P_1_2(R(m,n))
    /*4*/  + b_2)) / (TINY_Real + gB(m,n) * Lt(m,n)))
    /*1*/  + gAlp(m,n) * ((R(m,n) * b_1) / (TINY_Real
    /*3*/  + gB(m,n)) + (exp(2 * fconf(m,n)) * fA(m,n) * R(m,n)
    /*3*/  * ((-1 + 2 * pow2(Lt(m,n))) * P_1_2(R(m,n))
    /*4*/  + b_2)) / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * gB(m,n) * Lt(m,n))));

Jacobian[4][8]=
	gDAlp(m,n) / (TINY_Real + exp(4 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2));

Jacobian[4][9]=
	(-2 * gDAlp(m,n)) / (TINY_Real + exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2));

Jacobian[4][10]=
	k_g * ((exp(2 * fconf(m,n)) * fAlp(m,n)
    /*2*/  * P_1_2(R(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * gA(m,n)) + (exp(2 * fconf(m,n)) * gAlp(m,n) * (2
    /*3*/  * P_1_1(R(m,n)) + (1 - 2 * pow2(Lt(m,n)))
    /*3*/  * P_2_1(R(m,n)))) / (TINY_Real + 2 * exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * Lt(m,n)));

Jacobian[4][11]=
	k_g * (fAlp(m,n) * (-((exp(4 * fconf(m,n))
    /*4*/  * fA(m,n) * b_3) / (TINY_Real + exp(4 * gconf(m,n))
    /*4*/  * gA(m,n) * gB(m,n))) + (exp(2 * fconf(m,n))
    /*3*/  * (P_1_2(R(m,n)) - b_2)) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gB(m,n) * Lt(m,n))) + gAlp(m,n)
    /*1*/  * (-((exp(2 * fconf(m,n)) * b_1) / (TINY_Real
    /*4*/  + exp(2 * gconf(m,n)) * gB(m,n))) - (exp(4
    /*4*/  * fconf(m,n)) * fA(m,n) * ((-1 + 2 * pow2(Lt(m,n)))
    /*4*/  * P_1_2(R(m,n)) + b_2)) / (TINY_Real + exp(4
    /*4*/  * gconf(m,n)) * gA(m,n) * gB(m,n) * Lt(m,n))));

Jacobian[4][14]=
	2 * gA1(m,n) * gAlp(m,n);

Jacobian[4][25]=
	gtrK_r(m,n);

Jacobian[5][0]=
	k_f * (gAlp(m,n) * ((2 * exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * (P_2_0(R(m,n)) + b_0)) / (TINY_Real
    /*3*/  + exp(2 * fconf(m,n)) * fA(m,n) * Power(R(m,n),2))
    /*2*/  + (2 * (2 * P_1_1(R(m,n)) + b_1)) / (TINY_Real
    /*3*/  + Power(R(m,n),2) * Lt(m,n))) + fAlp(m,n) * ((2
    /*3*/  * (P_1_2(R(m,n)) + b_2)) / (TINY_Real
    /*3*/  + Power(R(m,n),2)) + (exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * ((2 - 4 * pow2(Lt(m,n))) * P_1_1(R(m,n)) + (3
    /*5*/  - 2 * pow2(Lt(m,n))) * P_2_1(R(m,n)) + 2 * b_1))
    /*2*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*3*/  * Power(R(m,n),2) * Lt(m,n))));

Jacobian[5][1]=
	(4 * fDAlp_r(m,n) * r(m,n) + fDAlp(m,n) * (8 - 4
    /*2*/  * fDA(m,n) * r(m,n) + 8 * fDB(m,n) * r(m,n) + 8
    /*2*/  * fDconf(m,n) * r(m,n))) / (TINY_Real + exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),2) * r(m,n)) + (2
    /*1*/  * fAlp(m,n) * p(m,n) * ftrK_r(m,n)) / (TINY_Real
    /*1*/  + exp(2 * fconf(m,n)) * fA(m,n) * Lt(m,n)) + k_f
    /*0*/  * (gAlp(m,n) * ((-2 * exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * (P_2_0(R(m,n)) + b_0)) / (TINY_Real + exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * Power(R(m,n),2)) - (2 * (2
    /*4*/  * P_1_1(R(m,n)) + b_1)) / (TINY_Real
    /*3*/  + Power(R(m,n),2) * Lt(m,n))) + fAlp(m,n) * ((-2
    /*3*/  * (P_1_2(R(m,n)) + b_2)) / (TINY_Real
    /*3*/  + Power(R(m,n),2)) + (exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * ((-2 + 4 * pow2(Lt(m,n))) * P_1_1(R(m,n)) + (-3
    /*5*/  + 2 * pow2(Lt(m,n))) * P_2_1(R(m,n)) - 2 * b_1))
    /*2*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*3*/  * Power(R(m,n),2) * Lt(m,n))));

Jacobian[5][3]=
	(-2 * fDAlp(m,n)) / (TINY_Real + exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),2));

Jacobian[5][5]=
	(2 * fAlp(m,n) * ftrK(m,n)) / (3 + TINY_Real);

Jacobian[5][6]=
	k_f * ((exp(2 * gconf(m,n)) * gAlp(m,n)
    /*2*/  * (-P_1_0(R(m,n)) + P_2_0(R(m,n)))) / (TINY_Real
    /*2*/  + exp(2 * fconf(m,n)) * fA(m,n) * Power(R(m,n),2))
    /*1*/  + (exp(2 * gconf(m,n)) * fAlp(m,n) * (-2
    /*3*/  * P_1_1(R(m,n)) + (3 - 2 * pow2(Lt(m,n)))
    /*3*/  * P_2_1(R(m,n)))) / (TINY_Real + 2 * exp(2
    /*3*/  * fconf(m,n)) * fA(m,n) * Power(R(m,n),2)
    /*2*/  * Lt(m,n)));

Jacobian[5][7]=
	k_f * (gAlp(m,n) * ((exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * (P_1_0(R(m,n)) + b_0)) / (TINY_Real + exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * gB(m,n) * Power(R(m,n),2))
    /*2*/  + (-b_1 - 2 * R(m,n) * b_2) / (TINY_Real + gB(m,n)
    /*3*/  * Power(R(m,n),2) * Lt(m,n))) + fAlp(m,n) * (-(b_3
    /*3*/  / (TINY_Real + gB(m,n) * R(m,n))) + (exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * (-2 * (-1 + pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) + b_1)) / (TINY_Real + exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * gB(m,n) * Power(R(m,n),2)
    /*3*/  * Lt(m,n))));

Jacobian[5][10]=
	(2 * fDAlp_r(m,n) * r(m,n) + fDAlp(m,n) * (4 - 2
    /*2*/  * fDA(m,n) * r(m,n) + 4 * fDB(m,n) * r(m,n) + 4
    /*2*/  * fDconf(m,n) * r(m,n))) / (TINY_Real + exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),3) * r(m,n))
    /*0*/  + (fAlp(m,n) * p(m,n) * ftrK_r(m,n)) / (TINY_Real
    /*1*/  + exp(2 * fconf(m,n)) * Power(fA(m,n),2) * Lt(m,n))
    /*0*/  + k_f * (-((exp(2 * gconf(m,n)) * gAlp(m,n)
    /*3*/  * gA(m,n) * P_1_1(R(m,n))) / (TINY_Real + exp(2
    /*4*/  * fconf(m,n)) * Power(fA(m,n),2) * R(m,n))) + (exp(2
    /*3*/  * gconf(m,n)) * fAlp(m,n) * gA(m,n) * (2
    /*3*/  * P_1_1(R(m,n)) + (-3 + 2 * pow2(Lt(m,n)))
    /*3*/  * P_2_1(R(m,n)))) / (TINY_Real + 2 * exp(2
    /*3*/  * fconf(m,n)) * Power(fA(m,n),2) * Power(R(m,n),2)
    /*2*/  * Lt(m,n)));

Jacobian[5][11]=
	k_f * (gAlp(m,n) * (-((gA(m,n) * (P_1_0(R(m,n))
    /*5*/  + b_0)) / (TINY_Real + fA(m,n) * gB(m,n)
    /*4*/  * Power(R(m,n),3))) + (exp(2 * fconf(m,n)) * (b_1
    /*4*/  + 2 * R(m,n) * b_2)) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gB(m,n) * Power(R(m,n),3)
    /*3*/  * Lt(m,n))) + fAlp(m,n) * ((exp(2 * fconf(m,n))
    /*3*/  * b_3) / (TINY_Real + exp(2 * gconf(m,n)) * gB(m,n)
    /*3*/  * Power(R(m,n),2)) + (gA(m,n) * (2 * (-1
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) - b_1))
    /*2*/  / (TINY_Real + fA(m,n) * gB(m,n) * Power(R(m,n),3)
    /*3*/  * Lt(m,n))));

Jacobian[5][12]=
	fDAlp(m,n) / (TINY_Real + exp(4 * fconf(m,n))
    /*1*/  * Power(fA(m,n),2));

Jacobian[5][13]=
	(-2 * fDAlp(m,n)) / (TINY_Real + exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),2));

Jacobian[5][15]=
	2 * fA1(m,n) * fAlp(m,n);

Jacobian[5][25]=
	ftrK_r(m,n);

Jacobian[6][0]=
	(-2 * gDA(m,n) * gAlp(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(2 * gconf(m,n)) * Lt(m,n));

Jacobian[6][6]=
	gDA(m,n) * gBet(m,n) + gBet_r(m,n) + gAlp(m,n)
    /*0*/  * (-gA1(m,n) - (gDA(m,n) * p(m,n)) / (TINY_Real
    /*2*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n)));

Jacobian[6][8]=
	gBet(m,n) * gA(m,n) + gAlp(m,n) * p(m,n) * (-(1
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n)) * Lt(m,n))) + 1
    /*1*/  / (TINY_Real + exp(2 * gconf(m,n)) * Lt(m,n)));

Jacobian[6][14]=
	-(gAlp(m,n) * gA(m,n));

Jacobian[6][25]=
	gDA(m,n) * gA(m,n);

Jacobian[7][0]=
	(-2 * gAlp(m,n) * gB(m,n) * p(m,n) * (1 + gDB(m,n)
    /*2*/  * r(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*1*/  * gA(m,n) * r(m,n) * Lt(m,n));

Jacobian[7][6]=
	-((gAlp(m,n) * gB(m,n) * p(m,n) * (1 + gDB(m,n)
    /*3*/  * r(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2) * r(m,n) * Lt(m,n)));

Jacobian[7][7]=
	gDB(m,n) * gBet(m,n) + gBetr(m,n) + gAlp(m,n)
    /*0*/  * (gA1(m,n) / (2 + TINY_Real) + p(m,n) * ((-1
    /*3*/  - gDB(m,n) * r(m,n)) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * Lt(m,n) * r(m,n)) + (1
    /*3*/  + gDB(m,n) * r(m,n)) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * r(m,n) * Lt(m,n))));

Jacobian[7][9]=
	gBet(m,n) * gB(m,n) + gAlp(m,n) * p(m,n)
    /*0*/  * (-(gB(m,n) / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * Lt(m,n))) + gB(m,n) / (TINY_Real + exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * Lt(m,n)));

Jacobian[7][25]=
	(gB(m,n) * (1 + gDB(m,n) * r(m,n))) / (TINY_Real
    /*1*/  + r(m,n));

Jacobian[8][0]=
	(-2 * gDA(m,n) * gDAlp(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n))
    /*0*/  + gAlp(m,n) * ((-2 * gDA(m,n) * p_r(m,n))
    /*1*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*2*/  * Lt(m,n)) + (2 * p(m,n) * (-(gDA_r(m,n) * Lt(m,n))
    /*3*/  + pow2(gDA(m,n)) * Lt(m,n) + gDA(m,n) * (2
    /*4*/  * gDconf(m,n) * Lt(m,n) + Lt_r(m,n)))) / (TINY_Real
    /*2*/  + exp(2 * gconf(m,n)) * gA(m,n)
    /*2*/  * Power(Lt(m,n),2)));

Jacobian[8][2]=
	(-2 * gDA(m,n) * gAlp(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n));

Jacobian[8][6]=
	-((gDA(m,n) * gDAlp(m,n) * p(m,n)) / (TINY_Real
    /*2*/  + exp(2 * gconf(m,n)) * Power(gA(m,n),2) * Lt(m,n)))
    /*0*/  + gAlp(m,n) * (-((gDA(m,n) * p_r(m,n))
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * Power(gA(m,n),2) * Lt(m,n))) + (p(m,n)
    /*2*/  * (-(gDA_r(m,n) * Lt(m,n)) + pow2(gDA(m,n))
    /*3*/  * Lt(m,n) + gDA(m,n) * (2 * gDconf(m,n) * Lt(m,n)
    /*4*/  + Lt_r(m,n)))) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2) * Power(Lt(m,n),2)));

Jacobian[8][8]=
	gBet_r(m,n) + p(m,n) * (-(gAlp_r(m,n) / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n)))
    /*1*/  + gDAlp(m,n) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * gA(m,n) * Lt(m,n))) + gAlp(m,n) * (-(p_r(m,n)
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * Lt(m,n))) + p_r(m,n) / (TINY_Real + exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * Lt(m,n)) + p(m,n) * (((2
    /*4*/  * gconf_r(m,n)) / (TINY_Real + exp(2 * gconf(m,n))
    /*4*/  * gA(m,n)) + gA_r(m,n) / (TINY_Real + exp(2
    /*5*/  * gconf(m,n)) * Power(gA(m,n),2))) / (TINY_Real
    /*3*/  + Lt(m,n)) + Lt_r(m,n) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * Power(Lt(m,n),2)) + (-2
    /*3*/  * gDA(m,n) * Lt(m,n) - 2 * gDconf(m,n) * Lt(m,n)
    /*3*/  - Lt_r(m,n)) / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * Power(Lt(m,n),2))));

Jacobian[8][14]=
	-gDAlp(m,n);

Jacobian[8][25]=
	gDA_r(m,n);

Jacobian[9][0]=
	(-2 * gDAlp(m,n) * p(m,n) * (1 + gDB(m,n)
    /*2*/  * r(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*1*/  * gA(m,n) * r(m,n) * Lt(m,n)) + gAlp(m,n) * ((-2
    /*2*/  * p_r(m,n) * (1 + gDB(m,n) * r(m,n))) / (TINY_Real
    /*2*/  + exp(2 * gconf(m,n)) * gA(m,n) * r(m,n) * Lt(m,n))
    /*1*/  + (2 * p(m,n) * ((1 - gDB_r(m,n) * pow2(r(m,n))
    /*4*/  + gDA(m,n) * r(m,n) * (1 + gDB(m,n) * r(m,n)) + 2
    /*4*/  * gDconf(m,n) * r(m,n) * (1 + gDB(m,n) * r(m,n)))
    /*3*/  * Lt(m,n) + r(m,n) * (1 + gDB(m,n) * r(m,n))
    /*3*/  * Lt_r(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * gA(m,n) * Power(r(m,n),2) * Power(Lt(m,n),2)));

Jacobian[9][2]=
	(-2 * gAlp(m,n) * p(m,n) * (1 + gDB(m,n)
    /*2*/  * r(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*1*/  * gA(m,n) * r(m,n) * Lt(m,n));

Jacobian[9][6]=
	-((gDAlp(m,n) * p(m,n) * (1 + gDB(m,n) * r(m,n)))
    /*1*/  / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2) * r(m,n) * Lt(m,n))) + gAlp(m,n)
    /*0*/  * (-((p_r(m,n) * (1 + gDB(m,n) * r(m,n)))
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * Power(gA(m,n),2) * r(m,n) * Lt(m,n))) + (p(m,n)
    /*2*/  * ((1 - gDB_r(m,n) * pow2(r(m,n)) + gDA(m,n)
    /*4*/  * r(m,n) * (1 + gDB(m,n) * r(m,n)) + 2 * gDconf(m,n)
    /*4*/  * r(m,n) * (1 + gDB(m,n) * r(m,n))) * Lt(m,n)
    /*3*/  + r(m,n) * (1 + gDB(m,n) * r(m,n)) * Lt_r(m,n)))
    /*1*/  / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2) * Power(r(m,n),2)
    /*2*/  * Power(Lt(m,n),2)));

Jacobian[9][8]=
	(gAlp(m,n) * p(m,n) * (-1 - gDB(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[9][9]=
	gBet_r(m,n) + p(m,n) * (-(gAlp_r(m,n) / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n)))
    /*1*/  + gDAlp(m,n) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * gA(m,n) * Lt(m,n))) + gAlp(m,n) * (-(p_r(m,n)
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * Lt(m,n))) + p_r(m,n) / (TINY_Real + exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * Lt(m,n)) + p(m,n) * (((2
    /*4*/  * gconf_r(m,n)) / (TINY_Real + exp(2 * gconf(m,n))
    /*4*/  * gA(m,n)) + gA_r(m,n) / (TINY_Real + exp(2
    /*5*/  * gconf(m,n)) * Power(gA(m,n),2))) / (TINY_Real
    /*3*/  + Lt(m,n)) + Lt_r(m,n) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * Power(Lt(m,n),2))
    /*2*/  + (-(gDA(m,n) * Lt(m,n)) - 2 * gDconf(m,n) * Lt(m,n)
    /*3*/  - Lt_r(m,n)) / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * Power(Lt(m,n),2))));

Jacobian[9][25]=
	gDB_r(m,n) - 1 / (TINY_Real + Power(r(m,n),2));

Jacobian[10][1]=
	(2 * fDA(m,n) * fAlp(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(2 * fconf(m,n)) * Lt(m,n));

Jacobian[10][10]=
	-(fA1(m,n) * fAlp(m,n)) + fBet_r(m,n) + fDA(m,n)
    /*0*/  * gBet(m,n) - (fDA(m,n) * gAlp(m,n) * p(m,n))
    /*0*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*1*/  * Lt(m,n));

Jacobian[10][12]=
	fA(m,n) * gBet(m,n) - (fA(m,n) * gAlp(m,n)
    /*1*/  * p(m,n)) / (TINY_Real + exp(2 * gconf(m,n))
    /*1*/  * gA(m,n) * Lt(m,n)) - (fAlp(m,n) * p(m,n))
    /*0*/  / (TINY_Real + exp(2 * fconf(m,n)) * Lt(m,n));

Jacobian[10][15]=
	-(fAlp(m,n) * fA(m,n));

Jacobian[10][25]=
	fDA(m,n) * fA(m,n);

Jacobian[11][1]=
	(2 * exp(2 * gconf(m,n)) * fAlp(m,n) * gB(m,n)
    /*1*/  * p(m,n) * (1 + fDB(m,n) * r(m,n)) * R(m,n))
    /*0*/  / (TINY_Real + exp(4 * fconf(m,n)) * fA(m,n)
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[11][10]=
	(exp(2 * gconf(m,n)) * fAlp(m,n) * gB(m,n)
    /*1*/  * p(m,n) * (1 + fDB(m,n) * r(m,n)) * R(m,n))
    /*0*/  / (TINY_Real + exp(4 * fconf(m,n))
    /*1*/  * Power(fA(m,n),2) * r(m,n) * Lt(m,n));

Jacobian[11][11]=
	fDB(m,n) * gBet(m,n) + gBetr(m,n) + (gAlp(m,n)
    /*1*/  * p(m,n) * (-1 - fDB(m,n) * r(m,n))) / (TINY_Real
    /*1*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n) * r(m,n))
    /*0*/  + fAlp(m,n) * (fA1(m,n) / (2 + TINY_Real)
    /*1*/  + (p(m,n) * (-1 - fDB(m,n) * r(m,n))) / (TINY_Real
    /*2*/  + exp(2 * fconf(m,n)) * fA(m,n) * r(m,n) * Lt(m,n)));

Jacobian[11][13]=
	(exp(2 * gconf(m,n)) * gBet(m,n) * gB(m,n)
    /*1*/  * R(m,n)) / (TINY_Real + exp(2 * fconf(m,n)))
    /*0*/  - (gAlp(m,n) * gB(m,n) * p(m,n) * R(m,n))
    /*0*/  / (TINY_Real + exp(2 * fconf(m,n)) * gA(m,n)
    /*1*/  * Lt(m,n)) - (exp(2 * gconf(m,n)) * fAlp(m,n)
    /*1*/  * gB(m,n) * p(m,n) * R(m,n)) / (TINY_Real + exp(4
    /*2*/  * fconf(m,n)) * fA(m,n) * Lt(m,n));

Jacobian[11][25]=
	(exp(2 * gconf(m,n)) * gB(m,n) * (1 + fDB(m,n)
    /*2*/  * r(m,n)) * R(m,n)) / (TINY_Real + exp(2
    /*2*/  * fconf(m,n)) * r(m,n));

Jacobian[12][1]=
	(2 * fDA(m,n) * fDAlp(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(2 * fconf(m,n)) * fA(m,n) * Lt(m,n))
    /*0*/  + fAlp(m,n) * ((2 * fDA(m,n) * p_r(m,n))
    /*1*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*2*/  * Lt(m,n)) - (2 * p(m,n) * (-(fDA_r(m,n) * Lt(m,n))
    /*3*/  + pow2(fDA(m,n)) * Lt(m,n) + fDA(m,n) * (2
    /*4*/  * fDconf(m,n) * Lt(m,n) + Lt_r(m,n)))) / (TINY_Real
    /*2*/  + exp(2 * fconf(m,n)) * fA(m,n)
    /*2*/  * Power(Lt(m,n),2)));

Jacobian[12][3]=
	(2 * fDA(m,n) * fAlp(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(2 * fconf(m,n)) * fA(m,n) * Lt(m,n));

Jacobian[12][10]=
	(fDA(m,n) * fDAlp(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(2 * fconf(m,n)) * Power(fA(m,n),2) * Lt(m,n))
    /*0*/  + fAlp(m,n) * ((fDA(m,n) * p_r(m,n)) / (TINY_Real
    /*2*/  + exp(2 * fconf(m,n)) * Power(fA(m,n),2)
    /*2*/  * Lt(m,n)) + (p(m,n) * (fDA_r(m,n) * Lt(m,n)
    /*3*/  - pow2(fDA(m,n)) * Lt(m,n) - fDA(m,n) * (2
    /*4*/  * fDconf(m,n) * Lt(m,n) + Lt_r(m,n)))) / (TINY_Real
    /*2*/  + exp(2 * fconf(m,n)) * Power(fA(m,n),2)
    /*2*/  * Power(Lt(m,n),2)));

Jacobian[12][12]=
	gBet_r(m,n) + gAlp(m,n) * ((((2 * gconf_r(m,n))
    /*3*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n))
    /*3*/  + gA_r(m,n) / (TINY_Real + exp(2 * gconf(m,n))
    /*4*/  * Power(gA(m,n),2))) / (TINY_Real + Lt(m,n))
    /*2*/  + Lt_r(m,n) / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * Power(Lt(m,n),2))) * p(m,n) - p_r(m,n)
    /*1*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*2*/  * Lt(m,n))) + p(m,n) * (-(gAlp_r(m,n) / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n)))
    /*1*/  - fDAlp(m,n) / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * Lt(m,n))) + fAlp(m,n) * (-(p_r(m,n)
    /*2*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*3*/  * Lt(m,n))) + (p(m,n) * (2 * fDA(m,n) * Lt(m,n) + 2
    /*3*/  * fDconf(m,n) * Lt(m,n) + Lt_r(m,n))) / (TINY_Real
    /*2*/  + exp(2 * fconf(m,n)) * fA(m,n)
    /*2*/  * Power(Lt(m,n),2)));

Jacobian[12][15]=
	-fDAlp(m,n);

Jacobian[12][25]=
	fDA_r(m,n);

Jacobian[13][1]=
	(2 * fDAlp(m,n) * p(m,n) * (1 + fDB(m,n)
    /*2*/  * r(m,n))) / (TINY_Real + exp(2 * fconf(m,n))
    /*1*/  * fA(m,n) * r(m,n) * Lt(m,n)) + fAlp(m,n) * ((2
    /*2*/  * p_r(m,n) * (1 + fDB(m,n) * r(m,n))) / (TINY_Real
    /*2*/  + exp(2 * fconf(m,n)) * fA(m,n) * r(m,n) * Lt(m,n))
    /*1*/  - (2 * p(m,n) * ((1 - fDB_r(m,n) * pow2(r(m,n))
    /*4*/  + fDA(m,n) * r(m,n) * (1 + fDB(m,n) * r(m,n)) + 2
    /*4*/  * fDconf(m,n) * r(m,n) * (1 + fDB(m,n) * r(m,n)))
    /*3*/  * Lt(m,n) + r(m,n) * (1 + fDB(m,n) * r(m,n))
    /*3*/  * Lt_r(m,n))) / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * Power(r(m,n),2) * Power(Lt(m,n),2)));

Jacobian[13][3]=
	(2 * fAlp(m,n) * p(m,n) * (1 + fDB(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[13][10]=
	(fDAlp(m,n) * p(m,n) * (1 + fDB(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + exp(2 * fconf(m,n))
    /*1*/  * Power(fA(m,n),2) * r(m,n) * Lt(m,n)) + fAlp(m,n)
    /*0*/  * ((p_r(m,n) * (1 + fDB(m,n) * r(m,n))) / (TINY_Real
    /*2*/  + exp(2 * fconf(m,n)) * Power(fA(m,n),2) * r(m,n)
    /*2*/  * Lt(m,n)) + (p(m,n) * (-((1 - fDB_r(m,n)
    /*5*/  * pow2(r(m,n)) + fDA(m,n) * r(m,n) * (1 + fDB(m,n)
    /*6*/  * r(m,n)) + 2 * fDconf(m,n) * r(m,n) * (1 + fDB(m,n)
    /*6*/  * r(m,n))) * Lt(m,n)) - r(m,n) * (1 + fDB(m,n)
    /*4*/  * r(m,n)) * Lt_r(m,n))) / (TINY_Real + exp(2
    /*3*/  * fconf(m,n)) * Power(fA(m,n),2) * Power(r(m,n),2)
    /*2*/  * Power(Lt(m,n),2)));

Jacobian[13][12]=
	(fAlp(m,n) * p(m,n) * (1 + fDB(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[13][13]=
	gBet_r(m,n) + gAlp(m,n) * ((((2 * gconf_r(m,n))
    /*3*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n))
    /*3*/  + gA_r(m,n) / (TINY_Real + exp(2 * gconf(m,n))
    /*4*/  * Power(gA(m,n),2))) / (TINY_Real + Lt(m,n))
    /*2*/  + Lt_r(m,n) / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * Power(Lt(m,n),2))) * p(m,n) - p_r(m,n)
    /*1*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*2*/  * Lt(m,n))) + p(m,n) * (-(gAlp_r(m,n) / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n)))
    /*1*/  - fDAlp(m,n) / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * Lt(m,n))) + fAlp(m,n) * (-(p_r(m,n)
    /*2*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*3*/  * Lt(m,n))) + (p(m,n) * (fDA(m,n) * Lt(m,n) + 2
    /*3*/  * fDconf(m,n) * Lt(m,n) + Lt_r(m,n))) / (TINY_Real
    /*2*/  + exp(2 * fconf(m,n)) * fA(m,n) * Power(Lt(m,n),2)));

Jacobian[13][25]=
	fDB_r(m,n) - 1 / (TINY_Real + Power(r(m,n),2));

Jacobian[14][0]=
	(-8 * (-(gDAlp_r(m,n) * r(m,n)) + gDAlp(m,n) * (1
    /*3*/  + gDA(m,n) * r(m,n) + gDB(m,n) * r(m,n) + 4
    /*3*/  * gDconf(m,n) * r(m,n)))) / (TINY_Real + 3 * exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2) * r(m,n))
    /*0*/  + gAlp(m,n) * ((-8 * (2 * gDB(m,n) + pow2(gB(m,n))
    /*3*/  * (-gL(m,n) + gL_r(m,n) * r(m,n)))) / (TINY_Real + 3
    /*2*/  * exp(4 * gconf(m,n)) * Power(gB(m,n),2) * r(m,n))
    /*1*/  + (8 * (-2 * pow2(gDA(m,n)) * r(m,n) - 4
    /*3*/  * pow2(gDconf(m,n)) * r(m,n) - 2 * gDconf(m,n) * (1
    /*4*/  + gDB(m,n) * r(m,n)) + gDA(m,n) * (2 + 3
    /*4*/  * gDB(m,n) * r(m,n) - 2 * gDconf(m,n) * r(m,n))
    /*3*/  + r(m,n) * (gDA_r(m,n) - gDB_r(m,n) + 2
    /*4*/  * gDconf_r(m,n) + gsig(m,n) * (3 + 2 * gDB(m,n)
    /*5*/  * r(m,n))))) / (TINY_Real + 3 * exp(4 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2) * r(m,n)) - (2 * gA1_r(m,n)
    /*2*/  * p(m,n)) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * gA(m,n) * Lt(m,n))) + k_g * (fAlp(m,n) * ((-4
    /*3*/  * exp(2 * fconf(m,n)) * fA(m,n) * (2 * P_1_2(R(m,n))
    /*4*/  + b_2)) / (TINY_Real + 3 * exp(2 * gconf(m,n))
    /*3*/  * gA(m,n)) - (4 * (3 * P_1_1(R(m,n)) - 2
    /*4*/  * P_2_1(R(m,n)) + b_1)) / (TINY_Real + 3 * Lt(m,n)))
    /*1*/  + gAlp(m,n) * ((-4 * (3 * P_1_0(R(m,n)) - 2
    /*4*/  * P_2_0(R(m,n)) + b_0)) / (3 + TINY_Real) + (4
    /*3*/  * exp(2 * fconf(m,n)) * fA(m,n) * (2 * (-2
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) - 3 * (-1
    /*5*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)) - b_1))
    /*2*/  / (TINY_Real + 3 * exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * Lt(m,n))));

Jacobian[14][1]=
	k_g * (fAlp(m,n) * ((4 * exp(2 * fconf(m,n))
    /*3*/  * fA(m,n) * (2 * P_1_2(R(m,n)) + b_2)) / (TINY_Real
    /*3*/  + 3 * exp(2 * gconf(m,n)) * gA(m,n)) + (4 * R(m,n)
    /*3*/  * (b_2 + 2 * R(m,n) * b_3)) / (TINY_Real + 3
    /*3*/  * Lt(m,n))) + gAlp(m,n) * ((4 * R(m,n) * (b_1 + 2
    /*4*/  * R(m,n) * b_2)) / (3 + TINY_Real) - (4 * exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * (2 * (-2 + pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) - 3 * (-1 + pow2(Lt(m,n)))
    /*4*/  * P_2_1(R(m,n)) - b_1)) / (TINY_Real + 3 * exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * Lt(m,n))));

Jacobian[14][2]=
	(8 * gDAlp(m,n)) / (TINY_Real + 3 * exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2)) + (4 * gAlp(m,n)
    /*1*/  * (1 + gDA(m,n) * r(m,n) + gDB(m,n) * r(m,n) + 4
    /*2*/  * gDconf(m,n) * r(m,n))) / (TINY_Real + 3 * exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2) * r(m,n));

Jacobian[14][4]=
	gA1(m,n) * gAlp(m,n);

Jacobian[14][6]=
	(-4 * (-(gDAlp_r(m,n) * r(m,n)) + gDAlp(m,n) * (1
    /*3*/  + gDA(m,n) * r(m,n) + gDB(m,n) * r(m,n) + 4
    /*3*/  * gDconf(m,n) * r(m,n)))) / (TINY_Real + 3 * exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),3) * r(m,n))
    /*0*/  + gAlp(m,n) * ((4 * (-2 * pow2(gDA(m,n)) * r(m,n)
    /*3*/  - 4 * pow2(gDconf(m,n)) * r(m,n) - 2 * gDconf(m,n)
    /*3*/  * (1 + gDB(m,n) * r(m,n)) + gDA(m,n) * (2 + 3
    /*4*/  * gDB(m,n) * r(m,n) - 2 * gDconf(m,n) * r(m,n))
    /*3*/  + r(m,n) * (gDA_r(m,n) - gDB_r(m,n) + 2
    /*4*/  * gDconf_r(m,n) + gsig(m,n) * (3 + 2 * gDB(m,n)
    /*5*/  * r(m,n))))) / (TINY_Real + 3 * exp(4 * gconf(m,n))
    /*2*/  * Power(gA(m,n),3) * r(m,n)) - (gA1_r(m,n)
    /*2*/  * p(m,n)) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2) * Lt(m,n))) + k_g * ((-2 * exp(2
    /*3*/  * fconf(m,n)) * fAlp(m,n) * fA(m,n)
    /*2*/  * P_1_2(R(m,n))) / (TINY_Real + 3 * exp(2
    /*3*/  * gconf(m,n)) * Power(gA(m,n),2)) - (2 * exp(2
    /*3*/  * fconf(m,n)) * fA(m,n) * gAlp(m,n) * (P_1_1(R(m,n))
    /*3*/  + (-1 + pow2(Lt(m,n))) * P_2_1(R(m,n))))
    /*1*/  / (TINY_Real + 3 * exp(2 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2) * Lt(m,n)));

Jacobian[14][7]=
	(-8 * gDB(m,n) * gAlp(m,n)) / (TINY_Real + 3
    /*1*/  * exp(4 * gconf(m,n)) * Power(gB(m,n),3) * r(m,n))
    /*0*/  + k_g * (fAlp(m,n) * ((2 * exp(2 * fconf(m,n))
    /*3*/  * fA(m,n) * R(m,n) * b_3) / (TINY_Real + 3 * exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * gB(m,n)) - (2 * R(m,n)
    /*3*/  * (b_2 + 2 * R(m,n) * b_3)) / (TINY_Real + 3
    /*3*/  * gB(m,n) * Lt(m,n))) + gAlp(m,n) * ((-2 * R(m,n)
    /*3*/  * (b_1 + 2 * R(m,n) * b_2)) / (TINY_Real + 3
    /*3*/  * gB(m,n)) - (2 * exp(2 * fconf(m,n)) * fA(m,n)
    /*3*/  * R(m,n) * (2 * (-1 + pow2(Lt(m,n))) * P_1_2(R(m,n))
    /*4*/  - b_2)) / (TINY_Real + 3 * exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * gB(m,n) * Lt(m,n))));

Jacobian[14][8]=
	(2 * gDAlp(m,n)) / (TINY_Real + 3 * exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2)) + (2 * gAlp(m,n)
    /*1*/  * (-2 + 4 * gDA(m,n) * r(m,n) - 3 * gDB(m,n)
    /*2*/  * r(m,n) + 2 * gDconf(m,n) * r(m,n))) / (TINY_Real
    /*1*/  + 3 * exp(4 * gconf(m,n)) * Power(gA(m,n),2)
    /*1*/  * r(m,n));

Jacobian[14][9]=
	(2 * gDAlp(m,n)) / (TINY_Real + 3 * exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2)) + gAlp(m,n) * ((-2
    /*2*/  * (3 * gDA(m,n) - 2 * gDconf(m,n) + 2 * gsig(m,n)
    /*3*/  * r(m,n))) / (TINY_Real + 3 * exp(4 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2)) + 4 / (TINY_Real + 3 * exp(4
    /*3*/  * gconf(m,n)) * Power(gB(m,n),2) * r(m,n)));

Jacobian[14][10]=
	k_g * ((2 * exp(2 * fconf(m,n)) * fAlp(m,n)
    /*2*/  * P_1_2(R(m,n))) / (TINY_Real + 3 * exp(2
    /*3*/  * gconf(m,n)) * gA(m,n)) + (2 * exp(2 * fconf(m,n))
    /*2*/  * gAlp(m,n) * (P_1_1(R(m,n)) + (-1
    /*4*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)))) / (TINY_Real + 3
    /*2*/  * exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n)));

Jacobian[14][11]=
	k_g * (fAlp(m,n) * ((-2 * exp(4 * fconf(m,n))
    /*3*/  * fA(m,n) * b_3) / (TINY_Real + 3 * exp(4
    /*4*/  * gconf(m,n)) * gA(m,n) * gB(m,n)) - (2 * exp(2
    /*4*/  * fconf(m,n)) * (2 * P_1_2(R(m,n)) + b_2))
    /*2*/  / (TINY_Real + 3 * exp(2 * gconf(m,n)) * gB(m,n)
    /*3*/  * Lt(m,n))) + gAlp(m,n) * ((-2 * exp(2 * fconf(m,n))
    /*3*/  * (2 * P_1_1(R(m,n)) + b_1)) / (TINY_Real + 3
    /*3*/  * exp(2 * gconf(m,n)) * gB(m,n)) + (2 * exp(4
    /*4*/  * fconf(m,n)) * fA(m,n) * (2 * (-1 + pow2(Lt(m,n)))
    /*4*/  * P_1_2(R(m,n)) - b_2)) / (TINY_Real + 3 * exp(4
    /*4*/  * gconf(m,n)) * gA(m,n) * gB(m,n) * Lt(m,n))));

Jacobian[14][14]=
	gAlp(m,n) * gtrK(m,n);

Jacobian[14][16]=
	(-2 * gAlp(m,n)) / (TINY_Real + 3 * exp(4
    /*2*/  * gconf(m,n)) * r(m,n));

Jacobian[14][18]=
	(-2 * gAlp(m,n) * (3 + 2 * gDB(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + 3 * exp(4 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2));

Jacobian[14][25]=
	gA1_r(m,n);

Jacobian[15][0]=
	k_f * (gAlp(m,n) * ((4 * exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * (P_2_0(R(m,n)) + b_0)) / (TINY_Real + 3
    /*3*/  * exp(2 * fconf(m,n)) * fA(m,n) * Power(R(m,n),2))
    /*2*/  - (4 * (P_1_1(R(m,n)) - b_1)) / (TINY_Real + 3
    /*3*/  * Power(R(m,n),2) * Lt(m,n))) + fAlp(m,n) * ((-4
    /*3*/  * (P_1_2(R(m,n)) - b_2)) / (TINY_Real + 3
    /*3*/  * Power(R(m,n),2)) + (4 * exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * (2 * (-1 + pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) + pow2(Lt(m,n)) * P_2_1(R(m,n))
    /*4*/  + b_1)) / (TINY_Real + 3 * exp(2 * fconf(m,n))
    /*3*/  * fA(m,n) * Power(R(m,n),2) * Lt(m,n))));

Jacobian[15][1]=
	(-8 * (-(fDAlp_r(m,n) * r(m,n)) + fDAlp(m,n) * (1
    /*3*/  + fDA(m,n) * r(m,n) + fDB(m,n) * r(m,n) + 4
    /*3*/  * fDconf(m,n) * r(m,n)))) / (TINY_Real + 3 * exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),2) * r(m,n))
    /*0*/  + fAlp(m,n) * ((8 * (-(exp(-4 * fconf(m,n))
    /*4*/  * fL_r(m,n)) + (exp(-4 * fconf(m,n)) * fL(m,n))
    /*3*/  / r(m,n) - (2 * exp(-4 * gconf(m,n)) * fDB(m,n))
    /*3*/  / (pow2(gB(m,n)) * pow2(R(m,n)) * r(m,n)))) / (3
    /*2*/  + TINY_Real) + (8 * (-2 * pow2(fDA(m,n)) * r(m,n)
    /*3*/  - 4 * pow2(fDconf(m,n)) * r(m,n) - 2 * fDconf(m,n)
    /*3*/  * (1 + fDB(m,n) * r(m,n)) + fDA(m,n) * (2 + 3
    /*4*/  * fDB(m,n) * r(m,n) - 2 * fDconf(m,n) * r(m,n))
    /*3*/  + r(m,n) * (fDA_r(m,n) - fDB_r(m,n) + 2
    /*4*/  * fDconf_r(m,n) + fsig(m,n) * (3 + 2 * fDB(m,n)
    /*5*/  * r(m,n))))) / (TINY_Real + 3 * exp(4 * fconf(m,n))
    /*2*/  * Power(fA(m,n),2) * r(m,n)) + (2 * fA1_r(m,n)
    /*2*/  * p(m,n)) / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * Lt(m,n))) + k_f * (gAlp(m,n) * ((-4
    /*3*/  * exp(2 * gconf(m,n)) * gA(m,n) * (P_2_0(R(m,n))
    /*4*/  + b_0)) / (TINY_Real + 3 * exp(2 * fconf(m,n))
    /*3*/  * fA(m,n) * Power(R(m,n),2)) + (4 * (P_1_1(R(m,n))
    /*4*/  - b_1)) / (TINY_Real + 3 * Power(R(m,n),2)
    /*3*/  * Lt(m,n))) + fAlp(m,n) * ((4 * (P_1_2(R(m,n))
    /*4*/  - b_2)) / (TINY_Real + 3 * Power(R(m,n),2)) - (4
    /*3*/  * exp(2 * gconf(m,n)) * gA(m,n) * (2 * (-1
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) + pow2(Lt(m,n))
    /*4*/  * P_2_1(R(m,n)) + b_1)) / (TINY_Real + 3 * exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * Power(R(m,n),2)
    /*3*/  * Lt(m,n))));

Jacobian[15][3]=
	(8 * fDAlp(m,n)) / (TINY_Real + 3 * exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),2)) + (4 * fAlp(m,n)
    /*1*/  * (1 + fDA(m,n) * r(m,n) + fDB(m,n) * r(m,n) + 4
    /*2*/  * fDconf(m,n) * r(m,n))) / (TINY_Real + 3 * exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),2) * r(m,n));

Jacobian[15][5]=
	fA1(m,n) * fAlp(m,n);

Jacobian[15][6]=
	k_f * ((-2 * exp(2 * gconf(m,n)) * gAlp(m,n)
    /*2*/  * (P_1_0(R(m,n)) - P_2_0(R(m,n)))) / (TINY_Real + 3
    /*2*/  * exp(2 * fconf(m,n)) * fA(m,n) * Power(R(m,n),2))
    /*1*/  + (2 * exp(2 * gconf(m,n)) * fAlp(m,n)
    /*2*/  * (-P_1_1(R(m,n)) + pow2(Lt(m,n)) * P_2_1(R(m,n))))
    /*1*/  / (TINY_Real + 3 * exp(2 * fconf(m,n)) * fA(m,n)
    /*2*/  * Power(R(m,n),2) * Lt(m,n)));

Jacobian[15][7]=
	k_f * (gAlp(m,n) * ((2 * exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * (P_1_0(R(m,n)) + b_0)) / (TINY_Real + 3
    /*3*/  * exp(2 * fconf(m,n)) * fA(m,n) * gB(m,n)
    /*3*/  * Power(R(m,n),2)) - (2 * (P_1_1(R(m,n)) - b_1))
    /*2*/  / (TINY_Real + 3 * gB(m,n) * Power(R(m,n),2)
    /*3*/  * Lt(m,n))) + fAlp(m,n) * ((-2 * (P_1_2(R(m,n))
    /*4*/  - b_2)) / (TINY_Real + 3 * gB(m,n)
    /*3*/  * Power(R(m,n),2)) + (2 * exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * ((-1 + 2 * pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) + b_1)) / (TINY_Real + 3 * exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * gB(m,n) * Power(R(m,n),2)
    /*3*/  * Lt(m,n))));

Jacobian[15][10]=
	(-4 * (-(fDAlp_r(m,n) * r(m,n)) + fDAlp(m,n) * (1
    /*3*/  + fDA(m,n) * r(m,n) + fDB(m,n) * r(m,n) + 4
    /*3*/  * fDconf(m,n) * r(m,n)))) / (TINY_Real + 3 * exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),3) * r(m,n))
    /*0*/  + fAlp(m,n) * ((4 * (-2 * pow2(fDA(m,n)) * r(m,n)
    /*3*/  - 4 * pow2(fDconf(m,n)) * r(m,n) - 2 * fDconf(m,n)
    /*3*/  * (1 + fDB(m,n) * r(m,n)) + fDA(m,n) * (2 + 3
    /*4*/  * fDB(m,n) * r(m,n) - 2 * fDconf(m,n) * r(m,n))
    /*3*/  + r(m,n) * (fDA_r(m,n) - fDB_r(m,n) + 2
    /*4*/  * fDconf_r(m,n) + fsig(m,n) * (3 + 2 * fDB(m,n)
    /*5*/  * r(m,n))))) / (TINY_Real + 3 * exp(4 * fconf(m,n))
    /*2*/  * Power(fA(m,n),3) * r(m,n)) + (fA1_r(m,n)
    /*2*/  * p(m,n)) / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * Power(fA(m,n),2) * Lt(m,n))) + k_f * ((-2 * exp(2
    /*3*/  * gconf(m,n)) * gAlp(m,n) * gA(m,n)
    /*2*/  * P_1_1(R(m,n))) / (TINY_Real + 3 * exp(2
    /*3*/  * fconf(m,n)) * Power(fA(m,n),2) * R(m,n)) + (2
    /*2*/  * exp(2 * gconf(m,n)) * fAlp(m,n) * gA(m,n)
    /*2*/  * (P_1_1(R(m,n)) - pow2(Lt(m,n)) * P_2_1(R(m,n))))
    /*1*/  / (TINY_Real + 3 * exp(2 * fconf(m,n))
    /*2*/  * Power(fA(m,n),2) * Power(R(m,n),2) * Lt(m,n)));

Jacobian[15][11]=
	(-8 * exp(2 * fconf(m,n)) * fDB(m,n) * fAlp(m,n))
    /*0*/  / (TINY_Real + 3 * exp(6 * gconf(m,n))
    /*1*/  * Power(gB(m,n),3) * r(m,n) * Power(R(m,n),3)) + k_f
    /*0*/  * (gAlp(m,n) * ((-2 * gA(m,n) * (P_1_0(R(m,n))
    /*4*/  + b_0)) / (TINY_Real + 3 * fA(m,n) * gB(m,n)
    /*3*/  * Power(R(m,n),3)) + (2 * exp(2 * fconf(m,n))
    /*3*/  * (P_1_1(R(m,n)) - b_1)) / (TINY_Real + 3 * exp(2
    /*4*/  * gconf(m,n)) * gB(m,n) * Power(R(m,n),3)
    /*3*/  * Lt(m,n))) + fAlp(m,n) * ((2 * exp(2 * fconf(m,n))
    /*3*/  * (P_1_2(R(m,n)) - b_2)) / (TINY_Real + 3 * exp(2
    /*4*/  * gconf(m,n)) * gB(m,n) * Power(R(m,n),3)) - (2
    /*3*/  * gA(m,n) * ((-1 + 2 * pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) + b_1)) / (TINY_Real + 3 * fA(m,n)
    /*3*/  * gB(m,n) * Power(R(m,n),3) * Lt(m,n))));

Jacobian[15][12]=
	(2 * fDAlp(m,n)) / (TINY_Real + 3 * exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),2)) + (2 * fAlp(m,n)
    /*1*/  * (-2 + 4 * fDA(m,n) * r(m,n) - 3 * fDB(m,n)
    /*2*/  * r(m,n) + 2 * fDconf(m,n) * r(m,n))) / (TINY_Real
    /*1*/  + 3 * exp(4 * fconf(m,n)) * Power(fA(m,n),2)
    /*1*/  * r(m,n));

Jacobian[15][13]=
	(2 * fDAlp(m,n)) / (TINY_Real + 3 * exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),2)) + fAlp(m,n) * ((-2
    /*2*/  * (3 * fDA(m,n) - 2 * fDconf(m,n) + 2 * fsig(m,n)
    /*3*/  * r(m,n))) / (TINY_Real + 3 * exp(4 * fconf(m,n))
    /*2*/  * Power(fA(m,n),2)) + 4 / (TINY_Real + 3 * exp(4
    /*3*/  * gconf(m,n)) * Power(gB(m,n),2) * r(m,n)
    /*2*/  * Power(R(m,n),2)));

Jacobian[15][15]=
	fAlp(m,n) * ftrK(m,n);

Jacobian[15][17]=
	(-2 * fAlp(m,n)) / (TINY_Real + 3 * exp(4
    /*2*/  * fconf(m,n)) * r(m,n));

Jacobian[15][19]=
	(-2 * fAlp(m,n) * (3 + 2 * fDB(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + 3 * exp(4 * fconf(m,n))
    /*1*/  * Power(fA(m,n),2));

Jacobian[15][25]=
	fA1_r(m,n);

Jacobian[16][0]=
	k_g * gAlp(m,n) * (-8 * exp(4 * gconf(m,n))
    /*1*/  * gj(m,n) - (8 * exp(2 * fconf(m,n)) * fA(m,n)
    /*2*/  * p(m,n) * P_1_2(R(m,n)) * R(m,n)) / (TINY_Real
    /*2*/  + Power(gA(m,n),2))) + (gAlp(m,n) * p(m,n) * (4 - 2
    /*2*/  * gL_r(m,n) * pow2(gB(m,n)) * pow2(r(m,n))))
    /*0*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*1*/  * Power(gB(m,n),2) * Power(r(m,n),2) * Lt(m,n));

Jacobian[16][1]=
	(-4 * k_g * exp(2 * fconf(m,n)) * fA(m,n)
    /*1*/  * gAlp(m,n) * p(m,n) * (2 * P_1_1(R(m,n)) - 3
    /*2*/  * P_2_1(R(m,n)))) / (TINY_Real + Power(gA(m,n),2));

Jacobian[16][2]=
	(8 * gAsig(m,n) * gAlp(m,n) * pow2(r(m,n)))
    /*0*/  / (TINY_Real + Power(gA(m,n),2));

Jacobian[16][6]=
	(12 * gA1(m,n) * gDAlp(m,n) - 6 * gBet_rr(m,n))
    /*0*/  / (TINY_Real + 3 * Power(gA(m,n),3)) - (4 * k_g
    /*1*/  * exp(2 * fconf(m,n)) * fA(m,n) * gAlp(m,n) * p(m,n)
    /*1*/  * P_2_1(R(m,n))) / (TINY_Real + Power(gA(m,n),3))
    /*0*/  + gAlp(m,n) * ((-8 * (gAsig(m,n) * pow2(r(m,n))
    /*3*/  * (gDA(m,n) + gDB(m,n) + 6 * gDconf(m,n) + gsig(m,n)
    /*4*/  * r(m,n)) - gtrK_r(m,n))) / (TINY_Real + 3
    /*2*/  * Power(gA(m,n),3)) + (p(m,n) * (2 - gL_r(m,n)
    /*3*/  * pow2(gB(m,n)) * pow2(r(m,n)))) / (TINY_Real
    /*2*/  + exp(2 * gconf(m,n)) * Power(gA(m,n),2)
    /*2*/  * Power(gB(m,n),2) * Power(r(m,n),2) * Lt(m,n)));

Jacobian[16][7]=
	(4 * (gBet(m,n) - gBet_r(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + Power(gB(m,n),3) * Power(r(m,n),2))
    /*0*/  - (4 * k_g * exp(2 * fconf(m,n)) * fA(m,n)
    /*1*/  * gAlp(m,n) * p(m,n) * P_1_2(R(m,n)) * R(m,n))
    /*0*/  / (TINY_Real + Power(gA(m,n),2) * gB(m,n))
    /*0*/  + gAlp(m,n) * p(m,n) * (-4 / (TINY_Real + exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * Power(gB(m,n),3) * Lt(m,n)
    /*2*/  * Power(r(m,n),2)) + 4 / (TINY_Real + exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * Power(gB(m,n),3)
    /*2*/  * Power(r(m,n),2) * Lt(m,n)));

Jacobian[16][8]=
	(4 * gAsig(m,n) * gAlp(m,n) * pow2(r(m,n)))
    /*0*/  / (TINY_Real + 3 * Power(gA(m,n),2));

Jacobian[16][9]=
	(4 * gAsig(m,n) * gAlp(m,n) * pow2(r(m,n)))
    /*0*/  / (TINY_Real + 3 * Power(gA(m,n),2));

Jacobian[16][10]=
	(2 * k_g * exp(2 * fconf(m,n)) * gAlp(m,n)
    /*1*/  * p(m,n) * P_2_1(R(m,n))) / (TINY_Real
    /*1*/  + Power(gA(m,n),2));

Jacobian[16][11]=
	(4 * k_g * exp(4 * fconf(m,n)) * fA(m,n)
    /*1*/  * gAlp(m,n) * p(m,n) * P_1_2(R(m,n))) / (TINY_Real
    /*1*/  + exp(2 * gconf(m,n)) * Power(gA(m,n),2) * gB(m,n));

Jacobian[16][14]=
	(-2 * gDAlp(m,n)) / (TINY_Real + Power(gA(m,n),2));

Jacobian[16][16]=
	-gBet_r(m,n);

Jacobian[16][18]=
	(4 * gAsig(m,n) * gAlp(m,n) * pow3(r(m,n)))
    /*0*/  / (TINY_Real + 3 * Power(gA(m,n),2));

Jacobian[16][20]=
	(4 * gAlp(m,n) * pow2(r(m,n)) * (gDA(m,n)
    /*2*/  + gDB(m,n) + 6 * gDconf(m,n) + gsig(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + 3 * Power(gA(m,n),2));

Jacobian[16][25]=
	gL_r(m,n) - 2 / (TINY_Real + Power(gB(m,n),2)
    /*1*/  * Power(r(m,n),2));

Jacobian[17][0]=
	(-4 * k_f * exp(2 * gconf(m,n)) * fAlp(m,n)
    /*1*/  * gA(m,n) * p(m,n) * (2 * P_1_1(R(m,n))
    /*2*/  + P_2_1(R(m,n)))) / (TINY_Real + Power(fA(m,n),2)
    /*1*/  * Power(R(m,n),2));

Jacobian[17][1]=
	k_f * fAlp(m,n) * (-8 * exp(4 * fconf(m,n))
    /*1*/  * fj(m,n) + (8 * exp(2 * gconf(m,n)) * gA(m,n)
    /*2*/  * p(m,n) * P_1_1(R(m,n))) / (TINY_Real
    /*2*/  + Power(fA(m,n),2) * Power(R(m,n),2))) + (2
    /*1*/  * fAlp(m,n) * p(m,n) * (fL_r(m,n) - (2 * exp(4
    /*4*/  * fconf(m,n) - 4 * gconf(m,n))) / (pow2(gB(m,n))
    /*3*/  * pow2(r(m,n)) * pow2(R(m,n))))) / (TINY_Real
    /*1*/  + exp(2 * fconf(m,n)) * fA(m,n) * Lt(m,n));

Jacobian[17][3]=
	(8 * fAsig(m,n) * fAlp(m,n) * pow2(r(m,n)))
    /*0*/  / (TINY_Real + Power(fA(m,n),2));

Jacobian[17][6]=
	(-2 * k_f * exp(2 * gconf(m,n)) * fAlp(m,n)
    /*1*/  * p(m,n) * P_2_1(R(m,n))) / (TINY_Real
    /*1*/  + Power(fA(m,n),2) * Power(R(m,n),2));

Jacobian[17][7]=
	(-4 * k_f * exp(2 * gconf(m,n)) * fAlp(m,n)
    /*1*/  * gA(m,n) * p(m,n) * P_1_1(R(m,n))) / (TINY_Real
    /*1*/  + Power(fA(m,n),2) * gB(m,n) * Power(R(m,n),2));

Jacobian[17][10]=
	(12 * fA1(m,n) * fDAlp(m,n) - 6 * fBet_rr(m,n))
    /*0*/  / (TINY_Real + 3 * Power(fA(m,n),3)) + (4 * k_f
    /*1*/  * exp(2 * gconf(m,n)) * fAlp(m,n) * gA(m,n) * p(m,n)
    /*1*/  * P_2_1(R(m,n))) / (TINY_Real + Power(fA(m,n),3)
    /*1*/  * Power(R(m,n),2)) + fAlp(m,n) * ((-8 * (fAsig(m,n)
    /*3*/  * pow2(r(m,n)) * (fDA(m,n) + fDB(m,n) + 6
    /*4*/  * fDconf(m,n) + fsig(m,n) * r(m,n)) - ftrK_r(m,n)))
    /*1*/  / (TINY_Real + 3 * Power(fA(m,n),3)) + (p(m,n)
    /*2*/  * (fL_r(m,n) - (2 * exp(4 * fconf(m,n)))
    /*3*/  / (TINY_Real + exp(4 * gconf(m,n))
    /*4*/  * Power(gB(m,n),2) * Power(r(m,n),2)
    /*4*/  * Power(R(m,n),2)))) / (TINY_Real + exp(2
    /*3*/  * fconf(m,n)) * Power(fA(m,n),2) * Lt(m,n)));

Jacobian[17][11]=
	(4 * k_f * exp(2 * fconf(m,n)) * fAlp(m,n)
    /*1*/  * gA(m,n) * p(m,n) * P_1_1(R(m,n))) / (TINY_Real
    /*1*/  + Power(fA(m,n),2) * gB(m,n) * Power(R(m,n),3)) + (4
    /*1*/  * exp(6 * fconf(m,n)) * (gBet(m,n) - fBet_r(m,n)
    /*2*/  * r(m,n))) / (TINY_Real + exp(6 * gconf(m,n))
    /*1*/  * Power(gB(m,n),3) * Power(r(m,n),2)
    /*1*/  * Power(R(m,n),3)) - (4 * exp(6 * fconf(m,n))
    /*1*/  * gAlp(m,n) * p(m,n)) / (TINY_Real + exp(8
    /*2*/  * gconf(m,n)) * gA(m,n) * Power(gB(m,n),3) * Lt(m,n)
    /*1*/  * Power(r(m,n),2) * Power(R(m,n),3)) - (4 * exp(4
    /*2*/  * fconf(m,n)) * fAlp(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(6 * gconf(m,n)) * fA(m,n) * Power(gB(m,n),3)
    /*1*/  * Power(r(m,n),2) * Power(R(m,n),3) * Lt(m,n));

Jacobian[17][12]=
	(4 * fAsig(m,n) * fAlp(m,n) * pow2(r(m,n)))
    /*0*/  / (TINY_Real + 3 * Power(fA(m,n),2));

Jacobian[17][13]=
	(4 * fAsig(m,n) * fAlp(m,n) * pow2(r(m,n)))
    /*0*/  / (TINY_Real + 3 * Power(fA(m,n),2));

Jacobian[17][15]=
	(-2 * fDAlp(m,n)) / (TINY_Real + Power(fA(m,n),2));

Jacobian[17][17]=
	-fBet_r(m,n);

Jacobian[17][19]=
	(4 * fAsig(m,n) * fAlp(m,n) * pow3(r(m,n)))
    /*0*/  / (TINY_Real + 3 * Power(fA(m,n),2));

Jacobian[17][21]=
	(4 * fAlp(m,n) * pow2(r(m,n)) * (fDA(m,n)
    /*2*/  + fDB(m,n) + 6 * fDconf(m,n) + fsig(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + 3 * Power(fA(m,n),2));

Jacobian[17][25]=
	fL_r(m,n) - (2 * exp(4 * fconf(m,n))) / (TINY_Real
    /*1*/  + exp(4 * gconf(m,n)) * Power(gB(m,n),2)
    /*1*/  * Power(r(m,n),2) * Power(R(m,n),2));

Jacobian[18][0]=
	gAlp(m,n) * p(m,n) * ((-2 * (2 * gsig(m,n)
    /*3*/  + gsig_r(m,n) * r(m,n))) / (TINY_Real + exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * r(m,n) * Lt(m,n)) - (4
    /*2*/  * gA(m,n)) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * Power(gB(m,n),2) * Power(r(m,n),3) * Lt(m,n)));

Jacobian[18][6]=
	(4 * gA(m,n) * (gBet(m,n) - gBet_r(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + Power(gB(m,n),2) * Power(r(m,n),3))
    /*0*/  + gAlp(m,n) * ((4 * gAsig(m,n) * gA(m,n))
    /*1*/  / (TINY_Real + Power(gB(m,n),2)) + p(m,n) * (-4
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * Power(gB(m,n),2) * Lt(m,n) * Power(r(m,n),3))
    /*2*/  + (-2 * gsig(m,n) - gsig_r(m,n) * r(m,n))
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * Power(gA(m,n),2) * r(m,n) * Lt(m,n)) + 2
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * Power(gB(m,n),2) * Power(r(m,n),3) * Lt(m,n))));

Jacobian[18][7]=
	(4 * pow2(gA(m,n)) * (-gBet(m,n) + gBet_r(m,n)
    /*2*/  * r(m,n))) / (TINY_Real + Power(gB(m,n),3)
    /*1*/  * Power(r(m,n),3)) + gAlp(m,n) * ((-4 * gAsig(m,n)
    /*2*/  * pow2(gA(m,n))) / (TINY_Real + Power(gB(m,n),3))
    /*1*/  + p(m,n) * ((4 * gA(m,n)) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * Power(gB(m,n),3) * Lt(m,n)
    /*3*/  * Power(r(m,n),3)) - (4 * gA(m,n)) / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * Power(gB(m,n),3)
    /*3*/  * Power(r(m,n),3) * Lt(m,n))));

Jacobian[18][18]=
	2 * gBetr(m,n) + gAlp(m,n) * p(m,n) * (-2
    /*1*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*2*/  * Lt(m,n) * r(m,n)) + 2 / (TINY_Real + exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * r(m,n) * Lt(m,n)));

Jacobian[18][20]=
	(2 * gAlp(m,n) * pow2(gA(m,n))) / (TINY_Real
    /*1*/  + Power(gB(m,n),2));

Jacobian[18][25]=
	gsig_r(m,n) + (2 * gsig(m,n)) / (TINY_Real
    /*1*/  + r(m,n)) + (2 * pow2(gA(m,n))) / (TINY_Real
    /*1*/  + Power(gB(m,n),2) * Power(r(m,n),3));

Jacobian[19][1]=
	fAlp(m,n) * p(m,n) * ((2 * (2 * fsig(m,n)
    /*3*/  + fsig_r(m,n) * r(m,n))) / (TINY_Real + exp(2
    /*3*/  * fconf(m,n)) * fA(m,n) * r(m,n) * Lt(m,n)) + (4
    /*2*/  * exp(2 * fconf(m,n)) * fA(m,n)) / (TINY_Real
    /*2*/  + exp(4 * gconf(m,n)) * Power(gB(m,n),2)
    /*2*/  * Power(r(m,n),3) * Power(R(m,n),2) * Lt(m,n)));

Jacobian[19][10]=
	(4 * exp(4 * fconf(m,n)) * fA(m,n) * (gBet(m,n)
    /*2*/  - fBet_r(m,n) * r(m,n))) / (TINY_Real + exp(4
    /*2*/  * gconf(m,n)) * Power(gB(m,n),2) * Power(r(m,n),3)
    /*1*/  * Power(R(m,n),2)) - (4 * exp(4 * fconf(m,n))
    /*1*/  * fA(m,n) * gAlp(m,n) * p(m,n)) / (TINY_Real + exp(6
    /*2*/  * gconf(m,n)) * gA(m,n) * Power(gB(m,n),2)
    /*1*/  * Lt(m,n) * Power(r(m,n),3) * Power(R(m,n),2))
    /*0*/  + fAlp(m,n) * ((4 * exp(4 * fconf(m,n)) * fAsig(m,n)
    /*2*/  * fA(m,n)) / (TINY_Real + exp(4 * gconf(m,n))
    /*2*/  * Power(gB(m,n),2) * Power(R(m,n),2)) + p(m,n) * ((2
    /*3*/  * fsig(m,n) + fsig_r(m,n) * r(m,n)) / (TINY_Real
    /*3*/  + exp(2 * fconf(m,n)) * Power(fA(m,n),2) * r(m,n)
    /*3*/  * Lt(m,n)) - (2 * exp(2 * fconf(m,n))) / (TINY_Real
    /*3*/  + exp(4 * gconf(m,n)) * Power(gB(m,n),2)
    /*3*/  * Power(r(m,n),3) * Power(R(m,n),2) * Lt(m,n))));

Jacobian[19][11]=
	(4 * exp(6 * fconf(m,n)) * pow2(fA(m,n))
    /*1*/  * (-gBet(m,n) + fBet_r(m,n) * r(m,n))) / (TINY_Real
    /*1*/  + exp(6 * gconf(m,n)) * Power(gB(m,n),3)
    /*1*/  * Power(r(m,n),3) * Power(R(m,n),3)) + (4 * exp(6
    /*2*/  * fconf(m,n)) * gAlp(m,n) * p(m,n) * pow2(fA(m,n)))
    /*0*/  / (TINY_Real + exp(8 * gconf(m,n)) * gA(m,n)
    /*1*/  * Power(gB(m,n),3) * Lt(m,n) * Power(r(m,n),3)
    /*1*/  * Power(R(m,n),3)) + fAlp(m,n) * ((-4 * exp(6
    /*3*/  * fconf(m,n)) * fAsig(m,n) * pow2(fA(m,n)))
    /*1*/  / (TINY_Real + exp(6 * gconf(m,n))
    /*2*/  * Power(gB(m,n),3) * Power(R(m,n),3)) + (4 * exp(4
    /*3*/  * fconf(m,n)) * fA(m,n) * p(m,n)) / (TINY_Real
    /*2*/  + exp(6 * gconf(m,n)) * Power(gB(m,n),3)
    /*2*/  * Power(r(m,n),3) * Power(R(m,n),3) * Lt(m,n)));

Jacobian[19][19]=
	2 * gBetr(m,n) - (2 * gAlp(m,n) * p(m,n))
    /*0*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*1*/  * Lt(m,n) * r(m,n)) - (2 * fAlp(m,n) * p(m,n))
    /*0*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[19][21]=
	(2 * exp(4 * fconf(m,n)) * fAlp(m,n)
    /*1*/  * pow2(fA(m,n))) / (TINY_Real + exp(4 * gconf(m,n))
    /*1*/  * Power(gB(m,n),2) * Power(R(m,n),2));

Jacobian[19][25]=
	fsig_r(m,n) + (2 * fsig(m,n)) / (TINY_Real
    /*1*/  + r(m,n)) + (2 * exp(4 * fconf(m,n))
    /*1*/  * pow2(fA(m,n))) / (TINY_Real + exp(4 * gconf(m,n))
    /*1*/  * Power(gB(m,n),2) * Power(r(m,n),3)
    /*1*/  * Power(R(m,n),2));

Jacobian[20][0]=
	(-4 * (-(gDAlp_r(m,n) * r(m,n)) + gDAlp(m,n) * (1
    /*3*/  + gDA(m,n) * r(m,n) + gDB(m,n) * r(m,n) + 4
    /*3*/  * gDconf(m,n) * r(m,n)))) / (TINY_Real + exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2) * Power(r(m,n),3))
    /*0*/  + gAlp(m,n) * ((-2 * pow2(gB(m,n)) * (4
    /*3*/  * gsig_r(m,n) + gsig_rr(m,n) * r(m,n) + 2
    /*3*/  * pow2(gsig(m,n)) * r(m,n))) / (TINY_Real + exp(4
    /*3*/  * gconf(m,n)) * Power(gA(m,n),4) * r(m,n)) + (4
    /*2*/  * (gL(m,n) - gL_r(m,n) * r(m,n))) / (TINY_Real
    /*2*/  + exp(4 * gconf(m,n)) * Power(r(m,n),3)) + (2 * (-6
    /*3*/  * pow2(gDA(m,n)) * r(m,n) - 8 * pow2(gDconf(m,n))
    /*3*/  * r(m,n) - 4 * gDconf(m,n) * (1 + gDB(m,n)
    /*4*/  * r(m,n)) + gDA(m,n) * (8 * gDB(m,n) * r(m,n) - 4
    /*4*/  * gDconf(m,n) * r(m,n)) + r(m,n) * (4
    /*4*/  * gDconf_r(m,n) + gsig_gL_r(m,n) * pow2(gB(m,n))
    /*4*/  * pow2(r(m,n)) + 8 * gDB(m,n) * gsig(m,n) * r(m,n)
    /*4*/  - 2 * gsig_r(m,n) * r(m,n) + 2 * gL(m,n) * gsig(m,n)
    /*4*/  * pow2(gB(m,n)) * r(m,n)))) / (TINY_Real + exp(4
    /*3*/  * gconf(m,n)) * Power(gA(m,n),2) * Power(r(m,n),3))
    /*1*/  - (2 * p(m,n) * (2 * gAsig(m,n) + gAsig_r(m,n)
    /*3*/  * r(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * gA(m,n) * r(m,n) * Lt(m,n))) + k_g * (fAlp(m,n)
    /*1*/  * ((-2 * exp(2 * fconf(m,n)) * fA(m,n) * (2
    /*4*/  * P_1_2(R(m,n)) + b_2)) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * Power(r(m,n),2)) - (2 * (3
    /*4*/  * P_1_1(R(m,n)) - 2 * P_2_1(R(m,n)) + b_1))
    /*2*/  / (TINY_Real + Power(r(m,n),2) * Lt(m,n)))
    /*1*/  + gAlp(m,n) * ((-2 * (3 * P_1_0(R(m,n)) - 2
    /*4*/  * P_2_0(R(m,n)) + b_0)) / (TINY_Real
    /*3*/  + Power(r(m,n),2)) + (2 * exp(2 * fconf(m,n))
    /*3*/  * fA(m,n) * (2 * (-2 + pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) - 3 * (-1 + pow2(Lt(m,n)))
    /*4*/  * P_2_1(R(m,n)) - b_1)) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * Power(r(m,n),2)
    /*3*/  * Lt(m,n))));

Jacobian[20][1]=
	k_g * (fAlp(m,n) * ((2 * exp(2 * fconf(m,n))
    /*3*/  * fA(m,n) * (2 * P_1_2(R(m,n)) + b_2)) / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * gA(m,n) * Power(r(m,n),2))
    /*2*/  + (2 * R(m,n) * (b_2 + 2 * R(m,n) * b_3))
    /*2*/  / (TINY_Real + Power(r(m,n),2) * Lt(m,n)))
    /*1*/  + gAlp(m,n) * ((2 * R(m,n) * (b_1 + 2 * R(m,n)
    /*4*/  * b_2)) / (TINY_Real + Power(r(m,n),2)) - (2 * exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * (2 * (-2
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) - 3 * (-1
    /*5*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)) - b_1))
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * Power(r(m,n),2) * Lt(m,n))));

Jacobian[20][2]=
	(4 * gDAlp(m,n)) / (TINY_Real + exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2) * Power(r(m,n),2))
    /*0*/  + (2 * gAlp(m,n) * (1 + gDA(m,n) * r(m,n)
    /*2*/  + gDB(m,n) * r(m,n) + 4 * gDconf(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + exp(4 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2) * Power(r(m,n),3));

Jacobian[20][4]=
	gAsig(m,n) * gAlp(m,n);

Jacobian[20][6]=
	(-2 * (-(gDAlp_r(m,n) * r(m,n)) + gDAlp(m,n) * (1
    /*3*/  + gDA(m,n) * r(m,n) + gDB(m,n) * r(m,n) + 4
    /*3*/  * gDconf(m,n) * r(m,n)))) / (TINY_Real + exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),3) * Power(r(m,n),3))
    /*0*/  + gAlp(m,n) * ((-2 * pow2(gB(m,n)) * (4
    /*3*/  * gsig_r(m,n) + gsig_rr(m,n) * r(m,n) + 2
    /*3*/  * pow2(gsig(m,n)) * r(m,n))) / (TINY_Real + exp(4
    /*3*/  * gconf(m,n)) * Power(gA(m,n),5) * r(m,n)) + (-6
    /*2*/  * pow2(gDA(m,n)) * r(m,n) - 8 * pow2(gDconf(m,n))
    /*2*/  * r(m,n) - 4 * gDconf(m,n) * (1 + gDB(m,n) * r(m,n))
    /*2*/  + gDA(m,n) * (8 * gDB(m,n) * r(m,n) - 4
    /*3*/  * gDconf(m,n) * r(m,n)) + r(m,n) * (4
    /*3*/  * gDconf_r(m,n) + gsig_gL_r(m,n) * pow2(gB(m,n))
    /*3*/  * pow2(r(m,n)) + 8 * gDB(m,n) * gsig(m,n) * r(m,n)
    /*3*/  - 2 * gsig_r(m,n) * r(m,n) + 2 * gL(m,n) * gsig(m,n)
    /*3*/  * pow2(gB(m,n)) * r(m,n))) / (TINY_Real + exp(4
    /*3*/  * gconf(m,n)) * Power(gA(m,n),3) * Power(r(m,n),3))
    /*1*/  + (p(m,n) * (-2 * gAsig(m,n) - gAsig_r(m,n)
    /*3*/  * r(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2) * r(m,n) * Lt(m,n))) + k_g
    /*0*/  * (-((exp(2 * fconf(m,n)) * fAlp(m,n) * fA(m,n)
    /*3*/  * P_1_2(R(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * Power(gA(m,n),2) * Power(r(m,n),2))) - (exp(2
    /*3*/  * fconf(m,n)) * fA(m,n) * gAlp(m,n) * (P_1_1(R(m,n))
    /*3*/  + (-1 + pow2(Lt(m,n))) * P_2_1(R(m,n))))
    /*1*/  / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2) * Power(r(m,n),2) * Lt(m,n)));

Jacobian[20][7]=
	gAlp(m,n) * (-((gB(m,n) * (2 * gL(m,n) * gsig(m,n)
    /*4*/  + gsig_gL_r(m,n) * r(m,n))) / (TINY_Real + exp(4
    /*4*/  * gconf(m,n)) * Power(gA(m,n),2) * r(m,n)))
    /*1*/  + (gB(m,n) * (4 * gsig_r(m,n) + gsig_rr(m,n)
    /*3*/  * r(m,n) + 2 * pow2(gsig(m,n)) * r(m,n)))
    /*1*/  / (TINY_Real + exp(4 * gconf(m,n))
    /*2*/  * Power(gA(m,n),4) * r(m,n))) + k_g * (fAlp(m,n)
    /*1*/  * ((exp(2 * fconf(m,n)) * fA(m,n) * R(m,n) * b_3)
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * gB(m,n) * Power(r(m,n),2)) - (R(m,n) * (b_2 + 2
    /*4*/  * R(m,n) * b_3)) / (TINY_Real + gB(m,n)
    /*3*/  * Power(r(m,n),2) * Lt(m,n))) + gAlp(m,n)
    /*1*/  * (-((R(m,n) * (b_1 + 2 * R(m,n) * b_2))
    /*3*/  / (TINY_Real + gB(m,n) * Power(r(m,n),2))) - (exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * R(m,n) * (2 * (-1
    /*5*/  + pow2(Lt(m,n))) * P_1_2(R(m,n)) - b_2))
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * gB(m,n) * Power(r(m,n),2) * Lt(m,n))));

Jacobian[20][8]=
	gDAlp(m,n) / (TINY_Real + exp(4 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2) * Power(r(m,n),2)) + (2 * (3
    /*2*/  * gDA(m,n) - 2 * gDB(m,n) + gDconf(m,n))
    /*1*/  * gAlp(m,n)) / (TINY_Real + exp(4 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2) * Power(r(m,n),2));

Jacobian[20][9]=
	gDAlp(m,n) / (TINY_Real + exp(4 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2) * Power(r(m,n),2)) + (gAlp(m,n)
    /*1*/  * (-4 * gDA(m,n) + 2 * gDconf(m,n) - 4 * gsig(m,n)
    /*2*/  * r(m,n))) / (TINY_Real + exp(4 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2) * Power(r(m,n),2));

Jacobian[20][10]=
	k_g * ((exp(2 * fconf(m,n)) * fAlp(m,n)
    /*2*/  * P_1_2(R(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*2*/  * gA(m,n) * Power(r(m,n),2)) + (exp(2
    /*3*/  * fconf(m,n)) * gAlp(m,n) * (P_1_1(R(m,n)) + (-1
    /*4*/  + pow2(Lt(m,n))) * P_2_1(R(m,n)))) / (TINY_Real
    /*2*/  + exp(2 * gconf(m,n)) * gA(m,n) * Power(r(m,n),2)
    /*2*/  * Lt(m,n)));

Jacobian[20][11]=
	k_g * (fAlp(m,n) * (-((exp(4 * fconf(m,n))
    /*4*/  * fA(m,n) * b_3) / (TINY_Real + exp(4 * gconf(m,n))
    /*4*/  * gA(m,n) * gB(m,n) * Power(r(m,n),2))) - (exp(2
    /*4*/  * fconf(m,n)) * (2 * P_1_2(R(m,n)) + b_2))
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n)) * gB(m,n)
    /*3*/  * Power(r(m,n),2) * Lt(m,n))) + gAlp(m,n)
    /*1*/  * (-((exp(2 * fconf(m,n)) * (2 * P_1_1(R(m,n))
    /*5*/  + b_1)) / (TINY_Real + exp(2 * gconf(m,n)) * gB(m,n)
    /*4*/  * Power(r(m,n),2))) + (exp(4 * fconf(m,n))
    /*3*/  * fA(m,n) * (2 * (-1 + pow2(Lt(m,n)))
    /*4*/  * P_1_2(R(m,n)) - b_2)) / (TINY_Real + exp(4
    /*4*/  * gconf(m,n)) * gA(m,n) * gB(m,n) * Power(r(m,n),2)
    /*3*/  * Lt(m,n))));

Jacobian[20][16]=
	gAlp(m,n) * (-((pow2(gB(m,n)) * (2 * gsig(m,n)
    /*4*/  + gsig_r(m,n) * r(m,n))) / (TINY_Real + 2 * exp(4
    /*4*/  * gconf(m,n)) * Power(gA(m,n),2) * r(m,n))) - 1
    /*1*/  / (TINY_Real + exp(4 * gconf(m,n))
    /*2*/  * Power(r(m,n),3)));

Jacobian[20][18]=
	gAlp(m,n) * ((2 * gsig(m,n) * pow2(gB(m,n)))
    /*1*/  / (TINY_Real + exp(4 * gconf(m,n))
    /*2*/  * Power(gA(m,n),4)) + (-4 * gDB(m,n) - gL(m,n)
    /*2*/  * pow2(gB(m,n))) / (TINY_Real + exp(4 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2) * r(m,n)));

Jacobian[20][20]=
	2 * gBetr(m,n) + gAlp(m,n) * (gtrK(m,n) + p(m,n)
    /*1*/  * (-2 / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * Lt(m,n) * r(m,n)) + 2 / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * r(m,n) * Lt(m,n))));

Jacobian[20][25]=
	gAsig_r(m,n) + (2 * gAsig(m,n)) / (TINY_Real
    /*1*/  + r(m,n));

Jacobian[21][0]=
	k_f * (gAlp(m,n) * ((2 * exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * (P_2_0(R(m,n)) + b_0)) / (TINY_Real
    /*3*/  + exp(2 * fconf(m,n)) * fA(m,n) * Power(r(m,n),2)
    /*3*/  * Power(R(m,n),2)) - (2 * (P_1_1(R(m,n)) - b_1))
    /*2*/  / (TINY_Real + Power(r(m,n),2) * Power(R(m,n),2)
    /*3*/  * Lt(m,n))) + fAlp(m,n) * ((-2 * (P_1_2(R(m,n))
    /*4*/  - b_2)) / (TINY_Real + Power(r(m,n),2)
    /*3*/  * Power(R(m,n),2)) + (2 * exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * (2 * (-1 + pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) + pow2(Lt(m,n)) * P_2_1(R(m,n))
    /*4*/  + b_1)) / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*3*/  * Power(r(m,n),2) * Power(R(m,n),2) * Lt(m,n))));

Jacobian[21][1]=
	(-4 * (-(fDAlp_r(m,n) * r(m,n)) + fDAlp(m,n) * (1
    /*3*/  + fDA(m,n) * r(m,n) + fDB(m,n) * r(m,n) + 4
    /*3*/  * fDconf(m,n) * r(m,n)))) / (TINY_Real + exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),2) * Power(r(m,n),3))
    /*0*/  + fAlp(m,n) * ((-2 * exp(4 * gconf(m,n))
    /*2*/  * pow2(gB(m,n)) * pow2(R(m,n)) * (4 * fsig_r(m,n)
    /*3*/  + fsig_rr(m,n) * r(m,n) + 2 * pow2(fsig(m,n))
    /*3*/  * r(m,n))) / (TINY_Real + exp(8 * fconf(m,n))
    /*2*/  * Power(fA(m,n),4) * r(m,n)) + (4 * (fL(m,n)
    /*3*/  - fL_r(m,n) * r(m,n))) / (TINY_Real + exp(4
    /*3*/  * fconf(m,n)) * Power(r(m,n),3)) + (2 * (4 * exp(4
    /*4*/  * fconf(m,n)) * fDA(m,n) * (2 * fDB(m,n)
    /*4*/  - fDconf(m,n)) * r(m,n) - 6 * exp(4 * fconf(m,n))
    /*3*/  * pow2(fDA(m,n)) * r(m,n) - 8 * exp(4 * fconf(m,n))
    /*3*/  * pow2(fDconf(m,n)) * r(m,n) - 4 * exp(4
    /*4*/  * fconf(m,n)) * fDconf(m,n) * (1 + fDB(m,n)
    /*4*/  * r(m,n)) + r(m,n) * (4 * exp(4 * fconf(m,n))
    /*4*/  * fDconf_r(m,n) + exp(4 * gconf(m,n))
    /*4*/  * fsig_fL_r(m,n) * pow2(gB(m,n)) * pow2(r(m,n))
    /*4*/  * pow2(R(m,n)) + 8 * exp(4 * fconf(m,n)) * fDB(m,n)
    /*4*/  * fsig(m,n) * r(m,n) - 2 * exp(4 * fconf(m,n))
    /*4*/  * fsig_r(m,n) * r(m,n) + 2 * exp(4 * gconf(m,n))
    /*4*/  * fL(m,n) * fsig(m,n) * pow2(gB(m,n)) * pow2(R(m,n))
    /*4*/  * r(m,n)))) / (TINY_Real + exp(8 * fconf(m,n))
    /*2*/  * Power(fA(m,n),2) * Power(r(m,n),3)) + (2 * p(m,n)
    /*2*/  * (2 * fAsig(m,n) + fAsig_r(m,n) * r(m,n)))
    /*1*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*2*/  * r(m,n) * Lt(m,n))) + k_f * (gAlp(m,n) * ((-2
    /*3*/  * exp(2 * gconf(m,n)) * gA(m,n) * (P_2_0(R(m,n))
    /*4*/  + b_0)) / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n)
    /*3*/  * Power(r(m,n),2) * Power(R(m,n),2)) + (2
    /*3*/  * (P_1_1(R(m,n)) - b_1)) / (TINY_Real
    /*3*/  + Power(r(m,n),2) * Power(R(m,n),2) * Lt(m,n)))
    /*1*/  + fAlp(m,n) * ((2 * (P_1_2(R(m,n)) - b_2))
    /*2*/  / (TINY_Real + Power(r(m,n),2) * Power(R(m,n),2))
    /*2*/  - (2 * exp(2 * gconf(m,n)) * gA(m,n) * (2 * (-1
    /*5*/  + pow2(Lt(m,n))) * P_1_1(R(m,n)) + pow2(Lt(m,n))
    /*4*/  * P_2_1(R(m,n)) + b_1)) / (TINY_Real + exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * Power(r(m,n),2)
    /*3*/  * Power(R(m,n),2) * Lt(m,n))));

Jacobian[21][3]=
	(4 * fDAlp(m,n)) / (TINY_Real + exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),2) * Power(r(m,n),2))
    /*0*/  + (2 * fAlp(m,n) * (1 + fDA(m,n) * r(m,n)
    /*2*/  + fDB(m,n) * r(m,n) + 4 * fDconf(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + exp(4 * fconf(m,n))
    /*1*/  * Power(fA(m,n),2) * Power(r(m,n),3));

Jacobian[21][5]=
	fAsig(m,n) * fAlp(m,n);

Jacobian[21][6]=
	k_f * ((exp(2 * gconf(m,n)) * gAlp(m,n)
    /*2*/  * (-P_1_0(R(m,n)) + P_2_0(R(m,n)))) / (TINY_Real
    /*2*/  + exp(2 * fconf(m,n)) * fA(m,n) * Power(r(m,n),2)
    /*2*/  * Power(R(m,n),2)) + (exp(2 * gconf(m,n))
    /*2*/  * fAlp(m,n) * (-P_1_1(R(m,n)) + pow2(Lt(m,n))
    /*3*/  * P_2_1(R(m,n)))) / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * Power(r(m,n),2) * Power(R(m,n),2)
    /*2*/  * Lt(m,n)));

Jacobian[21][7]=
	k_f * (gAlp(m,n) * ((exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * (P_1_0(R(m,n)) + b_0)) / (TINY_Real + exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * gB(m,n) * Power(r(m,n),2)
    /*3*/  * Power(R(m,n),2)) + (-P_1_1(R(m,n)) + b_1)
    /*2*/  / (TINY_Real + gB(m,n) * Power(r(m,n),2)
    /*3*/  * Power(R(m,n),2) * Lt(m,n))) + fAlp(m,n)
    /*1*/  * ((-P_1_2(R(m,n)) + b_2) / (TINY_Real + gB(m,n)
    /*3*/  * Power(r(m,n),2) * Power(R(m,n),2)) + (exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * ((-1 + 2 * pow2(Lt(m,n)))
    /*4*/  * P_1_1(R(m,n)) + b_1)) / (TINY_Real + exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * gB(m,n) * Power(r(m,n),2)
    /*3*/  * Power(R(m,n),2) * Lt(m,n))));

Jacobian[21][10]=
	(-2 * (-(fDAlp_r(m,n) * r(m,n)) + fDAlp(m,n) * (1
    /*3*/  + fDA(m,n) * r(m,n) + fDB(m,n) * r(m,n) + 4
    /*3*/  * fDconf(m,n) * r(m,n)))) / (TINY_Real + exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),3) * Power(r(m,n),3))
    /*0*/  + fAlp(m,n) * ((-2 * exp(4 * gconf(m,n))
    /*2*/  * pow2(gB(m,n)) * pow2(R(m,n)) * (4 * fsig_r(m,n)
    /*3*/  + fsig_rr(m,n) * r(m,n) + 2 * pow2(fsig(m,n))
    /*3*/  * r(m,n))) / (TINY_Real + exp(8 * fconf(m,n))
    /*2*/  * Power(fA(m,n),5) * r(m,n)) + (4 * exp(4
    /*3*/  * fconf(m,n)) * fDA(m,n) * (2 * fDB(m,n)
    /*3*/  - fDconf(m,n)) * r(m,n) - 6 * exp(4 * fconf(m,n))
    /*2*/  * pow2(fDA(m,n)) * r(m,n) - 8 * exp(4 * fconf(m,n))
    /*2*/  * pow2(fDconf(m,n)) * r(m,n) - 4 * exp(4
    /*3*/  * fconf(m,n)) * fDconf(m,n) * (1 + fDB(m,n)
    /*3*/  * r(m,n)) + r(m,n) * (4 * exp(4 * fconf(m,n))
    /*3*/  * fDconf_r(m,n) + exp(4 * gconf(m,n))
    /*3*/  * fsig_fL_r(m,n) * pow2(gB(m,n)) * pow2(r(m,n))
    /*3*/  * pow2(R(m,n)) + 8 * exp(4 * fconf(m,n)) * fDB(m,n)
    /*3*/  * fsig(m,n) * r(m,n) - 2 * exp(4 * fconf(m,n))
    /*3*/  * fsig_r(m,n) * r(m,n) + 2 * exp(4 * gconf(m,n))
    /*3*/  * fL(m,n) * fsig(m,n) * pow2(gB(m,n)) * pow2(R(m,n))
    /*3*/  * r(m,n))) / (TINY_Real + exp(8 * fconf(m,n))
    /*2*/  * Power(fA(m,n),3) * Power(r(m,n),3)) + (p(m,n) * (2
    /*3*/  * fAsig(m,n) + fAsig_r(m,n) * r(m,n)))
    /*1*/  / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * Power(fA(m,n),2) * r(m,n) * Lt(m,n))) + k_f
    /*0*/  * (-((exp(2 * gconf(m,n)) * gAlp(m,n) * gA(m,n)
    /*3*/  * P_1_1(R(m,n))) / (TINY_Real + exp(2 * fconf(m,n))
    /*3*/  * Power(fA(m,n),2) * Power(r(m,n),2) * R(m,n)))
    /*1*/  + (exp(2 * gconf(m,n)) * fAlp(m,n) * gA(m,n)
    /*2*/  * (P_1_1(R(m,n)) - pow2(Lt(m,n)) * P_2_1(R(m,n))))
    /*1*/  / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * Power(fA(m,n),2) * Power(r(m,n),2)
    /*2*/  * Power(R(m,n),2) * Lt(m,n)));

Jacobian[21][11]=
	fAlp(m,n) * (-((exp(2 * gconf(m,n)) * gB(m,n) * (2
    /*4*/  * fL(m,n) * fsig(m,n) + fsig_fL_r(m,n) * r(m,n))
    /*3*/  * R(m,n)) / (TINY_Real + exp(6 * fconf(m,n))
    /*3*/  * Power(fA(m,n),2) * r(m,n))) + (exp(2 * gconf(m,n))
    /*2*/  * gB(m,n) * (4 * fsig_r(m,n) + fsig_rr(m,n)
    /*3*/  * r(m,n) + 2 * pow2(fsig(m,n)) * r(m,n)) * R(m,n))
    /*1*/  / (TINY_Real + exp(6 * fconf(m,n))
    /*2*/  * Power(fA(m,n),4) * r(m,n))) + k_f * (gAlp(m,n)
    /*1*/  * (-((gA(m,n) * (P_1_0(R(m,n)) + b_0)) / (TINY_Real
    /*4*/  + fA(m,n) * gB(m,n) * Power(r(m,n),2)
    /*4*/  * Power(R(m,n),3))) + (exp(2 * fconf(m,n))
    /*3*/  * (P_1_1(R(m,n)) - b_1)) / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gB(m,n) * Power(r(m,n),2)
    /*3*/  * Power(R(m,n),3) * Lt(m,n))) + fAlp(m,n) * ((exp(2
    /*4*/  * fconf(m,n)) * (P_1_2(R(m,n)) - b_2))
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n)) * gB(m,n)
    /*3*/  * Power(r(m,n),2) * Power(R(m,n),3)) + (gA(m,n)
    /*3*/  * ((1 - 2 * pow2(Lt(m,n))) * P_1_1(R(m,n)) - b_1))
    /*2*/  / (TINY_Real + fA(m,n) * gB(m,n) * Power(r(m,n),2)
    /*3*/  * Power(R(m,n),3) * Lt(m,n))));

Jacobian[21][12]=
	fDAlp(m,n) / (TINY_Real + exp(4 * fconf(m,n))
    /*1*/  * Power(fA(m,n),2) * Power(r(m,n),2)) + (2 * (3
    /*2*/  * fDA(m,n) - 2 * fDB(m,n) + fDconf(m,n))
    /*1*/  * fAlp(m,n)) / (TINY_Real + exp(4 * fconf(m,n))
    /*1*/  * Power(fA(m,n),2) * Power(r(m,n),2));

Jacobian[21][13]=
	fDAlp(m,n) / (TINY_Real + exp(4 * fconf(m,n))
    /*1*/  * Power(fA(m,n),2) * Power(r(m,n),2)) + (fAlp(m,n)
    /*1*/  * (-4 * fDA(m,n) + 2 * fDconf(m,n) - 4 * fsig(m,n)
    /*2*/  * r(m,n))) / (TINY_Real + exp(4 * fconf(m,n))
    /*1*/  * Power(fA(m,n),2) * Power(r(m,n),2));

Jacobian[21][17]=
	fAlp(m,n) * (-((exp(4 * gconf(m,n))
    /*3*/  * pow2(gB(m,n)) * pow2(R(m,n)) * (2 * fsig(m,n)
    /*4*/  + fsig_r(m,n) * r(m,n))) / (TINY_Real + 2 * exp(8
    /*4*/  * fconf(m,n)) * Power(fA(m,n),2) * r(m,n))) - 1
    /*1*/  / (TINY_Real + exp(4 * fconf(m,n))
    /*2*/  * Power(r(m,n),3)));

Jacobian[21][19]=
	fAlp(m,n) * ((2 * exp(4 * gconf(m,n)) * fsig(m,n)
    /*2*/  * pow2(gB(m,n)) * pow2(R(m,n))) / (TINY_Real
    /*2*/  + exp(8 * fconf(m,n)) * Power(fA(m,n),4)) + ((-4
    /*3*/  * fDB(m,n)) / (TINY_Real + exp(4 * fconf(m,n)))
    /*2*/  - (exp(4 * gconf(m,n)) * fL(m,n) * pow2(gB(m,n))
    /*3*/  * pow2(R(m,n))) / (TINY_Real + exp(8 * fconf(m,n))))
    /*1*/  / (TINY_Real + Power(fA(m,n),2) * r(m,n)));

Jacobian[21][21]=
	2 * gBetr(m,n) - (2 * gAlp(m,n) * p(m,n))
    /*0*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*1*/  * Lt(m,n) * r(m,n)) + fAlp(m,n) * (ftrK(m,n) - (2
    /*2*/  * p(m,n)) / (TINY_Real + exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * r(m,n) * Lt(m,n)));

Jacobian[21][25]=
	fAsig_r(m,n) + (2 * fAsig(m,n)) / (TINY_Real
    /*1*/  + r(m,n));

Jacobian[22][0]=
	(-2 * gAlp(m,n) * p(m,n) * (2 * pfD(m,n)
    /*2*/  + pfD_r(m,n) * r(m,n))) / (TINY_Real + exp(2
    /*2*/  * gconf(m,n)) * gA(m,n) * r(m,n) * Lt(m,n));

Jacobian[22][6]=
	(gAlp(m,n) * p(m,n) * (-2 * pfD(m,n) - pfD_r(m,n)
    /*2*/  * r(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2) * r(m,n) * Lt(m,n));

Jacobian[22][22]=
	gBet_r(m,n) + 2 * gBetr(m,n) - gDAlp(m,n)
    /*0*/  * pfv(m,n) + gAlp(m,n) * (-pfv_r(m,n) - (2
    /*2*/  * pfv(m,n)) / (TINY_Real + r(m,n)) + p(m,n) * (-2
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * Lt(m,n) * r(m,n)) + 2 / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * r(m,n) * Lt(m,n))));

Jacobian[22][25]=
	pfD_r(m,n) + (2 * pfD(m,n)) / (TINY_Real + r(m,n));

Jacobian[23][0]=
	(-2 * gAlp(m,n) * p(m,n) * (2 * pfS(m,n)
    /*2*/  + pfS_r(m,n) * r(m,n))) / (TINY_Real + exp(2
    /*2*/  * gconf(m,n)) * gA(m,n) * r(m,n) * Lt(m,n));

Jacobian[23][2]=
	2 * gAlp(m,n) * pfS(m,n) * pfv(m,n);

Jacobian[23][6]=
	(gAlp(m,n) * p(m,n) * (-2 * pfS(m,n) - pfS_r(m,n)
    /*2*/  * r(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2) * r(m,n) * Lt(m,n));

Jacobian[23][8]=
	gAlp(m,n) * pfS(m,n) * pfv(m,n);

Jacobian[23][22]=
	-gDAlp(m,n);

Jacobian[23][23]=
	2 * gBet_r(m,n) + 2 * gBetr(m,n) - gDAlp(m,n)
    /*0*/  * pfv(m,n) + gAlp(m,n) * ((-(pfv_r(m,n) * r(m,n))
    /*2*/  + pfv(m,n) * (-2 + gDA(m,n) * r(m,n) + 2
    /*3*/  * gDconf(m,n) * r(m,n))) / (TINY_Real + r(m,n))
    /*1*/  + p(m,n) * (-2 / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * Lt(m,n) * r(m,n)) + 2 / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * gA(m,n) * r(m,n)
    /*3*/  * Lt(m,n))));

Jacobian[23][24]=
	-gDAlp(m,n);

Jacobian[23][25]=
	pfS_r(m,n) + (2 * pfS(m,n)) / (TINY_Real + r(m,n));

Jacobian[24][0]=
	(4 * gDAlp(m,n) * pfS(m,n)) / (TINY_Real + exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2)) - (2 * gAlp(m,n)
    /*1*/  * p(m,n) * (2 * pftau(m,n) + pftau_r(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*1*/  * r(m,n) * Lt(m,n));

Jacobian[24][6]=
	(2 * gDAlp(m,n) * pfS(m,n)) / (TINY_Real + exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),3)) + (gAlp(m,n)
    /*1*/  * p(m,n) * (-2 * pftau(m,n) - pftau_r(m,n)
    /*2*/  * r(m,n))) / (TINY_Real + exp(2 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2) * r(m,n) * Lt(m,n));

Jacobian[24][23]=
	gK1(m,n) * gAlp(m,n) * pfv(m,n) - gDAlp(m,n)
    /*0*/  / (TINY_Real + exp(4 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2));

Jacobian[24][24]=
	gBet_r(m,n) + 2 * gBetr(m,n) - gDAlp(m,n)
    /*0*/  * pfv(m,n) + gAlp(m,n) * (-pfv_r(m,n) - (2
    /*2*/  * pfv(m,n)) / (TINY_Real + r(m,n)) + p(m,n) * (-2
    /*2*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * Lt(m,n) * r(m,n)) + 2 / (TINY_Real + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * r(m,n) * Lt(m,n))));

Jacobian[24][25]=
	pftau_r(m,n) + (2 * pftau(m,n)) / (TINY_Real
    /*1*/  + r(m,n));

Jacobian[25][25]=
	gBet_r(m,n) - (gAlp_r(m,n) * p(m,n)) / (TINY_Real
    /*1*/  + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n))
    /*0*/  + gAlp(m,n) * ((((2 * gconf_r(m,n)) / (TINY_Real
    /*4*/  + exp(2 * gconf(m,n)) * gA(m,n)) + gA_r(m,n)
    /*3*/  / (TINY_Real + exp(2 * gconf(m,n))
    /*4*/  * Power(gA(m,n),2))) / (TINY_Real + Lt(m,n))
    /*2*/  + Lt_r(m,n) / (TINY_Real + exp(2 * gconf(m,n))
    /*3*/  * gA(m,n) * Power(Lt(m,n),2))) * p(m,n) - p_r(m,n)
    /*1*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)
    /*2*/  * Lt(m,n)));

Jacobian[25][26]=
	0.75;

Jacobian[26][0]=
	k_g * pow3(gAlp(m,n)) * (-8 * exp(4 * gconf(m,n))
    /*1*/  * gj(m,n) - (8 * exp(2 * fconf(m,n)) * fA(m,n)
    /*2*/  * p(m,n) * P_1_2(R(m,n)) * R(m,n)) / (TINY_Real
    /*2*/  + Power(gA(m,n),2))) + p(m,n) * pow3(gAlp(m,n))
    /*0*/  * ((-2 * (-2 + gDA_r(m,n) * pow2(r(m,n)) + 2
    /*3*/  * gDB_r(m,n) * pow2(r(m,n)))) / (TINY_Real + 3
    /*2*/  * exp(2 * gconf(m,n)) * Power(gA(m,n),3)
    /*2*/  * Power(r(m,n),2) * Lt(m,n)) - (2 * (-6 + 3
    /*3*/  * gL_r(m,n) * pow2(gB(m,n)) * pow2(r(m,n)) + 2
    /*3*/  * gL(m,n) * pow2(gB(m,n)) * r(m,n) * (2 + gDA(m,n)
    /*4*/  * r(m,n) + 2 * gDB(m,n) * r(m,n)))) / (TINY_Real + 3
    /*2*/  * exp(2 * gconf(m,n)) * gA(m,n) * Power(gB(m,n),2)
    /*2*/  * Power(r(m,n),2) * Lt(m,n)));

Jacobian[26][1]=
	(-4 * k_g * exp(2 * fconf(m,n)) * fA(m,n) * p(m,n)
    /*1*/  * pow3(gAlp(m,n)) * (2 * P_1_1(R(m,n)) - 3
    /*2*/  * P_2_1(R(m,n)))) / (TINY_Real + Power(gA(m,n),2));

Jacobian[26][2]=
	(8 * gAsig(m,n) * pow2(r(m,n)) * pow3(gAlp(m,n)))
    /*0*/  / (TINY_Real + Power(gA(m,n),2));

Jacobian[26][6]=
	(-4 * k_g * exp(2 * fconf(m,n)) * fA(m,n) * p(m,n)
    /*1*/  * pow3(gAlp(m,n)) * P_2_1(R(m,n))) / (TINY_Real
    /*1*/  + Power(gA(m,n),3)) - (2 * pow2(gAlp(m,n)) * (-2
    /*2*/  * gBet(m,n) + gBet_gDA_r(m,n) * pow2(r(m,n)) + 2
    /*2*/  * gBet_gDB_r(m,n) * pow2(r(m,n)) - 6 * gA1(m,n)
    /*2*/  * gDAlp(m,n) * pow2(r(m,n)) + gDA_r(m,n) * gBet(m,n)
    /*2*/  * pow2(r(m,n)) + 2 * gDB_r(m,n) * gBet(m,n)
    /*2*/  * pow2(r(m,n)) + 4 * gBet_rr(m,n) * pow2(r(m,n)) + 2
    /*2*/  * gBet_r(m,n) * r(m,n))) / (TINY_Real + 3
    /*1*/  * Power(gA(m,n),3) * Power(r(m,n),2))
    /*0*/  + pow3(gAlp(m,n)) * ((-8 * (gAsig(m,n)
    /*3*/  * pow2(r(m,n)) * (gDA(m,n) + gDB(m,n) + 6
    /*4*/  * gDconf(m,n) + gsig(m,n) * r(m,n)) - gtrK_r(m,n)))
    /*1*/  / (TINY_Real + 3 * Power(gA(m,n),3)) + p(m,n)
    /*1*/  * ((2 * (-2 + gDA_r(m,n) * pow2(r(m,n)) + 2
    /*4*/  * gDB_r(m,n) * pow2(r(m,n)))) / (TINY_Real + 3
    /*3*/  * exp(2 * gconf(m,n)) * Power(gA(m,n),4) * Lt(m,n)
    /*3*/  * Power(r(m,n),2)) + (2 - gDA_r(m,n) * pow2(r(m,n))
    /*3*/  - 2 * gDB_r(m,n) * pow2(r(m,n))) / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * Power(gA(m,n),4)
    /*3*/  * Power(r(m,n),2) * Lt(m,n)) + (6 - 3 * gL_r(m,n)
    /*3*/  * pow2(gB(m,n)) * pow2(r(m,n)) - 2 * gL(m,n)
    /*3*/  * pow2(gB(m,n)) * r(m,n) * (2 + gDA(m,n) * r(m,n)
    /*4*/  + 2 * gDB(m,n) * r(m,n))) / (TINY_Real + 3 * exp(2
    /*4*/  * gconf(m,n)) * Power(gA(m,n),2) * Power(gB(m,n),2)
    /*3*/  * Power(r(m,n),2) * Lt(m,n))));

Jacobian[26][7]=
	(4 * pow2(gAlp(m,n)) * (gBet(m,n) - gBet_r(m,n)
    /*2*/  * r(m,n))) / (TINY_Real + Power(gB(m,n),3)
    /*1*/  * Power(r(m,n),2)) - (4 * k_g * exp(2 * fconf(m,n))
    /*1*/  * fA(m,n) * p(m,n) * pow3(gAlp(m,n))
    /*1*/  * P_1_2(R(m,n)) * R(m,n)) / (TINY_Real
    /*1*/  + Power(gA(m,n),2) * gB(m,n)) + p(m,n)
    /*0*/  * pow3(gAlp(m,n)) * (-4 / (TINY_Real + exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * Power(gB(m,n),3) * Lt(m,n)
    /*2*/  * Power(r(m,n),2)) + 4 / (TINY_Real + exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * Power(gB(m,n),3)
    /*2*/  * Power(r(m,n),2) * Lt(m,n)));

Jacobian[26][8]=
	((2 * gBet(m,n) * gL(m,n)) / (3 + TINY_Real)
    /*1*/  + gBet_r(m,n) / (TINY_Real + 3 * Power(gA(m,n),2)))
    /*0*/  * pow2(gAlp(m,n)) + pow3(gAlp(m,n)) * ((4
    /*2*/  * gAsig(m,n) * pow2(r(m,n))) / (TINY_Real + 3
    /*2*/  * Power(gA(m,n),2)) + p(m,n) * ((-2 * gL(m,n))
    /*2*/  / (TINY_Real + 3 * exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * Lt(m,n)) + (2 * gL(m,n)) / (TINY_Real + 3 * exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * Lt(m,n))));

Jacobian[26][9]=
	((4 * gBet(m,n) * gL(m,n)) / (3 + TINY_Real) + (2
    /*2*/  * gBet_r(m,n)) / (TINY_Real + 3
    /*2*/  * Power(gA(m,n),2))) * pow2(gAlp(m,n))
    /*0*/  + pow3(gAlp(m,n)) * ((4 * gAsig(m,n) * pow2(r(m,n)))
    /*1*/  / (TINY_Real + 3 * Power(gA(m,n),2)) + p(m,n)
    /*1*/  * ((-4 * gL(m,n)) / (TINY_Real + 3 * exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * Lt(m,n)) + (4 * gL(m,n))
    /*2*/  / (TINY_Real + 3 * exp(2 * gconf(m,n)) * gA(m,n)
    /*3*/  * Lt(m,n))));

Jacobian[26][10]=
	(2 * k_g * exp(2 * fconf(m,n)) * p(m,n)
    /*1*/  * pow3(gAlp(m,n)) * P_2_1(R(m,n))) / (TINY_Real
    /*1*/  + Power(gA(m,n),2));

Jacobian[26][11]=
	(4 * k_g * exp(4 * fconf(m,n)) * fA(m,n) * p(m,n)
    /*1*/  * pow3(gAlp(m,n)) * P_1_2(R(m,n))) / (TINY_Real
    /*1*/  + exp(2 * gconf(m,n)) * Power(gA(m,n),2) * gB(m,n));

Jacobian[26][14]=
	(-4 * gDAlp(m,n) * pow2(gAlp(m,n))) / (TINY_Real
    /*1*/  + 3 * Power(gA(m,n),2));

Jacobian[26][16]=
	((2 * gDA(m,n) * gBet(m,n) + 4 * gDB(m,n)
    /*2*/  * gBet(m,n) - gBet_r(m,n) + 4 * gBetr(m,n))
    /*1*/  * pow2(gAlp(m,n))) / (3 + TINY_Real) + p(m,n)
    /*0*/  * pow3(gAlp(m,n)) * ((-2 * (2 + gDA(m,n) * r(m,n)
    /*3*/  + 2 * gDB(m,n) * r(m,n))) / (TINY_Real + 3 * exp(2
    /*3*/  * gconf(m,n)) * gA(m,n) * Lt(m,n) * r(m,n)) + (2
    /*2*/  * (2 + gDA(m,n) * r(m,n) + 2 * gDB(m,n) * r(m,n)))
    /*1*/  / (TINY_Real + 3 * exp(2 * gconf(m,n)) * gA(m,n)
    /*2*/  * r(m,n) * Lt(m,n)));

Jacobian[26][18]=
	(4 * gAsig(m,n) * pow3(gAlp(m,n)) * pow3(r(m,n)))
    /*0*/  / (TINY_Real + 3 * Power(gA(m,n),2));

Jacobian[26][20]=
	(4 * pow2(r(m,n)) * pow3(gAlp(m,n)) * (gDA(m,n)
    /*2*/  + gDB(m,n) + 6 * gDconf(m,n) + gsig(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + 3 * Power(gA(m,n),2));

Jacobian[26][25]=
	Bq_r(m,n) + pow2(gAlp(m,n)) * ((-6 + 2 * gL(m,n)
    /*2*/  * pow2(gB(m,n)) * r(m,n) * (2 + gDA(m,n) * r(m,n)
    /*3*/  + 2 * gDB(m,n) * r(m,n))) / (TINY_Real + 3
    /*2*/  * Power(gB(m,n),2) * Power(r(m,n),2)) + (gDA_r(m,n)
    /*2*/  + 2 * gDB_r(m,n) - 2 / (TINY_Real
    /*3*/  + Power(r(m,n),2))) / (TINY_Real + 3
    /*2*/  * Power(gA(m,n),2)));

Jacobian[26][26]=
	-eta;

