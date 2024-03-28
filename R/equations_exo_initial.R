# Title: Model Zezza

# Purpose : This is the model developed in the paper called
# U.S. growth, the housing market, and the distribution of
# income.
#
# Equations that appear are in the paper are includes with the
# same number they have in the paper. They are marked as such [xx].
# Equations that do not appear in the paper but are present in
# Gennaros code are marked as such [xxCx]. They belong to
# a equation that is in the paper they have the same number and
# only a C as a suffix. If they only appear in the code and are
# not helper equations they are marked as such [00xxC].

# Author: Alexander Toplitsch
# Contact details: alexander.toplitsch@s.wu.ac.at

# Date script created: Fr 02.12.2023 21:03 -------------
# Date script last modified:

# Check necessary packages for installation and call them

#############################################################
##########################Model##############################
#############################################################

# Set up the model equations
model_eqs <- sfcr_set(
  #----------------------------------#
  # Capitalists ("rich households")
  #----------------------------------#
  # [1] - Disposable income - Sum of wage income, distributed profits from firms and banks,
  # interest income from deposits and treasuries, and rents, net of direct taxes
  # paid to the government.
  Yc ~ wc * Nc + Rents + FD + rm[-1] * Mc[-1] + FB + rb[-1] * Bh[-1] - Tdc,
  # [2] Saving - augments the stock of wealth.
  Shc ~ Yc - Cc,
  # [3] Consumption inflator.
  Cc ~ cc * p,
  # [4] Consumption - depend on expected disposable income, the open stock of wealth  and
  # expected capital gains, all measured in real terms.
  cc ~ alpha_1c * yc_e + alpha_2c * vc[-1] + alpha_3c *
    (cge_e + cghc_e - p_e * vc[-1]/(1 + p_e)),
  # [5] Stock of Wealth - includes saving and capitals gains which are given from changes in
  # market price for equities [7] and homes [8].
  Vc ~ Vc[-1] + Shc + CGE + CGHc,
  # [5C1] Wealth deflator.
  vc ~ Vc / p,
  # [6] Disposable income deflator.
  yc ~ Yc / p,
  # [7] Changes in market price for equities -  Capital gains.
  CGE ~ (pe - pe[-1]) * E[-1],
  # [8] Changes in market price for homes - Capital gains.
  CGHc ~ (ph - ph[-1]) * Hc[-1],
  #-----Portfolio Choice-----#
  # [9] Cash - depend on current consumption.
  HPc ~ eta * Cc,
  # [10] Bank deposits.
  Mc ~ Vc - HPc - Bh - E * pe - ph * Hc,
  # [11] Bonds.
  Bh ~ Bh_port * (Vc_e - HPc),
  # [11C1] Bonds portfolio matrix.
  Bh_port ~ (lambda_10 - lambda_11 * rrm  - lambda_12 *
               rre_e - lambda_13 * (Yc_e/Vc_e) + lambda_14 *
               rrb - lambda_15 * rrh_e),
  # [12] Price of Equities.
  pe ~ (pe_port * (Vc_e - HPc) - xi * (I-FU)) / E[-1],
  # [12C1] Equities portfolio matrix.
  pe_port ~ (lambda_20 - lambda_21 * rrm + lambda_22 *
               rre_e - lambda_23 * (Yc_e/Vc_e) - lambda_24 *
               rrb - lambda_25 * rrh_e),
  # [13] Homes.
  Hc ~ Hc_port * (Vc_e - HPc) / ph,
  # [13C1] Homes portfolio matrix.
  Hc_port ~ (lambda_30 - lambda_31 * rrm - lambda_32 *
               rre_e - lambda_33 * (Yc_e/Vc_e) - lambda_34 *
               rrb + lambda_35 * rrh_e),
  # [14] Return on equities - Given by distributed profits and expected capital gains.
  re ~ (FD + CGE_e)/(pe[-1] * E[-1]),
  # [15] Return on Housing - Given by rents and expected capital gains.
  rh ~ (Rents + CGHc_e)/(ph[-1] * Hc[-1]),
  #----------------------------------#
  # Workers ("other households")
  #----------------------------------#
  # [16] Disposable income - Sum of wage income, rent to pay,
  # interest payments from bank deposits, interest payments for mortgages,
  # and rents, net of direct taxes paid to the government.
  Yo ~ wo * No - Rents + rm[-1] * Mo[-1] - rmo[-1] * MOo[-1] - Tdo,
  # [17] Saving
  Sho ~ Yo - Co,
  # [18] Consumption inflator.
  Co ~ co * p,
  # [19] Consumption - depends on expected real disposable income, past
  # real wealth, and expected real capital gains on homes minus past
  # wealth normalized by expected inflation.
  co ~ alpha_0o + alpha_1o * yo_e + alpha_2o * vo[-1] +
    alpha_3o * (cgho_e - p_e * vo[-1]/(1 + p_e)) - alpha_4o * morp * MOo[-1]/Yo[-1],
  # + im * No * (cc[-1]/Nc - co[-1]/No) Imitation parameter removed.
  # [20] Wealth - Past wealth plus saving plus capital gains from homes.
  Vo ~ Vo[-1] + Sho + CGHo,
  # [20C1] Wealth deflator.
  vo ~ Vo / p,
  # [21] Disposable income deflator.
  yo ~ Yo/p,
  # [22] Capital gains from homes.
  CGHo ~ (ph - ph[-1]) * Ho[-1],
  # [MEMO] Real income of households combined.
  yc_yo ~ yc + yo,
  # [MEMO] Nominal income of households combined.
  Yc_Yo ~ yc_yo * p,
  # [MEMO] Nominal consumption of households combined.
  C ~ Cc + Co,
  # [MEMO] Nominal income of households combined - check.
  Yc_Yo0 ~ (wo * No + wc * Nc) + FD + rm[-1] * Mo[-1] + FB + rb[-1] * Bh[-1] - Td,
  # [MEMO] Expected real income of households combined.
  yc_yo_e ~ yc_yo[-1] * (1 + y_e),
  #-----Portfolio Choice-----#
  # [24] Cash - depends on current consumption.
  HPo ~ eta * Co,
  # [25] Bank deposits - residual.
  Mo ~ Vo - HPo - ph * Ho + MOo,
  # [26] Demand for Homes - depends on population growth,
  # expected real income and lagged debt repayment ratios.
  Ho_d ~ Ho[-1] * ((No/No[-1])-1) + mu_1 * ((yo/yo[-1])-1) -
    mu_2 * delta_debt_rep,
  # [26C1] Number of Homes.
  Ho ~ Ho[-1] + Ho_d,
  # [0026C1] Debt repayment ratio.
  debt_rep ~ (rmo[-1] + morp[-1]) * MOo_1[-1] / Yo[-1],
  # [0026C2] Change in debt repayment.
  delta_debt_rep ~ debt_rep - debt_rep[-1],
  # [0026C3] Mo lag1
  MOo_1 ~ MOo[-1],
  # [27] Change Mortgages.
  MOo ~ MOo[-1] * (1 - morp) +
    MOo_cond * (ph * (Ho_d) - Sho),
  # [0027C1] Condition change Mortgages.
  MOo_cond ~ if(ph * Ho_d - Sho > cond) {1} else {0},
  # [28] Share of rented homes owned by capitalist.
  Rents ~ rent * Hc[-1],
  # [29] Rent increases
  rent ~ rent[-1] * (1 + y_e),
  #----------------------------------#
  # Nonfinancial firms
  #----------------------------------#
  # [30] Investment decision - Growth of k depends on actual profits, Tobin's q,
  # borrowing costs from banks and utilization rate. This the the growth rate.
  kgr ~ (iota_0 + iota_1 * FU[-1]/K_1[-1] -
           iota_2 * rll[-1] * lev[-1] +
           iota_3 * q[-1] +
           iota_4 * (ut[-1] - unorm)),
  # [30C1] Leverage.
  lev ~ L/K,
  # [30C2] Tobin's Q.
  q ~ pe * E/K,
  # [30C3] Real Capital.
  k ~ (1 + kgr) * k[-1],
  # [30C4] Capital in non real terms.
  K ~ k * p,
  # [30C5] Real Investment,
  i ~ k - k[-1],
  # [30C6] Nominal Investment.
  I ~ i * p,
  # [30C7] Real rate on loans.
  rll ~ (1 + rl) / (1 + p_e)-1,
  # [0030C1] Capital from one periods ago.
  K_1 ~ K[-1],
  # [31] Prices - are set with a mark-up on wages.
  p ~ (1 + rho) * w / (prod * (1 - tau)),
  # [32] Total Profits - are determined relative to the wage bill.
  FT ~ rho * WB,
  # [33] Distributed Profits - fixed share net of taxes and interest payments
  # is distributed to capitalists.
  FD ~ (1 - beta) * (FT - rl[-1] * L[-1] - TF),
  # [34] Mark-up - depends on relative strength of workers and capitalists.
  rho ~ ((rho_1 * (prodg - wo_g) + 1)*(1 + rho[-1])) - 1,
  # [35] Utilization rate - ratio of real sales to "normal" sales, which in
  # turn are in a fixed ratio (lambda) to the stock of real capital.
  ut ~ s/sfc,
  # [35C1] Normal Sales - fixed ratio of real capital.
  sfc ~ (lambda * k[-1]),
  # [36] Retained Profits.
  FU ~ FT - rl[-1] * L[-1] - FD - TF,
  # [37] New equities issued.
  E ~ E[-1] + xi * (I - FU) / pe,
  # [38] Loans - rewritten to display capital.
  L ~ L[-1] + I - FU - pe * (E - E[-1]),
  # # [Memo] Firms wealth.
  Vf ~ p * K + HU * ph - L - pe * E,
  #----------------------------------#
  # Banks
  #----------------------------------#
  # [39] Demand for Bonds.
  Bb ~ chi_1 * M,
  # [40] Reserve requirement.
  HPb ~ chi_2 * M,
  # [39C1] Deposits.
  M ~ Mo + Mc,
  # [41] Advancements from Central Bank - if internal funds are
  # not sufficient to cover for demand for loans the banks gets
  # advances from the central bank.
  A ~ L + HPb + Bb + MOo - M,
  # [Check].
  A0 ~ HP - Bcb,
  # [42] Interest rates on loans.
  rl ~ ra + spread_1,
  # [43] Interest rates on mortgages.
  rmo ~ ra + spread_2,
  # [44] Interest rates on deposits.
  rm ~ ra + spread_3,
  # [45] Banks profits, all distributed.
  FB ~ rl[-1] * L[-1] + rb[-1] * Bb[-1] + rmo[-1] * MOo[-1] -
    (rm[-1] * M[-1] + ra[-1] * A[-1]),
  #----------------------------------#
  # Central bank
  #----------------------------------#
  # [46] Accommodates demand for advances and buys bonds that are
  # not absorbed by Households and Banks.
  Bcb ~ B - Bh - Bb,
  # [47] Interest income is redistributed to the government.
  FC ~ ra[-1] * A[-1] + rb[-1] * Bcb[-1],
  # [47C1] Rate on Advances.
  ra ~ ( 1 + rra) * (1 + p_e) - 1,
  # [MEMO] Total High powered money.
  HP ~ HPb + HPo + HPc,
  #----------------------------------#
  # Government
  #----------------------------------#
  # [48] Deficit -  Collects taxes production, wages, and profits,
  # and any deficit is financed by issuing bonds.
  GD ~ (G + rb[-1] * Bh[-1] + rb[-1] * Bcb[-1] + rb[-1] * Bb[-1]) -
    (IT + Td + TF + FC),
  # [48C1] Rate on bonds.
  rb ~ (1 + rrb) * (1 + p_e)-1,
  # [49] Taxes on production.
  IT ~ tau * S,
  # [50] Tax on capitalists.
  Tdc ~ tau_d * wc * Nc,
  # [51] Tax on workers.
  Tdo ~ tau_d * wo * No,
  # [52] Total tax from wages.
  Td ~ Tdc + Tdo,
  # [53] Tax on profits.
  TF ~ tau_f * FT,
  # [54] New bonds issued.
  B ~ B[-1] + GD,
  # [55] Government spending inflator.
  G ~ g * p,
  # [56] Real Government spending.
  g ~ g[-1] * (1 + y_e),
  # [MEMO] Taxes.
  TT ~ IT + Td + TF,
  #----------------------------------#
  # The housing market
  #----------------------------------#
  # Demand for houses is laid down for capitalists and
  # workers in equations [13] and [26].

  # [57] Number of unsold homes - changes when the number of
  # newly built homes exceed the demand for homes. Can't be negative.
  HU ~ HU_cond * ((HU[-1] + HN - (Hc - Hc[-1])-(Ho - Ho[-1]))) + 0,
  # [0057C1] Condition helper.
  HU_cond ~ if(((HU[-1] + HN - (Hc - Hc[-1])-(Ho - Ho[-1]))) > cond) {1} else {0},
  # [0057C2] HU lag 1.
  HU_1 ~ HU[-1],
  # [58] Check for HND.
  HN ~ HN_cond * HND + 0,
  # [0058C1] Condition for HND.
  HN_cond ~ if(HND - HU[-1] > cond) {1} else{0},
  # [58C1] Supply of new homes - is a function of expected demand
  # and past capital gains.
  HND ~ HND_cond * (v_1 * (Hc[-1] * y_e + (Ho - Ho[-1])) +
                      v_2 * (ph_3 - ph_3[-1])) + 0,
  # [0058C1] Condition for Supply of new homes.
  HND_cond ~ if(v_1 * (Hc[-1] * y_e + (Ho - Ho[-1])) +
                v_2 * (ph_3 - ph_3[-1]) > 0) {1} else{0},
  # [0058C2] ph lag 1,
  ph_1 ~ ph[-1],
  # [0058C3] ph lag 2.
  ph_2 ~ ph_1[-1],
  # [0058C4] ph lag 3.
  ph_3 ~ ph_2[-1],
  # [59] Market price of Homes.
  ph ~ ph_cond1 * ph1 + ph_cond2 * ph2,
  # [0059C1] Condition for market price of Homes.
  ph_cond1 ~ if(HU < 10) {1} else {0},
  # [0059C2] Condition for market price of Homes.
  ph_cond2 ~ if(HU >= 10) {1} else {0},
  # [59C1] Total supply of Houses.
  HNS ~ HN + HU[-1],
  # [59C2] Helper for House prices ph2.
  ph2 ~ ph[-1] + (- v_3 * (HU_1 - HU_1[-1])) ,
  # [59C3] Helper for House prices ph1.
  ph1 ~ ph[-1] + ph[-1] * v_3 * (Hc[-1] * y_e + (Ho - Ho[-1]) - HNS[-1]),
  # [59C4] Lag 1 of Ho
  Ho_1 ~ Ho[-1],
  #----------------------------------#
  # Aggregate demand, un-/employment, wages
  #----------------------------------#
  # [60] Aggregate Demand -  Sales
  s ~ c + i + ih + g,
  # [60C1] Helper for Aggregate Demand - total consumption.
  c ~ cc + co,
  # [60C2] Real value of houses.
  ih ~ HN * p,
  # [60C3] Nominal value of houses
  IH ~ ih * p,
  # [61] Aggregate Demand inflator.
  S ~ s * p,
  # [62] Number of workers
  N ~ s/prod,
  # [63] Share of capitalists
  Nc ~ omega_c * N,
  # [64] Share of workers
  No ~ N - Nc,
  # [65] growth in real income
  y ~ (s-s[-1])/s[-1],
  # [66] Unemployment rate - Follows some sort of Okun's laws
  ur ~ ur_cond * ((ur[-1]) - (y - ny) / okun),
  # [0066C1] Condition Unemployment rate
  ur_cond ~ if((ur[-1] - (y - ny) / okun) > 0) {1} else {0},
  # [66C1] helper for ur.
  ur1 ~ ur1_cond1 * 0 + ur1_cond2 * 0.2 + ur1_cond3 * ur,
  # [0066C1] helper for ur.
  ur1_cond1 ~ if(ur <= cond) {1} else {0},
  # [0066C2] helper for ur.
  ur1_cond2 ~ if(ur >= 0.2) {1} else {0},
  # [0066C3] helper for ur.
  ur1_cond3 ~ if(ur > cond & ur < 0.2) {1} else {0},
  # [67] Wage Bill - All wages paid out.
  WB ~ w * N,
  # [67C1] Wage.
  w ~ (wc * omega_c + wo * (1 - omega_c)),
  # [68] Wage for capitalists.
  wc ~ wc[-1] * (1 + (wc_g)),
  # [69] Wage for workers.
  wo ~ wo[-1] * (1 + (wo_g)),
  # [70] Wage growth/inflation for capitalists.
  wc_g ~ p_e + omega * prodg_e + shockwc_g,
  # [71] Wage growth/inflation for workers.
  wo_g ~ p_e + omega * prodg_e + shockwo_g,
  # [72] Omega is the wage share,
  omega ~ o_0 - o_2 * sqrt(ur1/o_1),
  # [73] Productivity gains.
  prodg ~ pi_0 - pi_1 * dut + shockprodg,
  # [73C1] Production of Capitalists.
  prodc ~ prodc[-1] * (1 + prodg),
  # [73C1] Production of Workers.
  prodo ~ prodo[-1] * (1 + prodg),
  # [73C3] Total production.
  prod ~ omega_c * prodc + (1 - omega_c) * prodo,
  # [73C4] Change in utilization rate.
  dut ~ dut_cond * ((ut - unorm) / 100),
  # [0073C4] Change in utilization rate.
  dut_cond ~ if((ut - unorm) > 0) {1} else {0},
  # Accounting MEMO
  WBo ~ wo * No,
  WBc ~ WB - WBo,
  #----------------------------------#
  # Expectations X_e = X[-1] + sigma * (X_e[-1] - X[-1])
  #----------------------------------#
  # [74] Expected inflation
  p_e ~ pgr[-1] + sigma_p_e * (p_e[-1] - pgr[-1]),
  # [74C1] Inflation.
  pgr ~ (p - p[-1])/p[-1],
  # [75] Expected productivity growth
  prodg_e ~ prodg[-1] + sigma_pg * (prodg_e[-1] - prodg[-1]),
  # [76] Expected income growth
  y_e ~ y[-1] + sigma_yg * (y_e[-1] - y[-1]),
  # [77] Expected worker real income
  yo_e ~ yo_e[-1] * (1 + y_e),
  # [78] Expected capitalist real income.
  yc_e ~ yc_e[-1] * (1 + y_e),
  # [78] Expected capitalist nominal income.
  Yc_e ~ Yc[-1] * (1 + y_e) * (1 + p_e),
  # [79] Expected capital gains for homes (capitalists)
  CGHc_e ~ (ph_e - ph[-1]) * Hc[-1],
  # [80] Helper for Expected growth on capital gains for homes (capitalists).
  phg_e ~ ph[-1] / ph_1[-1] -1 + sigma_pe *
    (phg_e - (ph[-1] / ph_1[-1] -1)) + shockphg_e,
  # [81] Expected price of houses.
  ph_e ~ ph[-1] * (1 + phg_e),
  # [82] Expected real capital gains for homes (capitalists)
  cghc_e ~ (ph_e - ph[-1]) * Hc[-1] / p[-1] * (1 + p_e),
  # [83] Expected capital gains for homes (workers)
  CGHo_e ~ (ph_e - ph[-1]) * Ho[-1],
  # [84] Expected real capital gains for homes (workers)
  cgho_e ~ (ph_e - ph[-1]) * Ho[-1] / p[-1] * (1 + p_e),
  # [85] Expected capital gains on equities
  CGE_e ~ (pe_e - pe[-1]) * E[-1],
  # [86] Helper for Expected growth of capital gains on equities.
  peg_e ~ (pe[-1] - pe_1[-1]) / pe_1[-1] + sigma_pe *
    (peg_e[-1] - ((pe[-1] - pe_1[-1]) / pe_1[-1])) + shockpeg_e,
  # [87] pe lag 1,
  pe_1 ~ pe[-1],
  # [88] Expected price of equities.
  pe_e ~ pe[-1] * (1 + peg_e),
  # [89] Expected real capital gains on equities
  cge_e ~ ((pe_e - pe[-1]) * E[-1]) / (p[-1] * (1 + p_e)),
  # [90] Expected wealth for capitalists.
  Vc_e ~ Vc[-1] * (1 + y_e) - Cc + CGE_e + CGHc_e,
  # [91] Expected return on equities
  re_e ~ re[-1] + sigma * (re_e[-1] - re[-1]),
  # [92] Expected real return on equities
  rre_e ~ (1 + re_e) / (1 + p_e) -1,
  # [93] Expected return on houses.
  rh_e ~ rh[-1] + sigma * (rh_e[-1] - rh[-1]),
  # [94] Expected real return on houses.
  rrh_e ~ (1 + rh_e) / (1 + p_e) - 1,
)

# Set up parameter and initial conditions

model_ext <- sfcr_set(
  alpha_1c	~	0.7	, #	Consumption function	Parameters
  alpha_0o	~	0	, #	Consumption function	Parameters
  alpha_1o	~	0.8	, #	Consumption function	Parameters
  alpha_4o	~	10	, #	Consumption function	Parameters
  alpha_2c	~	0.025	, #	Consumption function	Parameters
  alpha_3c	~	0.08	, #	Consumption function	Parameters
  alpha_2o	~	0.025	, #	Consumption function	Parameters
  alpha_3o	~	0.08	, #	Consumption function	Parameters
  iota_0	~	-0.05	, #	Investment function	Parameters
  iota_1	~	2	, #	Investment function	Parameters
  iota_2	~	1	, #	Investment function	Parameters
  iota_3	~	0.2	, #	Investment function	Parameters
  iota_4	~	0.4	, #	Investment function	Parameters
  lambda	~	1.3	, #	Full capacity production	Parameters
  beta	~	0.1	, #	Retained profits	Parameters
  rho_1	~	0	, #	Mark-up	Parameters
  tau_d	~	0.2	, #	Tax rates	Parameters
  tau	~	0.1	, #	Tax rates	Parameters
  tau_f	~	0.4	, #	Tax rates	Parameters
  sigma_yg	~	0.75	, #	Expectations	Parameters
  sigma_p_e	~	0.75	, #	Expectations	Parameters
  sigma_pe	~	0.75	, #	Expectations	Parameters
  sigma_pg	~	0.75	, #	Expectations	Parameters
  sigma	~	0.75	, #	Expectations	Parameters
  eta	~	0.2	, #	Cash	Parameters
  chi_2	~	0.25	, #	Cash	Parameters
  xi	~	0.25	, #	Share of investment financed by issuing equities	Parameters
  lambda_20	~	0.18	, #	Parameters in asset demand function	Parameters
  lambda_21	~	0.3	, #	Parameters in asset demand function	Parameters
  lambda_22	~	0.25	, #	Parameters in asset demand function	Parameters
  lambda_23	~	0.01	, #	Parameters in asset demand function	Parameters
  lambda_24	~	0.1	, #	Parameters in asset demand function	Parameters
  lambda_25	~	0.05	, #	Parameters in asset demand function	Parameters
  lambda_10	~	0.5	, #	Parameters in asset demand function	Parameters
  lambda_11	~	0.45	, #	Parameters in asset demand function	Parameters
  lambda_12	~	0.25	, #	Parameters in asset demand function	Parameters
  lambda_13	~	0.01	, #	Parameters in asset demand function	Parameters
  lambda_14	~	0.25	, #	Parameters in asset demand function	Parameters
  lambda_15	~	0.05	, #	Parameters in asset demand function	Parameters
  lambda_30	~	0.18	, #	Parameters in asset demand function	Parameters
  lambda_31	~	0.45	, #	Parameters in asset demand function	Parameters
  lambda_32	~	0.25	, #	Parameters in asset demand function	Parameters
  lambda_33	~	0.01	, #	Parameters in asset demand function	Parameters
  lambda_34	~	0.25	, #	Parameters in asset demand function	Parameters
  lambda_35	~	0.1	, #	Parameters in asset demand function	Parameters
  okun	~	3	, #	Okun	Parameters
  ny	~	0.0338	, #	Normal growth	Parameters
  rra	~	0.025	, #	Real rate Advances	Parameters
  rrb	~	0.03	, #	Real rate Bonds	Parameters
  spread_1	~	0.005	, #	Banks spreads on central bank interest rate	Parameters
  spread_3	~	-0.025	, #	Banks spreads on central bank interest rate	Parameters
  spread_2	~	0.0025	, #	Banks spreads on central bank interest rate	Parameters
  omega_c	~	0.05	, #	Wage share capitalists	Parameters
  pi_0	~	0.02	, #	Production parameter	Parameters
  shockprodg	~	0	, #	Shock	Parameters
  shockwc_g	~	0	, #	Shock to wage growth capitalist	Parameters
  shockwo_g	~	0	, #	Shock to wage growth worker	Parameters
  chi_1	~	0.6	, #	Demand for bonds parameter	Parameters
  rrm	~	0.02	, #	Real rate on deposits	Parameters
  unorm	~	0.75	, #	Normal rate of utilization	Parameters
  pi_1	~	0	, #	Producivity gains	Parameters
  o_0	~	2	, #	Wage share parameter	Parameters
  o_1	~	0.05	, #	Wage share parameter	Parameters
  o_2	~	1	, #	Wage share parameter	Parameters
  v_1	~	0.5	, #	Parameter Number of new homes	Parameters
  v_2	~	1	, #	Parameter Number of new homes	Parameters
  v_3	~	0.0005	, #	House price parameter	Parameters
  shockphg_e	~	0	, #	Shock to expected growth of house price	Parameters
  im	~	0	, #	Impersonation parameter workers	Parameters
  shockpeg_e	~	0, #	Shock expected growth price of equities	Parameters
  morp	~	0.01	, #	Mortgages repayment rate	Parameters
  mu_1 ~ 0.001,
  mu_2 ~ 20,
  cond ~ 0
)

model_init <- sfcr_set(
  G ~ 1244.444, #	Nominal Government Debt	Initial value
  rho	~	0.4	, #	Mark-up	Initial value
  Bh_port	~	0.4	, #	Assets	Initial value
  pe_port	~	0.25	, #	Assets Initial value
  Hc_port	~	0.1	, #	Assets	Initial value
  omega	~	1	, #	Sum of wage share	Initial value
  rb	~	0.03	, #	Rate Bonds	Initial value
  ra	~	0.025	, #	Rate Advances (Leitzins)	Initial value
  rl	~	0.03	, #	Rate Loans	Initial value
  rm	~	0.02	, #	Rate Deposits	Initial value
  g	~	800, #	Government Debt	Initial value
  wo	~	0.833333333333333	, #	Wage workers	Initial value
  wc	~	4.16666666666667	, #	Wage capitalists	Initial value
  prodo	~	0.833333333333333	, #	Production workers	Initial value
  prodc	~	4.16666666666667	, #	Production capitalists	Initial value
  w	~	1	, #	Unit wage	Initial value
  prod	~	1	, #	Unit production	Initial value
  p	~	1.55555555555556	, #	Price	Initial value
  pe	~	5	, #	Equity price	Initial value
  pe_1	~	5	, #	Equity price	Initial value
  s	~	2128	, #	Real sales in aggregate demand	Initial value
  S	~	3310.22222222222	, #	Nominal sales in aggregate demand	Initial value
  c	~	1170.4	, #	Real total consumption	Initial value
  C	~	1820.62222222222	, #	Nominal total consumption	Initial value
  i	~	157.6	, #	Real investment	Initial value
  I	~	245.155555555555	, #	Nominal investment	Initial value
  yc_yo	~	1447.04	, #	Real total income	Initial value
  Yc_Yo	~	2250.95111111111	, #	Nominal total income	Initial value
  GD	~	193.581795555556	, #	Government deficit (Debt growth)	Initial value
  B	~	2383.36	, #	Bonds	Initial value
  Bh	~	66.2044444444444	, #	Households bonds	Initial value
  Bb	~	1125.47555555556	, #	Banks bonds	Initial value
  Bcb	~	1191.68	, #	CB Bonds (Residual)	Initial value
  E	~	410.467555555555	, #	Equities	Initial value
  L	~	2184.74666666667	, #	Loans	Initial value
  A	~	1721.31555555556	, #	CB Advances	Initial value
  M	~	4501.90222222222	, #	Total deposits	Initial value
  FU	~	165.511111111111	, #	Retained profits	Initial value
  y	~	0.02	, #	Real income growth	Initial value
  y_e	~	0.02	, #	Expected real income growth	Initial value
  pgr	~	0	, #	Inflation	Initial value
  p_e	~	0	, #	Expected inflation	Initial value
  pe_e	~	5	, #	Expected price of equities	Initial value
  K	~	3310.22222222222	, #	Capital	Initial value
  K_1 ~ 3310.22222222222,
  k	~	2128	, #	Real capital	Initial value
  lev	~	0.66	, #	Leverage ratio	Initial value
  q	~	0.62	, #	Tobin's Q	Initial value
  sfc	~	2766.4	, #	Normal Sales - fixed ratio of real capital	Initial value
  ut	~	0.769230769230769	, #	Utlization ratio	Initial value
  prodg	~	0.02	, #	Production growth	Initial value
  N	~	2128	, #	Number of People	Initial value
  ur	~	0.05	, #	Unemployment rate	Initial value
  ur1	~	0.05	, #	Unemployment rate	Initial value
  prodg_e	~	0.02	, #	Expected production growth	Initial value
  peg_e	~	0	, #	Expected growth of price of equities	Initial value
  wc_g	~	0.02	, #	Wage growth capitalist	Initial value
  wo_g	~	0.02	, #	Wage growth worker	Initial value
  wc_g	~	0.02	, #	Wage growth capitalist	Initial value
  wo_g	~	0.02	, #	Wage growth worker	Initial value
  WB	~	2128	, #	Wage bill	Initial value
  FT	~	851.2	, #	Total profits	Initial value
  TF	~	340.48	, #	Taxes on Profits	Initial value
  FD	~	400.65984	, #	Distributed profits	Initial value
  CGE	~	0	, #	Capital gains on equities	Initial value
  CGE_e	~	0	, #	Expected capital gain on equities	Initial value
  re	~	0.195221198156682	, #	Return on equities	Initial value
  re_e	~	0.195221198156682	, #	Expected return on equities	Initial value
  rre_e	~	0.195221198156682	, #	Expected real return on equities	Initial value
  rll	~	0.03	, #	Real rate on loans	Initial value
  dut	~	0	, #	Change in rate of utilization	Initial value
  Nc	~	106.4	, #	Number of capitalists	Initial value
  No	~	2021.6	, #	Number of woekers	Initial value
  Mc	~	1383	, #	Deposits capitalists	Initial value
  Mo	~	2000	, #	Deposits workers	Initial value
  Td	~	450.190222222222	, #	Total tax on labour	Initial value
  Tdc	~	450.190222222222	, #	Tax on labour capitalists	Initial value
  Tdo	~	0	, #	Tax on labour workers	Initial value
  rent	~	0.5	, #	Rent of house	Initial value
  Rents	~	375	, #	Total rent	Initial value
  yc	~	1218.80851008	, #	Real income capitalist	Initial value
  yo	~	228.23148992	, #	Real income worker	Initial value
  Vc	~	6620.44444444444	, #	Nominal wealth capitalists	Initial value
  Vo	~	0	, #	Nominal wealth workers	Initial value
  vo	~	0	, #	Real wealth workers	Initial value
  vc	~	4256	, #	Real wealth capitlists	Initial value
  Vc_e	~	6620.44444444444	, #	Expected nominal wealth capitalists	Initial value
  Yc	~	861.078705777778	, #	Nominal Income capitalists	Initial value
  Yo	~	355.026762097778	, #	Nominal Income workers	Initial value
  CGHc	~	0	, #	Nominal capital gains in houses capitalists	Initial value
  CGHo	~	0	, #	Nominal capital gais in houses workers	Initial value
  cgho_e	~	0	, #	Expected real capital gains in houses workers	Initial value
  ph	~	3	, #	Price houses	Initial value
  ph_1	~	3	, #	Price houses	Initial value lag1
  ph_2	~	3	, #	Price houses	Initial value lag2
  ph_3	~	3	, #	Price houses	Initial value lag3
  Hc	~	37	, #	Houses capitalits	Initial value
  Ho	~	1430	, #	Houses workers	Initial value
  Ho_1	~	1430	, #	Houses workers	Initial value lag1
  phg_e	~	0	, #	Expected housing price growth	Initial value
  MOo ~ 1000, # Mortgages
  MOo_1 ~ 1000, # Mortgages lag 1
  rmo	~	0.025	, #	Rate on Mortgages	Initial value
  rh	~	0.01875	, #	Return on houses	Initial value
  rh_e	~	0.01875	, #	Expected retur on houses	Initial value
  HU	~	100	, #	Number of unsold homes	Initial value
  HU_1	~	100	, #	Number of unsold homes	Initial value lag1
  HN	~	0	, #	Number of new homes	Initial value
  ih	~	0	, #	Real value of houses	Initial value
  IH	~	0	, #	Nominal value of houses	Initial value
  Ho_d	~	0	, #	Change in number of houses workers	Initial value
  HND	~	0	, #	Demand for new homes	Initial value
  HNS	~	0	, #	Supply of homes	Initial value
  cc	~	1170.4	, #	Real consumption capitalits	Initial value
  co	~	1170.4	, #	Real consumption workers	Initial value
  Cc ~ 1820.62, # Nominal consumption capitalits	Initial value
  Co ~ 1820.62, # Nominal consumption workers	Initial value
  debt_rep ~ 0.099, # Debt Repayment
  delta_debt_rep ~ 0, # Change in Debt Repayment
  HPc ~ 364.1244444, # High-powerd money capitalists
  HPo ~ 364.1244444, # High-powerd money workers
  Sho ~ -1465.595, # Savings workers
  Shc ~ -959.544, # Savings capitalists
  kgr ~ 0.161892308, # Investmentc(growth of capital stock)
  Hc_port ~ 0.447456565, # Housing portfolio capitalists
  Yc_e ~ 861.0787058, # Expected nominal income capitalists
  rrh_e ~ 0.01875, # Expected real return on houses
  FB ~ 13.61377778, # Banks Profit
  Yc_Yo0 ~ 2134.069529, # Total nominal income
  yc_yo_e ~ 1447.040, # Expected real total income
  Vf ~ 1212.150123, # Firms wealth
  HPb ~ 1125.475556, # Demand for High-powered money in Banks
  HP ~ 2979.2, # Sum of High-powered money
  A0 ~ 1787.520, # Advances check
  TT ~ 1121.692, # Total taxes
  IT ~ 331.02222, # Taxes on production
  ph1 ~ 3,
  ph2 ~ 3,
  WBo ~ 1684.666667, # Wage bill Workers
  WBc ~ 443.3333333, # Wage bill capitalist
  yo_e ~ 1218.809, # Real income workers
  yc_e ~ 228.231, # Real income capitalists
  HU_cond ~ 1,
  HN_cond ~ 1,
  HND_cond ~ 1,
  ph_cond1 ~ 0,
  ph_cond2 ~ 1,
  MOo_cond ~ 1,
  dut_cond ~ 1,
  ur1_cond1 ~ 0,
  ur1_cond2 ~ 0,
  ur1_cond3 ~ 1,
  ur_cond ~ 1,
  cond ~ 0
)

model_eqs_df <- unlist(as.character(model_eqs))
write.xlsx(model_eqs_df, "model_equations.xlsx")

model_init_df <- unlist(as.character(model_init))
write.xlsx(model_init_df, "model_initial.xlsx")

model_ext_df <- unlist(as.character(model_ext))
write.xlsx(model_ext_df, "model_exo.xlsx")

