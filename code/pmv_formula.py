# Implementation of PMV calculation according to Fanger equation for thermal comfort

import math

def calculate_PMV_value(Temp = 23.5, MeanR = 25.5, VEL = 0.2 , hum = 60, MET = 1.2, CLO = 0.5, EXW = 0.0, WV = 0.0):
    TA = Temp # Room Temperature
    TR = MeanR # Mean Radiation Temperature
    RH = hum # Relative Humidity

    WME = EXW # = 0.0 # External work, normally around 0
    PA = WV # = 0.0 # Water vapour pressure

    FNPS = math.exp(16.6536 - 4030.183 / (TA + 235)) # saturated vapour pressure, KPA

    if PA == 0:
        PA = RH * 10 * FNPS

    ICL = 0.155 * CLO # ICL = thermal insulation of the clothing
    M = MET * 58.15 # M = Metabolic rate in W/m^2
    W = WME * 58.15 # W = external work in W/m^2
    MW = M - W # internal heat production in the human body

    if ICL <= 0.078:
        FCL = 1 + 1.29 * ICL # FCL = clothing area factor
    else:
        FCL = 1.05 + 0.645 * ICL  # 0.645

    HCF = 12.1 * math.sqrt(VEL)

    TAA = TA + 273
    TRA = TR + 273

    # CALCULATE SURFACE TEMPERATURE OF CLOTHING BY ITERATION

    TCLA = TAA + (35.5 - TA) / 3.5 * (6.45 * (ICL + 0.1)) 
    # first guess for surface temperature of clothing

    P1 = ICL * FCL
    P2 = P1 * 3.96
    P3 = P1 * 100
    P4 = P1 * TAA
    P5 =308.7 - 0.028 * MW + P2 * pow((TRA / 100), 4) #35.7 ->308.7
    XN = TCLA / 100
    XF = XN
    N = 0  # number of iteration
    EPS = 0.00001

    # comeback
    while True:
        #print("goto function")
        XF = (XF + XN) / 2
        HCN = 2.38 * pow((abs(100 * XF - TAA)), 0.25)
        if HCF > HCN:
            HC = HCF
        else:
            HC = HCN
        XN = (P5 + P4 * HC - P2 * pow(XF, 4)) / (100 + P3 * HC) #FFFFFFNNN
        #print(XF, XN, HC)
        N = N + 1
        if N > 150:
            print("FEHLER mit pmv")
        if abs(XN - XF) <= EPS:
            break
    # print(N)

    TCL = 100 * XN - 273
    # HEAT LOSS COMPONENTS
    HL1 = 3.05 * 0.001 * (5733 - 6.99 - MW - PA)
    if MW > 58.15:
        HL2 = 0.42 * (MW - 58.15)
    else:
        HL2 = 0

    HL3 = 1.7 * 0.00001 * M * (5867 - PA)
    HL4 = 0.0014 * M * (34 - TA)
    HL5 = 3.96 * FCL * (pow(XN, 4) - pow((TRA / 100), 4))
    HL6 = FCL * HC * (TCL - TA)
    # CALCULATE PMV AND PPD
    TS = 0.303 * math.exp(-0.036 * M) + 0.028
    PMV = TS * (MW - HL1 - HL2 - HL3 - HL4 - HL5 - HL6)
    PPD = 100 - 95 * math.exp(-0.03353 * pow(PMV, 4) - 0.2179 * pow(PMV, 2))

    return PMV, PPD
