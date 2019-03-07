## Translated into mrgsolve by Miao Li based on PFOS desolve code from Wei-Chun Chou (Original code is from PFOS PBPK model, Worley and Fisher, 2015)

RatPBPK.code <- '
$PARAM @annotated
QCC                 :  14.000 : L/h/kg^0.75,         Cardiac output (Brown 1997)
QLC                 :   0.183 : Unitless,            Fraction blood flow to liver (Brown 1997)
QKC                 :   0.141	: Unitless,            Fraction blood flow to kidney (Brown 1997)
Htc                 :   0.46  : Unitless,            Hematocrit for Rat (Davies 1993)
BW                  :   0.3   : kg,                  Bodyweight (measrument data if available)
VLC                 :   0.035 : Unitless,            Fractional liver tissue (Brown 1997)
VKC                 :   0.0084: Unitless,            Fractional kidney tissue (Brown 1997)
VPlasC              :   0.0312: L/kg BW,             Fractional plasma (Davies 1993)
VfilC               : 0.00084 : L/kg BW,             Fraction vol. of filtrate; 10% of Kidney volume; (Worley and Fisher et al., 2015)
VPTCC               : 1.35e-4 : L/kg kidney,         Volume of proximal tubule cells (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)(Hsu et al., 2014)
FVBK                :   0.160 : Unitless,            Blood volume fraction of kidney (Brown, 1997)
PL                  :   3.72  : Unitless,            Liver/ plasma PC; (Loccisano et al., 2012) 
PK                  :   0.80  : Unitless,            Kidney/ plasma PC; (Loccisano et al., 2012)
PRest               :   0.20  : Unitless,            Restofbody/ plasma PC; (Loccisano et al., 2012)
MW                  : 500.126 : g/mol,               PFOS molecular mass 
MKC                 : 0.0084  : Unitless,            Fraction mass of kidney (percent of BW); Brown, 1997
Free                : 0.09    : Unitless,            Free fraction; (Worley and Fisher et al., 2015) 
Vmax_baso_invitro   : 393.45  : pmol/mg,             Protein/min, Vmax of basolateral transporter; averaged in vitro value of rOAT1 and rOAT3 (Nakagawa, 2007); initial value asumsed the same as PFOS (Worley and Fisher, 2015)
Km_baso             :  27.2   : mg/L,                Km of basolateral transpoter, averaged in vitro value of rOAT1 and rOAT3 (Nakagawa, 2007); initial value asumsed the same as PFOS (Worley and Fisher, 2015)
Vmax_apical_invitro :  9300   : pmol/mg protein/min, Vmax of apical transporter; averaged invitro value of Oatp1a1 (Weaver, 2010); initial value asumsed the same as PFOS (Worley and Fisher, 2015)
Km_apical           :  52.3   : mg/L,                Km of apical transpoter, in vitro value for Oatp1a1 (Weaver, 2010);initial value asumsed the same with PFOS (Worley and Fisher, 2015) 
RAFbaso             :  3.99   : Unitless             Relative activity factor, basolateral transpoters (male) (fit to data); initial value asumsed the same as PFOS (Worley and Fisher, 2015)
RAFapi              :  4.07   : Unitless             Relative acitivty factor, apical transpoters (fit to data); 0.001356 (female);initial value asumsed the same as PFOS (Worley and Fisher, 2015)
protein             :  2.0e-6 : mg protein/PTCs,     Amount of protein in proximal tubule cells (Addis et al., 1936)
GFRC                :  62.1   : L/hr/kg kiney,       Glomerular filtration rate (male) (Corley, 2005)
Kdif                :   0.001 : L/h,                 Diffusion rate from proximal tubule cells to kidney serum; Assumed the same as PFOS (Worley and Fisher, 2015)
Kabsc               :   2.12  : 1/(h*BW^-0.25),      Rate of absorption of PFOS from small intestine to liver (initial value assumed the same as PFOS (2.12) from Worley and Fisher, 2015 and then re-fitting)                        
KunabsC             : 7.05e-5 : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces (initial value assumed the same as PFOS (7.05e-5) from Worley and Fisher, 2015 and then re-fitting)
GEC                 :   0.540 : 1/(h*BW^0.25),       Gastric emptying time  (Yang et al., 2013)
K0C                 :   1.000 : 1/(h*BW^-0.25),      Rate of uptake from the stomach into the liver (initial value assumed the same as PFOS (1) from Worley and Fisher, 2015 and then re-fitting)
KeffluxC            :   2.49  : 1/(h*BW^-0.25),      Rate of clearance of PFOS from proximal tubule cells into blood (initial value assumed the same as PFOS (2.49) from Worley and Fisher, 2015 and then re-fitting)
KbileC              :  0.004  : 1/(h*BW^-0.25),      Biliary elimination rate (male); liver to feces storage (initial value assumed the same as PFOS (0.004) from Worley and Fisher, 2015 and then re-fitting)  
KurineC             :  1.60   : 1/(h*BW^-0.25),      Rate of urine elimination from urine storage (male) (initial value assumed from Worley and Fisher, 2015 and then re-fitting)  

$MAIN
double QC = QCC*pow(BW, 0.75)*(1-Htc);               // L/h, Cardiac output (adjusted for plasma)
double QK = QKC*QC;                                  // L/h, Plasma flow to kidney
double QL = QLC*QC;                                  // L/h, Plasma flow to liver
double QRest = QC-QK-QL;                             // L/h, Plasma flow to the rest of body

double VL = VLC*BW;                                  // L,   Volume of liver 
double VK = VKC*BW;                                  // L,   Volume of kidney 
double VPlas = VPlasC*BW;                            // L,   Volume of plasma  
double Vfil = VfilC*BW;                              // L,   Volume of filtrate  
double VPTC = VK*VPTCC;                              // L,   Volume of proximal tubule cells 
double VKb = VK*FVBK;                                // L,   Volume of blood in the kidney; fraction blood volume of kidney (0.16) from Brown, 1997
double VRest = (0.93*BW)-VL-VK-VPlas-VPTC-Vfil;      // L,   Rest of body; volume of remaining tissue (L); Revised original equation (VR = (0.93*BW) - VPlas - VPTC - Vfil - VL) from Worley and Fisher, 2015  
double ML = VL*1000;                                 // g,   Liver weight in gram
double MK = VK*1000;                                 // g,   kidney weight in gram

double PTC = MKC*6e7*1000;                           // cells/kg BW, Number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney); Revised original equation (PTC = MKC*6e7) from Worley and Fisher, 2015
double MPTC = VPTC*1000;                             // g,           mass of the proximal tubule cells (assuming density 1 kg/L)
double Vmax_basoC = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000);         // mg/h/kg BW^0.75, Vmax of basolateral transporters (average Oat1 and Oat3) equation from Worley and Fisher, 2015
double Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000);      // mg/h/kg BW^0.75, Vmax of apical transpoters in in vitro studies (Oatp1a1) equation from Worley and Fisher, 2015
double Vmax_baso = Vmax_basoC*pow(BW,0.75);          // mg/h, Vmax of basolateral transporters
double Vmax_apical = Vmax_apicalC*pow(BW,0.75);      // mg/h, Vmax of apical transpoters
double Kbile = KbileC*pow(BW,(-0.25));               // 1/h, Biliary elimination, liver to feces storage
double Kurine = KurineC*pow(BW,(-0.25));             // 1/h, Urinary elimination; 
double Kefflux = KeffluxC*pow(BW,(-0.25));           // 1/h, Efflux clearance rate from PTC to blood 
double GFR = GFRC*(MK/1000);                         // L/h, Glomerular filtration rate, scaled to mass of kidney 

//GI tract parameters
double Kabs = Kabsc*pow(BW,(-0.25));                 // 1/h, Rate of absorption of PFOS from small intestine to liver
double Kunabs = KunabsC*pow(BW,(-0.25));             // 1/h, Rate of unabsorbed dose to appear in feces
double GE = GEC*pow(BW,(-0.25));                     // 1/h, Gasric emptying time 
double K0 = K0C*pow(BW,(-0.25));                     // 1/h, Rate of uptake from the stomach into the liver

$CMT A_baso A_apical Adif Aefflux ACI APlas_free AUCCA_free APTC AUCCPTC AFil AUCfil Aurine AKb AUCKb ARest AUCCRest AST
AabsST ASI AabsSI Afeces AL AUCCL 

$ODE

// Concentrations in the tissues and in the venous palsma leaving each of the tissues
double CA_free = APlas_free/VPlas;                   // mg/L, Free PFOS concentration in the plasma
double CA = CA_free/Free;                            // mg/L, Concentration of total PFOS in the plasma

double CL = AL/VL;                                   // mg/L, Concentration of PFOS in the liver compartment
double CVL = CL/PL;                                  // mg/L, Concentration of PFOS in plasma leaving liver

double CKb = AKb/VKb;                                // mg/L, Concetraitons of PFOS in venous plasma leaving kidney 
double CVK = CKb;                                    // mg/L, Concentration of PFOS in plasma leaving kidney
double CK = CVK*PK;                                  // mg/L, Concetraitons of PFOS in Kidney compartment

double CRest = ARest/VRest;                          // mg/L, Concentration of PFOS in the rest of the body
double CVRest = ARest/(VRest*PRest);                 // mg/L, Concentration of PFOS in the venous palsma leaving the rest of body

// Kidney compartment plus 2 subcompartment (Proximal Tubule cells: PTCs, Filtrate: Fil)
// Concentration in kidney, PTCs and fil

double CPTC = APTC/VPTC;                             // mg/L, Concetraitons of PFOS in PTCs
double Cfil = AFil/Vfil;                             // mg/L, Concetraitons of PFOS in Fil

// Virtural compartment: 
// Basolateral (baso) 
// transport, Diffusion (dif), 
// Apical (apical) transport, and 
// efflux clearance (efflux)
// Clerance (CL) via glormerular filtration 

double RA_baso = (Vmax_baso*CKb)/(Km_baso+CKb);                                             // mg/h, Rate of basolateral transpoters 
dxdt_A_baso = RA_baso;                                                                      // mg,   Amount of basolateral transpoters
double RA_apical = (Vmax_apical*Cfil)/(Km_apical + Cfil);                                   // mg/h, Rate of apical transpoter 
dxdt_A_apical = RA_apical;                                                                  // mg,   Amount of apical transpoter
double Rdif = Kdif*(CKb - CPTC);                                                            // mg/h, Rate of diffusion from into the PTC
dxdt_Adif = Rdif;                                                                           // mg,   Amount moved via glomerular filtration
double RAefflux = Kefflux*APTC;                                                             // mg/h, Rate of efflux clearance rate from PTC to blood
dxdt_Aefflux = RAefflux;                                                                    // mg,   Amount of efflux clearance rate from PTC to blood
double RCI = CA*GFR*Free;                                                                   // mg/h, Rate of clerance (CL) to via glomerular filtration (GFR) 
dxdt_ACI = RCI;                                                                             // mg,   Amount of clearance via GFR

// {PFOS distribution in each compartment}
// PFOS in plasma
double RPlas_free = (QRest*CVRest*Free)+(QK*CVK*Free)+(QL*CVL*Free)-(QC*CA*Free)+RAefflux;  // mg/h,   Rate of change in the plasma
dxdt_APlas_free = RPlas_free;                                                               // mg,     Amount of free PFOS in the plasma
dxdt_AUCCA_free = CA_free;                                                                  // mg*h/L, Area under curve of free PFOS in plasma compartment

// Proximal Tubule Cells (PTCs)
double RPTC = Rdif + RA_apical + RA_baso - RAefflux;                                        // mg/h,   Rate of change in PTCs 
dxdt_APTC = RPTC;                                                                           // mg,     Amount in the PTCs
dxdt_AUCCPTC = CPTC;                                                                        // mg*h/L, Area under curve of PFOS in the compartment of PTCs

// Proximal Tubule Lumen/ Filtrate (Fil)
double Rfil = CA*GFR*Free - RA_apical - AFil*Kurine;                                        // mg/h,   Rate of change in the Fil
dxdt_AFil = Rfil;                                                                           // mg,     Amount in the Fil
dxdt_AUCfil = Cfil;                                                                         // mg*h/L, Area under curve of PFOS in the compartment of Fil

// Urine elimination
double Rurine = Kurine*AFil;                                                                // mg/h,   Rate of change in urine
dxdt_Aurine = Rurine;                                                                       // mg,     Amount in urine
//double percentOD_in_urine = (Aurine/Odose)*100;                                           // %,      Percent of oral dose in the urine  

// Kidney compartment
double RKb = QK*(CA-CVK)*Free - CA*GFR*Free - Rdif - RA_baso;                               // mg/h,   Rate of change in Kidney compartment
dxdt_AKb = RKb;                                                                             // mg,     Amount in kidney compartment
dxdt_AUCKb = CK;                                                                            // mg*h/L, Area under curve of PFOS in the Kidney compartment

// PFOS in the compartment of rest of body, flow-limited model
double RRest = QRest*(CA-CVRest)*Free;                                                      // mg/h,   Rate of change in rest of body
dxdt_ARest = RRest;                                                                         // mg,     Amount in rest of body 
dxdt_AUCCRest = CRest;                                                                      // mg*h/L, Area under curve of PFOS in the compartment of rest of body
    
// Gastrointestinal (GI) tract
// Stomach compartment
double RST = - K0*AST - GE*AST;                                                             // mg/h,   Rate of chagne in stomach caomprtment
dxdt_AST = RST; // mg, Amount in Stomach
double RabsST = K0*AST;                                                                     // mg/h,   Rate of absorption in the stomach
dxdt_AabsST = RabsST;                                                                       // mg,     Amount absorbed in the stomach

// Small intestine compartment
double RSI = GE*AST - Kabs*ASI - Kunabs*ASI;                                                // mg/h,   Rate of chagne in small intestine caomprtment
dxdt_ASI = RSI;                                                                             // mg,     Amount in small intestine
double RabsSI = Kabs*ASI;                                                                   // mg/h,   Rate of absorption in the small intestine
dxdt_AabsSI = RabsSI;                                                                       // mg,     Amount absorbed in the small intestine
double Total_oral_uptake = AabsSI + AabsST;                                                 // mg,     Total oral uptake in the GI

// Biliary excretion
double Abile = Kbile*AL;                                                                    // mg,     Amount of PFOS in bile excretion
double amount_per_gram_liver = (AL/ML)*1000;                                                // ug/g,   Amount of PFOS in liver per gram liver

// Feces compartment
double Rfeces = Kbile*AL + Kunabs*ASI;                                                      // mg/h,   Rate of change in feces compartment
dxdt_Afeces = Rfeces;                                                                       // mg,     Amount of the feces compartment
// double percentOD_in_feces = (Afeces/Odose)*100;                                          // %,      Percent of the oral dose in the feces

// PFOS in liver compartment, flow-limited model
double RL = QL*(CA-CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST;                                // mg/h,   Rate of chagne in liver caomprtment
dxdt_AL = RL;                                                                               // mg,     Amount in liver compartment
dxdt_AUCCL = CL;                                                                            // mg*h/L, Area under curve of PFOS in liver compartment

// {Mass balance equations}
double Qbal = QC-QL-QK-QRest;
double Tmass = APlas_free + ARest + AKb + AFil + APTC + AL + AST + ASI;
double Loss = Aurine + Afeces;
double Input = AabsSI + AabsST;
//double Bal = Input- Tmass - Loss;

$TABLE
capture Plasma = CA_free/Free;
capture Liver  = AL/VL;
capture Kidney = CK;
capture AUC_CA = AUCCA_free;
capture AUC_CL = AUCCL;
capture AUC_CK = AUCKb;
'
