##translated into mrgsolve by Miao Li based on PFOA desolve code from Wei-Chun Chou (revised Worley et al. 2015)

MonkeyPBPK.code <- '
$PARAM @annotated
QCC                 :   18.96 : L/h/kg^0.75, Cardiac output,male (Brown 1997)
QLC                 :   0.194 : Fraction blood flow to liver (Brown 1997)
QKC                 :   0.123	: Fraction blood flow to kidney (Brown 1997)
Htc                 :   0.42  : Hematocrit for Monkey (Choi et al., 2016, table2)
BW                  :   3.5   : kg, Bodyweight (Seacat, et al., 2002)
VLC                 :   0.026 : Fractional liver tissue (Brown 1997)
VKC                 :   0.004 : Fractional kidney tissue (Brown 1997)
VPlasC              :   0.0448: L/kg BW, Fractional plasma (Davies 1993)
VfilC               :   4e-4  : L/kg BW, Fraction vol. of filtrate; (10% of kidney volume)
VPTCC               : 1.35e-4 : L/kg kidney, Vol. of proximal tubule cells (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)(Hsu et al., 2014)
FVBK                :   0.160 : Blood volume fraction of kidney (Brown, 1997)
PL                  :   3.72  : Liver/plasma PC (Luccisanno et al., 2011)
PK                  :   0.8   : Kidney/plasma PC (Luccisanno et al., 2011)
PRest               :   0.2   : Restofbody/plasma PC (Luccisanno et al., 2011)
MW                  : 500.126 : g/mol, PFOS molecular mass 
MKC                 : 0.0084  : Fraction mass of kidny (percent of BW); Brown, 1997
Free                : 0.025   : Free fraction; (Luccisanno et al., 2011) 
Vmax_baso_invitro   : 439.2   : pmol/mg protein/min, Vmax of basolateral transporter; averaged in vitro value of rOAT1 and rOAT3 (Nakagawa, 2007);(initial value assumed the same with PFOA in rat from worley et al., 2017 and then re-fitting)
Km_baso             :  20.1   : mg/L, Km of basolateral transpoter, averaged in vitro value of rOAT1 and rOAT3 (Nakagawa, 2007); (initial value assumed the same with PFOA in rat from worley et al., 2015 and then re-fitting)
Vmax_apical_invitro :  37400  : pmol/mg protein/min, Vmax of apical transporter; averaged invitro value of Oatp1a1 (Weaver, 2010); (initial value assumed the same with PFOA in rat from worley et al., 2015 and then re-fitting)
Km_apical           :  77.5   : mg/L, Km of apical transpoter, in vitro value for Oatp1a1 (Weaver, 2010); (initial value assumed the same with PFOA in rat from worley et al., 2015 and then re-fitting) 
RAFbaso             :  1.0    : Relative activity factor, basolateral transpoters (male) (fit to data); 0.01356 (female)
RAFapi              :  7e-4   : Relative acitivty factor, apical transpoters (fit to data); 0.001356 (female); (initial value assumed the same with PFOA in rat from worley et al., 2015 and then re-fitting)
protein             :  2.0e-6 : mg protein/proximal tubuel cell, Amount of protein in proximal tubule cells (Addis et al., 1936)
GFRC                :  21.85  : L/hr/kg kiney, Glomerular filtration rate (male); (IWAMA, 2014)
Kdif                :   0.001 : L/h, Diffusion rate from proximal tubule cells; Assumed the same with PFOA Worley et al. (2015)
Kabsc               :   2.12  : 1/(h*BW^-0.25), Rate of absorption of chemical from small intestine to liver (initial value assumed the same with PFOA (2.12) from worley et al., 2015 and then re-fitting)                        
KunabsC             : 7.06e-5 : 1/(h*BW^-0.25), Rate of unabsorbed dose to appear in feces (initial value assumed the same with PFOA (7.05e-5) from worley et al., 2015 and then re-fitting)
GEC                 :   2.340 : 1/(h*BW^0.25),  Gastric emptying time; Fisher et al. (2011)  
K0C                 :   1.000 : 1/(h*BW^-0.25), Rate of uptake from the stomach into the liver (initial value assumed the same with PFOA (1) from worley et al., 2015 and then re-fitting)
KeffluxC            :   0.1   : 1/(h*BW^-0.25), Rate of clearance of PFOS from proximal tubule cells into blood (initial value assumed the same with PFOA (2.49) from worley et al., 2015 and then re-fitting)
KbileC              :  0.0001 : 1/(h*BW^-0.25), Biliary elimination rate (male); liver to feces storage (initial value assumed the same with PFOA (0.004) from worley et al., 2015 and then re-fitting)  
KurineC             :  0.062   : 1/(h*BW^-0.25), Rate of urine elimination from urine storage (male) (initial value from Loccisano et al., 2011 and then re-fitting)  

$MAIN
double QC = QCC*pow(BW, 0.75)*(1-Htc);   // L/h, Cardiac output (adjusted for plasma)
double QK = QKC*QC;                     // L/h, Plasma flow to kidney
double QL = QLC*QC;                    // L/h, Plasma flow to liver
double QRest = QC-QK-QL;              // L/h, Plasma flow to the rest of body

double VL = VLC*BW;     // volume of liver (L)
double VK = VKC*BW;     // volume of kidney (L)
double VPlas = VPlasC*BW;  // volume of plasma (L) 
double Vfil = VfilC*BW; // volume of filtrate (L) 
double VPTC = VK*VPTCC; // volume of proximal tubule cells (L)
double VKb = VK*FVBK;   // volume of blood in the kidney (L); fraction blood volume of kidney (0.16) from Brown, 1997
double VRest = (0.93*BW)-VL-VK-VPlas-VPTC-Vfil; //Rest of body; volume of remaining tissue (L); Revised original equation (VR = (0.93*BW) - VPlas - VPTC - Vfil - VL) from Worley et al. (2015)  
double ML = VL*1000;    // g, liver weight in gram
double MK = VK*1000;    // g, kidney weight in gram

double PTC = MKC*6e7*1000; // cells/kg BW, Number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney,Hsu et al., 2014); Revised original equation (PTC = MKC*6e7) from Worley et al. (2015)
double MPTC = VPTC*1000; // g, mass of the proximal tubule cells (assuming density 1 kg/L)
double Vmax_basoC = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000); // mg/h/kg BW^0.75, Vmax of basolateral transporters (average Oat1 and Oat3) equation from Worley et al. (2015)
double Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000); // mg/h/kg BW^0.75, Vmax of apical transpoters in in vitro studies (Oatp1a1) equation from Worley et al. (2015)
double Vmax_baso = Vmax_basoC*pow(BW,0.75); // mg/h
double Vmax_apical = Vmax_apicalC*pow(BW,0.75); // mg/h
double Kbile = KbileC*pow(BW,(-0.25)); // 1/h, Biliary elimination, liver to feces storage
double Kurine = KurineC*pow(BW,(-0.25)); // 1/h, Urinary elimination; from filtrate
double Kefflux = KeffluxC*pow(BW,(-0.25)); // 1/h, Efflux clearance rate from PTC to blood 
double GFR = GFRC*(MK/1000); // L/h, Glomerular filtration rate, scaled to mass of kidney 

//GI tract parameters
double Kabs = Kabsc*pow(BW,(-0.25));     // 1/h, rate of absorption of chemical from small intestine to liver
double Kunabs = KunabsC*pow(BW,(-0.25)); // 1/h, rate of unabsorbed dose to appear in feces
double GE = GEC*pow(BW,(-0.25));         // 1/h, Gasric emptying time 
double K0 = K0C*pow(BW,(-0.25));         // 1/h, Rate of uptake from the stomach into the liver

$CMT A_baso A_apical Adif Aefflux ACI APlas_free AUCCA_free APTC AUCCPTC AFil AUCfil Aurine AKb AUCKb ARest AUCCRest AST
AabsST ASI AabsSI Afeces AL AUCCL 

$ODE

// Concentrations in the tissues and in the capillary blood of the tissues
double CA_free = APlas_free/VPlas; // mg/L, Free PFOS oncentration in the plasma
double CA = CA_free/Free; // mg/L, Concentration of total PFOS in plasma

double CL = AL/VL; // mg/L, Concentration of parent drug in the tissue of liver
double CVL = CL/PL; // mg/L, Concentration of parent drug in the capillary blood of liver

double CKb = AKb/VKb; // mg/L, Concetraitons of PFOS in Kidney blood
double CVK = CKb; // mg/L, concentration in venous blood leaving kidney
double CK = CVK*PK; // mg/L, Concetraitons of PFOS in Kidney compartment

double CRest = ARest/VRest; // mg/L, Crest drug concentration in the rest of the body (mg/L)
double CVRest = ARest/(VRest*PRest); // mg/L, Concentration of parent drug in the capillary blood of the rest of body

// Kidney compartment plus 2 subcompartment (Proximal Tubule cells: PTC, Filtrate: Fil)
// Concentration in kidney, PTC and fil

double CPTC = APTC/VPTC; // mg/L, Concetraitons of PFOS in PTC blood
double Cfil = AFil/Vfil; // mg/L, Concetraitons of PFOS in FIL blood

// Virtural compartment: 
// Basolateral (baso) 
// transport, Diffusion (dif), 
// Apical (apical) transport, and 
// efflux clearance (efflux)
// Clerance (CL) via glormerular filtration 
double RA_baso = (Vmax_baso*CKb)/(Km_baso+CKb); // mg/h, Rate of basolateral transpoters 
dxdt_A_baso = RA_baso; // mg, Amount of basolateral transpoters
double RA_apical = (Vmax_apical*Cfil)/(Km_apical + Cfil); // mg/h, Rate of apical transpoter 
dxdt_A_apical = RA_apical; // mg, Amount of apical transpoter
double Rdif = Kdif*(CKb - CPTC); // mg/h, Rate of diffusion from into the PTC
dxdt_Adif = Rdif; // mg, Amount moved via glomarular filtration
double RAefflux = Kefflux*APTC; // mg/h, Rate of efflux clearance rate from PTC to blood
dxdt_Aefflux = RAefflux;  //  mg, Amount of efflux clearance rate from PTC to blood
double RCI = CA*GFR*Free; //  mg/h, Rate of clerance (CL) to via glormerular filtration (GFR) 
dxdt_ACI = RCI; // mg, Amount of clearance via GFR

// {PFOS distribution in each compartment}
// PFOS in plasma
double RPlas_free = (QRest*CVRest*Free)+(QK*CVK*Free)+(QL*CVL*Free)-(QC*CA*Free)+RAefflux; // mg/h, Rate of change in the plasma
dxdt_APlas_free = RPlas_free; // mg, Amount in the plasma
dxdt_AUCCA_free = CA_free; // mg*h/L, Area under curve of PFOS in liver compartment

// Proximal Tubule Cells (PTC)
double RPTC = Rdif + RA_apical + RA_baso - RAefflux; // mg/h, Rate of change in PTC 
dxdt_APTC = RPTC; // mg, Amount moved in PTC
dxdt_AUCCPTC = CPTC; // mg*h/L, Area under curve of PFOS in the compartment of PTC

// Proximal Tubule Lumen/ Filtrate (Fil)
double Rfil = CA*GFR*Free - RA_apical - AFil*Kurine; // mg/h, Rate of change in Fil
dxdt_AFil = Rfil; // mg, Amount moved in Fil
dxdt_AUCfil = Cfil; // mg*h/L, Area under curve of PFOS in the compartment of Fil

// Urine elimination
double Rurine = Kurine*AFil; // mg/h, Rate of change in urine
dxdt_Aurine = Rurine; // mg, Amount in urine
//double percentOD_in_urine = (Aurine/Odose)*100; //Percent of oral dose in the urine  

// Kidney compartment
double RKb = QK*(CA-CVK)*Free - CA*GFR*Free - Rdif - RA_baso; // mg/h, Rate of change in Kidney caomprtment
dxdt_AKb = RKb; // mg, Amount in kidney compartment
dxdt_AUCKb = CK; // mg*h/L, Area under curve of PFOS in the Kidney compartment

// PFOA in the compartment of rest of body, flow-limited model
double RRest = QRest*(CA-CVRest)*Free; // mg/h, Rate of change in rest of body
dxdt_ARest = RRest; // mg, Amount in rest of body 
dxdt_AUCCRest = CRest; // mg*h/L, Area under curve of PFOS in the compartment of rest of body
    
// Gastrointestinal (GI) tract
// Stomach compartment
double RST = - K0*AST - GE*AST; // mg/h, Rate of chagne in Stomach caomprtment
dxdt_AST = RST; // mg, Amount in Stomach
double RabsST = K0*AST; //  mg/h, Rate of absorption in the stomach
dxdt_AabsST = RabsST; // mg, Amount absorbed in the stomach

// Small intestine compartment
double RSI = GE*AST - Kabs*ASI - Kunabs*ASI; // mg/h, Rate of chagne in Small intestine caomprtment
dxdt_ASI = RSI; // mg, Amount in Small intestine
double RabsSI = Kabs*ASI; // mg/h, Rate of absorption in the Small intestine
dxdt_AabsSI = RabsSI; // mg, Amount absorbed in the Small intestine
double Total_oral_uptake = AabsSI + AabsST; // mg, Total oral uptake in the GI

// Biliary excretion
double Abile = Kbile*AL; // mg, Amount of PFOS in bile excretion
double amount_per_gram_liver = (AL/ML)*1000; // ug/g, Amount of PFOS in liver per gram liver

// Feces compartment
double Rfeces = Kbile*AL + Kunabs*ASI; // mg/h, Rate of change in feces compartment
dxdt_Afeces = Rfeces; // mg, Amount of the feces compartment
//double percentOD_in_feces = (Afeces/Odose)*100; Percent of the oral dose in the feces

// PFOS in liver compartment, flow-limited model
double RL = QL*(CA-CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST; // mg/h, Rate of chagne in liver caomprtment
dxdt_AL = RL; // mg, Amount in liver compartment
dxdt_AUCCL = CL; // mg*h/L, Area under curve of PFOS in liver compartment

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
