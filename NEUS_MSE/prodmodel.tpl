// SCHAEFER MODEL FOR MS_PROD OUTPUT
// SINGLE SPECIES
// GAVIN FAY
// 9 Feb 2012

DATA_SECTION

 init_int Nsp
 init_int rFase
 init_vector r_init(1,Nsp)
 init_int KFase
 init_vector K_init(1,Nsp)
 init_int zFase
 init_vector z_init(1,Nsp)
 init_int thetaFase
 init_vector theta_init(1,Nsp)
 init_int Fyear
 //init_int Syear
 //init_int yrBreak1; // first year of change in catchability
 //init_int yrBreak2; // second year of change in catchability
 init_int Lyear
 init_matrix ObsCat(Fyear,Lyear,1,Nsp)
 init_int NBio
 init_int Nsurvey
 !!int ncol=Nsp+2;
 init_matrix ObsBio(1,NBio,1,ncol)
 //init_matrix ObsCV(1,NBio,1,ncol)
 !!cout << Fyear << " " << ObsCat(Fyear) << endl;
 !!cout << ObsBio(1) << endl;

 //init_vector q1_init(1,Nsp)
 //init_int qAdjFase // q adjustment estimation phase
 //init_vector qAdj1_init(1,Nsp)
 //init_vector qAdj2_init(1,Nsp)

 matrix logcpue(1,NBio,1,Nsp)

 int i
 int j
 int iyr
 int iobs


PARAMETER_SECTION

 init_number dummy(-1);

 //init_bounded_vector log_r(1,Nsp,-10.1,0.4,rFase)
 init_vector log_r(1,Nsp,rFase) 
 //init_bounded_vector log_K(1,Nsp,-10,16.52356,KFase)
 init_vector log_K(1,Nsp,KFase)
 //init_bounded_vector log_z(1,Nsp,0.00001,10,zFase)
 init_vector log_z(1,Nsp,zFase)
// init_bounded_vector logit_theta(1,Nsp,-10,0.,thetaFase)
 init_vector logit_theta(1,Nsp,thetaFase)

// init_vector log_q1(1,Nsp)
// init_vector log_qAdj1(1,Nsp,qAdjFase)
// init_vector log_qAdj2(1,Nsp,qAdjFase)

 number eps

 sdreport_vector r(1,Nsp)
 sdreport_vector K(1,Nsp)
 vector z(1,Nsp)
 sdreport_vector theta(1,Nsp) // initial biomass adjustmant as scalar of carrying capacity
 sdreport_matrix q(1,Nsp,1,Nsurvey) // baseline catchability
 //matrix qAdj(Fyear,Lyear,1,Nsp)
 //sdreport_vector qAdj1(1,Nsp) // catchability adjustment 1
 //sdreport_vector qAdj2(1,Nsp) // catchability adjustment 2
 sdreport_vector sigma(1,Nsp) 
 number gamma
 number temp
 vector n(1,Nsurvey)

 vector logcpue_pred(Fyear,Lyear)
 vector resid(Fyear,Lyear)
 matrix bio_pred(Fyear,Lyear,1,Nsp)

 sdreport_matrix Bio(Fyear,Lyear+1,1,Nsp)

 // report derived variables
 sdreport_vector MSY(1,Nsp);
 sdreport_vector Bmsy(1,Nsp);
 sdreport_matrix ratioBmsy(Fyear,Lyear+1,1,Nsp);

 objective_function_value objfun

PRELIMINARY_CALCS_SECTION

 log_r = log(r_init);
 log_K = log(K_init);
 log_z = log(z_init-1.);
 //logit_theta = log(elem_div(theta_init,(1.-theta_init)));
 logit_theta = log(theta_init);

// log_q1 = log(q1_init);
// //cout << log_q1 << endl;
// log_qAdj1 = log(qAdj1_init);
// log_qAdj2 = log(qAdj2_init);


PROCEDURE_SECTION

  // constant added to cpue predictions
  eps = 1.e-07;

 objfun += square(dummy);

 //parameter values
 r = mfexp(log_r);
 K = mfexp(log_K);
 z = 1.+eps+mfexp(log_z);
 theta = mfexp(logit_theta); //elem_div(mfexp(logit_theta),(1.+mfexp(logit_theta)));

 //q1 = mfexp(log_q1);
 //cout << q1 << endl;
 //qAdj1 = mfexp(log_qAdj1);
 //qAdj2 = mfexp(log_qAdj2);
 //cout << r << endl;
 //cout << K << endl;
 //cout << theta << endl;

 //get_biomass
 
 for (i=1;i<=Nsp;i++)
  {
   gamma = pow(z(i),(z(i)/(z(i)-1)))/(z(i)-1);

   Bio(Fyear,i) = theta(i)*K(i);
   //qAdj(1,i) = q1(i);

 for (iyr=Fyear+1;iyr<=Lyear+1;iyr++)
  {
   dvariable fpen1=0.;
   //dvariable temp = gamma*y(i)*Bio(iyr-1,i)-gamma*y(i)*K(i)*pow(Bio(iyr-1,i)/K(i),z(i));
   //Bio(iyr,i) = posfun(Bio(iyr-1,i)+temp,1.,fpen1);
   Bio(iyr,i) = posfun(Bio(iyr-1,i)*(1.+r(i)*(1.-pow((Bio(iyr-1,i)/K(i)),z(i)-1))/(z(i)-1)),1.,fpen1);
   dvariable sr = 1.-ObsCat(iyr-1,i)/Bio(iyr,i);
   dvariable kcat=ObsCat(iyr-1,i);
   objfun+=1000*fpen1;
   if(sr< 0.001)
    {
     dvariable fpen=0.;
     kcat=Bio(iyr,i)*(1.-posfun(sr,0.001,fpen));
     objfun+=10000*fpen;
     // cout << " kludge "<<iy <<" "<<kcat<<" "<<cat(iy)<<" "<<fpen<<endl;
    }
    Bio(iyr,i)-=kcat;
    bio_pred(iyr-1,i) = 0.5*(Bio(iyr-1,i)+Bio(iyr,i));
    //cout << i << " " << iyr << " " << Bio(iyr,i) << endl;

//   // vector of q adjustments 
//   if(iyr < (yrBreak1 - Fyear + 2)){
//        qAdj(iyr-1,i) = q1(i);
//   } else if (iyr > (yrBreak2 - Fyear + 1)){
//        qAdj(iyr-1,i) = q1(i)*qAdj1(i);
//   } else
//        qAdj(iyr-1,i) = q1(i)*qAdj2(i);
   
  }

  
  //cout << "bio_pred" << endl;
  //cout << bio_pred << endl;
  //cout << "qAdj" << endl;
  //cout << qAdj << endl;
 }
 

 //objective function value
 for (i=1;i<=Nsp;i++)
  { 
  q = 0.;
  n = 0;
  // MLE for catchability (Polacheck et al 1993 equation 11)
  for (j=1;j<=NBio;j++) {
   q(i,ObsBio(j,2)) +=  log(ObsBio(j,3)/(eps+bio_pred(ObsBio(j,1),i)));
     n(ObsBio(j,2)) += 1;
  }
  for (j=1;j<=Nsurvey;j++) q(i,j) = mfexp(q(i,j)/n(j));

  //q(i) = mfexp(sum(log(elem_div(column(ObsBio,2),eps+column(bio_pred,i)))/NBio));

  // predictions
  for (j=1;j<=NBio;j++)
   logcpue_pred(ObsBio(j,1)) = log(eps+bio_pred(ObsBio(j,1),i)) + log(q(i,ObsBio(j,2)));

  //cout << "test logcpue_pred" << endl;

  resid = 0.;
  // calculate model residuals
  for(j=1;j<=NBio;j++)
    resid(ObsBio(j,1)) = log(ObsBio(j,3))-logcpue_pred(ObsBio(j,1));


  //cout << "test resid" << endl;
  //cout << resid << endl;
  
  // MLE for observation error standard deviation (Polacheck et al 1993 equation 10c)
  sigma(i) = sqrt(norm2(resid)/NBio);
  cout << i << " " << q(i) << " " << sigma(i) << endl;
  objfun += NBio*log(sigma(i)) + NBio/2;

  }
  cout << objfun << endl;

  for (i=1;i<=Nsp;i++) {
//----------- Calculate MSY reference points (Gompertz only) ------------------
  MSY(i) = (r(i)*K(i))/(mfexp(1)*log_K(i));
  Bmsy(i) = K(i)/mfexp(1);

  for (iyr=Fyear;iyr<=Lyear+1;iyr++){
    ratioBmsy(iyr,i) = (eps+Bio(iyr,i))/Bmsy(i);
  }
 }
 
  
REPORT_SECTION


 report << objfun << endl; 
 report << r << endl;
 report << K << endl;
 report << z << endl;
 report << theta << endl;
 report << q << endl;
 report << sigma << endl;
 for (iyr=Fyear;iyr<=Lyear+1;iyr++)
  report << iyr << " " << Bio(iyr) << " " << ratioBmsy(iyr,1) << endl;
// for (iyr=Syear;iyr<=Lyear;iyr++)
//  report << iyr << " " << resid(iyr-Syear+Fyear) << endl;
 for (j=1;j<=NBio;j++)
  report << ObsBio(j,1) << " " << log(ObsBio(j,3)) << " " << logcpue_pred(ObsBio(j,1)) << endl;

 