clear all;
diary paper_CWM_Temp

% INITIAL CONDITIONS
moisto = 0.0; %moisture content
ash=5*(1-moisto)/100;	%0.05;
dafo=1-moisto-ash;  %0.855;	%0.95;
rhodry = 750;
rhowet = rhodry/(1-moisto); 	%833.33   %750;
rhoo = rhowet*1.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
width =0.225; %Only oven length is considered, such that WIDTH = 0.45/2 = 0.225
nseg = 21; 
dx = width/nseg; %0.010714286 
delt = 600;       %GOOD solution
         t=13.8333*3600;        %total time
nt = t/delt;      %number of timesteps
doo= dx*dx/delt;
v0 = dx/delt;

wall = 0.05; %0.1
M = 1.5+wall/dx;  M=fix(M);
N = M+nseg-1;     N=fix(N);

Wcond = 2.3; Wrho = 2800; Wcp = 1200;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inlet_temp = 1150*ones(N,1); %1050; %deg

CH4	= 0.06943; C2H6	= 0.01166; CO	= 0.02272; CO2	= 0.01065; TAR	= 0.109658059;
H2O	= 0.048970557; NH3	= 0.012348159; H2S	= 0.015721205; Hydrogen	= 0.01440518;
car=0.838;  hy =0.053;  ox =0.071;  nit=0.018;  sul=0.02;

moist=moisto*ones(N,1);
daf=dafo*ones(N,1);   %daf=dafo;
atmassfrmla=1/(car/12+hy+ox/16+nit/14+sul/32);
atmass=atmassfrmla*ones(N,1);
carr=car*ones(N,1);
hyy=hy*ones(N,1);
oxx=ox*ones(N,1);
nitt=nit*ones(N,1);
sull=sul*ones(N,1);

H2Ozz = CO2zz = C2H6zz = NH3zz = H2Szz = TARzz = CH4zz = COzz = H2zz = zeros(N,1);
e0 = 1 - rhoo * (0.000668 * dafo + ash / 3000);

for i = 1:M-1
    Q(i) = 0;
    K(i) = Wcond;
    CVOL(i) = Wrho * Wcp;
end

##for i = 1:M-1   
##    %T(i) = inlet_temp(1);		%inlet temp = 1323K
##    T(i) = 1150;
##end

for i = M:N      %for(i = M; i <= N; i++) 
    T(i) = 20;      %ATOmospheric temp = 20deg
end

XB = 0;
Lboil = 2.257E6;
Qinput = 0;
time = 0;

%% Outside time loop
probe_locations = [0.05; 0.10625; 0.1625; 0.21875; 0.275];
%desired_temperatures = 0:10:100;
%desired_temperatures = 0:0.5:101;
%desired_temperatures = 300:50:1300;
%desired_temperatures = 1050:1:1300;
desired_temperatures = 0:100:900;

recorded_times = zeros(length(desired_temperatures), length(probe_locations));
for z=1:nt 
  for i = M:N      %for(i = M; i <= N; i++) 
    TK=T(i)+273;  E0 = 183000;  RT=8.3144*TK;  
    %VOLATILE EVOLUTION MODEL  
    if T(i)>350
      FE_H2O = @(EA) exp(-((((EA)-E0)/48000).^8));    FE_CO2 = @(EA) exp(-((((EA)-E0)/78000).^4));
      FE_C2H6 = @(EA) exp(-((((EA)-E0)/61000).^4));   FE_NH3 = @(EA) exp(-((((EA)-E0)/106000).^4));
      FE_H2S = @(EA) exp(-((((EA)-E0)/114000).^4));   FE_TAR = @(EA) exp(-((((EA)-E0)/48000).^8));
      FE_CH4 = @(EA)exp(-((((EA)-E0)/110000).^2));    FE_CO =  @(EA)exp(-((((EA)-E0)/93000).^4));
      FE_H2 =  @(EA)exp(-((((EA)-E0)/165000).^4)); 
      
      A_PVM = @(EA) 1-exp(-1.3E13.*delt.*exp(-(EA)./RT));
      A_PVM_dmDt = @(EA) 1.3E13*exp((-EA./RT)-(1.3E13.*delt.*exp(-EA./RT)));     
            
      H2Ox = @(EA) -H2O*(-8/(48000^8)).*(((EA)-E0).^7).*FE_H2O(EA).*A_PVM(EA);    
      CO2x = @(EA) -CO2*(-4/(78000^4)).*(((EA)-E0).^3).*FE_CO2(EA).*A_PVM(EA);
      C2H6x = @(EA) -C2H6*(-4/(61000^4)).*(((EA)-E0).^3).*FE_C2H6(EA).*A_PVM(EA);
      NH3x = @(EA) -NH3*(-4/(106000^4)).*(((EA)-E0).^3).*FE_NH3(EA).*A_PVM(EA);
      H2Sx = @(EA) -H2S*(-4/(114000^4)).*(((EA)-E0).^3).*FE_H2S(EA).*A_PVM(EA);
      TARx = @(EA) -TAR*(-8/(48000^8)).*(((EA)-E0).^7).*FE_TAR(EA).*A_PVM(EA);
      CH4x = @(EA) -CH4*(-2/(110000^2)).*(((EA)-E0).^1).*FE_CH4(EA).*A_PVM(EA);
      COx = @(EA) -CO*(-4/(93000^4)).*(((EA)-E0).^3).*FE_CO(EA).*A_PVM(EA);
      H2x = @(EA) -Hydrogen*(-4/(165000^4)).*(((EA)-E0).^3).*FE_H2(EA).*A_PVM(EA);
    
      H2O_dm = @(EA) H2O*(-8/(48000^8)).*(((EA)-E0).^7).*FE_H2O(EA).*A_PVM_dmDt(EA);    
      CO2_dm = @(EA) CO2*(-4/(78000^4)).*(((EA)-E0).^3).*FE_CO2(EA).*A_PVM_dmDt(EA);
      C2H6_dm = @(EA) C2H6*(-4/(61000^4)).*(((EA)-E0).^3).*FE_C2H6(EA).*A_PVM_dmDt(EA);
      NH3_dm = @(EA) NH3*(-4/(106000^4)).*(((EA)-E0).^3).*FE_NH3(EA).*A_PVM_dmDt(EA);
      H2S_dm = @(EA) H2S*(-4/(114000^4)).*(((EA)-E0).^3).*FE_H2S(EA).*A_PVM_dmDt(EA);
      TAR_dm = @(EA) TAR*(-8/(48000^8)).*(((EA)-E0).^7).*FE_TAR(EA).*A_PVM_dmDt(EA);
      CH4_dm = @(EA) CH4*(-2/(110000^2)).*(((EA)-E0).^1).*FE_CH4(EA).*A_PVM_dmDt(EA);
      CO_dm = @(EA) CO*(-4/(93000^4)).*(((EA)-E0).^3).*FE_CO(EA).*A_PVM_dmDt(EA);
      H2_dm = @(EA) Hydrogen*(-4/(165000^4)).*(((EA)-E0).^3).*FE_H2(EA).*A_PVM_dmDt(EA);
      
      H2Oz = integral(H2Ox, E0, Inf);    CO2z = integral(CO2x, E0, Inf);    C2H6z = integral(C2H6x, E0, Inf);
      NH3z = integral(NH3x, E0, Inf);    H2Sz = integral(H2Sx, E0, Inf);    TARz = integral(TARx, E0, Inf);
      CH4z = integral(CH4x, E0, Inf);    COz = integral(COx, E0, Inf);    H2z = integral(H2x, E0, Inf);
      totalz = H2Oz + CO2z + C2H6z + NH3z + H2Sz + TARz + CH4z + COz + H2z;
      
      H2Ozz(i) = integral(H2O_dm, E0, Inf);    CO2zz(i) = integral(CO2_dm, E0, Inf);    C2H6zz(i) = integral(C2H6_dm, E0, Inf);
      NH3zz(i) = integral(NH3_dm, E0, Inf);    H2Szz(i) = integral(H2S_dm, E0, Inf);    TARzz(i) = integral(TAR_dm, E0, Inf);
      CH4zz(i) = integral(CH4_dm, E0, Inf);    COzz(i) = integral(CO_dm, E0, Inf);    H2zz(i) = integral(H2_dm, E0, Inf);
      Total_dmDt = H2Ozz(i) + CO2zz(i) + C2H6zz(i) + NH3zz(i) + H2Szz(i) + TARzz(i) + CH4zz(i) + COzz(i) + H2zz(i); 
      
      H2Oy = H2O - H2Oz;    CO2y = CO2 - CO2z;    C2H6y = C2H6 - C2H6z;
      NH3y = NH3 - NH3z;    H2Sy = H2S - H2Sz;    TARy = TAR - TARz;
      CH4y = CH4 - CH4z;    COy = CO - COz;       H2y = Hydrogen - H2z;

      Y=H2Oy+CO2y+C2H6y+NH3y+H2Sy+TARy+CH4y+COy+H2y+0.68443684;
      VM = 1- Y;    
      
      carr(i)=(0.6707481032+0.75*CH4y+0.8*C2H6y+0.4286*COy+0.2727*CO2y+0.85*TARy)/Y;
      hyy(i)= (0.00136887368+0.25*CH4y+0.2*C2H6y+0.082*TARy+H2y+0.1111*H2Oy+0.1765*NH3y+0.0588*H2Sy)/Y;
      oxx(i)= (0.00136887368+0.5714*COy+0.7273*CO2y+0.049*TARy+0.8889*H2Oy)/Y;
      nitt(i)=(0.0068443684+0.009*TARy+0.8235*NH3y)/Y;
      sull(i)=(0.00410662104+0.01*TARy+0.9412*H2Sy)/Y;    
      daf(i)=Y*dafo;   
      atmass(i)=1/((carr(i)/12)+(hyy(i))+(oxx(i)/16)+(nitt(i)/14)+(sull(i)/32));    
        %Y=real(Y);      car=real(car);      hy=real(hy);      ox=real(ox);      nit=real(nit);      sul=real(sul);
    end      
    
      %PHYSICAL PROPERTIES MODEL
      mass=daf(i)+ash+moist(i);
      rho=rhoo*mass;
      Csteam_equ(i) = 1810 + 0.66*T(i); %Csteam = 2000;        
      densdaf = 1/(((0.0053*carr(i))/12)+(hyy(i)*0.00577)+((oxx(i)*0.00346)/16)+((nitt(i)*0.0669)/14)+((sull(i)*0.0384)/32));
      dens=mass/(daf(i)/densdaf+ash/3000+moist(i)/1000);
      e=1-rho./dens;    %eshr=max(0, 1-(1-e) /(1-e0));
      effectivepor = 1- ((1-e).^(1/3));
      F1=exp(380/TK);
      E1=F1*(380/TK/(F1-1))^2;
      F2=exp(1800/TK);
      E2=F2*(1800/TK/(F2-1))^2;
      Cdaf=8314.4*(E1+2*E2)/atmass(i);
      Cash=754+.586*T(i);
      C=(daf(i)*Cdaf+ash*Cash+4186.8*moist(i))/mass;    
      CVOL(i) = C.*rho;        %J/kg.K * kg/m3 = J/K.m3

      %HEAT Calculation

      H2Oxhs =  -13468000+ (1810*T(i)) + (0.33*T(i).^2);     CO2xhs =  -8964000 + (900*T(i))  + (0.22*T(i).^2);	
      C2H6xhs =	-2866000 + (1920*T(i)) + (1.46*T(i).^2);     NH3xhs =	-2760000 + (2080*T(i)) + (0.84*T(i).^2); 
      H2Sxhs =  -630000  + (980*T(i))  + (0.24*T(i).^2);     TARxhs =  -1936000 + (1390*T(i)) + (1.83*T(i).^2);
      CH4xhs =  -4719000 + (2210*T(i)) + (1.56*T(i).^2);     COxhs =   -3971000 + (1020*T(i)) + (0.11*T(i).^2);
      H2xhs =   -0       + (14250*T(i))+ (0.53*T(i).^2); 
      
      %hsH2O = -daf*H2Oxhs*H2Ozz;    hsCO2 = -daf*CO2xhs*CO2zz;    hsC2H6 = -daf*C2H6xhs*C2H6zz;  
      %hsNH3 = -daf*NH3xhs*NH3zz;    hsH2S = -daf*H2Sxhs*H2Szz;    hsTAR = -daf*TARxhs*TARzz;
      %hsCH4 = -daf*CH4xhs*CH4zz;    hsCO = -daf*COxhs*COzz ;      hsH2 = -daf*H2xhs*H2zz;    
      %totalhs = hsH2O + hsCO2 + hsC2H6 + hsNH3 + hsH2S + hsTAR + hsCH4 + hsCO;      

      valH_C = 33620000;    valphi_C = 32760000;
      valH_H = 141930000;   valphi_H = 141770000;
      valH_O = -14530000;   valphi_O = 0;
      valH_N = 0;           valphi_N = 0;
      valH_S = 9420000;     valphi_S = 9260000;

      mu_C = 12; mu_H = 1; mu_O = 16; mu_N = 14; mu_S = 32;
      
      G1=exp(380/TK);    H1=380*(1/(G1-1));
      G2=exp(1800/TK);   H2=3600*(1/(G2-1));    
      f_T = (8314.4)*(H1+H2-156);          %J/kg

      hf_H2O_H = 0.1111*(valH_H - valphi_H + f_T/mu_H);    hf_H2O_O = 0.8889*(valH_O - valphi_O + f_T/mu_O);
      hf_CO2_C = 0.2727*(valH_C - valphi_C + f_T/mu_C);    hf_CO2_O = 0.7273*(valH_O - valphi_O + f_T/mu_O);
      hf_C2H6_C = 0.8*(valH_C - valphi_C + f_T/mu_C);      hf_C2H6_H = 0.2*(valH_H - valphi_H + f_T/mu_H);
      hf_NH3_N = 0.1765*(valH_N - valphi_N + f_T/mu_N);    hf_NH3_H = 0.8235*(valH_H - valphi_H + f_T/mu_H);
      hf_H2S_H = 0.0588*(valH_H - valphi_H + f_T/mu_H);    hf_H2S_S = 0.9412*(valH_S - valphi_S + f_T/mu_S);
      hf_TAR_C = 0.85*(valH_C - valphi_C + f_T/mu_C);      hf_TAR_H = 0.082*(valH_H - valphi_H + f_T/mu_H);    
      hf_TAR_O = 0.049*(valH_O - valphi_O + f_T/mu_O);     hf_TAR_N = 0.009*(valH_N - valphi_N + f_T/mu_N);   
      hf_TAR_S = 0.01*(valH_S - valphi_S + f_T/mu_S);
      hf_CH4_C = 0.75*(valH_C - valphi_C + f_T/mu_C);      hf_CH4_H = 0.25*(valH_H - valphi_H + f_T/mu_H);
      hf_CO_C = 0.4286*(valH_C - valphi_C + f_T/mu_C);     hf_CO_O = 0.5714*(valH_O - valphi_O + f_T/mu_O);
      hf_H2_H = 1.0*(valH_H - valphi_H + f_T/mu_H);
      
      hf_H2O = hf_H2O_H + hf_H2O_O; 
      hf_CO2 = hf_CO2_C + hf_CO2_O; 
      hf_C2H6 = hf_C2H6_C + hf_C2H6_H;      
      hf_NH3 = hf_NH3_N + hf_NH3_H;      
      hf_H2S = hf_H2S_H + hf_H2S_S;                 
      hf_TAR = hf_TAR_C + hf_TAR_H + hf_TAR_O + hf_TAR_N + hf_TAR_S;
      hf_CH4 = hf_CH4_C + hf_CH4_H;
      hf_CO = hf_CO_C + hf_CO_O;      
      hf_H2 = hf_H2_H; 
         
      dq_mj_H2O = daf(i)*(hf_H2O - H2Oxhs);    dq_mj_CO2 = daf(i)*(hf_CO2 - CO2xhs);    dq_mj_C2H6 = daf(i)*(hf_C2H6 - C2H6xhs);
      dq_mj_NH3 = daf(i)*(hf_NH3 - NH3xhs);    dq_mj_H2S = daf(i)*(hf_H2S - H2Sxhs);    dq_mj_TAR = daf(i)*(hf_TAR - TARxhs);
      dq_mj_CH4 = daf(i)*(hf_CH4 - CH4xhs);    dq_mj_CO = daf(i)*(hf_CO - COxhs);       dq_mj_H2 = daf(i)*(hf_H2 - H2xhs);
      
      %TOTAL HEAT dq_dt for overall COAL COMBUSTION //
      dq_dt_H2O = -dq_mj_H2O*H2Ozz(i);    dq_dt_CO2 = -dq_mj_CO2*CO2zz(i);    dq_dt_C2H6 = -dq_mj_C2H6*C2H6zz(i);
      dq_dt_NH3 = -dq_mj_NH3*NH3zz(i);    dq_dt_H2S = -dq_mj_H2S*H2Szz(i);    dq_dt_TAR = -dq_mj_TAR*TARzz(i);
      dq_dt_CH4 = -dq_mj_CH4*CH4zz(i);    dq_dt_CO = -dq_mj_CO*COzz(i);       dq_dt_H2 = -dq_mj_H2*H2zz(i);    
      
      %Integrating above equation and calculating mass of individual specie releasing heat //
      %dq_dt_H2O = dq_mj_H2O*H2Oz;    dq_dt_CO2 = dq_mj_CO2*CO2z;    dq_dt_C2H6 = dq_mj_C2H6*C2H6z;
      %dq_dt_NH3 = dq_mj_NH3*NH3z;    dq_dt_H2S = dq_mj_H2S*H2Sz;    dq_dt_TAR = dq_mj_TAR*TARz;
      %dq_dt_CH4 = dq_mj_CH4*CH4z;    dq_dt_CO = dq_mj_CO*COz;       dq_dt_H2 = dq_mj_H2*H2z;
      
      %W/kg
      Total_dq_dt = dq_dt_H2O + dq_dt_CO2 + dq_dt_C2H6 + dq_dt_NH3 + dq_dt_H2S + dq_dt_TAR + dq_dt_CH4 + dq_dt_CO + dq_dt_H2;
      Q(i) = Total_dq_dt* dafo * rhoo;  %W/m3
      % Q=real(Q);   
    
      K0 = ((densdaf/4511)^3.5)*sqrt(TK); 
      %if (e < 0.5)	%if (T < 773) %623
      if (e < e0)
        KparGas = (7.45E-5*TK);
        KparRad = (2.28E-10*(TK^3));
        Kadd = KparGas + KparRad;
        K1 = effectivepor/Kadd;
        K2 = (1-effectivepor)/K0;		
        K(i) = (0.6*moist(i)) + ((1-moist(i))/ (K1 + K2));

      else  %(totalpor(T) >= pR_)
        pR = 0.5;
        rad_Int = (pR*(1-e))/(1-pR);
        rad_Ext = 1-((1-e)/(1-pR));
        K3= (1-e)*K0;
        K4= e*(4.96E-4*TK);
        K5= rad_Int*(3.42E-11*(TK^3));
        K6= rad_Ext*(2.28E-9*(TK^3));
        K(i) =	K3 + K4 + K5 + K6;
      end
    end
    %endfor

    % TEMPERATURE HISTORY MODEL   
    % TDMA 
        L(1) = 0;
        G(1) = 1;    
        U(1) = 0;
        Z(1) = 1050;  %TF Inlet condition

      for i=2:N-1
          KL = 2/(1/K(i)+1/K(i-1));      %KL=real(KL); %K(i)=real(K(i)); K(i-1)=real(K(i-1)); K(i+1)=real(K(i+1));
          KU = 2/(1/K(i+1)+1/K(i));      %KU=real(KU);
          DIFFL = KL/CVOL(i); %K/CVOL(i);
          DIFFU = KU/CVOL(i); %K/CVOL(i);
          F_DIFFL = DIFFL/doo;
          F_DIFFU = DIFFU/doo;
          REACT = Q(i)/CVOL(i); %Q/CVOL;    %(W/m3)/(J/K.m3)   = (J/s.m3)/(J/K.m3) = K/s
          F_REACT = REACT*delt;
          L(i) = -F_DIFFL;
          G(i) = 1+F_DIFFL+F_DIFFU;    
          U(i) = -F_DIFFU;
          Z(i) = T(i)+F_REACT;    %Z(i) = T(i);
      end

    % Outlet condition - SYMMETRY  
        L(N) = -F_DIFFL;             %L(N) = -F_DIFFU;
        G(N) = 1+F_DIFFL;            %G(N) = 1+F_DIFFU;
        U(N) = 0;
        Z(N) = T(N)+ Q(N)*delt/CVOL(N); %F_REACT(N-1);
        %Z(N) = T(N-1)+ Q(N-1)*dt/CVOL(N-1);

      IB = IC = N+1;
      for i = N:-1:M
         if (T(i)==100 && moist(i)>0)
           IB = i;
         end
      end 

      for i = N:-1:IB+1 %IB+1 = IC
         if (T(i)<100)
           IC = i;
         end
      end 

      for i=IB:IC-1   %IB = IC-1
          L(i) = 0;
          G(i) = 1;    
          U(i) = 0;
          Z(i) = 100;
      end 

      % Forward  elimination
      for i=2:N
        FACTOR = L(i)/G(i-1);
        G(i)=G(i)-FACTOR*U(i-1);
        Z(i)=Z(i)-FACTOR*Z(i-1);
      end  
      % Backward Substitution
      for i = N-1:-1:1
         Z(i)=Z(i)-Z(i+1)*U(i)/G(i+1);
      end
          
      for i = 1:N %temperature calculation at the diagonal
        T(i)=Z(i)/G(i);
      end
      for i=M:N
        if (T(i)>100 && moist(i)>0)
          T(i)=100;
        end
      end
      
      if IB<N+1
        QB=K(IB-1)*(T(IB-1)-100)/dx;            
        QBMAX=Lboil*(dx-XB)*rhoo*moist(IB)/delt;  
        if QB>QBMAX
           QB=QBMAX;
           XB=0;
           moist(IB)=0;
        %end
        else XB=XB+(QB*delt)/(rhoo*moist(IB)*Lboil);
        end
        if IC<N+1
          QC=K(IC)*(100-T(IC))/dx;              
          QERR=QB-QC;                          
          for i=IC:N
             if QERR > 0  %if (!(QERR > 0)) break;
                QSENS=CVOL(i)*(100-T(i))*v0;
             if QSENS>QERR
               T(i)=T(i)+QERR/(CVOL(i)*v0);
               moist(i)=moist(i)+QERR/(Lboil*rhoo*v0);
               QERR=0;
             %end
             else
               T(i)=100;
               moist(i)=moist(i)+QSENS/(Lboil*rhoo*v0);
               QERR=QERR-QSENS;
             end
          end
        end
        else
          V=Csteam_equ(i)*QB/Lboil;
          for i = IB-1:-1:M
            %if (!(eshr)) break ;
            %if(eshr(i)==0)   
            if((T(i)>0)&&(T(i)<483))   
              %fprintf('This loop works\n');   
              VEL=V/CVOL(i);
              VELL=VEL/v0;
              T(i)=T(i)-VELL*(T(i-1)-T(i+1))/2;
            end  
          end
        end
      end

    if M > 2  
      Qinput = Qinput + max(0, (Wcond * (T(1) - T(M - 1)) / ((M - 2) * v0)));      
    end
    
    Heat_req = Qinput * 1E-6 / (rhoo * width);  %MJ/kg
    fprintf('%0.0f\n', [z*delt]); 
    %% Inside time loop
##    for i = M : N
##        x_this = (i - M) * dx + 0.05;
##        x_next = x_this + dx;
##       
##        % x_M --- x_M+1 --(p1)-- x_M+2 --- x_M+3 --(p2)-- x_M+4 ...
##        for j = 1 : length(probe_locations)
##            if x_this <= probe_locations(j) && x_next >= probe_locations(j)    %if x_this <= probe_locations(j) && x_next >= probe_locations(j)
##                for k = 1 : length(desired_temperatures)
##                    if (recorded_times(k, j) == 0)                       
##                      % constant interpolation    
##                      temperature_this = T(i);    %moist(i);                        
##                      if (temperature_this >= desired_temperatures(k))
##                          %fprintf('z\n');
##                          %fprintf('%0.5f\n',[z/3600]);                       
##                          recorded_times(k, j) = (z/3600)*delt;     %not works
##                      end
##                    
##                    end                
##                end
##            end
##        end
##    end

%plot(T)
      
end
diary off