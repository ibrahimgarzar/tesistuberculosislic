close all
clear all
%%

% Fß5 (k10) y Fß7 (k14)

%umbral=([0.01,0.4,0.8,1,5.8,16.6,22.6,28.6,34.6,40.6,46.6,52.6,58.6,64.6,69.4,79]);
%umbral2=([1.00E-11,1.00E-10,2.20E-10,3.10E-10,4.30E-10,7.00E-10,8.50E-10,1.00E-09,3.00E-09,5.50E-09,8.00E-09,1.00E-08,3.50E-08,6.00E-08,8.00E-08,1.00E-07]);

%10000 iteraciones
%umbral=([0.01:0.01:0.1,0.1:0.1:1,1:1:80]); %ptf
%umbral2=([1e-11:1e-11:1e-10,1e-10:2e-11:1e-9,1e-9:4e-10:1e-8,1e-8:4e-9:0.9e-7]); %delta

%prueba de sanidad con delta fijo
%umbral=([0.00001:1:19,19.1664743109897,20:1:80]);
%umbral2=(1.68164470296576e-09);

%10000 iteraciones con el mismo step
umbral=(0.00119:1.7000e-04:0.0181); % ß1
umbral2=(0.001:(0.37-0.001)/100:0.37); % ß5

matriz= NaN(length(umbral2),length(umbral));

for ii = 1:length(umbral2)

    for i = 1:length(umbral)

k2=1.68164470296576e-09;
k4=1677.71529516028;
k5=3377.90445527686;
k6=8881221.80265334;
k9=4.89887206065021e-10;
k10=umbral2(ii); %0.0968962836928783;
k13=0.000208989059062582;
k14=168.675283030697;
k1=umbral(i); %0.00217688531461890;
k3=0.0116821323710160;
k8=2.67341167302623;
k12=1.46921454819336;
k15=23979612.7030956;
k16=19.1664743109897;
k17=9.23387740720957;


  syms M_t Mf_t T_t Tf_t

dydt1= 0==Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;
dydt2= 0==M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
dydt3= 0==k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
dydt4= 0==k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;

equations = [dydt1 dydt2 dydt3 dydt4];
vars=[M_t Mf_t T_t Tf_t];

range = [NaN NaN; NaN NaN;NaN NaN; NaN NaN];
sol = vpasolve(equations, vars, range);


dM=Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;
dMf=M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)  -Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
dT=k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
dTf=k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;

J=jacobian([dM dMf dT dTf], vars);

stable_positive_real_solution_matrix=[];
counter_stable_positive_real_solution=0;


for sol_num=1:1:length(sol.M_t);

if isreal([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)]);
     if sum(double(([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])>=0))==4;
            Jeval=subs(J, vars, [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)]);
        eigenvals=eig(Jeval);  
           if (sum(double(eigenvals)<0))==4;
            stable_positive_real_solution_matrix=[stable_positive_real_solution_matrix; 
          double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
           counter_stable_positive_real_solution=counter_stable_positive_real_solution+1;
            
     end
     end
     end
       
           [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)];
         
end

if counter_stable_positive_real_solution > 0
 MtotS=sum(stable_positive_real_solution_matrix(:,[1,2]),2);
 TtotS=sum(stable_positive_real_solution_matrix(:,[3,4]),2);
end

 if counter_stable_positive_real_solution == 2
matriz(ii,i)=2;
 end

 if counter_stable_positive_real_solution == 1
if MtotS>TtotS  
   matriz(ii,i)=3;
else
   matriz(ii,i)=1;
end
 end

if counter_stable_positive_real_solution == 0
matriz(ii,i)=4;
 end

    end
end
