function [Estado_final,Xx] = Solve(indicador,indicador_funcion)
% Total population
N_total = 1;
% Initial conditions for compartments
S_0 = 0.70 * N_total;
E_0 = 0.25 * N_total;
I_0 = 0.02 * N_total;
A_0 = 0.03 * N_total;
R_0 = 0.0000 * N_total;
W_0 = 1 / 100000;
% Control bounds
b_1 = 0.2;
b_2 = 0.2;
a_1 = 0.0;
a_2 = 0.0;
% Discretization parameters
N_res = 100;
M_res = 100;
% Cost function parameters
c1 = 1;
c2 = 1;
% Model parameters
m_p = 0.0018; % Mortality rate
beta_p = 0.1; % Transmission rate of infected individuals
kappa = 0.5; % Asymptomatic proportion relative to infected
beta_W = 0.8; % Proportion of asymptomatic contribution
n_p = m_p;
delta_p = 0.5; % Infection proportion in asymptomatic
w_p = 0.1923; % Incubation period
w_pprima = 0.1923; % Latent period (less than incubation)
gamma_p = 0.1724; % Recovery rate for infected
gamma_pprima = 0.1724; % Recovery rate for asymptomatic
mu_p = 0.1; % Shedding coefficient of infected individuals
c = 0.5; % Coefficient for asymptomatic shedding
epsilon = 0.1; % Virus lifespan in the market
% Transmission coefficients
b_p = beta_p * N_total;
b_W = (mu_p * beta_W * N_total) / epsilon;
% Time discretization
T = 10;
h = T / N_res;
Xx = 0:h:T;
% Initialize state variables
for i = 1:M_res+1
    S(i) = S_0;
    E(i) = E_0;
    I(i) = I_0;
    A(i) = A_0;
    R(i) = R_0;
    W(i) = W_0;
    u1(i) = 0;
    u2(i) = 0;
end
% Initialize adjoint variables
for i=M_res+1+N_res:M_res+1+N_res+M_res
lamda_1(i)=1;
lamda_2(i)=0;
lamda_3(i)=1;
lamda_4(i)=0;
lamda_5(i)=-1;
lamda_6(i)=0;
end
% Forward simulation of state equations and backward simulation of adjoint equations
for j=M_res+1:M_res+1+N_res-1
% Forward equations for state variables
S(j+1)=S(j)+h*(n_p-m_p*S(j)-b_p*S(j)*(I(j)+kappa*A(j))-b_W*S(j)*W(j)-u1(j));
E(j+1)=E(j)+h*(b_p*S(j)*(I(j)+kappa*A(j))+b_W*S(j)*W(j)-(1-delta_p)*w_p*E(j)-delta_p*w_pprima*E(j)-m_p*E(j));
I(j+1)=I(j)+h*((1-delta_p)*w_p*E(j)-(gamma_p+m_p)*I(j)-u2(j));
A(j+1)=A(j)+h*(delta_p*w_pprima*E(j)-(gamma_pprima+m_p)*A(j));
R(j+1)=R(j)+h*(gamma_p*I(j)+gamma_pprima*A(j)-m_p*R(j));
W(j+1)=W(j)+h*(epsilon*(I(j)+c*A(j)-W(j)));
% Backward equations for costate variables
lamda_1(2*M_res+2+N_res-j-1)=lamda_1(2*M_res+2+N_res-j)-h*(-1+m_p*lamda_1(2*M_res+2+N_res-j)+b_p*lamda_1(2*M_res+2+N_res-j)*(I(j+1)+kappa*A(j+1))+b_W*lamda_1(2*M_res+2+N_res-j)*W(j+1)+u1(j)*lamda_1(2*M_res+2+N_res-j)-lamda_2(2*M_res+2+N_res-j)*b_p*(I(j+1)+kappa*A(j+1))-lamda_2(2*M_res+2+N_res-j)*b_W*W(j+1)  );
lamda_2(2*M_res+2+N_res-j-1)=lamda_2(2*M_res+2+N_res-j)-h*(-lamda_2(2*M_res+2+N_res-j)*(delta_p-1)*w_p+lamda_2(2*M_res+2+N_res-j)*delta_p*w_pprima+m_p*lamda_2(2*M_res+2+N_res-j)-lamda_3(2*M_res+2+N_res-j)*(1-delta_p)*w_p-lamda_4(2*M_res+2+N_res-j)*delta_p*w_pprima );
lamda_3(2*M_res+2+N_res-j-1)=lamda_3(2*M_res+2+N_res-j)-h*(-1+lamda_1(2*M_res+2+N_res-j)*b_p*S(j+1)-lamda_2(2*M_res+2+N_res-j)*b_p*S(j+1)+lamda_3(2*M_res+2+N_res-j)*(gamma_p+m_p)+lamda_3(2*M_res+2+N_res-j)*u2(j)-lamda_5(2*M_res+2+N_res-j)*gamma_p-epsilon*lamda_6(2*M_res+2+N_res-j));
lamda_4(2*M_res+2+N_res-j-1)=lamda_4(2*M_res+2+N_res-j)-h*(lamda_1(2*M_res+2+N_res-j)*b_p*S(j+1)*kappa-lamda_2(2*M_res+2+N_res-j)*b_p*S(j+1)*kappa+lamda_4(2*M_res+2+N_res-j)*(gamma_pprima+m_p)-lamda_5(2*M_res+2+N_res-j)*gamma_pprima-lamda_6(2*M_res+2+N_res-j)*epsilon*c);
lamda_5(2*M_res+2+N_res-j-1)=lamda_5(2*M_res+2+N_res-j)-h*(1+lamda_5(2*M_res+2+N_res-j)*m_p);
lamda_6(2*M_res+2+N_res-j-1)=lamda_6(2*M_res+2+N_res-j)-h*(lamda_1(2*M_res+2+N_res-j)*b_W*S(j+1)-lamda_2(2*M_res+2+N_res-j)*b_W*S(j+1)+lamda_6(2*M_res+2+N_res-j)*epsilon);
% Control characterization
if indicador==1
u1(j+1)=0;
u2(j+1)=0;
end
if indicador==2
u1(j+1)=0;
u2(j+1)=min(max(a_2,(I(j+1)*lamda_3(2*M_res+2+N_res-j))/(2*c2)),b_2);
end
if indicador==3
u1(j+1)=min(max(a_1,(S(j+1)*lamda_1(2*M_res+2+N_res-j))/(2*c1)),b_1);
u2(j+1)=0;
end
if indicador==4
u1(j+1)=min(max(a_1,(S(j+1)*lamda_1(2*M_res+2+N_res-j))/(2*c1)),b_1);
u2(j+1)=min(max(a_2,(I(j+1)*lamda_3(2*M_res+2+N_res-j))/(2*c2)),b_2);
end
end
S_f=S(M_res+1:M_res+2+N_res-1);
E_f=E(M_res+1:M_res+2+N_res-1);
I_f=I(M_res+1:M_res+2+N_res-1);
A_f=A(M_res+1:M_res+2+N_res-1);
R_f=R(M_res+1:M_res+2+N_res-1);
W_f=W(M_res+1:M_res+2+N_res-1);
u1_f=u1(M_res+1:M_res+2+N_res-1);
u2_f=u2(M_res+1:M_res+2+N_res-1);
%Save and plot state and control functions
datos = [S_f, E_f, I_f, A_f, R_f, W_f, u1_f, u2_f];
if indicador_funcion==1
    Estado_final=S_f;
end
if indicador_funcion==2
    Estado_final=E_f;
end
if indicador_funcion==3
    Estado_final=I_f;
end
if indicador_funcion==4
    Estado_final=A_f;
end
if indicador_funcion==5
    Estado_final=R_f;
end
if indicador_funcion==6
    Estado_final=W_f;
end
if indicador_funcion==7
    Estado_final=u1_f;
end
if indicador_funcion==8
    Estado_final=u2_f;
end
end