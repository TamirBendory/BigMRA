% merging few data files
clear; 

load 'data_exp_gpu_1500';

Temp_M1 = M1;
Temp_M2 = M2;
Temp_M3 = M3;
Temp_m = m_eff(m_eff>0);

load 'data_exp_gpu_900';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];


load 'data_exp_gpu_2000';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];

load 'data_exp_gpu_2000_2';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];


load 'data_exp_gpu_600';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];


load 'data_exp_chai_2800';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];

load 'data_exp_chai';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];


load 'data_exp_gpu';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];

load 'data_exp_latte';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];

load 'data_exp_latte_7900';
M1 = M1 + Temp_M1;
M2 = M2 + Temp_M2;
M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
m_eff = [m_eff; Temp_m];
