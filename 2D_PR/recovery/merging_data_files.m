% merging two data file
clear; 

load 'data_exp_gpu_first_round';

Temp_M1 = M1;
Temp_M2 = M2;
Temp_M3 = M3;
Temp_m = m_eff(m_eff>0);

load 'data_exp_gpu_second_round';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];


load 'data_exp_gpu_third_round';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];

load 'data_exp_gpu_fourth_round';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];


load 'data_exp_gpu_fifth_round';
M1 = M1 + Temp_M1;
M2 = M2 + Temp_M2;
M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
m_eff = [m_eff; Temp_m];

