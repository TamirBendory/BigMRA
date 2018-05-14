% merging few data files
clear; 

load 'data_files/data_exp_9200_latte';

Temp_M1 = M1;
Temp_M2 = M2;
Temp_M3 = M3;
Temp_m = m_eff(m_eff>0);

load 'data_files/data_exp_2200_polar';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];

load 'data_files/data_exp_2400_polar_v2';
Temp_M1 = M1 + Temp_M1;
Temp_M2 = M2 + Temp_M2;
Temp_M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
Temp_m = [m_eff; Temp_m];


load 'data_files/data_exp_4200_chai';
M1 = M1 + Temp_M1;
M2 = M2 + Temp_M2;
M3 = M3 + Temp_M3;
m_eff = m_eff(m_eff>0);
m_eff = [m_eff; Temp_m];