V1 = ohmFlux(:,1); h1 = 'V1';
V3 = ohmFlux(:,3); h3 = 'V3';
V4 = ohmFlux(:,4); h4 = 'V4';
V5 = ohmFlux(:,5); h5 = 'V5';
V6 = ohmFlux(:,6); h6  = 'V6';
V7 = ohmFlux(:,7); h7 = 'V7';
h8 = 'V_source';
e = 1.60217662E-19; % C
hbar = 1.054571800E-34; % Js
B_exp = B * hbar/e * 10^16; % T


fid = fopen('data.txt','w');
fprintf(fid, [ 'B (T)' ' ' h1 ' ' h3 ' ' h4 ' ' h5 ' ' h6 ' ' h7 ' ' h8 '\r\n']);
fprintf(fid, '%f %f %f %f %f %f %f %f \r\n',[B_exp', V1, V3, V4, V5, V6, V7, V_source]');
fclose(fid);true