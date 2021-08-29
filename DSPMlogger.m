function DSPMlogger(fname)
%%% DSPM Logger %%%
load([fname, '.mat']); 
fid = fopen([fname '.txt'],'wt'); 

brk = '***********************';

fprintf(fid,'%s \n',brk);
fprintf(fid,'\n'); % print blank line
fprintf(fid,'%s \n',fname);
date = sprintf('date & time = %s',datestr(clock,21));
fprintf(fid,'%s  \n',date);
fprintf(fid,'%s \n',brk); % print break
fprintf(fid,'\n'); % print blank line

fprintf(fid,'%s \n','ACOUSTIC PARAMETERS:');
fprintf(fid,'%s %d %s \n','frequency f = ',f,'[Hz]');
fprintf(fid,'%s %d %s \n','wave number k = ',k,'[m^-1]');
fprintf(fid,'%s %d %s \n','angular frequency = ',omega,'[s^-1]');
fprintf(fid,'%s %d %s \n','wavelength = ',wl,'[m]');
fprintf(fid,'%s %d %s \n','transducer velocity amplitude = ',v0,'[m/s]');
fprintf(fid,'\n'); % print blank line

fprintf(fid,'%s  \n','MEDIUM PARAMETERS:');
fprintf(fid,'%s %d %s  \n','density = ',rho,'[kg.m^-3]');
fprintf(fid,'%s %d %s  \n','sound velocity c = ',c,'[m/s]');
fprintf(fid,'\n'); % print blank line

if dropOn
fprintf(fid,'%s  \n','DROPLET PARAMETERS:');
fprintf(fid,'%s %d %s  \n','droplet density = ',rho_f,'[kg.m^-3]');
fprintf(fid,'%s %d %s  \n','droplet sound velocity c = ',c_f,'[m/s]');
fprintf(fid,'%s %d %s  \n','droplet radius rd = ',rd,'[m]');
fprintf(fid,'%s %d %s  \n','droplet position = ',z0,'[m]');
fprintf(fid,'%s %d  \n','Is Rigid? = ',isRigid);
fprintf(fid,'\n'); % print blank line
end

fprintf(fid,'%s  \n','GEOMETRY PARAMETERS:');
fprintf(fid,'%s %d %s  \n','source radius = ',rs,'[m]');
fprintf(fid,'%s %d %s  \n','pitch distance pS = ',pS,'[m]');
fprintf(fid,'%s %d %s  \n','transducer radius = ',R0,'[m]');
fprintf(fid,'%s %d %s  \n','transducer position = ',L0,'[m]');
fprintf(fid,'%s %d %s  \n','reflector radius = ',R_ref,'[m]');
fprintf(fid,'%s %d %s  \n','reflector position = ',R_ref,'[m]');
fprintf(fid,'\n'); % print blank line

fprintf(fid,'%s \n',brk); % print break
fprintf(fid,'\n'); % print blank line
fprintf(fid,'%s  \n','DESCRITIZATION:');
fprintf(fid,'%s %s \n','descritization type = ',des_typ);
fprintf(fid,'%s %d \n','number of transducer points = ',numTR);
fprintf(fid,'%s %d \n','number of reflector points = ',numR);
fprintf(fid,'%s %d \n','number of drop points = ',numD);

fprintf(fid,'%s \n',brk); % print break
fprintf(fid,'\n'); % print blank line
fprintf(fid,'%s  \n','DSPM RESULTS:');
fprintf(fid,'%s %d \n','scale factor = ',scale_factor);
fprintf(fid,'%s %d  %s\n','antinodal pressure = ',abs(antiNodalPressure),'[Pa]');
fprintf(fid,'%s %d  %s\n','droplet BC residue = ',err,'[m/s]');

fprintf(fid,'%s \n',brk); % print break
fprintf(fid,'\n'); % print blank line
fprintf(fid,'%s  \n','RADIATION PRESSURE RESULTS:');
fprintf(fid,'%s %d  %s\n','max radiation pressure = ',max(p_rad),'[Pa]');
fprintf(fid,'%s %d  %s\n','min radiation pressure = ',min(p_rad),'[Pa]');
fprintf(fid,'%s %d  %s\n','X-Projection = ',mean(PX),'[Pa]');
fprintf(fid,'%s %d  %s\n','Y-Projection = ',mean(PY),'[Pa]');
fprintf(fid,'%s %d  %s\n','Z-Projection = ',mean(PZ),'[Pa]');
fprintf(fid,'%s  \n','RADIATION PRESSURE RESULTS (ABSOLUTE VALUES):');
fprintf(fid,'%s %d  %s\n','max radiation pressure = ',abs(max(p_rad)),'[Pa]');
fprintf(fid,'%s %d  %s\n','min radiation pressure = ',abs(min(p_rad)),'[Pa]');
fprintf(fid,'%s %d  %s\n','X-Projection = ',abs(mean(PX)),'[Pa]');
fprintf(fid,'%s %d  %s\n','Y-Projection = ',abs(mean(PY)),'[Pa]');
fprintf(fid,'%s %d  %s\n','Z-Projection = ',abs(mean(PZ)),'[Pa]');

fprintf(fid,'%s \n',brk); % print break
fprintf(fid,'%s \n',brk); % print break
fprintf(fid,'%s \n',brk); % print break
fclose(fid); 

end
