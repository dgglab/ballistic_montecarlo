function runBsweepBS(savepath,fnameout,frmgrp,BSp,BSn,B,N_inject,p_ifbounce_then_scatter,...
    p_bounce_off_ohmic,L_scatter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% runBsweepBS.m
%
% Parsing function used for calling man script of simulation in a logically
% paralized fasion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of parallel pools
n_pools=2;
% Number of simulations per block. code is excecuted in blocks so that in
% the case of early termination, there are consecutive simulations and not
% random elements selected from the list.
N_perblock=18;


%array containing the radius of curvature
save([savepath, fnameout,'/frameData.mat'],'frmgrp','N_inject','B');

rng('shuffle') % Seed rng by time
r = rng;
seed = r.Seed;
save([savepath, fnameout,'/rngData.mat'],'seed');

B_N=length(B);

% % Start parallel pools for simulations
% if parpool('size')==n_pools
% else
%     if parpool('size')~=0
%         parpool close
%     end
%     eval(['matlabpool open local ',num2str(n_pools)]);
% end

%ppool = parpool('local',n_pools);
   
for r=1:N_perblock:B_N
    tic
    parfor n=r:min((r+N_perblock-1),B_N)
        rng(seed+n)
        if ~exist([savepath,fnameout,'/',fnameout,'_',num2str(n,'%05d'),'.mat'],'file')
            t=tic;
            B_cur=B(n);
            
            [BS,edgenum,arcdata]=caustics_simulation_function_BS(frmgrp,BSp,BSn,B_cur,N_inject,...
                p_ifbounce_then_scatter,p_bounce_off_ohmic,L_scatter);

            parsave([savepath,fnameout,'/',fnameout,'_',num2str(n,'%05d'),'classes.mat'],edgenum,frmgrp,BS);
            

            fid=fopen([savepath,fnameout,'/',fnameout,'_',num2str(n,'%05d'),'arcdata.txt'],'w+');
            fwrite(fid,arcdata,'double');
            fclose(fid);
            t2=toc(t);
            time = datestr(now,'HH:MM:SS.FFF');
            fprintf('\n%d/%d:  %0.1f  %s\n\n',n,B_N,t2, time)
        end
    end
    toc
end

%delete(ppool)
 
