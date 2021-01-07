function output(st,t,population,deads,deathcauses,sex_index,offpsringproduced,parents,mattime_index)
    filenamepop=strcat('output/pop_means',num2str(st),'.txt');
    filenamesd=strcat('output/pop_sds',num2str(st),'.txt');
    filenamepop_m=strcat('output/pop_means_m',num2str(st),'.txt');
    filenamepop_f=strcat('output/pop_means_f',num2str(st),'.txt');
    filenamedeaths=strcat('output/d',num2str(st),'.txt');
    filenamedeaths_m=strcat('output/dm',num2str(st),'.txt');
    filenamedeaths_f=strcat('output/df',num2str(st),'.txt');
    filenameinfo=strcat('output/info',num2str(st),'.txt');
    filenameparents=strcat('output/parents',num2str(st),'.txt');
    filenamematuremean_m=strcat('output/pop_means_adult_m',num2str(st),'.txt');
    filenamematuremean_f=strcat('output/pop_means_adult_f',num2str(st),'.txt');

    % means of all the individuals who are alive
    pop_full=population(population(1:end,1)~=0,1:end);
    means=mean(pop_full,1);
    sds=std(pop_full,1);
    data=horzcat(t,means);
    dlmwrite(filenamepop,data,'-append');
    data=horzcat(t,sds);
    dlmwrite(filenamesd,data,'-append');
        
    % by sex
    males=pop_full(pop_full(:,sex_index)==0,:);
    females=pop_full(pop_full(:,sex_index)==1,:);
    means_m=mean(males,1);
    means_f=mean(females,1);
    data_m=horzcat(t,means_m);
    dlmwrite(filenamepop_m,data_m,'-append');
    data_f=horzcat(t,means_f);
    dlmwrite(filenamepop_f,data_f,'-append');
    
    means_mat_m=mean(males(males(:,mattime_index)>0,:),1);
    means_mat_f=mean(females(females(:,mattime_index)>0,:),1);
    data_mat_m=horzcat(t,means_mat_m);
    dlmwrite(filenamematuremean_m,data_mat_m,'-append');
    data_mat_f=horzcat(t,means_mat_f);
    dlmwrite(filenamematuremean_f,data_mat_f,'-append');

    % means of all those who are dead, also by sex
%     means_d=mean(deads,1);
%     means_d_m=mean(deads(deads(:,sex_index)==0,:),1);
%     means_d_f=mean(deads(deads(:,sex_index)==1,:),1);
%     data2=horzcat(t,means_d,deathcauses);
%     dlmwrite(filenamedeaths,data2,'-append');
%     data2_m=horzcat(t,means_d_m);
%     dlmwrite(filenamedeaths_m,data2_m,'-append');
%     data2_f=horzcat(t,means_d_f);
%     dlmwrite(filenamedeaths_f,data2_f,'-append');
      timestep=t.*ones(size(deads,1),1);
      data_death=horzcat(timestep,deads);
      dlmwrite(filenamedeaths,data_death,'-append');
    if isempty(deads)
        info_temp= horzcat(t,size(pop_full,1),size(pop_full(pop_full(:,sex_index)==0,:),1),...
            size(pop_full(pop_full(:,sex_index)==1,:),1),0,0,...
            0,offpsringproduced);
    else
        info_temp= horzcat(t,size(pop_full,1),size(pop_full(pop_full(:,sex_index)==0,:),1),...
            size(pop_full(pop_full(:,sex_index)==1,:),1),size(deads,1),size(deads(deads(:,sex_index)==0,:),1),...
            size(deads(deads(:,sex_index)==1,:),1),offpsringproduced);
    end
    dlmwrite(filenameinfo,info_temp,'-append');
    dlmwrite(filenameparents,[ones(size(parents,1),1).*t,parents],'-append');
end