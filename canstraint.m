function canstraint(curr_f,curr_m,targetsize,starttime,endtime,matFile,ind,st,evolveTraits,wbs,wcancer)

    % gets the optimised strategies
    % initiates a population with 50% males and 50% females
    % with some variation around the mean trait values 
    % evaluate individuals each time step, records deaths
    % reproduction to fill the empty spots

    % curr_f = the currency that the females use
    % curr_m = the currency that the males use
    % targetsize = target mature body size
    % starttime = 1 if a new population, defined differently if not
    % endtime = how many time steps
    % ind = index for the seed
    % st = which strategy
    % evolveTraits = 0, ontogeny-related traits will not evolve / 1, ontogeny-related traits will evolve
    % w_bs and w_cancer = weight for phenotype (how much shared vs. sex-specific)
    % between 0 and 1

    % random seed
    t1 = datetime('now','Format','dd-MMM-yyyy HH:mm:ss.SSS');
    s=second(t1,'secondofday');
    seed=s*ind; 
    rng(seed);
    mkdir('output');
    
    w_bs =[wbs, 1-wbs];
    w_cancer = [wcancer, 1-wcancer];
    
    if starttime > 1
        warning('off', 'MATLAB:MKDIR:DirectoryExists');
        mkdir('output');
        endtimeNew=endtime;
        load(matFile);    
        endtime=endtimeNew; % overriding the previous endtime
        clear endtimeNew;
    else    
        % try different values of these for robustness
        % init_var = variation for body size and other traits around the mean
        % that comes from the input 
        % also relevant for the mutation steps
        populationsize=1000;init_var=0.05;

        A = importdata('strategies_0.01_1-4_300.txt');
        % the strategies optimised 300 rounds at extmort 0.01 for 4 replicates
        % previously
        % strategy has 7 traits, each with 6 alleles
        % body size has 6 alleles
        % strategy: P-Q-H-A-S-T-X

        % HARDCODED, need to change these in case it is necessary
        extmort= A(1,4);cellmort=A(1,2);cancer_danger=A(1,3);
        mort=[cellmort,extmort];
        celldeath=1;withextmort=1;
        nofonco=4;

        strategies_opt=A(A(1:end,5)==targetsize,1:end);

        % which, out of the four optimised strategies, will be used
        initial_pop_str=strategies_opt(st,6:end);

        if evolveTraits
            alt=nofonco-1;
            strategy_bounds=[0.0, 0.0, 1.0, 1.0;... % prob asym
                  0.0, 0.0, 1.0, 1.0;... % prob dif
                  5, 5, 500, 500;... % tel length
                  0,0,alt,alt;... % apop_thr
                  0.0, 0.0, 1.0, 1.0;... % apop_percent
                  3,3,50,50;... % nof layer
                  0.01,0.01,200,200]; % div_prop
            index=1;indexh=1;
            for i=1:size(initial_pop_str,2)
                strbounds2(index:index+5,:)=strategy_bounds(i,2);
                strbounds3(index:index+5,:)=strategy_bounds(i,3);
                strbounds2h(indexh:indexh+2,:)=strategy_bounds(i,2);
                strbounds3h(indexh:indexh+2,:)=strategy_bounds(i,3);
                index=index+6;
                indexh=indexh+3;
                evol_bounds_diploid=[strbounds2,strbounds3];                
                evol_bounds_haploid=[strbounds2h,strbounds3h];
            end
        else
            evol_bounds_haploid=[];
            evol_bounds_diploid=[];
        end

        % initializes the population with the inputs above, creates variation
        % around the mean trait values
        population=createpop(initial_pop_str,populationsize,targetsize,init_var,evol_bounds_diploid,w_bs,w_cancer);

        % population data structure and indices are as follows:
        % GENOTYPE STRATEGY (1:42) male-female-shared // male-female-shared
        % GENOTYPE TARGET BODYSIZE (43:48)
        % PHENOTYPE STRATEGY-BODYSIZE (49:56)
        % SEX (57), CURRENTSIZE (58), AGE (59), BUDGET (60), 
        % MATURATION TIME (61), OFFSPRING COUNT (62)
        str_index=49:55;targetsize_index=56;sex_index=57;
        currentsize_index=58;age_index=59;budget_index=60;
        mattime_index=61;offcount_index=62;

        Ns{populationsize}=[]; % each individual's body with cells placed in state space

        for i=1:populationsize
            focal_str=population(i,str_index);
            % these three things always need to be integers even when they
            % evolve, because they determine the dimensions of a matrix, which
            % is our 'individual'
            H=focal_str(3); % i.e. how many times a cell can divide (Hayflick limit)             
            K=nofonco; % how many steps until cancer
            T=focal_str(6); % differentiation levels
            N_temp=zeros([H K T]);N_temp(1,1,1)=1; % starting with one cell
            Ns{i}=N_temp; 
        end

        t=1;
    end
    
    while t <= endtime
        % we will keep track of individuals who die
        numdeads_t=0;deathcauses=zeros(1,8);deads=[];
        for i=1:populationsize
            focal_str=population(i,str_index);
            t_bodysize=population(i,targetsize_index);
            age=population(i,age_index);
            mat_time=population(i,mattime_index);
            N=Ns{i};
            if focal_str(1) > 0   % if not an empty spot
            [deathcause,mat_time,addedbudget,new_bs,N]=onelife_per_t(focal_str,mort,...
                cancer_danger,t_bodysize,withextmort,celldeath,N,mat_time,age,nofonco);
                if deathcause~=0
                    % here things related to death
                    numdeads_t=numdeads_t+1;
                    deads(numdeads_t,:)=horzcat(population(i,:),deathcause);
                    deathcauses(deathcause)=deathcauses(deathcause)+1;
                    population(i,1:end)=zeros(1,size(population,2)); % empty the spot
                else
                    % update the indiviual in the population
                    population(i,age_index)=age+1; 
                    population(i,mattime_index)=mat_time; 
                    population(i,currentsize_index)=new_bs;
                    population(i,budget_index)=population(i,budget_index)+addedbudget;
                    Ns{i}=N;
                end
            end
        end
        % make new offspring if there are empty spaces
        [population,Ns,offspringproduced,parents]=reproduce(population,curr_f,curr_m,Ns,nofonco,init_var,...
            evol_bounds_haploid,sex_index,budget_index,currentsize_index,mattime_index,...
            offcount_index,age_index,w_bs,w_cancer);
        
        output(st,t,population,deads,deathcauses,sex_index,offspringproduced,parents,mattime_index);
        t=t+1;    
        if sum(population(1:end,1))==0
            break; % nobody is alive
        end
    end
    filenamews=strcat('output/w',num2str(st),'_',num2str(endtime)); 
    % output the workspace in case we want to continue
    save(filenamews);
end