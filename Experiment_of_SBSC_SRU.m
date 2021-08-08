clc
clear all
close all
addpath('failure_probability_evaluation','infill_criterion','test_function_system_reliability','Stopping_Criteria');
%% import the basci information of the test function
% four_branches_system multi_modal_function  parallel_three_modes AK_SYSi_parallel_three_modes
for fun_name={'multi_modal_function','four_branches_system','parallel_three_modes','AK_SYSi_parallel_three_modes'}
    test_function=char(fun_name);
    [num_vari,num_obj,num_initial_sample,mu,sigma,design_space,type,stopping_thresholds]=test_function_system_reliability(test_function);
    % define the infill criterion
    for error=stopping_thresholds
        % define the times of repeat
        num_trials=30;
        sample_x=cell(num_trials,num_obj);
        sample_y=cell(num_trials,num_obj);
        record.result=zeros(num_trials,4+num_obj);
        for run=1:num_trials
            % define the number of the kriging models
            kriging_obj = cell(1,num_obj);
            fprintf('-------------Test_function=%s--------------------- error threhold=%d th ----------------------\n',test_function,error)
            %% main iteration process for the components
            % design sapce of the constrain optimization problem  define as
            sample_x_initial= repmat(design_space(1,:),num_initial_sample,1)+repmat(design_space(2,:)-design_space(1,:),num_initial_sample,1).*...
                lhsdesign(num_initial_sample,num_vari,'criterion','maximin','iterations',1000);
            sample_y_initial=feval(test_function,sample_x_initial);
            %  constructing the initial Kriging model
            for ii=1:num_obj
                sample_x{run,ii}=sample_x_initial;
                sample_y{run,ii}=sample_y_initial(:,ii);
                kriging_obj{ii} = dacefit(  sample_x{run,ii},   sample_y{run,ii},'regpoly0','corrgauss',10*ones(1,num_vari),0.1*ones(1,num_vari),30*ones(1,num_vari));
            end
            clear sample_x_initial sample_y_initial
            gen=1;
            num_search=10^4*gen;
            search_x=MCS_Population_Generation(mu,sigma,num_search);
            [~,~,cov_estimate,~]=system_reliabiliy_evaluation_single_fidelity(search_x,num_obj,kriging_obj,test_function,type);
            iteration=0;
            % big circle
            while cov_estimate>0.05 || iteration<=2
                search_x=MCS_Population_Generation(mu,sigma,num_search);
                stop_value=bootstrap_stop_system_single_fidelity(search_x,kriging_obj,0.05,error,type);
                while (stop_value>error || iteration<=2) && iteration<=300
                    % determine the candidate sample and its component index
                    search_x=MCS_Population_Generation(mu,sigma,num_search);
                    %   [U_value,U_index]=Infill_SYSi_RU_without_distance_threshold(search_x,kriging_obj,type);
                    [U_value,U_index]=Infill_SYSi_RU(search_x,search_x,sample_x,run,kriging_obj,type);
                    % [U_value,U_index]=Infill_Weighted_VB(search_x,search_x,kriging_obj,type);
                    [bestobj,Index]=min(U_value);
                    xselected=search_x(Index,:);
                    Index_comp=U_index(Index,:);
                    % update corresponding kriging model
                    if ~ismember( sample_x{run,Index_comp} ,xselected)
                        fselected_temp=feval(test_function,xselected);
                        fselected=fselected_temp(:,Index_comp);
                        sample_x{run,Index_comp} =[sample_x{run,Index_comp} ;xselected];
                        sample_y{run,Index_comp}=[ sample_y{run,Index_comp};fselected];
                        kriging_obj{Index_comp} = dacefit( sample_x{run,Index_comp} ,sample_y{run,Index_comp},'regpoly0','corrgauss',10*ones(1,num_vari),0.1*ones(1,num_vari),30*ones(1,num_vari));
                    else
                        Mu=repmat(mu,num_search,1);
                        Sigma=repmat(sigma,num_search,1);
                        search_x=normrnd(Mu,Sigma);
                    end
                    iteration=iteration+1;
                    %           if num_vari==2
                    %                                 figure_plot_system_single_fidelity(design_space,sigma,test_function,kriging_obj,search_x,sample_x,run,iteration);
                    %           end
                    stop_value=bootstrap_stop_system_single_fidelity(search_x,kriging_obj,0.05,error,type);
                    % print the current information to the screen
                    [pf_estimate,pf_real,~,real_relative_error]=system_reliabiliy_evaluation_single_fidelity(search_x,num_obj,kriging_obj,test_function,type);
                    fprintf(' run=%d;component=%d th;iteration=%d;pf_estimated=%f; pf_real=%f; stop_value=%f real_error=%f\n', run,Index_comp,iteration,pf_estimate,pf_real,stop_value,real_relative_error);
                end
                num_search=2*num_search;
                
                %% calc      lulate the failure probility
                [~,~,cov_estimate,~]=system_reliabiliy_evaluation_single_fidelity(search_x,num_obj,kriging_obj,test_function,type);
                gen=gen+1;
            end
            
            search_x_system=MCS_Population_Generation(mu,sigma,5*10^6);
            [pf_estimate,pf_real,cov_estimate,real_relative_error]=system_reliabiliy_evaluation_single_fidelity(search_x_system,num_obj,kriging_obj,test_function,type);
            fprintf(' run=%d;gen=%d;pf_estimated=%f; pf_real=%f; stop_value=%f real_error=%f\n', run,gen,pf_estimate,pf_real,stop_value,real_relative_error);
            for jj=1:num_obj
                num_sample_cost(1,jj)=size(sample_x{run,jj},1);
            end
            record.result(run,:)=[pf_estimate,pf_real,cov_estimate,real_relative_error,num_sample_cost];
            save(strcat('Results/',test_function,'/',mfilename,'_',num2str(1000*error),'.mat'),'record','sample_x','sample_y');
        end
    end
    
end
