function [obj,p_safe_fail,p_fail_safe]=eek_stopping_criterion(search_x,kriging_model,type)
switch type
    case 'parallel'
        %% First stage
        % calculate the predicted value
        num_obj = length(kriging_model);
        u_g = zeros(size(search_x,1), num_obj);
        mse_g = zeros(size(search_x,1), num_obj);
        for ii = 1: num_obj
            [u_g(:, ii), mse_g(:, ii)] = predictor(search_x, kriging_model{ii});
        end
        mse_g=sqrt(max(0,mse_g));
        index_fail=find(sum(u_g<=0,2)==num_obj);
        index_safe=find(sum(u_g<=0,2)< num_obj);
       %% Second stage: identify the samples located in the fail and failure areas respectively
        %% find the samples located in the failure regions to calculate their wrong sign expectation in safety regions
        Fail_y=u_g(index_fail,:);
        Fail_mse=mse_g(index_fail,:);
        % calculate the  wrong sign prediction of each component
        for ii=1:num_obj
           p_fail_safe_Comp(:,ii)=Gaussian_CDF(-abs(Fail_y(:,ii)./Fail_mse(:,ii)));
        end
        p_fail_safe=1-prod((1-p_fail_safe_Comp),2);
        p_fail_fail=prod((1-p_fail_safe_Comp),2);
        % P_fail=normcdf(0,Safe_y,Safe_mse);
          mu_safe=size(Fail_y,1)*mean(p_fail_safe);
          sigma_safe=size(Fail_y,1)*mean(p_fail_safe.*p_fail_fail);
          F_con_interval_up=poissinv(0.975,mu_safe);
          F_con_interval_low=poissinv(0.025,mu_safe);    
          CI_safe_max=max(F_con_interval_up,F_con_interval_low);
      
        %% Third stage: identify the samples located in the safe and failure areas respectively
        Safe_y=u_g(index_safe,:);
        Safe_mse=mse_g(index_safe,:);
        % calculate the  wrong sign prediction of each component
        for ii=1:num_obj
            p_safe_fail_Comp(:,ii)=Gaussian_CDF(-abs(Safe_y(:,ii)./Safe_mse(:,ii)));
        end
        for jj=1:size(Safe_y,1)
            flag_fail=find(Safe_y(jj,:)<0);
            flag_safe=find(Safe_y(jj,:)>=0);
            p_safe_fail(jj,:)=prod(1-p_safe_fail_Comp(jj,flag_fail)).*prod(p_safe_fail_Comp(jj,flag_safe));
        end
           p_safe_safe=1-p_safe_fail;
           mu_fail=size(Safe_y,1)*mean(p_safe_fail);
           std_fail=size(Safe_y,1)*mean(p_safe_fail.*p_safe_safe);
         CI_fail_max_up=norminv(0.975,mu_fail,std_fail);%Kriging预测的失效域内错误分类的候选点个数的上界
         CI_fail_max_low=norminv(0.025,mu_fail,std_fail);
         CI_fail_max=max(CI_fail_max_up,CI_fail_max_low);              
        %% Fourth stage: calculate the final value of maximum estimated relative error
        Predict_fail=find(sum(u_g<=0,2)==num_obj);
        N_f=size(Predict_fail,1);
        error_1=abs((N_f/(N_f-CI_safe_max))-1);
        error_2=abs((N_f/(N_f+CI_fail_max))-1);
        obj=max(error_1,error_2);
        
    case 'series'
        %% First stage: identify the highly-uncertain samples among the MCS population
        % calculate the predicted value
        num_obj = length(kriging_model);
        u_g = zeros(size(search_x,1), num_obj);
        mse_g = zeros(size(search_x,1), num_obj);
   
        for ii = 1: num_obj
            [u_g(:, ii), mse_g(:, ii)] = predictor(search_x, kriging_model{ii});
        end
        mse_g=sqrt(max(0,mse_g));
        index_fail=find(sum(u_g>=0,2)  <num_obj);
        index_safe=find(sum(u_g>=0,2)== num_obj);
        
        %% Second stage: identify the samples located in the safe and failure areas respectively
        % find the samples located in the Safe regions to calculate their wrong sign expectation in failure regions
        Safe_y=u_g(index_safe,:);
        Safe_mse=mse_g(index_safe,:);
        % calculate the  wrong sign prediction of each component
        for ii=1:num_obj
            p_safe_fail_Comp(:,ii)=Gaussian_CDF(-abs(Safe_y(:,ii)./Safe_mse(:,ii)));
        end
        p_safe_fail=1-prod((1-p_safe_fail_Comp),2);
        p_safe_safe=prod((1-p_safe_fail_Comp),2);
        mu_fail=size(Safe_y,1)*mean(p_safe_fail);
        std_fail=size(Safe_y,1)*mean(p_safe_fail.*p_safe_safe);
       CI_fail_max_up=norminv(0.975,mu_fail,std_fail);
       CI_fail_max_low=norminv(0.025,mu_fail,std_fail);    
       CI_fail_max=max(CI_fail_max_up,CI_fail_max_low);
        % estimate the bootstrapt confidence interval
    
        %% Third stage: identify the samples located in the fail but actual safe
        Fail_y=u_g(index_fail,:);
        Fail_mse=mse_g(index_fail,:);
        % calculate the  wrong sign prediction of each component
          for ii=1:num_obj
            p_fail_safe_Comp(:,ii)=Gaussian_CDF(-abs(Fail_y(:,ii)./Fail_mse(:,ii)));
          end
         for jj=1:size(Fail_y,1)
            flag_fail=find(Fail_y(jj,:)<0);
            flag_safe=find(Fail_y(jj,:)>=0);
            p_fail_safe(jj,:)=prod(p_fail_safe_Comp(jj,flag_fail)).*prod(1-p_fail_safe_Comp(jj,flag_safe));
         end
        p_fail_fail=1-p_fail_safe;
        mu_safe=size(Fail_y,1)*mean(p_fail_safe);
        std_safe=size(Fail_y,1)*mean(p_fail_safe.*p_fail_fail);
        CI_safe_max_up=poissinv(0.975,mu_safe);
        CI_safe_max_low=poissinv(0.025,mu_safe);  
        CI_safe_max=max(CI_safe_max_up,CI_safe_max_low);
                
        %% Fourth stage: calculate the final value of maximum estimated relative error
        Predict_fail=find(sum(u_g>=0,2)<num_obj);
        N_f=size(Predict_fail,1);
        error_1=abs((N_f/(N_f-CI_safe_max))-1);
        error_2=abs((N_f/(N_f+CI_fail_max))-1);
        obj=max(error_1,error_2);
        
end
end
