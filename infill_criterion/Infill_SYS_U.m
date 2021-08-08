function [obj,index]=Infill_SYS_U(x,kriging_model,type)

switch type
    
    case 'parallel'
        % calculate the values of x
        num_obj = length(kriging_model);
        u_g = zeros(size(x,1), num_obj);
        mse_g = zeros(size(x,1), num_obj);
        for ii = 1: num_obj
            [u_g(:, ii), mse_g(:, ii)] = predictor(x, kriging_model{ii});
        end
        s_g=sqrt(max(0,mse_g));
        [y,index]=max(u_g,[],2);
        for jj=1:size(u_g,1)
            s(jj,:)=s_g(jj,index(jj,:));
        end
        obj=abs(y)./s;
        
    case'series'
        % calculate the values of x
        num_obj = length(kriging_model);
        u_g = zeros(size(x,1), num_obj);
        mse_g = zeros(size(x,1), num_obj);
        for ii = 1: num_obj
            [u_g(:, ii), mse_g(:, ii)] = predictor(x, kriging_model{ii});
        end
        s_g=sqrt(max(0,mse_g));
        [y,index]=min(u_g,[],2);
        for jj=1:size(u_g,1)
            s(jj,:)=s_g(jj,index(jj,:));
        end
        obj=abs(y)./s;
    otherwise
        fprintf('This type of system has not been defined!!!\n')
end



end
