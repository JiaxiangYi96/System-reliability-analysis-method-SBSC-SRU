function [obj,index]=Infill_SYSi_RU(x,search_x,sample_x,run,kriging_model,type)

switch type
    
    case 'parallel'
        % calculate the values of x
        num_obj = length(kriging_model);
        u_g = zeros(size(x,1), num_obj);
        mse_g = zeros(size(x,1), num_obj);
        for ii = 1: num_obj
            [u_g(:, ii), mse_g(:, ii)] = predictor(x, kriging_model{ii});
            s_g=sqrt(max(0,mse_g));
            %% Compute the distance threshold
            d0=zeros(size(sample_x{run,ii},1),size(sample_x{run,ii},1));
            for jj=1:size(sample_x{run,ii},1)
                x_temp=sample_x{run,ii}(jj,:);
                x_temp=repmat(x_temp,size(sample_x{run,ii},1),1);
                d0(jj,:)=sqrt(sum((x_temp-sample_x{run,ii}).^2,2))';
            end
            d0(d0==0)=Inf;
            d_thete=min(min(d0));
            d1=zeros(size(x,1),1);
            for kk = 1:size(x,1)
                d2=sqrt(sum((repmat(x(kk,:),size(sample_x{run,ii},1),1)-sample_x{run,ii}).^2,2));
                d2= d_thete/2-min(d2);
                d1(kk,:)=max(d2,0);
            end
            %
            U_comp(:,ii)=abs(u_g(:,ii))./(Gaussian_PDF(abs(u_g(:,ii))./s_g(:,ii)).*s_g(:,ii));
            penalty(:,ii)=10^10*d1;
            %             Max_error_Comp(:,ii)=bootstrap_stop_single_fidelity(search_x,kriging_model{ii},0.05);
        end
        %          if  ~isempty(find(Max_error_Comp==Inf))
        %             Max_error_Comp(find(Max_error_Comp==Inf))=100;
        %             weight=Max_error_Comp./sum(Max_error_Comp);
        %         else
        %             weight=Max_error_Comp./sum(Max_error_Comp);
        %         end
        %
        weight=1;
        %% identify the two situations
        index_fail=find(sum(u_g<=0,2)==num_obj);
        index_safe=find(sum(u_g<=0,2)< num_obj);
        
        if size(index_fail,1)~=0
            [obj(index_fail,:),index(index_fail,:)]=min(U_comp(index_fail,:)./weight,[],2);
        end
        if size(index_safe,1)~=0
            flag=u_g(index_safe,:)>0;
            [obj(index_safe,:),index(index_safe,:)]=max(flag.*U_comp(index_safe,:).* weight,[],2);
        end
        for nnn=1:size(x,1)
            obj(nnn,:)=obj(nnn,:)+penalty(nnn,index(nnn,:));
        end
        
    case'series'
        % calculate the values of x
        num_obj = length(kriging_model);
        u_g = zeros(size(x,1), num_obj);
        mse_g = zeros(size(x,1), num_obj);
        for ii = 1: num_obj
            [u_g(:, ii), mse_g(:, ii)] = predictor(x, kriging_model{ii});
            s_g=sqrt(max(0,mse_g));
            %% Compute the distance threshold
            d0=zeros(size(sample_x{run,ii},1),size(sample_x{run,ii},1));
            for jj=1:size(sample_x{run,ii},1)
                x_temp=sample_x{run,ii}(jj,:);
                x_temp=repmat(x_temp,size(sample_x{run,ii},1),1);
                d0(jj,:)=sqrt(sum((x_temp-sample_x{run,ii}).^2,2))';
            end
            d0(d0==0)=Inf;
            d_thete=min(min(d0));
            d1=zeros(size(x,1),1);
            for kk = 1:size(x,1)
                d2=sqrt(sum((repmat(x(kk,:),size(sample_x{run,ii},1),1)-sample_x{run,ii}).^2,2));
                d2= d_thete/2-min(d2);
                d1(kk,:)=max(d2,0);
            end
            %
            U_comp(:,ii)=abs(u_g(:,ii))./(Gaussian_PDF(abs(u_g(:,ii))./s_g(:,ii)).*s_g(:,ii));
            penalty(:,ii)=10^10*d1;
%                         Max_error_Comp(:,ii)=bootstrap_stop_single_fidelity(search_x,kriging_model{ii},0.05);
        end
%                 if  ~isempty(find(Max_error_Comp==Inf))
%                     Max_error_Comp(find(Max_error_Comp==Inf))=10;
%                     weight=Max_error_Comp./sum(Max_error_Comp);
%                 else
%                     weight=Max_error_Comp./sum(Max_error_Comp);
%                 end
        weight=1;
        %%
        index_safe=find(sum(u_g>0,2)==num_obj);
        index_fail=find(sum(u_g>0,2)< num_obj);
        if sum(index_safe)>0
            [obj(index_safe,:),index(index_safe,:)]=min(U_comp(index_safe,:)./weight,[],2);
        end
        if  sum(index_fail)>0
            flag=u_g(index_fail,:)>0;
            [obj(index_fail,:),index(index_fail,:)]=max(flag.*U_comp(index_fail,:).*weight ,[],2);
        end
        
        for nnn=1:size(x,1)
            obj(nnn,:)=obj(nnn,:)+penalty(nnn,index(nnn,:));
        end
        
end

end
