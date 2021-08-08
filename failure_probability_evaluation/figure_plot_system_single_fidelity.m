function obj=figure_plot_system_single_fidelity(design_space,sigma,test_function,kriging_obj,search_x,sample_x,run,iteration)

num_test=100;
x1_plot=linspace(design_space(1,1)-2*sigma(1),design_space(2,1)+2*sigma(1),num_test);
x2_plot=linspace(design_space(1,2)-2*sigma(2),design_space(2,2)+2*sigma(2),num_test);
% calculate the real responses
for ii=1:1:num_test
    for jj=1:1:num_test
        X_temp=[x1_plot(ii),x2_plot(jj)];
        Y_temp=feval(test_function,X_temp);
   
        for kk = 1: 3
        Y_pre_temp(:, kk) = predictor(X_temp, kriging_obj{kk});
        end
            
        y_real_1(jj,ii)=Y_temp(:,1);
        y_real_2(jj,ii)=Y_temp(:,2);
        y_real_3(jj,ii)=Y_temp(:,3);
        y_pre_1(jj,ii)=Y_pre_temp(:,1);
        y_pre_2(jj,ii)=Y_pre_temp(:,2);
        y_pre_3(jj,ii)=Y_pre_temp(:,3);
    end
end
figure (iteration)
v=[-1e6,0,1e6];
plot(search_x(:,1),search_x(:,2),'y.','MarkerSize',2);
hold on
contour(x1_plot,x2_plot,y_real_1,v,'k-','LineWidth',2,'ShowText','on');
hold on
contour(x1_plot,x2_plot,y_real_2,v,'k:','LineWidth',2,'ShowText','on');
hold on
contour(x1_plot,x2_plot,y_real_3,v,'k--','LineWidth',2,'ShowText','on');
hold on
contour(x1_plot,x2_plot,y_pre_1,v,'m--','LineWidth',2,'ShowText','on');
hold on
contour(x1_plot,x2_plot,y_pre_2,v,'b--','LineWidth',2,'ShowText','on');
hold on
contour(x1_plot,x2_plot,y_pre_3,v,'c--','LineWidth',2,'ShowText','on');
hold on
plot(sample_x{run,1}(1:10,1),sample_x{run,1}(1:10,2),'ko','MarkerFaceColor','k','MarkerSize',6,'LineWidth',1.5);
hold on
plot(sample_x{run,1}(11:end,1),sample_x{run,1}(11:end,2),'ms','MarkerSize',10,'LineWidth',1.5);
hold on
plot(sample_x{run,2}(1:10,1),sample_x{run,2}(1:10,2),'ko','MarkerFaceColor','k','MarkerSize',6,'LineWidth',1.5);
hold on
plot(sample_x{run,2}(11:end,1),sample_x{run,2}(11:end,2),'bo','MarkerSize',10,'LineWidth',1.5);
hold on
plot(sample_x{run,3}(1:10,1),sample_x{run,3}(1:10,2),'ko','MarkerFaceColor','k','MarkerSize',6,'LineWidth',1.5);
hold on
plot(sample_x{run,3}(11:end,1),sample_x{run,3}(11:end,2),'c^','MarkerSize',10,'LineWidth',1.5);
hold on

set(gca,'XTick',design_space(1,1):(design_space(2,1)-design_space(1,1))/5:design_space(2,1))
set(gca,'YTick',design_space(1,2):(design_space(2,2)-design_space(1,2))/5:design_space(2,2))
set(gca,'fontname','Times New Roman','LineWidth',1.5,'fontsize',18)
xlabel('\it x_1','fontname','Times New Roman','fontsize',16)
ylabel('\it x_2','fontname','Times New Roman','fontsize',16)
drawnow
obj=[];
end