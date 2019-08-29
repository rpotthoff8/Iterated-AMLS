
function plotter_2D(data,f_hat,err_hist,c,R,method,fignum)
%Function plots the progress of noisy function approximation

i = size(err_hist.train,1)-1; %Track of number of RBFs placed (iteration number)
h = figure(fignum); %Set figure
% set(h,'Position',[0 138 1120 840]);
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

if isfield(data,'noise_std')
    noize = 2*data.noise_std*ones(size(data.x,1),1); %Noise level for display on residual plot
end

subplot(221) %Plot function and approximation
if size(data.x,2)==2
    plot3(data.x(:,1),data.x(:,2),data.f0,'.');
elseif size(data.x,2)==1
    plot(data.x,data.f0,'.')
end
full_legends={'Train Data','Test Data','Approximation'};
include=boolean([1,0,1]);
hold on
if isfield(data,'x_test')
    include(2)=1;
    if size(data.x,2)==2
        plot3(data.x_test(:,1),data.x_test(:,2),data.f0_test,'.');
    elseif size(data.x,2)==1
        plot(data.x_test,data.f0_test,'.')
    end
end

if size(data.x,2)==2
    plot3(data.x_full(:,1),data.x_full(:,2),f_hat,'.');
elseif size(data.x,2)==1
    plot(data.x_full,f_hat,'.')
end

legend(char(full_legends(include)),'Location','best')
title(['Iteration ', num2str(i)]);
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
grid on
hold off

subplot(222) %Error history plotting
semilogy(0:i,err_hist.train);
full_legends={'Training','Testing','True','Data Error'};
include=boolean([1,0,0,0]);
hold on
if isfield(err_hist,'test')
    semilogy(0:i,err_hist.test);
    include(2)=1;
end
if isfield(err_hist,'true')
    semilogy(0:i,err_hist.true);
    include(3)=1;
end
if isfield(data,'err_actual')
    semilogy(0:i,repelem(data.err_actual,size(err_hist.true,1)),'k');
    include(4)=1;
end
legend(char(full_legends(include)),'Location','best')
title(['Current Training RMS Error: ', num2str(err_hist.train(end))]);
xlabel('Iteration'); ylabel('RMS Error');
grid on
hold off

subplot(223) %Individual residual values plotting'
full_legends={'Training','Testing'};
include=boolean([1,0]);
if size(data.x,2)==2
    plot3(data.x(:,1),data.x(:,2),R.train,'.');
    xlabel('x'); ylabel('y'); zlabel('$$f(x) - \hat{f}(x)$$','Interpreter','Latex');
elseif size(data.x,2)==1
    plot(data.x,R.train,'.');
    xlabel('x'); ylabel('$$f(x) - \hat{f}(x)$$','Interpreter','Latex');
end

hold on
if isfield(R,'test')
    if size(data.x,2)==2
        plot3(data.x_test(:,1),data.x_test(:,2),R.test,'.',data.x(:,1),data.x(:,2),zeros(length(data.x),1),'k.');
    elseif size(data.x,2)==1
        plot(data.x_test,R.test,'.');
    end
    include(2)=1;
end

if size(data.x,2)==2
   plot3(data.x(:,1),data.x(:,2),zeros(length(data.x),1),'k.'); 
end
% if isfield(data,'noise_std')
%     plot(data.x,noize,'-');
%     if isfield(R,'test')
%         plot(data.x_test,abs(R.test),'.')
%         legend('Training','Noise 2*sigma','Testing')
%     else
%         legend('Training','Noise 2*sigma')
%     end
% elseif isfield(R,'test')
%     plot(data.x_test,abs(R.test),'.')
%     legend('Training','Testing')
% else
%     legend('Training')
% end
title('Residuals')
legend(char(full_legends(include)),'Location','best')
grid on
hold off

subplot(224) %Coefficient value plotting
semilogy([1:length(c)],abs(c),'-*','MarkerSize',10); %axis([1 100 10e-4 10e0]);
title('Absolute Coefficients')
xlabel('Basis Function #'); ylabel('$$|c|$$','Interpreter','Latex');
grid on

sgtitle(strcat('Method', " ", num2str(fignum),':'," ",method),'Interpreter','none');
% pause(eps);
end