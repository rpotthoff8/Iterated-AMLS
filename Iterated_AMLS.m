% Last Updated: 6/04/19

function [xc,eps,mult,h,err,err_hist]=Iterated_AMLS(iterations,data,options)
%
%INPUT VARIABLES:
%numRBF: number of RBFs to be used
%data: data to use in function, includes following:
%       data.x: x positions, COLUMN vector, MUST BE PROVIDED
%       data.f0: y values for x positions, COLUMN vector, MUST BE PROVIDED
%       data.x_test: test x positions, COLUMN vector, OPTIONAL
%       data.f0a_test: test y values, COLUMN vector, OPTIONAL
%       data.noise_std: guess for standard deviation of noise, used for
%       self termination, Single value, OPTIONAL
%       data.f_true: actual y values (without noise), COLUMN vector, OPTIONAL

%options: options for function, includes the following:
%       data.plot: if plotter should be used: 1 for yes, 0 for no
%       (default is to not plot)
%       data.terminate_percent: if data.noise_std is provided, code will
%       terminate when this percentage of residuals are <= 2*noise_std,
%       OPTIONAL

%Variables for plotter function
method=options.label;

%Check if test data is provided, set 'test_data=1' if it is
if isfield(data,'x_test')
    test_data=1 ;
    [~, I_xOrder]=sortrows([data.x;data.x_test]); %Sort x in order, save indices of order
else
    test_data=0;
end

%Check if true data is provide, set 'true_data=1' if it is
if isfield(data,'f_true')
    true_data=1;
else
    true_data=0;
end

%Check if plotter option is provide, if not, default is to turn plotter off
if isfield(options,'plot')
    plot=options.plot;
    fignum=options.fignum;
else
    plot=0;
end

if isfield(options,'self_terminate')
    self_terminate=options.self_terminate;
else
    self_terminate=0;
end

% %variables determined from x data provided
% Ndata = size(data.x,1);
% xmin = min(data.x);
% xmax = max(data.x);
% della = (xmax-xmin)/(Ndata-1); %Uniform distribution fill-distance
% h=max(della); %Used fill distance based on farthest spacing
h=0;
eps=0.4; %SET THIS/ OPTIMIZE THIS. %From 'Iterated Approximate Moving Least Squares Approximation', Fasshauer and Zhang, Page 235
mult=cell(iterations,1);
e_norm=zeros(iterations,1);
f_range=max(data.f0)-min(data.f0);
e_norm_norm=zeros(iterations,1);
% END initialize section

%Initial AMLS quasi-interpolant: %From 'Iterated Approximate Moving Least Squares Approximation', Fasshauer
%and Zhang, Equation 12
xc=data.x;
mult{1}=data.f0;

[phi.train,h] = basisFunction_AMLS(data.x,xc,eps,h);
Q.train=phi.train*mult{1};

R.train = data.f0 - Q.train;
err.train = sqrt(mean(R.train.^2));
%The error history and the current approximation are plotted together every
%iteration.
err_hist.train = [err.train];

%LOOCV: From: 'Scattered Data Approximation of Noisy Data via Iterated
%Moving Least Squares', Fasshauer and Zhang, page 8
sub=eye(size(phi.train,1))-phi.train;
inverse_interp=sub^0;
rippa_top=inverse_interp*data.f0;
e=rippa_top./diag(inverse_interp);
e_norm(1)=norm(e);
e_norm_norm(1)=e_norm(1)/f_range;


if test_data==1 %If test data provided, determine approximation/residuals for test datapoints
    phi.test = basisFunction_AMLS(data.x_test,xc,eps,h);
    Q.test =phi.test*mult{1};
    
    R.test = data.f0_test - Q.test;
    err.test = sqrt(mean(R.test.^2));
    err_hist.test = [err.test];
    Q.full=[Q.train;Q.test];
    Q.full=Q.full(I_xOrder,:);
else
    Q.full=Q.train;
end

if true_data==1 %If true data is provided, determine actual residual/error
    R.true=(data.f_true)-Q.full;
    err.true=sqrt(mean(R.true.^2));
    err_hist.true = [err.true];
end
c=1;
if plot==1 %Call plotter function if turned on
    plotter_2D(data,Q.full,err_hist,c,R,method,fignum)
end


%NOTE: For the time being, the code is run a predetermined number of
%iterations, the convergence criteria has not been defined.

for i=2:iterations %iterations
    
    %The aproximation is updated.
    mult{i}=R.train;
    
    Q.train = Q.train+phi.train*mult{i};
    R.train = (data.f0) - Q.train;
    err.train = sqrt(mean(R.train.^2));
    err_hist.train = [err_hist.train;err.train];
    
    %LOOCV: From: 'Scattered Data Approximation of Noisy Data via Iterated
    %Moving Least Squares', Fasshauer and Zhang, page 8
    inverse_interp=inverse_interp+sub^(i-1);
    rippa_top=inverse_interp*data.f0;
    e=rippa_top./diag(inverse_interp);
    e_norm(i)=norm(e);
    e_norm_norm(i)=e_norm(i)/f_range;

    if test_data==1 %If test data provided, determine approximation/residuals for test datapoints
        %The aproximation is updated.
        Q.test = Q.test+phi.test*mult{i};
        R.test = (data.f0_test) - Q.test;
        err.test = sqrt(mean(R.test.^2));
        err_hist.test = [err_hist.test;err.test];
        Q.full=[Q.train;Q.test];
        Q.full=Q.full(I_xOrder,:);
        
    else
        Q.full=Q.train;
    end
    
    if true_data==1 %If true data is provided, determine actual residual/error
        R.true=(data.f_true)-Q.full;
        err.true=sqrt(mean(R.true.^2));
        err_hist.true = [err_hist.true;err.true];
    end
    
    if plot==1 %Call plotter function if turned on
        plotter_2D(data,Q.full,err_hist,e_norm,R,method,fignum)
    end
    
    if self_terminate==1
        tol=-0.01;
        if e_norm_norm(i)-e_norm_norm(i-1)>=tol
            mult(i+1:end)=[];
            disp(strcat('Reached LOOCV termination criteria, # Iterations=',string(i)))

            break          
        end
    end
    
    e_norm_norm(i)-e_norm_norm(i-1)
end

%Set coefficeints to zero for RBFs past min training error
if test_data==1 && options.best_test==1
    [~,minItest]=min(err_hist.test);
    if minItest<i
        mult(minItest+1:end)=[];
        
        %The aproximation is updated.
        Q.train=Iterated_AMLS_construct(data.x,xc,eps,h,mult);
        R.train = (data.f0) - Q.train;
        err.train = sqrt(mean(R.train.^2));
        
        if test_data==1 %If test data provided, determine approximation/residuals for test datapoints
            %The aproximation is updated.
            Q.test = Iterated_AMLS_construct(data.x_test,xc,eps,h,mult);
            R.test = (data.f0_test) - Q.test;
            err.test = sqrt(mean(R.test.^2));
            
            Q.full=[Q.train;Q.test];
            Q.full=Q.full(I_xOrder,:);
        else
            Q.full=Q.train;
        end
        
        if true_data==1 %If true data is provided, determine actual residual/error
            R.true=(data.f_true)-Q.full;
            err.true=sqrt(mean(R.true.^2));
        end
    end
end

end

%% subroutine functions past this point.
%
%************************


%*********************