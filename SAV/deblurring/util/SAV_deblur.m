function [HR,error,Energy_iter] = SAV_deblur(f,kernel)
% parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%model parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.ep=1e-4;
par.eps=1e-4;
par.b=0.1;
par.lambda=0.9;
par.const=1;
%%%schme parameters%%%%%%%%%%%%%%%%%%%
par.tau_min=1e-4;
par.tau_max=0.1;
par.rho=0.7;
par.tol=0.08;
Time=20;
relative_tol=1e-4;

% Initialize%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi=zeros(size(f));
afa=1+par.b*(abs(DisCuNew(phi)));
du_x=dxf(phi);
du_y=dyf(phi);
grad_u=sqrt(du_x.^2+du_y.^2+par.eps);

data_term=(f-imfilter(phi,kernel,'circular','conv')).^2;
inter_data=sum(sum(data_term));
rn=sqrt(sum(sum(afa.*grad_u))+par.lambda/2*inter_data+par.const);
%%%g=imfilter(f,w,filter_mode,boundary_options,size_options)
Energy=sum(sum(0.5*par.ep*grad_u.^2))+rn^2;
t=2*par.tau_min;
Energy_iter=zeros(round(Time/par.tau_min)+1,1);
error=zeros(round(Time/par.tau_min)+1,1);
iter=1;
Energy_iter(iter,1)=Energy;
[l1,l2]=size(phi);
DtD = abs(psf2otf([1,-1],[l1, l2])).^2 + abs(psf2otf([1;-1],[l1, l2])).^2;
tau_n=par.tau_max;

while t<Time
    phi_old=phi;    
 
    [phi,tau_n,Energy_new,afa,du_x,du_y,grad_u,inter_data]=AdaptiveTimestep(par,afa,du_x,du_y,grad_u,phi,inter_data,f,DtD,tau_n,phi_old,rn,Energy_iter(iter,1),kernel);    
    iter=iter+1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    Energy_iter(iter,1) =Energy_new;
    error(iter)=norm( Energy_iter(iter,1)- Energy_iter(iter-1,1))/norm( Energy_iter(iter,1));
    details = [num2str(iter), 'th iteration: ', 'error: ', num2str(error(iter)), 'time', num2str(tau_n),'total: ', num2str(t), 'energy: ', num2str(Energy_new)];
    disp(details)
    if error(iter) < relative_tol
        break;
    end
    if iter>300
        break
    end
    t=t+tau_n;
 end
HR=phi;
end


