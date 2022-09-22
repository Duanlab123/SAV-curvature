function [ phi,tau_n,Energy_new,afa,du_x,du_y,grad_u,inter_data ] = AdaptiveTimestep( par,afa,du_x,du_y,grad_u,phi,inter_data,LR1,DtD,tau_n,phi_old,rn,Energy,upscaling,h,kernel)
   [phi]=new_order1(par,afa,du_x,du_y,grad_u,phi,inter_data,LR1,DtD,tau_n,rn,upscaling,h,kernel);
   %step1: compute u^{n+1}
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   data_term=(LR1-downsa(imfilter(phi,h,'circular'),upscaling,kernel)).^2;
   inter_data=upscaling^2*sum(sum(data_term));
   afa=1+par.b*(abs(DisCuNew(phi)));
   du_x=dxf(phi);
   du_y=dyf(phi);
   grad_u=sqrt(du_x.^2+du_y.^2+sqrt(par.ep));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   E1=sum(sum(afa.*grad_u))+par.lambda/2*inter_data;   
   r_new= sqrt(E1 +par.const);
   Energy_new=sum(sum(0.5*par.ep*grad_u.^2))+r_new^2;       
    
    err =norm(Energy-Energy_new)/norm( Energy_new);
    if err<par.tol
       tau_n = max(par.tau_min,min(par.tau_max,par.rho*sqrt(par.tol/err)*tau_n));   
    else
        Energy=Energy_new;
        tau_n= max(par.tau_min,min(par.tau_max,par.rho*sqrt(par.tol/err)*tau_n));
    %    [phi,tau_n,Energy_new]=AdaptiveTimestep(par,afa,du_x,du_y,grad_u,phi,LR1,DtD,tau_n,phi_old,rn,Energy,upscaling,h,kernel);
       [ phi,tau_n,Energy_new,afa,du_x,du_y,grad_u,inter_data ] = AdaptiveTimestep( par,afa,du_x,du_y,grad_u,phi,inter_data,LR1,DtD,tau_n,phi_old,rn,Energy,upscaling,h,kernel);
    end
end

