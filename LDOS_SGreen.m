% SSH Model
% H=[0 u*cos(kx)+1j*v*sin(kx);u*cos(kx)-1j*v*sin(kx) 0];
%
clear;
sigmax=[0,1;1,0];sigmay=1j*[0,-1;1,0];
L_x=101;
n_DOF=2;


% parameter of system
v=1;
w_range=v*(0:0.02:2);
Ene_list=zeros(numel(w_range),n_DOF*L_x);

% parameter of Green function
eta=1e-2;
nw=numel(w_range);nE=311;n_mat=n_DOF*L_x;
E_min=-1.5;E_max=1.5;
E_list=linspace(E_min,E_max,nE);

for i_w=1:numel(w_range)
    w=w_range(i_w);
    Ham_OPC=zeros(L_x*n_DOF,L_x*n_DOF);
    i_x=0;
    Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,i_x*n_DOF+1:(i_x+1)*n_DOF)=v*sigmax;
    Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,(i_x+1)*n_DOF+1:(i_x+2)*n_DOF)=[0,0;w,0];
    for i_x=1:L_x-2
        Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,(i_x-1)*n_DOF+1:i_x*n_DOF)=Ham_OPC((i_x-1)*n_DOF+1:i_x*n_DOF,i_x*n_DOF+1:(i_x+1)*n_DOF)';
        Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,i_x*n_DOF+1:(i_x+1)*n_DOF)=v*sigmax;
        Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,(i_x+1)*n_DOF+1:(i_x+2)*n_DOF)=[0,0;w,0];
    end
    i_x=L_x-1;
    Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,(i_x-1)*n_DOF+1:i_x*n_DOF)=Ham_OPC((i_x-1)*n_DOF+1:i_x*n_DOF,i_x*n_DOF+1:(i_x+1)*n_DOF)';
    Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,i_x*n_DOF+1:(i_x+1)*n_DOF)=v*sigmax;
    Ene_list(i_w,:)=eig(Ham_OPC);
end

figure(1)
subplot(1,3,1)
plot(w_range,Ene_list,'.b');
xlabel('w/v');ylabel('E');
ylim([E_min,E_max])

% Green Function


GR=zeros(nE,nw);
GA=zeros(nE,nw);
A=zeros(nE,nw);

for i_E=1:nE
    for i_w=1:nw
        w=w_range(i_w);
        
        Ham_OPC=zeros(L_x*n_DOF,L_x*n_DOF);
        i_x=0;
        Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,i_x*n_DOF+1:(i_x+1)*n_DOF)=v*sigmax;
        Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,(i_x+1)*n_DOF+1:(i_x+2)*n_DOF)=[0,0;w,0];
        for i_x=1:L_x-2
            Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,(i_x-1)*n_DOF+1:i_x*n_DOF)=Ham_OPC((i_x-1)*n_DOF+1:i_x*n_DOF,i_x*n_DOF+1:(i_x+1)*n_DOF)';
            Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,i_x*n_DOF+1:(i_x+1)*n_DOF)=v*sigmax;
            Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,(i_x+1)*n_DOF+1:(i_x+2)*n_DOF)=[0,0;w,0];
        end
        i_x=L_x-1;
        Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,(i_x-1)*n_DOF+1:i_x*n_DOF)=Ham_OPC((i_x-1)*n_DOF+1:i_x*n_DOF,i_x*n_DOF+1:(i_x+1)*n_DOF)';
        Ham_OPC(i_x*n_DOF+1:(i_x+1)*n_DOF,i_x*n_DOF+1:(i_x+1)*n_DOF)=v*sigmax;

        GR(i_E,i_w)=trace(1/2/pi*inv((E_list(i_E)+1j*eta)*eye(n_mat)-Ham_OPC));
        GA(i_E,i_w)=trace(1/2/pi*inv((E_list(i_E)-1j*eta)*eye(n_mat)-Ham_OPC));
    end
end

subplot(1,3,2)
A=1j*(GR-GA);
surf(w_range,E_list,A/2/pi);shading('interp');view(0,90);
xlabel('w/v');ylabel('E');
colorbar;
caxis([0 30]);

% Surface Green Function

eta=1e-2;
L_surface=5;L_iter=5;%L_surface
nw=numel(w_range);nE=311;n_mat=n_DOF*L_surface;

E_R_list=linspace(E_min,E_max,nE)+1j*eta;
E_A_list=linspace(E_min,E_max,nE)-1j*eta;
GR_S=zeros(nE,nw);
GA_S=zeros(nE,nw);
A_S=zeros(nE,nw);
for i_E=1:nE
    for i_w=1:nw
        w=w_range(i_w);
        Ham_S=zeros(n_mat,n_mat);
   % L_surface>=3
        i_surf=0;
        Ham_S(i_surf*n_DOF+1:(i_surf+1)*n_DOF,i_surf*n_DOF+1:(i_surf+1)*n_DOF)...
                =v*sigmax;
        Ham_S(i_surf*n_DOF+1:(i_surf+1)*n_DOF,(i_surf+1)*n_DOF+1:(i_surf+2)*n_DOF)...
                =[0,0;w,0];
            
        Ham_R0(i_surf*n_DOF+1:(i_surf+1)*n_DOF,i_surf*n_DOF+1:(i_surf+1)*n_DOF)...
                =v*sigmax;
        Ham_R0(i_surf*n_DOF+1:(i_surf+1)*n_DOF,(i_surf+1)*n_DOF+1:(i_surf+2)*n_DOF)...
                =[0,0;w,0];
        for i_surf=1:L_surface-2
            Ham_S(i_surf*n_DOF+1:(i_surf+1)*n_DOF,(i_surf-1)*n_DOF+1:i_surf*n_DOF)...
                =Ham_S((i_surf-1)*n_DOF+1:i_surf*n_DOF,i_surf*n_DOF+1:(i_surf+1)*n_DOF)';
            Ham_S(i_surf*n_DOF+1:(i_surf+1)*n_DOF,i_surf*n_DOF+1:(i_surf+1)*n_DOF)...
                =v*sigmax;
            Ham_S(i_surf*n_DOF+1:(i_surf+1)*n_DOF,(i_surf+1)*n_DOF+1:(i_surf+2)*n_DOF)...
                =[0,0;w,0];
            
            Ham_R0(i_surf*n_DOF+1:(i_surf+1)*n_DOF,(i_surf-1)*n_DOF+1:i_surf*n_DOF)...
                =Ham_R0((i_surf-1)*n_DOF+1:i_surf*n_DOF,i_surf*n_DOF+1:(i_surf+1)*n_DOF)';
            Ham_R0(i_surf*n_DOF+1:(i_surf+1)*n_DOF,i_surf*n_DOF+1:(i_surf+1)*n_DOF)...
                =v*sigmax;
            Ham_R0(i_surf*n_DOF+1:(i_surf+1)*n_DOF,(i_surf+1)*n_DOF+1:(i_surf+2)*n_DOF)...
                =[0,0;w,0];
        end
        i_surf=L_surface-1;
        Ham_S(i_surf*n_DOF+1:(i_surf+1)*n_DOF,(i_surf-1)*n_DOF+1:i_surf*n_DOF)...
            =Ham_S((i_surf-1)*n_DOF+1:i_surf*n_DOF,i_surf*n_DOF+1:(i_surf+1)*n_DOF)';
        Ham_S(i_surf*n_DOF+1:(i_surf+1)*n_DOF,i_surf*n_DOF+1:(i_surf+1)*n_DOF)...
            =v*sigmax;
        
        Ham_R0(i_surf*n_DOF+1:(i_surf+1)*n_DOF,(i_surf-1)*n_DOF+1:i_surf*n_DOF)...
            =Ham_R0((i_surf-1)*n_DOF+1:i_surf*n_DOF,i_surf*n_DOF+1:(i_surf+1)*n_DOF)';
        Ham_R0(i_surf*n_DOF+1:(i_surf+1)*n_DOF,i_surf*n_DOF+1:(i_surf+1)*n_DOF)...
            =v*sigmax;
        
        Ham_SR0=zeros(L_surface*n_DOF);
        Ham_R0R1=zeros(L_surface*n_DOF);
        Ham_SR0(L_surface*n_DOF,1)=w;
        Ham_R0R1(L_surface*n_DOF,1)=w;
        
% L_surface=1
%         Ham_S=v*sigmax;
%         Ham_R0=v*sigmax;
%         Ham_SR0=[0,0;w,0];
%         Ham_R0R1=[0,0;w,0];

% L_surface=2
%         Ham_S=[0,v,0,0;v,0,w,0;0,w,0,v;0,0,v,0];
%         Ham_R0=[0,v,0,0;v,0,w,0;0,w,0,v;0,0,v,0];
%         Ham_SR0=[0,0,0,0;0,0,0,0;0,0,0,0;w,0,0,0];
%         Ham_R0R1=[0,0,0,0;0,0,0,0;0,0,0,0;w,0,0,0];

        
        t0=inv(E_R_list(i_E)*eye(n_mat)-Ham_R0)*Ham_R0R1';
        t0_tilted=inv(E_R_list(i_E)*eye(n_mat)-Ham_R0)*Ham_R0R1;
        ti=t0;ti_tilted=t0_tilted;
        T=ti;T_i=ti_tilted;
        for i_iter=1:L_iter
            t_temp=pinv(eye(n_mat)-ti*ti_tilted-ti_tilted*ti)*ti*ti;
            t_temp_iter=pinv(eye(n_mat)-ti*ti_tilted-ti_tilted*ti)*ti_tilted*ti_tilted;
            ti=t_temp;
            ti_iter=t_temp_iter;
            T=T+T_i*ti;
            T_i=T_i*ti_tilted;
        end
        g_R0=inv(E_R_list(i_E)*eye(n_mat)-Ham_R0-Ham_R0R1*T);
        sum_R=Ham_SR0*g_R0*Ham_SR0';
        GR_S(i_E,i_w)=trace(inv(E_R_list(i_E)*eye(n_mat)-Ham_S-sum_R));
        
        t0=inv(E_A_list(i_E)*eye(n_mat)-Ham_R0)*Ham_R0R1';
        t0_tilted=inv(E_A_list(i_E)*eye(n_mat)-Ham_R0)*Ham_R0R1;
        ti=t0;ti_tilted=t0_tilted;
        T=ti;T_i=ti_tilted;
        for i_iter=1:L_iter
            t_temp=pinv(eye(n_mat)-ti*ti_tilted-ti_tilted*ti)*ti*ti;
            t_temp_iter=pinv(eye(n_mat)-ti*ti_tilted-ti_tilted*ti)*ti_tilted*ti_tilted;
            ti=t_temp;
            ti_iter=t_temp_iter;
            T=T+T_i*ti;
            T_i=T_i*ti_tilted;
        end
        g_R0=inv(E_A_list(i_E)*eye(n_mat)-Ham_R0-Ham_R0R1*T);
        sum_R=Ham_SR0*g_R0*Ham_SR0';
        GA_S(i_E,i_w)=trace(inv(E_A_list(i_E)*eye(n_mat)-Ham_S-sum_R));

    end
end

subplot(1,3,3)
%surf(w_range,E_list,-imag(GR_S)/pi);shading('interp');view(0,90);
A_S=1j*(GR_S-GA_S);
surf(w_range,E_list,A_S/2/pi);shading('interp');view(0,90);
xlabel('w/v');ylabel('E');
colorbar;
caxis([0 30]);