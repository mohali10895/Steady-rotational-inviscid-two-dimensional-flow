clc
clear 
close all

%% Inputs
d=10;
div=2;                                      %no. of divisions per unit
delta_x=1/div;
delta_y=1/div;
delta_x_star=delta_x/d;
delta_y_star=delta_y/d;
NP=9;
N_cell=4;
zeta=[-sqrt(0.6),0,sqrt(0.6),sqrt(0.6),0,-sqrt(0.6),-sqrt(0.6),0,sqrt(0.6)];
eta=[-sqrt(0.6),-sqrt(0.6),-sqrt(0.6),0,0,0,sqrt(0.6),sqrt(0.6),sqrt(0.6)];
X_star_vec=0:(1/(div*d)):3;
Y_star_vec=0:(1/(div*d)):1;
c_epsi=-((pi^2)/4);

%% Meshing
% part 1
L_x_P1=2*d;
L_y_P1=d;
NE_x_P1=div*L_x_P1;
NE_y_P1=div*L_y_P1;
NE_P1=NE_x_P1*NE_y_P1;
NN_x_P1=NE_x_P1+1;
NN_y_P1=NE_y_P1+1;
NN_P1=NN_x_P1*NN_y_P1;

X_P1=nan(NN_x_P1,NN_y_P1);
Y_P1=nan(NN_x_P1,NN_y_P1);
for i=1:NN_x_P1
    for j=1:NN_y_P1
        X_P1(i,j)=(i-1)*delta_x;
        Y_P1(i,j)=(j-1)*delta_y;
    end
end

NN_T_P1=zeros(NN_x_P1,NN_y_P1);
for j=1:NN_y_P1
    for i=1:NN_x_P1
        NN_T_P1(i,j)=i+(j-1)*NN_x_P1;
    end
end


% part 2
L_x_P2=d;
L_y_P2=d/10;
NE_x_P2=div*L_x_P2;
NE_y_P2=div*L_y_P2;
NE_P2=NE_x_P2*NE_y_P2;
NN_x_P2=NE_x_P2+1;
NN_y_P2=NE_y_P2+1;
NN_P2=NN_x_P2*NN_y_P2;

X_P2=nan(NN_x_P2,NN_y_P2);
Y_P2=nan(NN_x_P2,NN_y_P2);
for i=1:NN_x_P2
    for j=1:NN_y_P2
        X_P2(i,j)=(i-1)*delta_x+L_x_P1;
        Y_P2(i,j)=(j-1)*delta_y;
    end
end

NN_T_P2=zeros(NE_x_P2,NE_y_P2);
for j=1:NN_y_P2
    for i=1:NE_x_P2
        NN_T_P2(i,j)=i+(j-1)*NE_x_P2+max(max(NN_T_P1));
    end
end


%Total
X=nan(size(X_P1,1)+size(X_P2,1)-1,size(Y_P1,2));
Y=nan(size(X_P1,1)+size(X_P2,1)-1,size(Y_P1,2));
NN_tot=nan(size(X_P1,1)+size(X_P2,1)-1,size(Y_P1,2));
for i=1:size(X_P1,1)
    for j=1:size(Y_P1,2)
        X(i,j)=X_P1(i,j);
        Y(i,j)=Y_P1(i,j);
        NN_tot(i,j)=NN_T_P1(i,j);
    end
end
for i=(size(X_P1,1)+1):(size(X_P1,1)+size(X_P2,1)-1)
    for j=1:size(Y_P2,2)
        X(i,j)=X_P2(i-size(X_P1,1)+1,j);
        Y(i,j)=Y_P2(i-size(X_P1,1)+1,j);
        NN_tot(i,j)=NN_T_P2(i-size(X_P1,1),j);
    end
end

X_star=X./d;
Y_star=Y./d;

N_symm=NN_tot(:,1);
N_H_P1=NN_T_P1(:,size(NN_T_P1,2));
N_H_P2=NN_T_P2(:,size(NN_T_P2,2));
N_inlet=NN_T_P1(1,:).';
N_outlet=NN_T_P2(size(NN_T_P2,1),:).';
N_ver=zeros(size(NN_T_P1,2)-NN_y_P2,1);
for i=NN_y_P2:size(NN_T_P1,2)
    N_ver(i+1-NN_y_P2,1)=NN_T_P1(size(NN_T_P1,1),i);
end

%% Solver For Case a
K_global_a=zeros(max(max(NN_tot)),max(max(NN_tot)));
RHS_a=zeros(max(max(NN_tot)),1);
for n=1:NE_y_P1
    if n<=NE_y_P2
        var_N=NE_x_P1+NE_x_P2;
    elseif n>NE_y_P2
        var_N=NE_x_P1;
    end
    for m=1:var_N
        % Weights Calculation
        W_i=zeros(1,NP);
        for i=1:1:NP
            if zeta(i)~=0
                W_i(i)=5/9;
            elseif zeta(i)==0
                W_i(i)=8/9;
            end
        end
        W_j=zeros(1,NP);
        for j=1:1:NP
            if eta(j)~=0
                W_j(j)=5/9;
            elseif eta(j)==0
                W_j(j)=8/9;
            end
        end

        % K Matrix
        K=zeros(N_cell,N_cell);
        Func=zeros(1,NP);
        node_1=NN_tot(m,n);
        node_2=NN_tot(m+1,n);
        node_3=NN_tot(m+1,n+1);
        node_4=NN_tot(m,n+1);
        for i=1:1:N_cell
            for j=1:1:N_cell
                for N=1:1:NP
                    x_zeta=0.25*((1-eta(N))*delta_x_star+(1+eta(N))*delta_x_star);
                    y_zeta=0;
                    x_eta=0;
                    y_eta=0.25*((1-zeta(N))*delta_y_star+(1+zeta(N))*delta_y_star);
                    Jaco=abs(x_zeta*y_eta-y_zeta*x_eta);
                    if i==1
                        a=-1; b=-1;
                    elseif i==2
                        a=1; b=-1;
                    elseif i==3
                        a=1; b=1;
                    elseif i==4
                        a=-1; b=1;
                    end
                    N_i_zeta=0.25*a*(1+eta(N)*b);
                    N_i_eta=0.25*b*(1+zeta(N)*a);
                    if j==1
                        a=-1; b=-1;
                    elseif j==2
                        a=1; b=-1;
                    elseif j==3
                        a=1; b=1;
                    elseif j==4
                        a=-1; b=1;
                    end
                    N_j_zeta=0.25*a*(1+eta(N)*b);
                    N_j_eta=0.25*b*(1+zeta(N)*a);
                    N_i_x=(1/Jaco)*(y_eta*N_i_zeta-y_zeta*N_i_eta);
                    N_j_x=(1/Jaco)*(y_eta*N_j_zeta-y_zeta*N_j_eta);
                    N_i_y=(1/Jaco)*(x_eta*N_i_zeta-x_zeta*N_i_eta);
                    N_j_y=(1/Jaco)*(x_eta*N_j_zeta-x_zeta*N_j_eta);
                    Func(N)=Jaco*(N_i_x*N_j_x+N_i_y*N_j_y);
                end
            K(i,j)=sum(W_i.*W_j.*Func);
            end
        end
        K_global_a(node_1,node_1)=K_global_a(node_1,node_1)+K(1,1);
        K_global_a(node_1,node_2)=K_global_a(node_1,node_2)+K(1,2);
        K_global_a(node_1,node_3)=K_global_a(node_1,node_3)+K(1,3);
        K_global_a(node_1,node_4)=K_global_a(node_1,node_4)+K(1,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_a(node_2,node_1)=K_global_a(node_2,node_1)+K(2,1);
        K_global_a(node_2,node_2)=K_global_a(node_2,node_2)+K(2,2);
        K_global_a(node_2,node_3)=K_global_a(node_2,node_3)+K(2,3);
        K_global_a(node_2,node_4)=K_global_a(node_2,node_4)+K(2,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_a(node_3,node_1)=K_global_a(node_3,node_1)+K(3,1);
        K_global_a(node_3,node_2)=K_global_a(node_3,node_2)+K(3,2);
        K_global_a(node_3,node_3)=K_global_a(node_3,node_3)+K(3,3);
        K_global_a(node_3,node_4)=K_global_a(node_3,node_4)+K(3,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_a(node_4,node_1)=K_global_a(node_4,node_1)+K(4,1);
        K_global_a(node_4,node_2)=K_global_a(node_4,node_2)+K(4,2);
        K_global_a(node_4,node_3)=K_global_a(node_4,node_3)+K(4,3);
        K_global_a(node_4,node_4)=K_global_a(node_4,node_4)+K(4,4);
    end
end

% Boundery Conditions
for i=1:size(N_symm)
    K_trans=K_global_a(N_symm(i),N_symm(i));
    K_global_a(N_symm(i),:)=0;
    K_global_a(:,N_symm(i))=0;
    K_global_a(N_symm(i),N_symm(i))=K_trans;
end
for i=1:size(N_inlet)
    K_trans=K_global_a(N_inlet(i),N_inlet(i));
    K_global_a(N_inlet(i),:)=0;
    RHS_a(:,1)=RHS_a(:,1)-K_global_a(:,N_inlet(i))*Y_star_vec(i);
    K_global_a(:,N_inlet(i))=0;
    K_global_a(N_inlet(i),N_inlet(i))=K_trans;
    RHS_a(N_inlet(i),1)=K_trans*Y_star_vec(i);
end
for i=1:size(N_H_P1)
    K_trans=K_global_a(N_H_P1(i),N_H_P1(i));
    K_global_a(N_H_P1(i),:)=0;
    RHS_a(:,1)=RHS_a(:,1)-K_global_a(:,N_H_P1(i));
    K_global_a(:,N_H_P1(i))=0;
    K_global_a(N_H_P1(i),N_H_P1(i))=K_trans;
    RHS_a(N_H_P1(i),1)=K_trans;
end
for i=1:size(N_H_P2)
    K_trans=K_global_a(N_H_P2(i),N_H_P2(i));
    K_global_a(N_H_P2(i),:)=0;
    RHS_a(:,1)=RHS_a(:,1)-K_global_a(:,N_H_P2(i));
    K_global_a(:,N_H_P2(i))=0;
    K_global_a(N_H_P2(i),N_H_P2(i))=K_trans;
    RHS_a(N_H_P2(i),1)=K_trans;
end
for i=1:size(N_ver)
    K_trans=K_global_a(N_ver(i),N_ver(i));
    K_global_a(N_ver(i),:)=0;
    RHS_a(:,1)=RHS_a(:,1)-K_global_a(:,N_ver(i));
    K_global_a(:,N_ver(i))=0;
    K_global_a(N_ver(i),N_ver(i))=K_trans;
    RHS_a(N_ver(i),1)=K_trans;
end
epsi_t_a=linsolve(K_global_a,RHS_a);
epsi_a=NN_tot;
for i=1:length(epsi_t_a)
    [row,col] = find(epsi_a==i);
    epsi_a(row,col)=epsi_t_a(i,1);
end



%% Solver For Case b
K_global_b=zeros(max(max(NN_tot)),max(max(NN_tot)));
RHS_b=zeros(max(max(NN_tot)),1);
for n=1:NE_y_P1
    if n<=NE_y_P2
        var_N=NE_x_P1+NE_x_P2;
    elseif n>NE_y_P2
        var_N=NE_x_P1;
    end
    for m=1:var_N
        % Weights Calculation
        W_i=zeros(1,NP);
        for i=1:1:NP
            if zeta(i)~=0
                W_i(i)=5/9;
            elseif zeta(i)==0
                W_i(i)=8/9;
            end
        end
        W_j=zeros(1,NP);
        for j=1:1:NP
            if eta(j)~=0
                W_j(j)=5/9;
            elseif eta(j)==0
                W_j(j)=8/9;
            end
        end

        % K Matrix
        K=zeros(N_cell,N_cell);
        Func=zeros(1,NP);
        node_1=NN_tot(m,n);
        node_2=NN_tot(m+1,n);
        node_3=NN_tot(m+1,n+1);
        node_4=NN_tot(m,n+1);
        for i=1:1:N_cell
            for j=1:1:N_cell
                for N=1:1:NP
                    x_zeta=0.25*((1-eta(N))*delta_x_star+(1+eta(N))*delta_x_star);
                    y_zeta=0;
                    x_eta=0;
                    y_eta=0.25*((1-zeta(N))*delta_y_star+(1+zeta(N))*delta_y_star);
                    Jaco=abs(x_zeta*y_eta-y_zeta*x_eta);
                    if i==1
                        a=-1; b=-1;
                    elseif i==2
                        a=1; b=-1;
                    elseif i==3
                        a=1; b=1;
                    elseif i==4
                        a=-1; b=1;
                    end
                    N_i_zeta=0.25*a*(1+eta(N)*b);
                    N_i_eta=0.25*b*(1+zeta(N)*a);
                    N_i=0.25*(1+a*zeta(N))*(1+b*eta(N));
                    if j==1
                        a=-1; b=-1;
                    elseif j==2
                        a=1; b=-1;
                    elseif j==3
                        a=1; b=1;
                    elseif j==4
                        a=-1; b=1;
                    end
                    N_j_zeta=0.25*a*(1+eta(N)*b);
                    N_j_eta=0.25*b*(1+zeta(N)*a);
                    N_j=0.25*(1+a*zeta(N))*(1+b*eta(N));
                    N_i_x=(1/Jaco)*(y_eta*N_i_zeta-y_zeta*N_i_eta);
                    N_j_x=(1/Jaco)*(y_eta*N_j_zeta-y_zeta*N_j_eta);
                    N_i_y=(1/Jaco)*(x_eta*N_i_zeta-x_zeta*N_i_eta);
                    N_j_y=(1/Jaco)*(x_eta*N_j_zeta-x_zeta*N_j_eta);
                    Func(N)=Jaco*(N_i_x*N_j_x+N_i_y*N_j_y+c_epsi*N_i*N_j);
                end
            K(i,j)=sum(W_i.*W_j.*Func);
            end
        end
        K_global_b(node_1,node_1)=K_global_b(node_1,node_1)+K(1,1);
        K_global_b(node_1,node_2)=K_global_b(node_1,node_2)+K(1,2);
        K_global_b(node_1,node_3)=K_global_b(node_1,node_3)+K(1,3);
        K_global_b(node_1,node_4)=K_global_b(node_1,node_4)+K(1,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_b(node_2,node_1)=K_global_b(node_2,node_1)+K(2,1);
        K_global_b(node_2,node_2)=K_global_b(node_2,node_2)+K(2,2);
        K_global_b(node_2,node_3)=K_global_b(node_2,node_3)+K(2,3);
        K_global_b(node_2,node_4)=K_global_b(node_2,node_4)+K(2,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_b(node_3,node_1)=K_global_b(node_3,node_1)+K(3,1);
        K_global_b(node_3,node_2)=K_global_b(node_3,node_2)+K(3,2);
        K_global_b(node_3,node_3)=K_global_b(node_3,node_3)+K(3,3);
        K_global_b(node_3,node_4)=K_global_b(node_3,node_4)+K(3,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_b(node_4,node_1)=K_global_b(node_4,node_1)+K(4,1);
        K_global_b(node_4,node_2)=K_global_b(node_4,node_2)+K(4,2);
        K_global_b(node_4,node_3)=K_global_b(node_4,node_3)+K(4,3);
        K_global_b(node_4,node_4)=K_global_b(node_4,node_4)+K(4,4);
    end
end

% Boundery Conditions
for i=1:size(N_symm)
    K_trans=K_global_b(N_symm(i),N_symm(i));
    K_global_b(N_symm(i),:)=0;
    K_global_b(:,N_symm(i))=0;
    K_global_b(N_symm(i),N_symm(i))=K_trans;
end
for i=1:size(N_inlet)
    K_trans=K_global_b(N_inlet(i),N_inlet(i));
    K_global_b(N_inlet(i),:)=0;
    RHS_b(:,1)=RHS_b(:,1)-K_global_b(:,N_inlet(i))*sin(pi*Y_star_vec(i)/2);
    K_global_b(:,N_inlet(i))=0;
    K_global_b(N_inlet(i),N_inlet(i))=K_trans;
    RHS_b(N_inlet(i),1)=K_trans*sin(pi*Y_star_vec(i)/2);
end
for i=1:size(N_H_P1)
    K_trans=K_global_b(N_H_P1(i),N_H_P1(i));
    K_global_b(N_H_P1(i),:)=0;
    RHS_b(:,1)=RHS_b(:,1)-K_global_b(:,N_H_P1(i));
    K_global_b(:,N_H_P1(i))=0;
    K_global_b(N_H_P1(i),N_H_P1(i))=K_trans;
    RHS_b(N_H_P1(i),1)=K_trans;
end
for i=1:size(N_H_P2)
    K_trans=K_global_b(N_H_P2(i),N_H_P2(i));
    K_global_b(N_H_P2(i),:)=0;
    RHS_b(:,1)=RHS_b(:,1)-K_global_b(:,N_H_P2(i));
    K_global_b(:,N_H_P2(i))=0;
    K_global_b(N_H_P2(i),N_H_P2(i))=K_trans;
    RHS_b(N_H_P2(i),1)=K_trans;
end
for i=1:size(N_ver)
    K_trans=K_global_b(N_ver(i),N_ver(i));
    K_global_b(N_ver(i),:)=0;
    RHS_b(:,1)=RHS_b(:,1)-K_global_b(:,N_ver(i));
    K_global_b(:,N_ver(i))=0;
    K_global_b(N_ver(i),N_ver(i))=K_trans;
    RHS_b(N_ver(i),1)=K_trans;
end
epsi_t_b=linsolve(K_global_b,RHS_b);
epsi_b=NN_tot;
for i=1:length(epsi_t_b)
    [row,col] = find(epsi_b==i);
    epsi_b(row,col)=epsi_t_b(i,1);
end

%% Velocity
%Case a
%u
K_global_epsi_u_a=zeros(max(max(NN_tot)),max(max(NN_tot)));
for n=1:NE_y_P1
    if n<=NE_y_P2
        var_N=NE_x_P1+NE_x_P2;
    elseif n>NE_y_P2
        var_N=NE_x_P1;
    end
    for m=1:var_N
        % Weights Calculation
        W_i=zeros(1,NP);
        for i=1:1:NP
            if zeta(i)~=0
                W_i(i)=5/9;
            elseif zeta(i)==0
                W_i(i)=8/9;
            end
        end
        W_j=zeros(1,NP);
        for j=1:1:NP
            if eta(j)~=0
                W_j(j)=5/9;
            elseif eta(j)==0
                W_j(j)=8/9;
            end
        end

        % K Matrix
        K=zeros(N_cell,N_cell);
        Func=zeros(1,NP);
        node_1=NN_tot(m,n);
        node_2=NN_tot(m+1,n);
        node_3=NN_tot(m+1,n+1);
        node_4=NN_tot(m,n+1);
        for i=1:1:N_cell
            for j=1:1:N_cell
                for N=1:1:NP
                    x_zeta=0.25*((1-eta(N))*delta_x_star+(1+eta(N))*delta_x_star);
                    y_zeta=0;
                    x_eta=0;
                    y_eta=0.25*((1-zeta(N))*delta_y_star+(1+zeta(N))*delta_y_star);
                    Jaco=abs(x_zeta*y_eta-y_zeta*x_eta);
                    if i==1
                        a=-1; b=-1;
                    elseif i==2
                        a=1; b=-1;
                    elseif i==3
                        a=1; b=1;
                    elseif i==4
                        a=-1; b=1;
                    end
                    N_i=0.25*(1+a*zeta(N))*(1+b*eta(N));
                    if j==1
                        a=-1; b=-1;
                    elseif j==2
                        a=1; b=-1;
                    elseif j==3
                        a=1; b=1;
                    elseif j==4
                        a=-1; b=1;
                    end
                    N_j_zeta=0.25*a*(1+eta(N)*b);
                    N_j_eta=0.25*b*(1+zeta(N)*a);
                    N_j_y=(1/Jaco)*(x_eta*N_j_zeta-x_zeta*N_j_eta);
                    Func(N)=Jaco*(N_i*N_j_y);
                end
            K(i,j)=sum(W_i.*W_j.*Func);
            end
        end
        K_global_epsi_u_a(node_1,node_1)=K_global_epsi_u_a(node_1,node_1)+K(1,1);
        K_global_epsi_u_a(node_1,node_2)=K_global_epsi_u_a(node_1,node_2)+K(1,2);
        K_global_epsi_u_a(node_1,node_3)=K_global_epsi_u_a(node_1,node_3)+K(1,3);
        K_global_epsi_u_a(node_1,node_4)=K_global_epsi_u_a(node_1,node_4)+K(1,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_epsi_u_a(node_2,node_1)=K_global_epsi_u_a(node_2,node_1)+K(2,1);
        K_global_epsi_u_a(node_2,node_2)=K_global_epsi_u_a(node_2,node_2)+K(2,2);
        K_global_epsi_u_a(node_2,node_3)=K_global_epsi_u_a(node_2,node_3)+K(2,3);
        K_global_epsi_u_a(node_2,node_4)=K_global_epsi_u_a(node_2,node_4)+K(2,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_epsi_u_a(node_3,node_1)=K_global_epsi_u_a(node_3,node_1)+K(3,1);
        K_global_epsi_u_a(node_3,node_2)=K_global_epsi_u_a(node_3,node_2)+K(3,2);
        K_global_epsi_u_a(node_3,node_3)=K_global_epsi_u_a(node_3,node_3)+K(3,3);
        K_global_epsi_u_a(node_3,node_4)=K_global_epsi_u_a(node_3,node_4)+K(3,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_epsi_u_a(node_4,node_1)=K_global_epsi_u_a(node_4,node_1)+K(4,1);
        K_global_epsi_u_a(node_4,node_2)=K_global_epsi_u_a(node_4,node_2)+K(4,2);
        K_global_epsi_u_a(node_4,node_3)=K_global_epsi_u_a(node_4,node_3)+K(4,3);
        K_global_epsi_u_a(node_4,node_4)=K_global_epsi_u_a(node_4,node_4)+K(4,4);
    end
end
K_global_u=zeros(max(max(NN_tot)),max(max(NN_tot)));
for n=1:NE_y_P1
    if n<=NE_y_P2
        var_N=NE_x_P1+NE_x_P2;
    elseif n>NE_y_P2
        var_N=NE_x_P1;
    end
    for m=1:var_N
        % Weights Calculation
        W_i=zeros(1,NP);
        for i=1:1:NP
            if zeta(i)~=0
                W_i(i)=5/9;
            elseif zeta(i)==0
                W_i(i)=8/9;
            end
        end
        W_j=zeros(1,NP);
        for j=1:1:NP
            if eta(j)~=0
                W_j(j)=5/9;
            elseif eta(j)==0
                W_j(j)=8/9;
            end
        end

        % K Matrix
        K=zeros(N_cell,N_cell);
        Func=zeros(1,NP);
        node_1=NN_tot(m,n);
        node_2=NN_tot(m+1,n);
        node_3=NN_tot(m+1,n+1);
        node_4=NN_tot(m,n+1);
        for i=1:1:N_cell
            for j=1:1:N_cell
                for N=1:1:NP
                    x_zeta=0.25*((1-eta(N))*delta_x_star+(1+eta(N))*delta_x_star);
                    y_zeta=0;
                    x_eta=0;
                    y_eta=0.25*((1-zeta(N))*delta_y_star+(1+zeta(N))*delta_y_star);
                    Jaco=abs(x_zeta*y_eta-y_zeta*x_eta);
                    if i==1
                        a=-1; b=-1;
                    elseif i==2
                        a=1; b=-1;
                    elseif i==3
                        a=1; b=1;
                    elseif i==4
                        a=-1; b=1;
                    end
                    N_i=0.25*(1+a*zeta(N))*(1+b*eta(N));
                    if j==1
                        a=-1; b=-1;
                    elseif j==2
                        a=1; b=-1;
                    elseif j==3
                        a=1; b=1;
                    elseif j==4
                        a=-1; b=1;
                    end
                    N_j=0.25*(1+a*zeta(N))*(1+b*eta(N));
                    Func(N)=Jaco*(N_i*N_j);
                end
            K(i,j)=sum(W_i.*W_j.*Func);
            end
        end
        K_global_u(node_1,node_1)=K_global_u(node_1,node_1)+K(1,1);
        K_global_u(node_1,node_2)=K_global_u(node_1,node_2)+K(1,2);
        K_global_u(node_1,node_3)=K_global_u(node_1,node_3)+K(1,3);
        K_global_u(node_1,node_4)=K_global_u(node_1,node_4)+K(1,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_u(node_2,node_1)=K_global_u(node_2,node_1)+K(2,1);
        K_global_u(node_2,node_2)=K_global_u(node_2,node_2)+K(2,2);
        K_global_u(node_2,node_3)=K_global_u(node_2,node_3)+K(2,3);
        K_global_u(node_2,node_4)=K_global_u(node_2,node_4)+K(2,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_u(node_3,node_1)=K_global_u(node_3,node_1)+K(3,1);
        K_global_u(node_3,node_2)=K_global_u(node_3,node_2)+K(3,2);
        K_global_u(node_3,node_3)=K_global_u(node_3,node_3)+K(3,3);
        K_global_u(node_3,node_4)=K_global_u(node_3,node_4)+K(3,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_u(node_4,node_1)=K_global_u(node_4,node_1)+K(4,1);
        K_global_u(node_4,node_2)=K_global_u(node_4,node_2)+K(4,2);
        K_global_u(node_4,node_3)=K_global_u(node_4,node_3)+K(4,3);
        K_global_u(node_4,node_4)=K_global_u(node_4,node_4)+K(4,4);
    end
end
RHS_u=-1*K_global_epsi_u_a*epsi_t_a;
comp_u=linsolve(K_global_u,RHS_u);
U_a=NN_tot;
for i=1:length(comp_u)
    [row,col]=find(U_a==i);
    U_a(row,col)=comp_u(i,1);
end

%V
K_global_epsi_v=zeros(max(max(NN_tot)),max(max(NN_tot)));
for n=1:NE_y_P1
    if n<=NE_y_P2
        var_N=NE_x_P1+NE_x_P2;
    elseif n>NE_y_P2
        var_N=NE_x_P1;
    end
    for m=1:var_N
        % Weights Calculation
        W_i=zeros(1,NP);
        for i=1:1:NP
            if zeta(i)~=0
                W_i(i)=5/9;
            elseif zeta(i)==0
                W_i(i)=8/9;
            end
        end
        W_j=zeros(1,NP);
        for j=1:1:NP
            if eta(j)~=0
                W_j(j)=5/9;
            elseif eta(j)==0
                W_j(j)=8/9;
            end
        end

        % K Matrix
        K=zeros(N_cell,N_cell);
        Func=zeros(1,NP);
        node_1=NN_tot(m,n);
        node_2=NN_tot(m+1,n);
        node_3=NN_tot(m+1,n+1);
        node_4=NN_tot(m,n+1);
        for i=1:1:N_cell
            for j=1:1:N_cell
                for N=1:1:NP
                    x_zeta=0.25*((1-eta(N))*delta_x_star+(1+eta(N))*delta_x_star);
                    y_zeta=0;
                    x_eta=0;
                    y_eta=0.25*((1-zeta(N))*delta_y_star+(1+zeta(N))*delta_y_star);
                    Jaco=abs(x_zeta*y_eta-y_zeta*x_eta);
                    if i==1
                        a=-1; b=-1;
                    elseif i==2
                        a=1; b=-1;
                    elseif i==3
                        a=1; b=1;
                    elseif i==4
                        a=-1; b=1;
                    end
                    N_i=0.25*(1+a*zeta(N))*(1+b*eta(N));
                    if j==1
                        a=-1; b=-1;
                    elseif j==2
                        a=1; b=-1;
                    elseif j==3
                        a=1; b=1;
                    elseif j==4
                        a=-1; b=1;
                    end
                    N_j_zeta=0.25*a*(1+eta(N)*b);
                    N_j_eta=0.25*b*(1+zeta(N)*a);
                    N_j_x=(1/Jaco)*(y_eta*N_j_zeta-y_zeta*N_j_eta);
                    Func(N)=Jaco*(N_i*N_j_x);
                end
            K(i,j)=sum(W_i.*W_j.*Func);
            end
        end
        K_global_epsi_v(node_1,node_1)=K_global_epsi_v(node_1,node_1)+K(1,1);
        K_global_epsi_v(node_1,node_2)=K_global_epsi_v(node_1,node_2)+K(1,2);
        K_global_epsi_v(node_1,node_3)=K_global_epsi_v(node_1,node_3)+K(1,3);
        K_global_epsi_v(node_1,node_4)=K_global_epsi_v(node_1,node_4)+K(1,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_epsi_v(node_2,node_1)=K_global_epsi_v(node_2,node_1)+K(2,1);
        K_global_epsi_v(node_2,node_2)=K_global_epsi_v(node_2,node_2)+K(2,2);
        K_global_epsi_v(node_2,node_3)=K_global_epsi_v(node_2,node_3)+K(2,3);
        K_global_epsi_v(node_2,node_4)=K_global_epsi_v(node_2,node_4)+K(2,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_epsi_v(node_3,node_1)=K_global_epsi_v(node_3,node_1)+K(3,1);
        K_global_epsi_v(node_3,node_2)=K_global_epsi_v(node_3,node_2)+K(3,2);
        K_global_epsi_v(node_3,node_3)=K_global_epsi_v(node_3,node_3)+K(3,3);
        K_global_epsi_v(node_3,node_4)=K_global_epsi_v(node_3,node_4)+K(3,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_epsi_v(node_4,node_1)=K_global_epsi_v(node_4,node_1)+K(4,1);
        K_global_epsi_v(node_4,node_2)=K_global_epsi_v(node_4,node_2)+K(4,2);
        K_global_epsi_v(node_4,node_3)=K_global_epsi_v(node_4,node_3)+K(4,3);
        K_global_epsi_v(node_4,node_4)=K_global_epsi_v(node_4,node_4)+K(4,4);
    end
end

RHS_v=K_global_epsi_v*epsi_t_a;
comp_v=-1*linsolve(K_global_u,RHS_v);
V_a=NN_tot;
for i=1:length(comp_v)
    [row,col]=find(V_a==i);
    V_a(row,col)=comp_v(i,1);
end

%Case b
%u
K_global_epsi_u_b=zeros(max(max(NN_tot)),max(max(NN_tot)));
for n=1:NE_y_P1
    if n<=NE_y_P2
        var_N=NE_x_P1+NE_x_P2;
    elseif n>NE_y_P2
        var_N=NE_x_P1;
    end
    for m=1:var_N
        % Weights Calculation
        W_i=zeros(1,NP);
        for i=1:1:NP
            if zeta(i)~=0
                W_i(i)=5/9;
            elseif zeta(i)==0
                W_i(i)=8/9;
            end
        end
        W_j=zeros(1,NP);
        for j=1:1:NP
            if eta(j)~=0
                W_j(j)=5/9;
            elseif eta(j)==0
                W_j(j)=8/9;
            end
        end

        % K Matrix
        K=zeros(N_cell,N_cell);
        Func=zeros(1,NP);
        node_1=NN_tot(m,n);
        node_2=NN_tot(m+1,n);
        node_3=NN_tot(m+1,n+1);
        node_4=NN_tot(m,n+1);
        for i=1:1:N_cell
            for j=1:1:N_cell
                for N=1:1:NP
                    x_zeta=0.25*((1-eta(N))*delta_x_star+(1+eta(N))*delta_x_star);
                    y_zeta=0;
                    x_eta=0;
                    y_eta=0.25*((1-zeta(N))*delta_y_star+(1+zeta(N))*delta_y_star);
                    Jaco=abs(x_zeta*y_eta-y_zeta*x_eta);
                    if i==1
                        a=-1; b=-1;
                    elseif i==2
                        a=1; b=-1;
                    elseif i==3
                        a=1; b=1;
                    elseif i==4
                        a=-1; b=1;
                    end
                    N_i=0.25*(1+a*zeta(N))*(1+b*eta(N));
                    if j==1
                        a=-1; b=-1;
                    elseif j==2
                        a=1; b=-1;
                    elseif j==3
                        a=1; b=1;
                    elseif j==4
                        a=-1; b=1;
                    end
                    N_j_zeta=0.25*a*(1+eta(N)*b);
                    N_j_eta=0.25*b*(1+zeta(N)*a);
                    N_j_y=(1/Jaco)*(x_eta*N_j_zeta-x_zeta*N_j_eta);
                    Func(N)=Jaco*(N_i*N_j_y);
                end
            K(i,j)=sum(W_i.*W_j.*Func);
            end
        end
        K_global_epsi_u_b(node_1,node_1)=K_global_epsi_u_b(node_1,node_1)+K(1,1);
        K_global_epsi_u_b(node_1,node_2)=K_global_epsi_u_b(node_1,node_2)+K(1,2);
        K_global_epsi_u_b(node_1,node_3)=K_global_epsi_u_b(node_1,node_3)+K(1,3);
        K_global_epsi_u_b(node_1,node_4)=K_global_epsi_u_b(node_1,node_4)+K(1,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_epsi_u_b(node_2,node_1)=K_global_epsi_u_b(node_2,node_1)+K(2,1);
        K_global_epsi_u_b(node_2,node_2)=K_global_epsi_u_b(node_2,node_2)+K(2,2);
        K_global_epsi_u_b(node_2,node_3)=K_global_epsi_u_b(node_2,node_3)+K(2,3);
        K_global_epsi_u_b(node_2,node_4)=K_global_epsi_u_b(node_2,node_4)+K(2,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_epsi_u_b(node_3,node_1)=K_global_epsi_u_b(node_3,node_1)+K(3,1);
        K_global_epsi_u_b(node_3,node_2)=K_global_epsi_u_b(node_3,node_2)+K(3,2);
        K_global_epsi_u_b(node_3,node_3)=K_global_epsi_u_b(node_3,node_3)+K(3,3);
        K_global_epsi_u_b(node_3,node_4)=K_global_epsi_u_b(node_3,node_4)+K(3,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_epsi_u_b(node_4,node_1)=K_global_epsi_u_b(node_4,node_1)+K(4,1);
        K_global_epsi_u_b(node_4,node_2)=K_global_epsi_u_b(node_4,node_2)+K(4,2);
        K_global_epsi_u_b(node_4,node_3)=K_global_epsi_u_b(node_4,node_3)+K(4,3);
        K_global_epsi_u_b(node_4,node_4)=K_global_epsi_u_b(node_4,node_4)+K(4,4);
    end
end
K_global_u_b=zeros(max(max(NN_tot)),max(max(NN_tot)));
for n=1:NE_y_P1
    if n<=NE_y_P2
        var_N=NE_x_P1+NE_x_P2;
    elseif n>NE_y_P2
        var_N=NE_x_P1;
    end
    for m=1:var_N
        % Weights Calculation
        W_i=zeros(1,NP);
        for i=1:1:NP
            if zeta(i)~=0
                W_i(i)=5/9;
            elseif zeta(i)==0
                W_i(i)=8/9;
            end
        end
        W_j=zeros(1,NP);
        for j=1:1:NP
            if eta(j)~=0
                W_j(j)=5/9;
            elseif eta(j)==0
                W_j(j)=8/9;
            end
        end

        % K Matrix
        K=zeros(N_cell,N_cell);
        Func=zeros(1,NP);
        node_1=NN_tot(m,n);
        node_2=NN_tot(m+1,n);
        node_3=NN_tot(m+1,n+1);
        node_4=NN_tot(m,n+1);
        for i=1:1:N_cell
            for j=1:1:N_cell
                for N=1:1:NP
                    x_zeta=0.25*((1-eta(N))*delta_x_star+(1+eta(N))*delta_x_star);
                    y_zeta=0;
                    x_eta=0;
                    y_eta=0.25*((1-zeta(N))*delta_y_star+(1+zeta(N))*delta_y_star);
                    Jaco=abs(x_zeta*y_eta-y_zeta*x_eta);
                    if i==1
                        a=-1; b=-1;
                    elseif i==2
                        a=1; b=-1;
                    elseif i==3
                        a=1; b=1;
                    elseif i==4
                        a=-1; b=1;
                    end
                    N_i=0.25*(1+a*zeta(N))*(1+b*eta(N));
                    if j==1
                        a=-1; b=-1;
                    elseif j==2
                        a=1; b=-1;
                    elseif j==3
                        a=1; b=1;
                    elseif j==4
                        a=-1; b=1;
                    end
                    N_j=0.25*(1+a*zeta(N))*(1+b*eta(N));
                    Func(N)=Jaco*(N_i*N_j);
                end
            K(i,j)=sum(W_i.*W_j.*Func);
            end
        end
        K_global_u_b(node_1,node_1)=K_global_u_b(node_1,node_1)+K(1,1);
        K_global_u_b(node_1,node_2)=K_global_u_b(node_1,node_2)+K(1,2);
        K_global_u_b(node_1,node_3)=K_global_u_b(node_1,node_3)+K(1,3);
        K_global_u_b(node_1,node_4)=K_global_u_b(node_1,node_4)+K(1,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_u_b(node_2,node_1)=K_global_u_b(node_2,node_1)+K(2,1);
        K_global_u_b(node_2,node_2)=K_global_u_b(node_2,node_2)+K(2,2);
        K_global_u_b(node_2,node_3)=K_global_u_b(node_2,node_3)+K(2,3);
        K_global_u_b(node_2,node_4)=K_global_u_b(node_2,node_4)+K(2,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_u_b(node_3,node_1)=K_global_u_b(node_3,node_1)+K(3,1);
        K_global_u_b(node_3,node_2)=K_global_u_b(node_3,node_2)+K(3,2);
        K_global_u_b(node_3,node_3)=K_global_u_b(node_3,node_3)+K(3,3);
        K_global_u_b(node_3,node_4)=K_global_u_b(node_3,node_4)+K(3,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_u_b(node_4,node_1)=K_global_u_b(node_4,node_1)+K(4,1);
        K_global_u_b(node_4,node_2)=K_global_u_b(node_4,node_2)+K(4,2);
        K_global_u_b(node_4,node_3)=K_global_u_b(node_4,node_3)+K(4,3);
        K_global_u_b(node_4,node_4)=K_global_u_b(node_4,node_4)+K(4,4);
    end
end
RHS_u_b=-1*K_global_epsi_u_b*epsi_t_b;
comp_u_b=linsolve(K_global_u_b,RHS_u_b);
U_b=NN_tot;
for i=1:length(comp_u_b)
    [row,col]=find(U_b==i);
    U_b(row,col)=comp_u_b(i,1);
end

%V
K_global_epsi_v_b=zeros(max(max(NN_tot)),max(max(NN_tot)));
for n=1:NE_y_P1
    if n<=NE_y_P2
        var_N=NE_x_P1+NE_x_P2;
    elseif n>NE_y_P2
        var_N=NE_x_P1;
    end
    for m=1:var_N
        % Weights Calculation
        W_i=zeros(1,NP);
        for i=1:1:NP
            if zeta(i)~=0
                W_i(i)=5/9;
            elseif zeta(i)==0
                W_i(i)=8/9;
            end
        end
        W_j=zeros(1,NP);
        for j=1:1:NP
            if eta(j)~=0
                W_j(j)=5/9;
            elseif eta(j)==0
                W_j(j)=8/9;
            end
        end

        % K Matrix
        K=zeros(N_cell,N_cell);
        Func=zeros(1,NP);
        node_1=NN_tot(m,n);
        node_2=NN_tot(m+1,n);
        node_3=NN_tot(m+1,n+1);
        node_4=NN_tot(m,n+1);
        for i=1:1:N_cell
            for j=1:1:N_cell
                for N=1:1:NP
                    x_zeta=0.25*((1-eta(N))*delta_x_star+(1+eta(N))*delta_x_star);
                    y_zeta=0;
                    x_eta=0;
                    y_eta=0.25*((1-zeta(N))*delta_y_star+(1+zeta(N))*delta_y_star);
                    Jaco=abs(x_zeta*y_eta-y_zeta*x_eta);
                    if i==1
                        a=-1; b=-1;
                    elseif i==2
                        a=1; b=-1;
                    elseif i==3
                        a=1; b=1;
                    elseif i==4
                        a=-1; b=1;
                    end
                    N_i=0.25*(1+a*zeta(N))*(1+b*eta(N));
                    if j==1
                        a=-1; b=-1;
                    elseif j==2
                        a=1; b=-1;
                    elseif j==3
                        a=1; b=1;
                    elseif j==4
                        a=-1; b=1;
                    end
                    N_j_zeta=0.25*a*(1+eta(N)*b);
                    N_j_eta=0.25*b*(1+zeta(N)*a);
                    N_j_x=(1/Jaco)*(y_eta*N_j_zeta-y_zeta*N_j_eta);
                    Func(N)=Jaco*(N_i*N_j_x);
                end
            K(i,j)=sum(W_i.*W_j.*Func);
            end
        end
        K_global_epsi_v_b(node_1,node_1)=K_global_epsi_v_b(node_1,node_1)+K(1,1);
        K_global_epsi_v_b(node_1,node_2)=K_global_epsi_v_b(node_1,node_2)+K(1,2);
        K_global_epsi_v_b(node_1,node_3)=K_global_epsi_v_b(node_1,node_3)+K(1,3);
        K_global_epsi_v_b(node_1,node_4)=K_global_epsi_v_b(node_1,node_4)+K(1,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_epsi_v_b(node_2,node_1)=K_global_epsi_v_b(node_2,node_1)+K(2,1);
        K_global_epsi_v_b(node_2,node_2)=K_global_epsi_v_b(node_2,node_2)+K(2,2);
        K_global_epsi_v_b(node_2,node_3)=K_global_epsi_v_b(node_2,node_3)+K(2,3);
        K_global_epsi_v_b(node_2,node_4)=K_global_epsi_v_b(node_2,node_4)+K(2,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_epsi_v_b(node_3,node_1)=K_global_epsi_v_b(node_3,node_1)+K(3,1);
        K_global_epsi_v_b(node_3,node_2)=K_global_epsi_v_b(node_3,node_2)+K(3,2);
        K_global_epsi_v_b(node_3,node_3)=K_global_epsi_v_b(node_3,node_3)+K(3,3);
        K_global_epsi_v_b(node_3,node_4)=K_global_epsi_v_b(node_3,node_4)+K(3,4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K_global_epsi_v_b(node_4,node_1)=K_global_epsi_v_b(node_4,node_1)+K(4,1);
        K_global_epsi_v_b(node_4,node_2)=K_global_epsi_v_b(node_4,node_2)+K(4,2);
        K_global_epsi_v_b(node_4,node_3)=K_global_epsi_v_b(node_4,node_3)+K(4,3);
        K_global_epsi_v_b(node_4,node_4)=K_global_epsi_v_b(node_4,node_4)+K(4,4);
    end
end

RHS_v_b=K_global_epsi_v_b*epsi_t_b;
comp_v_b=-1*linsolve(K_global_u_b,RHS_v_b);
V_b=NN_tot;
for i=1:length(comp_v_b)
    [row,col]=find(V_b==i);
    V_b(row,col)=comp_v_b(i,1);
end

%% Meshing Plot
X=flip(X.');
Y=flip(Y.');
X_P1=flip(X_P1.');
Y_P1=flip(Y_P1.');
X_P2=flip(X_P2.');
Y_P2=flip(Y_P2.');
[X_1,Y_1]=meshgrid(X_P1,Y_P1);
[X_2,Y_2]=meshgrid(X_P2,Y_P2);
figure
hold all
xlabel('X')
ylabel('Y')
title('Mesh Generation')
for i=1:size(X_1,1)
plot(X_1(i,:),Y_1(i,:),'k')
end
for i=1:size(X_1,2)
plot(X_1(:,i),Y_1(:,i),'k')
end
for i=1:size(X_2,1)
plot(X_2(i,:),Y_2(i,:),'k')
end
for i=1:size(X_2,2)
plot(X_2(:,i),Y_2(:,i),'k')
end

%% Results For Case a
X_star=flip(X_star.');
Y_star=flip(Y_star.');
epsi_a=flip(epsi_a.');
U_a=flip(U_a.');
V_a=flip(V_a.');
figure
contour(X_star,Y_star,epsi_a,50)
colorbar
xlabel('X')
ylabel('Y')
title('Streamlines for 2D Flow in a Pipe - Case a) Short Inlet Duct Length')
grid on

figure
U_a=U_a./sqrt(U_a.^2+V_a.^2);
V_a=V_a./sqrt(U_a.^2+V_a.^2);
quiver(X_star,Y_star,U_a,V_a,1)
grid on
ylim([0 1])
xlim([0 3])
xlabel('X')
ylabel('Y')
title('Velocity for 2D Flow in a Pipe - Case a) Short Inlet Duct Length')

figure
U_a=U_a./sqrt(U_a.^2+V_a.^2);
V_a=V_a./sqrt(U_a.^2+V_a.^2);
hold on
quiver(X_star,Y_star,U_a,V_a,1)
contour(X_star,Y_star,epsi_a,50)
grid on
ylim([0 1])
xlim([0 3])
xlabel('X')
ylabel('Y')
title('Velocity with Streamlines - Case a) Short Inlet Duct Length')
colorbar

%% Results For Case b
epsi_b=flip(epsi_b.');
U_b=flip(U_b.');
V_b=flip(V_b.');
figure
contour(X_star,Y_star,epsi_b,100)
colorbar
xlabel('X')
ylabel('Y')
title('Streamlines for 2D Flow in a Pipe - Case b) Long Inlet Duct Length')
grid on

figure
U_b=U_b./sqrt(U_b.^2+V_b.^2);
V_b=V_b./sqrt(U_b.^2+V_b.^2);
quiver(X_star,Y_star,U_b,V_b,1)
grid on
ylim([0 1])
xlim([0 3])
xlabel('X')
ylabel('Y')
title('Velocity for 2D Flow in a Pipe - Case b) Long Inlet Duct Length')

figure
U_b=U_b./sqrt(U_b.^2+V_b.^2);
V_b=V_b./sqrt(U_b.^2+V_b.^2);
hold on
quiver(X_star,Y_star,U_b,V_b,1)
contour(X_star,Y_star,epsi_b,50)
grid on
ylim([0 1])
xlim([0 3])
xlabel('X')
ylabel('Y')
title('Velocity with Streamlines - Case b) Long Inlet Duct Length')
colorbar