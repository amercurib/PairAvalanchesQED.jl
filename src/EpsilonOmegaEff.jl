module EpsilonOmegaEff

import LinearAlgebra as linalg
using Einsum

#computation of invariants
export computing_invariants

function computing_invariants(EB_T)
    nx,ny,nz = size(EB_T)[1:3]
    G_inv = zeros((nx,ny,nz))
    F_inv = zeros((nx,ny,nz))
    for i = 1:nx, j = 1:ny, k = 1:nz, m1 = 1:3
        F_inv[i,j,k] += (EB_T[i,j,k,m1]*EB_T[i,j,k,m1]- EB_T[i,j,k,m1+3]*EB_T[i,j,k,m1+3])/2
        G_inv[i,j,k] += EB_T[i,j,k,m1]*EB_T[i,j,k,m1+3]
    end
    Eps_inv = sqrt.(sqrt.(F_inv.*F_inv .+ G_inv.*G_inv) .+ F_inv)
    return (F_inv,Eps_inv)
end


#Make field strength tensor
export make_TF

function make_TF(EB_T)
    nx,ny,nz = size(EB_T)[1:3]
    TF = zeros((nx,ny,nz,4,4))
    for i = 1:nx, j = 1:ny, k = 1:nz
        TF[i,j,k,1,2:4] = EB_T[i,j,k,1:3]
        TF[i,j,k,2:4,1] = EB_T[i,j,k,1:3]
        TF[i,j,k,2,3] = EB_T[i,j,k,6]
        TF[i,j,k,3,2] = -EB_T[i,j,k,6]
        TF[i,j,k,2,4] = -EB_T[i,j,k,5]
        TF[i,j,k,4,2] = EB_T[i,j,k,5]
        TF[i,j,k,3,4] = EB_T[i,j,k,4]
        TF[i,j,k,4,3] = -EB_T[i,j,k,4]
    end
    return TF
end


#Make the square of the field strength tensor
export make_square_field_tensor

function make_square_field_tensor(TF)
    nx,ny,nz = size(TF)[1:3]
    TF2 = zeros((nx,ny,nz,4,4))
    for i = 1:nx,j = 1:ny, k = 1:nz, m1 = 1:4, m2 = 1:4, msum = 1:4
        TF2[i,j,k,m1,m2] = TF[i,j,k,m1,msum]*TF[i,j,k,msum,m2]
    end
    return TF2
end


#Compute the derivaties from the field
export computing_derivatives_EB

function computing_derivatives_EB(EB_T,dx,dy,dz)
    nx,ny,nz = size(EB_T)[1:3]
    DEB_T = zeros((nx,ny,nz,6,4))
    
    for i = 2:nx-1
        DEB_T[i,:,:,:,2] = (EB_T[i+1,:,:,:] .- EB_T[i-1,:,:,:])/(2*dx)
    end
    DEB_T[1,:,:,:,2] = (EB_T[2,:,:,:] .- EB_T[1,:,:,:])/dx
    DEB_T[nx,:,:,:,2] = (EB_T[nx,:,:,:] .- EB_T[nx-1,:,:,:])/dx

    for j = 2:ny-1
    	DEB_T[:,j,:,:,3] = (EB_T[:,j+1,:,:] .- EB_T[:,j-1,:,:])/(2*dy)
    end
    DEB_T[:,1,:,:,3] = (EB_T[:,2,:,:] .- EB_T[:,1,:,:])/dy
    DEB_T[:,ny,:,:,3] = (EB_T[:,ny,:,:] .- EB_T[:,ny-1,:,:])/dy

    for k = 2:nz-1
    	DEB_T[:,:,k,:,4] = (EB_T[:,:,k+1,:] .- EB_T[:,:,k-1,:])/(2*dz)
    end
    DEB_T[:,:,1,:,4] = (EB_T[:,:,2,:] .- EB_T[:,:,1,:])/dz
    DEB_T[:,:,nz,:,4] = (EB_T[:,:,nz,:] .- EB_T[:,:,nz-1,:])/dz

    #E time derivative
    DEB_T[:,:,:,1,1] = (DEB_T[:,:,:,6,3] .- DEB_T[:,:,:,5,4]) # index for derivative are not the same than for components of field 1 = Ex  , derivative 1 = d/dt
    DEB_T[:,:,:,2,1] = (DEB_T[:,:,:,4,4] .- DEB_T[:,:,:,6,2])
    DEB_T[:,:,:,3,1] = (DEB_T[:,:,:,5,2] .- DEB_T[:,:,:,4,3])
    #B time derivative
    DEB_T[:,:,:,4,1] =  .-(DEB_T[:,:,:,3,3] .- DEB_T[:,:,:,2,4])#Bx is index 4 ; By 5;  Bz 6
    DEB_T[:,:,:,5,1] =  .-(DEB_T[:,:,:,1,4] .- DEB_T[:,:,:,3,2])
    DEB_T[:,:,:,6,1] =  .-(DEB_T[:,:,:,2,2] .- DEB_T[:,:,:,1,3])
    return DEB_T
end


#Make the tensor with the derivatives
export make_derivative_tensor

function make_derivative_tensor(DEB_T)
    nx,ny,nz = size(DEB_T)[1:3]
    DTF = zeros((nx,ny,nz,4,4,4))
    for i = 1:nx, j = 1:ny, k = 1:nz, i_der = 1:4
        DTF[i,j,k,1,2:4,i_der] = DEB_T[i,j,k,1:3,i_der]
        DTF[i,j,k,2:4,1,i_der] = DEB_T[i,j,k,1:3,i_der]
        DTF[i,j,k,2,3,i_der] = DEB_T[i,j,k,6,i_der]
        DTF[i,j,k,3,2,i_der] = -DEB_T[i,j,k,6,i_der]
        DTF[i,j,k,2,4,i_der] = -DEB_T[i,j,k,5,i_der]
        DTF[i,j,k,4,2,i_der] = DEB_T[i,j,k,5,i_der]
        DTF[i,j,k,3,4,i_der] = DEB_T[i,j,k,4,i_der]
        DTF[i,j,k,4,3,i_der] = -DEB_T[i,j,k,4,i_der]
    end
    return DTF
end


#Compute the J matrix 
export make_J_matrix

function make_J_matrix(Eps_inv,TF2)
    nx,ny,nz = size(Eps_inv)[1:3]
    J = zeros((nx,ny,nz,4,4))
    for i = 1:nx, j = 1:ny, k = 1:nz
        J[i,j,k,:,:] = 4*Eps_inv[i,j,k]*Eps_inv[i,j,k]*linalg.I - TF2[i,j,k,:,:]
    end
    return J
end


#Invert the J matrix 
export invert_J_matrix

function invert_J_matrix(J,F_inv)
    nx,ny,nz = size(F_inv)[1:3]
    Jm1 = zeros((nx,ny,nz,4,4))
    for i = 1:nx, j = 1:ny, k = 1:nz
        Jloc =  J[i,j,k,:,:]
        if linalg.det(Jloc) !=0 
            Jm1[i,j,k,:,:] = linalg.inv(linalg.factorize(Jloc))
        end
    end
    return Jm1
end


#compute the f1 eigen 4-vector
export computing_eigenvector_f1

function computing_eigenvector_f1(Eps_inv,EB_T)
    nx,ny,nz = size(Eps_inv)[1:3]
    f1 = zeros((nx,ny,nz,4))
    f1[:,:,:,1] .= 1
    mask = Eps_inv .>0
    epsilon = Eps_inv[mask]
    Ex = EB_T[mask,1]
    Ex, Ey, Ez, Bx, By, Bz = EB_T[mask,1],EB_T[mask,2],EB_T[mask,3],EB_T[mask,4],EB_T[mask,5],EB_T[mask,6]
    alpha = Bz ./ epsilon
    beta = By ./ epsilon
    A = (alpha.*By .- Bx)./(epsilon.*(1 .+alpha.*alpha))
    f1[mask,4] .= (A.*(alpha.*Ex.-Ey) .-(beta.*Ex.+Ez)) ./ (A.*(Bx .+ alpha.*By) .- epsilon.*(1 .+beta.*beta))
    f1[mask,3] .= ((Bx .+ alpha.*By).*f1[mask,4] .+ Ey .- alpha.*Ex )./(epsilon.*(1 .+alpha.*alpha))
    f1[mask,2] .= (Bz.*f1[mask,3] - By.*f1[mask,4] .+ Ex)./epsilon
    return f1
end


#lower f1 4 vector index
export lower_f1_index

function lower_f1_index(f1)
    nx,ny,nz = size(f1)[1:3]
    f1_low = ones((nx,ny,nz,4))
    f1_low[:,:,:,2:4] = -f1[:,:,:,2:4]
    return f1_low
end


#Compute weff from previous quantity
export compute_weff

function compute_weff(F_inv,DTF,Jm1,f1_low,f1)
    nx,ny,nz = size(F_inv)[1:3]
    weff = zeros((nx,ny,nz))
    @einsum weff[i,j,k] = DTF[i,j,k,m5,m4,m6]*f1[i,j,k,m6] * f1_low[i,j,k,m5] * Jm1[i,j,k,m4,m1]* DTF[i,j,k,m1,m2,m3]* f1[i,j,k,m3] * f1[i,j,k,m2]
    weff = sqrt.(weff)
    return weff
end


#compute epsilon and omega eff from field array and resolution
export weff_epsilon_from_field

function weff_epsilon_from_field(EB_T,dx,dy,dz)
    F_inv, Eps_inv = computing_invariants(EB_T)
    TF2 = make_square_field_tensor(make_TF(EB_T))
    Jm1 = invert_J_matrix(make_J_matrix(Eps_inv,TF2),F_inv)
    DTF = make_derivative_tensor(computing_derivatives_EB(EB_T,dx,dy,dz))
    f1 = computing_eigenvector_f1(Eps_inv,EB_T)
    f1_low = lower_f1_index(f1)
    weff = compute_weff(F_inv,DTF,Jm1,f1_low,f1)
    eps2_weff = Eps_inv .* Eps_inv .* weff
    return weff, Eps_inv, eps2_weff
end


end