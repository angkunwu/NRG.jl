#using LinearAlgebra
#using SparseArrays

"""
eye(n::Int)

Create identity matrix (sparse)
"""
function eye(n::Int)
    mat = zeros(n,n)
    for k = 1:n
        mat[k, k] = 1.0
    end
    return sparse(mat) # return a sparse matrix
end
"""
replace non-zero sparse matrix element into ones
"""
function spones(A::SparseMatrixCSC)
    rows, cols, vals = findnz(A)
    m = size(vals, 1)
    vals = ones(m)
    return sparse(rows,cols,vals, A.m, A.n)
end
function spones(A::SparseVector)
    rows, vals = findnz(A)
    m = size(vals, 1)
    vals = ones(m)
    return sparsevec(rows,vals, A.n)
end

"""
block diagalize a matrix, returns the reorganized matrix and the
permutation and a vector of block sizes for the reorganized matrix
"""
function blockdiagalization(A::SparseMatrixCSC)
    # Breadth first search using sparse matrix
    m = size(A, 1)
    G = spones(A) # sparse(A);
    flag = zeros(Int64, m) # Record block size of this site
    permVec = zeros(Int64,m) # permutation order
    # Perform breadth-first search. nnz nonzeros
    l = 1
    for k = 1:m
        if flag[k] == 0
            flag[k] = 1 # mark visited node
            nzrow = findnz(G[k,:])[1]
            nzcol = findnz(G[:,k])[1]
            block = union(nzrow, nzcol) # find non zero indexes from row and col index
            block = union(block, [k])
            vec = sparsevec(block, ones(size(block,1)), m) # cols, vals, size of vector; bfs using sparse matrix multiply
            vec = G * vec # first layer search
            vec = spones(vec) # make non-zero elements ones
            blocknext = findnz(vec)[1]
            if ~isempty(blocknext) # exclude empty blocks
                # whether the data in blocknext is found in block.%~isequal(blocknext,block)
                blocknext = union(blocknext, block)
                while blocknext != block  
                    block = blocknext
                    vec = sparsevec(block, ones(size(block,1)), m)
                    vec = G * vec
                    vec = spones(vec)
                    blocknext = findnz(vec)[1]
                    blocknext = union(blocknext, block)
                end
            end
            M = size(block)[1] # size of this block
            permVec[l:(l+M-1)] = block
            l=l+M
            for j=1:M
                flag[block[j]] = M # mark visited nodes
            end
        end      
    end
    BlockMat = A[permVec,permVec]
    return BlockMat, permVec, flag[permVec]
end

"""
diagonalize the matrix organized in block
- 'H0': is already block diagonalized
- 'Bsizes': record the block size for each site it belongs to
diagonalization in blocks should not change the Q,Sz quantum numbers
# Outputs
eigenvalues (Vector), eigenvectors(sparse matrix)
"""
function DiagInBlocks(H0, Bsizes)
    Ns = size(H0, 1)
    Energy = zeros(Ns)
    Vecs = spzeros(Ns,Ns)
    klow = 1 
    # identify the diagonal blocks
    while klow <= Ns
        khigh = klow + Bsizes[klow] - 1
        Es, Vs = eigen(Matrix(H0[klow:khigh,klow:khigh])); 
        Energy[klow:khigh] = Es
        Vecs[klow:khigh,klow:khigh] = Vs
        klow = khigh+1
    end
    return Energy, Vecs
end


"""
[GS,GSQ,GSSz,Mfimp1,Szimp]=NRGtrunBlock(N,Lam,epsilon,t,U,ef,V,Ns;h=0,tol)
# Arguments
- 'N': the number of NRG iteration (Wilson chain length)
- 'Lam': the discretization coefficient
- 'epsilon': Wilson chain onsite electron energies
- 't': Wilson chain hopping parameters
- 'U': Coulomb repulsion at the impurity site
- 'V': hybridization strength
- 'Ns': NRG spectrum truncation limit
## Implicit arguments
- 'h': local magnetic filed at the impurity site
- 'tol': same energy level limit

# Outputs
- 'Gs': Many-body energy levels
- 'GSQ': Total charge Q for each many-body level
- 'GSSz': Total S^z for each many-body level
- 'Mfimp1': diagonal part of impurity annilation operator
- 'Szimp': diagonal part of Sz matrix at impurity site for chi_loc

tol 10^(-8) for no magnetic field, otherwise 10^(-12)
"""
function runNRG(N::Int64,Lam::Number,epsilon::Vector{<:Number},t::Vector{<:Number},U::Number,ef::Number,V::Number,Ns::Int64; h=0.0,tol=10^(-12))
    GS = Vector{Float64}[]
    GSQ = Vector{Float64}[] # record Q at each step for each state
    GSSz = Vector{Float64}[] # record Sz at each step for each state

    #QSztemp = zeros(Ns*4, 2) # temp of QSz for higher dimension
    #Mfimp1=zeros(Ns,Ns,N)  # keep track of impurity annilation operator
    #Mfimp1 = zeros(Ns, N)  # keep track of diagonal part of impurity annilation operator
    #Mfimp2=zeros(Ns,Ns,N) # 1 for spin up; 2 for spin down
    #Szimp=zeros(Ns,Ns,N) # try to obtain Sz impurity matrice for \chi_loc
    #Szimp = zeros(Ns, N) # obtain diagonal part of Sz impurity matrice for \chi_loc
    Sz = zeros(4, 4) 
    Sz[2, 2], Sz[3, 3] = 1/2, 1/2

    # Iterative diagonalization
    H0 = zeros(16, 16) # initial Hamiltonian, the basis chosen as 1 (no electron)
    # 2 (one up), 3 (one down), 4 (double occupancy) with 4*4 impurity sites *
    # first conduction degree of freedom, -1 for anticommutation relation
    cup = zeros(4, 4);
    cup[1, 2], cup[3, 4] = 1, 1 # annihilation operator for spin up
    cdn=zeros(4,4)
    cdn[1, 3], cdn[2, 4] = 1, 1 #f2(2,4)=-1 # annihilation operator for spin down
    IM=zeros(4,4) # Matrix for antisymmetry
    IM[1, 1], IM[4, 4], IM[3, 3], IM[2, 2] = 1, 1, 1, 1 # IM[2, 2], IM[3, 3] = -1, -1
    IM = sparse(IM)
    Himp = ef*(cup'*cup+cdn'*cdn)+U*(cup'*cup*cdn'*cdn)+h*Sz # impurity site hamiltonian
    QSz = zeros(4,2)
    QSz[:, 1] = [-1, 0, 0, 1]
    QSz[:, 2] = [0, 0.5, -0.5, 0]

    function NRGaddSite(Es, QSz, cup_old, cdn_old, hop, onsite, iteration; tol=1E-12)
        m = size(Es, 1)
        Hold_HD = kron(sparse(diagm(Es)), eye(4)) # old Hamiltonian from last energies
        Honsite = onsite * (cup'*cup + cdn'*cdn)
        Honsite = kron(eye(m), Honsite)

        cup_old_HD = kron(cup_old, IM)
        cdn_old_HD = kron(cdn_old, IM)
        HD_cup = kron(eye(m), cup)
        HD_cdn = kron(eye(m), cdn)
        # ' include complex conjugate
        Hhop = hop *(cup_old_HD'*HD_cup+HD_cup'*cup_old_HD+cdn_old_HD'*HD_cdn+HD_cdn'*cdn_old_HD)
        if iteration > 0
            Hnew = sqrt(Lam)*Hold_HD + Lam^((iteration-1)/2)*(Honsite + Hhop)
        else
            Hnew = Hold_HD + Honsite + Hhop
        end
        # update charge Q and Sz index for current Hilbert space
        NewQSz=zeros(4*m,2)
        for k=1:m
            NewQSz[4*k-3,1]=QSz[k, 1]-1
            NewQSz[4*k-2,1]=QSz[k, 1]
            NewQSz[4*k-1,1]=QSz[k, 1]
            NewQSz[4*k,1]=QSz[k, 1]+1
            NewQSz[4*k-3,2]=QSz[k, 2]
            NewQSz[4*k-2,2]=QSz[k, 2]+0.5
            NewQSz[4*k-1,2]=QSz[k, 2]-0.5
            NewQSz[4*k,2]=QSz[k, 2]
        end
        Hnew,permu, Bsizes = blockdiagalization(Hnew) # block diagalize H0, permutation sequence is ind
        NewQSz = NewQSz[permu, :]
        #display(Bsizes[1:16])
        #display(NewQSz[1:16,:])
        Egs, Vecs = DiagInBlocks(Hnew, Bsizes)
        ind = sortperm(Egs)
        NewEs = Egs[ind]
        NewEs = NewEs.- NewEs[1] # off set GS to zero
        NewVesc = Vecs[:, ind]
        NewQSz = NewQSz[ind,:] # block diagonalization should not change good quantum number

        if 4*m > Ns # reach space limit
            Nsdyn = Ns # dynamically take more states around space limit Ns
            while Nsdyn < 4*m && NewEs[Nsdyn+1] - NewEs[Nsdyn] < tol
                Nsdyn += 1
            end
            # truncate space
            NewEs = NewEs[1:Nsdyn]
            NewVesc = NewVesc[:,1:Nsdyn]
        end

        mnew = size(NewEs, 1)
        # set boundary for zero to avoid numerical divergence
        for k = 1:mnew-1
            if NewEs[k+1] - NewEs[k] < tol
                NewEs[k+1] = NewEs[k]
            end
        end

        return NewEs, NewQSz[1:mnew,:], permu, NewVesc
    end

    # tranform the last operator
    function transformC(cop, permu, Vs)
        mnew = size(Vs, 1)
        cop_old = kron(eye(div(mnew,4)), cop);
        cop_old = cop_old[permu, permu]
        cop_old = Vs' * cop_old * Vs
        return cop_old
    end

    Es = diag(Himp)
    cup_old, cdn_old= cup, cdn
    for iteration = 0:N-1
        if iteration > 0
            hop = t[iteration]
        else
            hop = V
        end
        @time NewEs, NewQSz, permu, NewVesc = NRGaddSite(Es, QSz, cup_old, cdn_old, hop, epsilon[iteration+1], iteration)
	
	if iteration > 0 # record results only from the NRG site 2
        	push!(GS, NewEs)
        	push!(GSQ, NewQSz[:,1])
        	push!(GSSz, NewQSz[:,2])
	end
        cup_old = transformC(cup, permu, NewVesc)
        cdn_old = transformC(cdn, permu, NewVesc)

        Es = NewEs
        QSz = NewQSz
        @show iteration
    end

    return GS, GSQ, GSSz
end
