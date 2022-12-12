using SymPy
using LinearAlgebra
using Latexify
using LaTeXStrings
using SymbolicUtils
using Kronecker

#sympy.interactive.printing.init_printing(use_unicode=False, wrap_line=False)


function sympDiag(coV)
    Nv = size(coV)
    Nv = Nv[1]
    N = Nv/2
    N = convert(Int64,N)
    
    Ω1 = [0 1; -1 0]
    idN = Matrix(1I,N,N)
    
    Ω = Ω1⊗idN
    B = 1im*V*Ω
    #latexify(B,render=true)
    P,D = B.diagonalize() #This does not return the eigenvalues in the order that we need
    latexify(D,render=true)
    latexify(P,render=true)
    
    # This will calculate the eigenvalues and eigenvectors, placing them in a data structure
    mat = B.eigenvects()
    #print(mat)

    # Extracting the symplectic eigenvectors from the data structure and normalizing
    P = zeros(1,length(mat))
    P = P'
    for i in collect(1:length(mat))
        if i == 1 # For the first step, we have to change the first row of the P matrix
            if mat[i][2] == 1 #Checking if there is degeneracy
                evec = mat[i][3]
                evec = evec[1]
                
                #Normalizing
                j = 0
                for ii in collect(1:length(evec))
                    j = j + abs(evec[ii])^2
                end
                evec = (1/sqrt(j))*evec
                P = evec
            else
                for jj in collect(1:length(mat[i][3]))
                    
                    if jj == 1
                        evec = mat[i][3]
                        evec = evec[jj]
                        
                        #Normalizing
                        j = 0
                        for ii in collect(1:length(evec))
                            j = j + abs(evec[ii])^2
                        end
                        evec = (1/sqrt(j))*evec
                        P = evec
                    else
                        evec = mat[i][3]
                        evec = evec[jj]

                        j = 0
                        for tt in collect(1:length(evec))
                            j = j + abs(evec[tt])^2
                        end
                        evec = (1/sqrt(j))*evec

                        P = [P evec]
                    end
                end
            end
        else # Now, performing for when we are not working with the first collumn
            if mat[i][2] == 1
                evec = mat[i][3]
                evec = evec[1]
                
                #Normalizing
                j = 0
                for tt in collect(1:length(evec))
                    j = j + abs(evec[tt])^2
                end
                evec = (1/sqrt(j))*evec

                P = [P evec]
            else
                for jj in collect(1:mat[i][2])
                    evec = mat[i][3]
                    evec = evec[jj]
                    #print(evec)
                    #print(jj)

                    #Normalizing
                    j = 0
                    for tt in collect(1:length(evec))
                        j = j + abs(evec[tt])^2
                    end
                    evec = (1/sqrt(j))*evec
                    #print(evec)
                    P = [P evec]
                    
                end
            end
        end
    end

    # Collecting the symplectic eigenvalues in a julia vector
    D = []
    for i in collect(1:length(mat))
        if mat[i][2] == 1
            eval = mat[i][1]
            append!(D, Sym[eval])
        else
            eval = mat[i][1]
            for ii in collect(1:mat[i][2])
                append!(D, Sym[eval])
            end
        end
    end
    
    #Now, getting the order that they need to be in for getting our desired answer
    di = []
    for i = collect(1:length(D))
        b = i+1
        for j = collect(b:length(D))
            if sympy.simplify(D[i]+D[j]) == 0
                if length(di) == 0
                    append!(di,[i j])
                    break
                else
                    if j != last(di)
                        append!(di,[i j])
                        break
                    end
                end
            end
        end
    end
    print(di)
    D1 = D
    P1 = P

    # Now, rearranging P and -D⊕D to be in the correct order
    if length(di) == 2
        mDD = sympy.diag(D[di[1]],D[di[2]]) 
        P = [P.col(di[1] - 1) P.col(di[2] - 1)]
    elseif length(di) == 4
        mDD = sympy.diag(D[di[1]],D[di[3]],D[di[2]],D[di[4]])
        P = [P.col(di[1] - 1) P.col(di[3] - 1) P.col(di[2] - 1) P.col(di[4] - 1)]
    end

    #print(D[di[2]])
    latexify(mDD,render=true)
    latexify(P,render=true)
    #latexify(P*mDD*inv(P),render=true)
    #latexify(transpose(P)*Ω*P,render=true)

    #=
    sz = size(D)
    DpD = zeros(sz[1],sz[1])
    sz2 = sz[1]/2
    sz2 = convert(Int64,sz2)
    for ii in collect(1:sz2)
        DpD[ii,ii] = D[ii+sz2,ii+sz2]
        DpD[ii+sz2,ii+sz2] = D[ii+sz2,ii+sz2]
    end

    DpD
    =#

    # Now, getting D⊕D
    if length(di) == 2
        DD = sympy.diag(D1[di[2]],D1[di[2]])
    elseif length(di) == 4
        DD = sympy.diag(D1[di[2]],D1[di[4]],D1[di[2]],D1[di[4]])
    end
    
    # Can be used to check the result
    #latexify(DD,render=true)
    #latexify(P,render=true)

    # Now calculating our S matrix
    U2 = (1/sqrt(2))*[1 1; 1im -1im ]
    U2T = (1/sqrt(2))*[1 -1im; 1 1im ] #Shouldnt this be the transpose conjugate? 
    S = P*(U2T⊗idN)
    #We can verify our S matrix is actually symplectic
    # S*Ω*transpose(S)

    # Finally, calculating our covariance matrix to confirm our results
    Vf = S*DD*transpose(S)
    Vf = sympy.simplify(Vf) #Using sympy to simplify our expression
    
    # I have found that this exports it to the terminal in the nicest form, even though the latex makes no sense
    sympy.latex(Vf)
    latexify(Vf,render=true)


end

# We begin by defining our variables
@vars λ γ

# Covariance matrix for a displaced vacuum state, in the end we just need to replace λ with 1
#V = (1/2)*[λ 0;0 λ]

# Covariance matrix for a Single mode thermal state
#V = (1/2)*[(1 + λ^2)/(1 - λ^2) 0 0 0; 0 (1 + λ^2)/(1 - λ^2) 0 0; 0 0 (1 + λ^2)/(1 - λ^2) 0; 0 0 0 (1 + λ^2)/(1 - λ^2)]

# Covariance matrix for a 2 mode thermal state
#V = (1/2)*[(1 + λ^2)/(1 - λ^2) 0 0 0; 0 (1 + γ^2)/(1 - γ^2) 0 0; 0 0 (1 + λ^2)/(1 - λ^2) 0; 0 0 0 (1 + γ^2)/(1 - γ^2)]

# Now, using our function to perform the symplectic diagonalization
sympDiag(V)
#Latexify(mat,render=true)
#latexify(B1,render=true)


#directsum(A,B) = [A zeros(size(A,1), size(B,2)); zeros(size(B,1), size(A,2)) B]
#directsum(D,D)

#latexify(sympDiag(V),render=true)



#D = diagm(val)
#S = P*(U2T⊗idN)
#S*DirectSum.⊕(D,D)*transpose(S)
