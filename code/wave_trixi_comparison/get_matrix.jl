A,b = linear_structure(semi)
A = Matrix(A)
B = zeros(size(A))
C = zeros(size(A))
for i in range(1,size(A)[1])
    if i % 2 == 0
        shift1 = size(A)[1]/2
        B[Int(shift1 + i/2),:] = A[i,:]
    else
        B[Int((i+1)/2),:] = A[i,:]
    end
end
for i in range(1,size(A)[1])
    if i % 2 == 0
        shift1 = size(A)[1]/2
        C[:,Int(shift1 + i/2)] = B[:,i]
    else
        C[:,Int((i+1)/2)] = B[:,i]
    end
end
show(stdout, "text/plain",C)