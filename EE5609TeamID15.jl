

function rankconsistencyTeamID15(A::Matrix{Int64}, b::Vector{Int64})
    consistent::Bool = false
    AugmentedMat = Augment(A, b)
    consistent,Rank,tag = UpperEchelonNew(AugmentedMat)
    Rank = convert(Int64, Rank)
    #Now we need to "de-Augment AugmentedMat to find U"
    U = AugmentedMat[begin:end, begin:end-1]
    # println("~~~~~~~~Row reduction process complete~~~~~")
    # display(U)
    return U,Rank,consistent
end 

function Augment(A, b)
    return [A b]
end

function AddGf(i, j)
Add_matrix = [  0 1 2 3  ;  1  0  3  2 ; 2 3 0 1 ; 3 2 1 0 ]
    return Add_matrix[i+1, j+1]
end

## Multiplication table
function MulGf(i, j)
    Mul_matrix = [ 0 0 0 0 ; 0  1  2  3 ; 0  2  3  1 ; 0 3 1 2 ]  
    return Mul_matrix[i+1, j+1]
end


Inv_vec = [1, 3, 2]
function InvGf(i)
    return Inv_vec[i]
end


"""
Precondition : Matrix has Non zero element in pivot_row and pivot_col. (depends on the OP of the function PartialPivotNew)
Args : Pivot column and Pivot Row ( Postion of the pivot and then it makes zero , beneath it) 
Returns: Matirx with Reduced Row j.
"""

function RowSubtractNew(A, pivot_row, pivot_col, j)
    # @printf("Subtracting Row %d from Row %d \n", pivot_row , j)
    factor = MulGf(InvGf(A[pivot_row , pivot_col]), A[j, pivot_col])
    # @printf("factor of subtraction = %d \n",factor)
    for index in pivot_col:size(A, 2)
        A[j , index] = AddGf(A[j , index], MulGf( A[pivot_row, index] ,factor))
    end
end

"""
This function creates a Pivot by finding the maximum element in the particular pivot column
and then performs an exchange of the row
Precondition: The previous positions are all zero
Args: Matrix and Position of the ***Supposed*** Next Pivot
Returns: Matrix and Boolean. 
Boolean = True if Pivot was found
        = False if No Pivot was found 
"""
function PartialPivotNew(A, sup_pivot_row, sup_pivot_col)
    maxIndex = sup_pivot_row 
    maxElement = A[sup_pivot_row, sup_pivot_col]
    isZeroCol = true 
    for i in sup_pivot_row:size(A, 1) 
        if(A[i, sup_pivot_col] != 0)
            isZeroCol = false
        end

        if(A[maxIndex, sup_pivot_col] < A[i, sup_pivot_col])
            maxIndex = i
            maxElement = A[i, sup_pivot_col]
        else
            continue
        end
    end
    #some optimizations could be done here ?
    if(sup_pivot_row != maxIndex)
    # @printf("Exchanging row %d woth row %d", sup_pivot_row, maxIndex)
    temp = A[sup_pivot_row,:]
    A[sup_pivot_row, :] = A[maxIndex, :]
    A[maxIndex , :] = temp
    end

    return A,isZeroCol

end

""" 
Precondition: This function assumes that this row has alread a pivot at pivot_row and pivot_col
Args : The pivot row and pivot column
Returns : The Matrix after Finding the Pivot on this Particular row
"""

function PivotSubtractNew(A, pivot_row, pivot_col)
    for i in pivot_row+1:size(A, 1)
        if(A[i, pivot_col] == 0)
            continue
        end
        RowSubtractNew(A, pivot_row, pivot_col, i)
    end
end

"""
This function creates a row reduced Echelon form of A  using RSNew, PSNew, PPNew.
"""
function UpperEchelonNew(A)
    tag = zeros(Int64, size(A, 2))
    col = 1
    for row in 1:size(A, 1) - 1
        isZeroCol = true
        while(isZeroCol && col <= size(A,2))
            A, isZeroCol = PartialPivotNew(A, row, col)
            if(isZeroCol)
                # @printf("The column %d is a Zero Column \n", col)
                col = col + 1
            end
        end

        if(!isZeroCol)
            # display(A)
            # @printf("(%d, %d) is a pivot position\n", row, col)
            tag[col] = 1
            pivot_row = row
            pivot_col = col
            PivotSubtractNew(A, pivot_row, pivot_col)
            #start next time with the next column, this column has been analyzed3 now
        end

        if(col > size(A, 2))
            #Now no more rows need to be checked since there are no spaces for pivots left( all zeroes beneath me anyway)
            # @printf("Ending the Process prematurely at row = %d\n",  row)
            break
        end

        if(row == size(A, 1) - 1)
            #then we need to check the position in which the last row has a pivot(if it has one )
            # println("Inside the check for the final Row Pivot")
            for column in col+1:size(A, 2)
                if(A[row+1, column] > 0)
                    tag[column] = tag[column] + 1
                    break #crucial
                end
            end
        end
        
        col = col + 1
    end

    consistent = false
    if(tag[size(A, 2)] == 0)
        consistent = true
    end

    Rank = sum(tag) - tag[size(A, 2)]

    # println("Rank = ",Rank)
    # println("Consistent = ",consistent)
    # println(tag)

    return consistent,Rank,tag

end

