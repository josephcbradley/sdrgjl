#=
Functions for counting relatinoships in spin chains. 
Trying to keep everything as Ints for speed. =#
export region_test, complement_A_count, two_intervals_count

function region_test(i, k, l)
    i >= k && i <= l
end

function complement_A_count(singlets, A)
    #=
    Given singlets and a sub-region A, count how many singlets shared connections 
    between A and B there are.

    This returns the integer count, which is quick and can be multiplied by ln(2) later.

    N.B this cannot be mutating as Ints are immutable.
    =#

    i, j = A #unpack region coordinates

    l = j - i  #calculate interval width

    #Iterate through the singlets - if both spins in the singlet are in A, 
    #add one to the count.

    count = 0
    # A singlet is 'between' A and B if either:
        # a is in A and b is not, or 
        # b is in A and a is no
    for singlet in singlets 
        a, b = singlet
        if region_test(a, i, j) && !region_test(b, i, j)
            count += 1
        elseif !region_test(a, i, j) && region_test(b, i, j)
            count += 1
        end
    end

    count
end

function two_intervals_count(singlets, i, l, r)
    #= 
    For two intervals of equal length l and separation r, 
    with the first interval having 1D coordinates i, j, 
    calculate the number of shared singlets.

    N.B this cannot be mutating as Ints are immutable. 
    =#

    # let j be the end coordinates of A 

    j = i + l - 1

    # Let x,y be the coordiantes of the second interval 
    x = i + l + r
    y = j + l + r

    #Iterate through the singlets - if one is in A₁ and the other in 
    # A₂, or vice versa, += 2

    count = 0 
    # A singlet is in A₁:A₂ if either:
        # a is in A₁ and b is in A₂, or 
        # b is in A₂ and a is in A₁
    for singlet in singlets 
        a, b = singlet
        if region_test(a, i, j) && region_test(b, x, y)
            count += 1
        elseif region_test(b, i, j) && region_test(a, x, y)
            count += 1
        end
    end

    count
end