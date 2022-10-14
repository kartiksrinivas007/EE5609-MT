using JLD2
using BenchmarkTools

function eval_the_prog()

    @load "EvaluationSet_for_Prog_1.jld2"
    pass_vector = Vector{Bool}(undef, length(Samples))

    for indx in 1:length(Samples)

        A = Samples[indx].A
        b = Samples[indx].b

        U, r, consistent = rankconsistencyTeamID15(A,b); # call your function here
        pass_vector[indx] = (Samples[indx].U == U) && (Samples[indx].r == r) && (Samples[indx].consistent == consistent)
        println("Sample index $indx correctness: $(pass_vector[indx])")
    end

    # timetaken = @belapsed rankconsistencyTeamID15($Samples[25].A,$Samples[25].b) evals=1 samples=5 seconds=10

    return pass_vector# ,timetaken

end

include("EE5609TeamID15.jl") # include your file here
pass_vector = eval_the_prog();
println("\n------Score------")
display(sum(pass_vector))
# display(timetaken)
println("-----------------")

