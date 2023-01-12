for i in 1:N
    for b in 1:B
        println(db²[i,permute[i,b],b],"    $(db²[lp[i,o(b,B)],lp[i,b],b])     $i $b: $(permute[i,b]) ")
        if b == 1
            println("$i db²[$(lp[i,B]),$(lp[i,1]),$b]      $(db²[lp[i,B],lp[i,1],b])")
        end
    end
end