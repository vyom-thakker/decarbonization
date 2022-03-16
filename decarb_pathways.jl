using JuMP
Fusing CSV
using LinearAlgebra
using DataFrames
using BARON


A_db = CSV.File("./A_mat.csv",header=1,delim=",") |> DataFrame
B_db = CSV.File("./B_mat.csv",header=1,delim=",") |> DataFrame
f_db = CSV.File("./f_val.csv",header=1,delim=",") |> DataFrame
f_val=f_db[!,2]
A_mat=Matrix(A_db[:,2:size(A_db,2)])
B_mat=Matrix(B_db[:,2:size(B_db,2)])
nodes_val=names(A_db)[2:45]
product_val=A_db[!,1]


function opt_decarb()
    dcrb =  Model(BARON.Optimizer)
    @variable(dcrb,s[j=1:size(A_mat,2)] ≥ 0)
    @variable(dcrb,f[i=1:size(A_mat,1)] ≥ 0 )
    @variable(dcrb,g[k=1:size(B_mat,1)])
    @constraint(dcrb,[i=1:size(A_mat,1)],f[i] ≥ f_val[i]/2)
    @constraint(dcrb,[i=1:size(A_mat,1)], f[i]==sum(A_mat[i,j]*s[j] for j in 1:size(A_mat,2)))
    @constraint(dcrb,[k=1:size(B_mat,1)],g[k]==sum(B_mat[k,j]*s[j] for j in 1:size(B_mat,2)))
    @variable(dcrb,obj)
    @constraint(dcrb,obj==sum(g[k] for k in 3:4)+g[8])
    @objective(dcrb,Min,obj)
    optimize!(dcrb)
    return JuMP.value.(s),JuMP.value.(f),JuMP.value.(obj)
end


s,f,obj=opt_decarb()

res_f=DataFrame(pdt=product_val, final_demand=f)
res_s=DataFrame(node=nodes_val, scaling_factor=s)


CSV.write("./design_scaling.csv",res_s)
