
#fitted parameters from datasets  (see Fig10 from Lorenz-Spreen)
# get parameters from captions, translate every parameter combination to (Tx,w)
# N=100 where does that fall on the N=300 phase diag?
using CSV,DataFrames

save_datatable = true

#Google Books
rs_books = [8.0,12.0,15.0,17.0,18.0,20.0]
N_books = 300 
c_books = 2.4
α_books = 0.003
dataset_names_books = ["Books_1870-1890","Books_1900-1920","Books_1930-1950","Books_1950-1970","Books_1970-1990","Books_1990-2010"]

#Movies
rs_mov = [18.0, 20.0, 22.0, 25.0, 30.0, 33.0, 35.0, 37.0]
N_mov = 100
c_mov = 2.4
α_mov = 0.001
dataset_names_mov = ["Movies_1980-1985","Movies_1985-1990","Movies_1990-1995","Movies_1995-2000","Movies_2000-2005","Movies_2005-2010","Movies_2010-2015","Movies_2015-2018"]

#Google Trends
rs_trends = [15.0, 17.0, 19.0, 21.0, 22.0, 26.0, 27.0, 30.0]
N_trends = 100
c_trends = 2.0
α_trends = 0.0005
dataset_names_trends = ["Google_2010","Google_2011","Google_2012","Google_2013","Google_2014","Google_2015","Google_2016","Google_2017"]

#Reddit 
rs_reddit = [25.0, 29.0, 30.0, 33.0, 34.0, 35.0]
N_reddit = 100
c_reddit = 5.4
α_reddit = 0.001
dataset_names_reddit = ["Reddit_2010","Reddit_2011","Reddit_2012","Reddit_2013","Reddit_2014","Reddit_2015"]

#publications 
rs_pub = [6.0,6.0,6.0,7.0,7.0]
N_pub = 300
c_pub = 2.4
α_pub = 0.003
dataset_names_pub = ["Publications_1990-2000","Publications_1995-2005","Publications_2000-2010","Publications_2005-2015","Publications_2010-2020"]

#Twitter
rs_twitter = [9.0,10.0,11.0,12.0]
N_twitter = 300 
c_twitter = 2.4
α_twitter = 0.005
dataset_names_twitter = ["Twitter_2013","Twitter_2014","Twitter_2015","Twitter_2016"]

K = 1.0
dataset_names = [dataset_names_books,dataset_names_mov,dataset_names_trends,dataset_names_reddit,dataset_names_pub,dataset_names_twitter]
Ns_vals = [N_books,N_mov,N_trends,N_reddit,N_pub,N_twitter]
r_vals = [rs_books,rs_mov,rs_trends,rs_reddit,rs_pub,rs_twitter]
c_vals = [c_books,c_mov,c_trends,c_reddit,c_pub,c_twitter]
α_vals = [α_books,α_mov,α_trends,α_reddit,α_pub,α_twitter]

#generate table with corresponding (y,w) values 
df = DataFrame(Dataset=String[],r=Float64[],c=Float64[],alpha=Float64[],N=Float64[],w=Float64[],y=Float64[])
for i in 1:length(Ns_vals)
    dataset_nameset = dataset_names[i]
    rs = r_vals[i]
    c = c_vals[i]
    α =  α_vals[i]
    N = Ns_vals[i]
    resc_params = rescaled_params.(rs, rs, K, c, α, N)
    ws = [resc_params[i][1] for i in 1:length(resc_params)]
    ys = [resc_params[i][2] for i in 1:length(resc_params)] ./ (Ns_vals[i] - 1)

    #add to dataframes
    for j in 1:length(ws)
        push!(df,(dataset_nameset[j],rs[j],c,α,N,ws[j],ys[j]))
    end

end

save_datatable && CSV.write("../../data/lyapcurve_datasets_parameters/datasets_parameters_precise.txt",df)




#=
#wiki 
rs_wiki = [45.0, 48.0, 50.0]
N_wiki = 100
resc_params_wiki = rescaled_params.(rs_wiki, rs_wiki, K, 2.4, 0.0001, N_wiki)
ws_wiki = [resc_params_wiki[i][1] for i in 1:length(resc_params_wiki)]
norm_tscale_wiki = [resc_params_wiki[i][2] for i in 1:length(resc_params_wiki)] ./ (N_wiki - 1)
=#