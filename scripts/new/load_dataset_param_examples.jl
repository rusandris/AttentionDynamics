using AttentionDynamics

#system sizes 
N_books = 300
N_mov = 100
N_trends = 100
N_reddit = 100
N_pub = 300
N_wiki = 100
N_twitter = 300 

#constant parameter values for each dataset 
p_const_mov = [25,25,1.0,2.4,0.001,N_mov]
p_const_trends = [22,22,1.0,2.0,0.0005,N_trends]
p_const_reddit = [30,30,1.0, 5.4,0.001,N_reddit]
p_const_pub = [6,6,1.0, 2.4,0.003,N_pub]
p_const_wiki = [48,48,1.0, 2.4,0.0001,N_wiki]
p_const_books = [15,15,1.0,2.4,0.003,N_books]

#(w,g) values for each dataset 
w_mov,g_mov = rescaled_params(p_const_mov...)
w_tr,g_tr = rescaled_params(p_const_trends...)
w_red,g_red = rescaled_params(p_const_reddit...)
w_pub,g_pub = rescaled_params(p_const_pub...)
w_wiki,g_wiki = rescaled_params(p_const_wiki...)
w_books,g_books = rescaled_params(p_const_books...)

#yvalues for each dataset (g_rel_levels)
y_mov = g_mov / (N_mov-1)
y_tr = g_tr / (N_trends-1)
y_red = g_red / (N_reddit-1)
y_pub = g_pub / (N_pub-1)
y_wiki = g_wiki / (N_wiki-1)
y_books = g_books / (N_books-1)