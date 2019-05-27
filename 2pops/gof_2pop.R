#!/shared/home/croux/.conda/envs/R_env/bin/Rscript
# #!/usr/bin/Rscript
for(i in commandArgs()){
        tmp = strsplit(i, '=')
        if(tmp[[1]][1] == 'nameA'){ nameA = tmp[[1]][2] }
        if(tmp[[1]][1] == 'nameB'){ nameB = tmp[[1]][2] }
}

# simulations
x = read.table(paste('ABC_', nameA, '_', nameB, '/gof/simulations.txt', sep=''), h=T)

# observation
y = read.table(paste('ABC_', nameA, '_', nameB, '/ABCstat_global.txt', sep=''), h=T)

# outfile
outfile = paste("ABC_", nameA, "_", nameB, "/gof/goodness_of_fit_test.txt", sep='')


# function to compute the pval
pvalue = function(distribution, obs){
	median_x = median(distribution)
	if(obs==median_x){
		pval=0.5
	}else if(as.numeric(obs)>median_x){
		pval = length(which(distribution>as.numeric(obs)))/length(distribution)
	}else{
		pval = length(which(distribution<as.numeric(obs)))/length(distribution)
	}
	return(pval)
}

# measure the pval over all statistics
ss = c(2:31, 40:48)

stats = NULL
pvals = NULL
mean_exp = NULL
mean_obs = NULL
for(i in ss){
	stats = c(stats, colnames(x)[i])
	pvals = c(pvals, round(pvalue(x[,i], y[i]), 5))
	mean_exp = c(mean_exp, round(mean(x[,i]), 5))
	mean_obs = c(mean_obs, round(as.numeric(y[i]), 5))
}

pvals_fdr_corrected = round(p.adjust(pvals, "fdr"), 5)

res = data.frame(stats, mean_exp, mean_obs, pvals_fdr_corrected)

write.table(x=res, file=outfile, quote=FALSE, sep='\t', col.names=T, row.names=F)


