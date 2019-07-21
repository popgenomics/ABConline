#!/shared/home/croux/.conda/envs/R_env/bin/Rscript
# #!/usr/bin/Rscript
for(i in commandArgs()){
        tmp = strsplit(i, '=')
        if(tmp[[1]][1] == 'timeStamp'){ timeStamp = tmp[[1]][2] }
}

### Summary Stats
# simulations
x = read.table(paste(timeStamp, '/gof/simulations.txt', sep=''), h=T)

# observation
y = read.table(paste(timeStamp, '/ABCstat_global.txt', sep=''), h=T)


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

# outfile
outfile = paste(timeStamp, "/gof/goodness_of_fit_test.txt", sep='')
write.table(x=res, file=outfile, quote=FALSE, sep='\t', col.names=T, row.names=F)


### jSFS
# expected sfs
exp_sfs = read.table(paste(timeStamp, '/gof/simulations_jsfs.txt', sep=''), h=T)
exp_sfs_2 = apply(exp_sfs, MARGIN=2, FUN="median")

# observed sfs
obs_sfs = read.table(paste(timeStamp, '/ABCjsfs.txt', sep=''), h=T)

# compute the pvalue
tested_sfs = NULL
for(i in 1:length(exp_sfs)){
	tested_sfs = c(tested_sfs, pvalue(exp_sfs[,i], obs_sfs[i]))
}

tested_sfs = round(p.adjust(tested_sfs, "fdr"), 5)

# all matrixes 
## obs | exp | exp-obs | pval
sfs = rbind(obs_sfs, exp_sfs_2, exp_sfs_2-obs_sfs, tested_sfs)

outfile_sfs = paste(timeStamp, "/gof/gof_sfs.txt", sep='')
write.table(x=sfs, file=outfile_sfs, quote=FALSE, sep='\t', col.names=T, row.names=F)
 
