#!/shared/home/croux/.conda/envs/R_env/bin/Rscript
# #!/usr/bin/Rscript
for(i in commandArgs()){
        tmp = strsplit(i, '=')
        if(tmp[[1]][1] == 'timeStamp'){ timeStamp = tmp[[1]][2] }
        if(tmp[[1]][1] == 'sub_dir'){ sub_dir = tmp[[1]][2] }
        if(tmp[[1]][1] == 'nIterations_gof'){ nIterations_gof = as.integer(tmp[[1]][2]) }
        if(tmp[[1]][1] == 'writeDistribution'){ writeDistribution = as.logical(tmp[[1]][2]) }
}

### Summary Stats
# simulations
x = read.table(paste(timeStamp, '/', sub_dir, '/simulations.txt', sep=''), h=T)

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

if( writeDistribution==TRUE){

	### Summary Stats
	prior_ss = gof1_ss = gof2_ss = NULL
	prior_sfs = gof1_sfs = gof2_sfs = NULL
	# simulations
	for(i in 0:(nIterations_gof-1)){
		prior_ss = rbind(prior_ss, read.table(paste(timeStamp, '/best_model/best_model_', i, '/ABCstat.txt', sep=''), h=T))
		gof1_ss = rbind(gof1_ss, read.table(paste(timeStamp, '/gof/gof_', i, '/ABCstat.txt', sep=''), h=T))
		gof2_ss = rbind(gof2_ss, read.table(paste(timeStamp, '/gof_2/gof_', i, '/ABCstat.txt', sep=''), h=T))
		prior_sfs = rbind(prior_sfs, read.table(paste(timeStamp, '/best_model/best_model_', i, '/ABCjsfs.txt', sep=''), h=T))
		gof1_sfs = rbind(gof1_sfs, read.table(paste(timeStamp, '/gof/gof_', i, '/ABCjsfs.txt', sep=''), h=T))
		gof2_sfs = rbind(gof2_sfs, read.table(paste(timeStamp, '/gof_2/gof_', i, '/ABCjsfs.txt', sep=''), h=T))

	}

	# observation
	obs_ss = read.table(paste(timeStamp, '/ABCstat_global.txt', sep=''), h=T)
	obs_sfs = read.table(paste( timeStamp, '/ABCjsfs.txt', sep=''), h=T)

	# output for the web interface
	origin = c('observed dataset', rep('prior', nrow(prior_ss)), rep('posterior', nrow(gof1_ss)), rep('optimized posterior', nrow(gof2_ss)))
	PCA = rbind(cbind(obs_ss, obs_sfs), cbind(prior_ss, prior_sfs), cbind(gof1_ss, gof1_sfs), cbind(gof2_ss, gof2_sfs))
	PCA = cbind(PCA, origin)
	write.table(PCA, paste(timeStamp, '/distribution_PCA.txt', sep=''), col.names=T, row.names=F, quote=F, sep='\t')
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
outfile = paste(timeStamp, "/", sub_dir, "/goodness_of_fit_test.txt", sep='')
write.table(x=res, file=outfile, quote=FALSE, sep='\t', col.names=T, row.names=F)


### jSFS
# expected sfs
exp_sfs = read.table(paste(timeStamp, '/', sub_dir, '/simulations_jsfs.txt', sep=''), h=T)
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

outfile_sfs = paste(timeStamp, "/", sub_dir, "/gof_sfs.txt", sep='')
write.table(x=sfs, file=outfile_sfs, quote=FALSE, sep='\t', col.names=T, row.names=F)
 