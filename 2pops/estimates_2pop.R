#!/shared/home/croux/.conda/envs/R_env/bin/Rscript
# #!/usr/bin/Rscript
for(i in commandArgs()){
	tmp = strsplit(i, '=')
	if(tmp[[1]][1] == 'nameA'){ nameA = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nameB'){ nameB = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nMin'){ nMin = as.integer(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'sub_dir_sim'){ sub_dir_sim = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nSubdir'){ nSubdir = as.integer(tmp[[1]][2]) } # number of subdirectories where simulations were ran
	if(tmp[[1]][1] == 'ncores'){ ncores = as.integer(tmp[[1]][2]) } # number of cores for the random forest
	if(tmp[[1]][1] == 'ntree'){ ntree = as.integer(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'outgroup'){ outgroup = as.integer(tmp[[1]][2]) } # 0: no outgroup, no SFS used. 1: outgroup, SFS used
}

outfile = paste('ABC_', nameA, '_', nameB, '/', sub_dir_sim, '/report_', nameA, '_', nameB, '.txt', sep='')

# colors
coul = c('#ffffcc', '#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#0c2c84')
coul = colorRampPalette(coul)

# observed data
obs_ss = read.table(paste('ABC_', nameA, '_', nameB, '/ABCstat_global.txt', sep=''), h=T)
obs_ss = obs_ss[, -grep('min', colnames(obs_ss))]
obs_ss = obs_ss[, -grep('max', colnames(obs_ss))]

######################################################
# parameters of the best model, IM_2M_2N and SI_2N
source('/shared/home/croux/softwares/ABConline/2pops/get_parameters.R')
#source("/home/croux/Documents/ABConline/2pops/get_parameters.R")

# IM_2M_2N
model_tmp = 'IM_2M_2N'
write(paste('\n#####\n\nparameters of model: ', model_tmp, sep=''), outfile, append=T)
posterior = get_posterior(nameA=nameA, nameB=nameB, nSubdir=nSubdir, sub_dir_sim=sub_dir_sim, model=model_tmp)
write('param\tHPD2.5%\tmedian\tHPD%97.5', outfile, append=T)
for(i in 1:ncol(posterior)){
	write(paste(colnames(posterior)[i], as.numeric(quantile(posterior[,i], 0.025)), as.numeric(quantile(posterior[,i], 0.5)), as.numeric(quantile(posterior[,i], 0.975)), sep='\t'), outfile, append=T)
}



# SI_2N
model_tmp = 'SI_2N'
write(paste('\n#####\n\nparameters of model: ', model_tmp, sep=''), outfile, append=T)
posterior = get_posterior(nameA=nameA, nameB=nameB, nSubdir=nSubdir, sub_dir_sim=sub_dir_sim, model=model_tmp)
write('param\tHPD2.5%\tmedian\tHPD%97.5', outfile, append=T)
for(i in 1:ncol(posterior)){
	write(paste(colnames(posterior)[i], as.numeric(quantile(posterior[,i], 0.025)), as.numeric(quantile(posterior[,i], 0.5)), as.numeric(quantile(posterior[,i], 0.975)), sep='\t'), outfile, append=T)
}


# AM_2M_2N
model_tmp = 'AM_2M_2N'
write(paste('\n#####\n\nparameters of model: ', model_tmp, sep=''), outfile, append=T)
posterior = get_posterior(nameA=nameA, nameB=nameB, nSubdir=nSubdir, sub_dir_sim=sub_dir_sim, model=model_tmp)
write('param\tHPD2.5%\tmedian\tHPD%97.5', outfile, append=T)
for(i in 1:ncol(posterior)){
	write(paste(colnames(posterior)[i], as.numeric(quantile(posterior[,i], 0.025)), as.numeric(quantile(posterior[,i], 0.5)), as.numeric(quantile(posterior[,i], 0.975)), sep='\t'), outfile, append=T)
}


# SC_2M_2N
model_tmp = 'SC_2M_2N'
write(paste('\n#####\n\nparameters of model: ', model_tmp, sep=''), outfile, append=T)
posterior = get_posterior(nameA=nameA, nameB=nameB, nSubdir=nSubdir, sub_dir_sim=sub_dir_sim, model=model_tmp)
write('param\tHPD2.5%\tmedian\tHPD%97.5', outfile, append=T)
for(i in 1:ncol(posterior)){
	write(paste(colnames(posterior)[i], as.numeric(quantile(posterior[,i], 0.025)), as.numeric(quantile(posterior[,i], 0.5)), as.numeric(quantile(posterior[,i], 0.975)), sep='\t'), outfile, append=T)
}

