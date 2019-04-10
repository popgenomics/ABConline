#!/usr/bin/Rscript
library('abcrf')
# model_comp_2pop.R nameA=txn nameB=ama nreps=20 Nref=100000 ntree=1000 ncores=6
for(i in commandArgs()){
	tmp = strsplit(i, '=')
	if(tmp[[1]][1] == 'nameA'){ nameA = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nameB'){ nameB = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nreps'){ nreps = as.integer(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'Nref'){ Nref = as.integer(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'ntree'){ ntree = as.integer(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'ncores'){ ncores = as.integer(tmp[[1]][2]) }
}
#nameA = 'txn'
#nameB = 'ama'
#nreps = 20 # number of sub-directory where simulations were run for each model
#Nref = 100000
#ntree = 1000
#ncores = 6
nsims_monolocus = 10000 # number of monolocus simulations

outfile = paste('ABC_', nameA, '_', nameB, '/report_', nameA, '_', nameB, '.txt', sep='')

# observed data
obs_ss = read.table(paste('ABC_', nameA, '_', nameB, '/ABCstat_global.txt', sep=''), h=T)
obs_ss = obs_ss[, -grep('min', colnames(obs_ss))]
obs_ss = obs_ss[, -grep('max', colnames(obs_ss))]
obs_sfs = read.table(paste('ABC_', nameA, '_', nameB, '/ABCjsfs.txt', sep=''), h=T)

ss_obs = cbind(obs_ss, obs_sfs)

# simulated data
models = c('SC_1M_1N', 'SC_1M_2N', 'SC_2M_1N', 'SC_2M_2N', 'AM_1M_1N', 'AM_1M_2N', 'AM_2M_1N', 'AM_2M_2N', 'IM_1M_1N', 'IM_1M_2N', 'IM_2M_1N', 'IM_2M_2N', 'SI_1N', 'SI_2N')
migration = c('migration', 'migration', 'migration', 'migration', 'isolation', 'isolation', 'isolation', 'isolation', 'migration', 'migration', 'migration', 'migration', 'isolation', 'isolation')
ss_sim = list()
params_sim = list()
all_models_sim = NULL

for(m in models){
	ss_sim_tmp = NULL
	params_sim_tmp = NULL
	for(rep in seq(0, nreps-1, 1)){
		# statistics
		tmp_ss = read.table(paste('ABC_', nameA, '_', nameB, '/', m, '_', rep, '_beta/ABCstat.txt', sep=''), h=T)
		tmp_ss = tmp_ss[, -grep('min', colnames(tmp_ss))]
		tmp_ss = tmp_ss[, -grep('max', colnames(tmp_ss))]
		tmp_sfs = read.table(paste('ABC_', nameA, '_', nameB, '/', m, '_', rep, '_beta/ABCjsfs.txt', sep=''), h=T)
		tmp = cbind(tmp_ss, tmp_sfs)
		ss_sim_tmp = rbind(ss_sim_tmp, tmp)
		
		# params
		tmp_params = read.table(paste('ABC_', nameA, '_', nameB, '/', m, '_', rep, '_beta/priorfile.txt', sep=''), h=T)
		params_sim_tmp = rbind(params_sim_tmp, tmp_params)
	}
	# statistics
	ss_sim[[m]] = ss_sim_tmp 
	all_models_sim = rbind(all_models_sim, ss_sim_tmp)
	
	# params
	params_sim[[m]] = params_sim_tmp
}


# remove uninformative statistics: those with no variation
ss_2_remove = c(1)
for(m in models){
	for(i in 2:ncol(ss_obs)){
		if( sd(ss_sim[[m]][, i])==0 ){
			ss_2_remove = c(ss_2_remove, i)
		}
	}
}
ss_2_remove = unique(ss_2_remove)

# model comparison #1 --> all models
modIndexes = NULL
for(m in models){
	modIndexes = c(modIndexes, rep(m, nrow(ss_sim[[m]])))
}

mod = abcrf(modIndexes~., data = data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
predicted_model = predict(mod, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)

write('MODEL COMPARISON #1: 14 models', outfile, append=F)
write('#confusion matrix:', outfile, append=T)
write.table(mod$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
write(paste('\n#best model among 14 models: ', predicted_model$allocation, sep=''), outfile, append=T)
write(paste('#proba best model among 14 models: ', predicted_model$post.prob, sep=''), outfile, append=T)
write('\n#votes:', outfile, append=T)
write.table(t(as.matrix(predicted_model$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)

# model comparison #2 --> two models: isolation versus migration
modIndexes = NULL
for(i in 1:length(models)){
	modIndexes = c(modIndexes, rep(migration[i], nrow(ss_sim[[models[i]]])))
}

mod_iso_mig = abcrf(modIndexes~., data = data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
predicted_model_iso_mig = predict(mod_iso_mig, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)

write('\n#####\n\nMODEL COMPARISON #2: 2 models', outfile, append=T)
write('#confusion matrix:', outfile, append=T)
write.table(mod_iso_mig$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
write(paste('\n#best model between migration and isolation: ', predicted_model_iso_mig$allocation, sep=''), outfile, append=T)
write(paste('#proba best model between migration and isolation: ', predicted_model_iso_mig$post.prob, sep=''), outfile, append=T)
write('\n#votes:', outfile, append=T)
write.table(t(as.matrix(predicted_model_iso_mig$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)


######################################################

# parameters of the best model, IM_2M_2N and SI_2N
best_model = predicted_model$allocation
write(paste('\n#####\n\nParameters of the best model which is: ', best_model, sep=''), outfile, append=T)
write(paste('To convert in demographic units (in generations), times have to be multiplied by ', 4*Nref, sep=''), outfile, append=T)
write(paste('To convert in demographic units (in number of individuals), sizes (N) have to be multiplied by ', Nref, '\n', sep=''), outfile, append=T)
test_IM = 0
test_SI = 0
if(best_model != 'IM_2M_2N'){
	test_IM = 1
	parameters_IM = NULL
}

if(best_model != 'SI_2N'){
	test_SI = 1
	parameters_SI = NULL
}

parameters_best = NULL
for(rep in seq(0, nreps-1, 1)){
		parameters_best = rbind(parameters_best, read.table(paste('ABC_', nameA, '_', nameB, '/', best_model, '_', rep, '_beta/priorfile.txt', sep=''), h=T))
		if( test_IM == 1){
			parameters_IM = rbind(parameters_IM, read.table(paste('ABC_', nameA, '_', nameB, '/IM_2M_2N_', rep, '_beta/priorfile.txt', sep=''), h=T))
		}
		
		if( test_SI == 1){
			parameters_SI = rbind(parameters_SI, read.table(paste('ABC_', nameA, '_', nameB, '/SI_2N_', rep, '_beta/priorfile.txt', sep=''), h=T))
		}
}

# best model
param_id = NULL
expectation = NULL
median = NULL
variance = NULL
variance_cdf = NULL
quantile_0025 = NULL
quantile_0975 = NULL

prior_min = NULL
prior_median = NULL
prior_max = NULL
prior_variance = NULL
prior_quantile_0025 = NULL
prior_quantile_0975 = NULL

pdf(paste('ABC_', nameA, '_', nameB, '/parameters_', best_model, '_', nameA, '_', nameB, '.pdf', sep=''), bg='white')
for(i in 1:ncol(parameters_best)){
	parameter_id = colnames(parameters_best)[i]
	param = parameters_best[,i]
	data = data.frame(param, ss_sim[[best_model]][, -ss_2_remove])
	model = regAbcrf(param~., data, ntree = ntree, ncores=ncores)
	predict_param = predict(model, ss_obs[, -ss_2_remove], data)
	
	param_id = c(param_id, parameter_id) 
	expectation = c(expectation, predict_param$expectation)
	median = c(median, predict_param$med)
	variance = c(variance, predict_param$variance)
	variance_cdf = c(variance_cdf, predict_param$variance.cdf)
	quantile_0025 = c(quantile_0025, predict_param$quantiles[1])
	quantile_0975 = c(quantile_0975, predict_param$quantiles[2])

	prior_min = c(prior_min, min(param))
	prior_median = c(prior_median, median(param))
	prior_max = c(prior_max, max(param))
	prior_variance = c(prior_variance, var(param))
	prior_quantile_0025 = c(prior_quantile_0025, quantile(param, 0.025))
	prior_quantile_0975 = c(prior_quantile_0975, quantile(param, 0.975))
	plot(model$model.rf$predictions, param, main=colnames(parameters_best)[i], xlab='estimated value', ylab='real value', cex.axis=1.2, cex.lab=1.2, pch=16, cex=1.2, col=rgb(0,0,0,0.25), xlim=range(param)); abline(v=predict_param$med, col='red', lwd=2)
}
dev.off()

expectation = round(expectation,5)
median = round(median,5)
variance = round(variance, 5)
variance_cdf = round(variance_cdf, 5)
quantile_0025 = round(quantile_0025, 5)
quantile_0975 = round(quantile_0975, 5)
prior_min = round(prior_min, 5)
prior_median = round(prior_median, 5)
prior_max = round(prior_max, 5)
prior_variance = round(prior_variance, 5)
prior_quantile_0025 = round(prior_quantile_0025, 5)
prior_quantile_0975 = round(prior_quantile_0975, 5)

res_best = data.frame(param_id, expectation, median, prior_median, variance, variance_cdf, prior_variance, quantile_0025, prior_quantile_0025, quantile_0975, prior_quantile_0975, prior_min, prior_max)

write.table(res_best, outfile, append=T, col.names=T, row.names=F, sep='\t', quote=F)


if( test_IM==1 ){
	write('\n#####\n\nParameters of the IM_2M_2N model', outfile, append=T)
	write(paste('To convert in demographic units (in generations), times have to be multiplied by ', 4*Nref, sep=''), outfile, append=T)
	write(paste('To convert in demographic units (in number of individuals), sizes (N) have to be multiplied by ', Nref, '\n', sep=''), outfile, append=T)
	
	# IM model
	param_id = NULL
	expectation = NULL
	median = NULL
	variance = NULL
	variance_cdf = NULL
	quantile_0025 = NULL
	quantile_0975 = NULL

	prior_min = NULL
	prior_median = NULL
	prior_max = NULL
	prior_variance = NULL
	prior_quantile_0025 = NULL
	prior_quantile_0975 = NULL

	pdf(paste('ABC_', nameA, '_', nameB, '/parameters_IM_2M_2N_', nameA, '_', nameB, '.pdf', sep=''), bg='white')
	for(i in 1:ncol(parameters_IM)){
		parameter_id = colnames(parameters_IM)[i]
		param = parameters_IM[,i]
		data = data.frame(param, ss_sim[['IM_2M_2N']][, -ss_2_remove])
		model = regAbcrf(param~., data, ntree = ntree, ncores=ncores)
		predict_param = predict(model, ss_obs[, -ss_2_remove], data)
		
		param_id = c(param_id, parameter_id) 
		expectation = c(expectation, predict_param$expectation)
		median = c(median, predict_param$med)
		variance = c(variance, predict_param$variance)
		variance_cdf = c(variance_cdf, predict_param$variance.cdf)
		quantile_0025 = c(quantile_0025, predict_param$quantiles[1])
		quantile_0975 = c(quantile_0975, predict_param$quantiles[2])

		prior_min = c(prior_min, min(param))
		prior_median = c(prior_median, median(param))
		prior_max = c(prior_max, max(param))
		prior_variance = c(prior_variance, var(param))
		prior_quantile_0025 = c(prior_quantile_0025, quantile(param, 0.025))
		prior_quantile_0975 = c(prior_quantile_0975, quantile(param, 0.975))
		
		plot(model$model.rf$predictions, param, main=colnames(parameters_IM)[i], xlab='estimated value', ylab='real value', cex.axis=1.2, cex.lab=1.2, pch=16, cex=1.2, col=rgb(0,0,0,0.25), xlim=range(param)); abline(v=predict_param$med, col='red', lwd=2)
	}
	dev.off()
	
	expectation = round(expectation,5)
	median = round(median,5)
	variance = round(variance, 5)
	variance_cdf = round(variance_cdf, 5)
	quantile_0025 = round(quantile_0025, 5)
	quantile_0975 = round(quantile_0975, 5)
	prior_min = round(prior_min, 5)
	prior_median = round(prior_median, 5)
	prior_max = round(prior_max, 5)
	prior_variance = round(prior_variance, 5)
	prior_quantile_0025 = round(prior_quantile_0025, 5)
	prior_quantile_0975 = round(prior_quantile_0975, 5)

	res_IM = data.frame(param_id, expectation, median, prior_median, variance, variance_cdf, prior_variance, quantile_0025, prior_quantile_0025, quantile_0975, prior_quantile_0975, prior_min, prior_max)

	write.table(res_IM, outfile, append=T, col.names=T, row.names=F, sep='\t', quote=F)
}

if( test_SI==1 ){
	write('\n#####\n\nParameters of the SI_2N model', outfile, append=T)
	write(paste('To convert in demographic units (in generations), times have to be multiplied by ', 4*Nref, sep=''), outfile, append=T)
	write(paste('To convert in demographic units (in number of individuals), sizes (N) have to be multiplied by ', Nref, '\n', sep=''), outfile, append=T)
	
	# SI model
	param_id = NULL
	expectation = NULL
	median = NULL
	variance = NULL
	variance_cdf = NULL
	quantile_0025 = NULL
	quantile_0975 = NULL

	prior_min = NULL
	prior_median = NULL
	prior_max = NULL
	prior_variance = NULL
	prior_quantile_0025 = NULL
	prior_quantile_0975 = NULL

	pdf(paste('ABC_', nameA, '_', nameB, '/parameters_SI_2N_', nameA, '_', nameB, '.pdf', sep=''), bg='white')
	for(i in 1:ncol(parameters_SI)){
		parameter_id = colnames(parameters_SI)[i]
		param = parameters_SI[,i]
		data = data.frame(param, ss_sim[['SI_2N']][, -ss_2_remove])
		model = regAbcrf(param~., data, ntree = ntree, ncores=ncores)
		predict_param = predict(model, ss_obs[, -ss_2_remove], data)
		
		param_id = c(param_id, parameter_id) 
		expectation = c(expectation, predict_param$expectation)
		median = c(median, predict_param$med)
		variance = c(variance, predict_param$variance)
		variance_cdf = c(variance_cdf, predict_param$variance.cdf)
		quantile_0025 = c(quantile_0025, predict_param$quantiles[1])
		quantile_0975 = c(quantile_0975, predict_param$quantiles[2])

		prior_min = c(prior_min, min(param))
		prior_median = c(prior_median, median(param))
		prior_max = c(prior_max, max(param))
		prior_variance = c(prior_variance, var(param))
		prior_quantile_0025 = c(prior_quantile_0025, quantile(param, 0.025))
		prior_quantile_0975 = c(prior_quantile_0975, quantile(param, 0.975))
		
		plot(model$model.rf$predictions, param, main=colnames(parameters_SI)[i], xlab='estimated value', ylab='real value', cex.axis=1.2, cex.lab=1.2, pch=16, cex=1.2, col=rgb(0,0,0,0.25), xlim=range(param)); abline(v=predict_param$med, col='red', lwd=2)

	}
	dev.off()
	
	expectation = round(expectation,5)
	median = round(median,5)
	variance = round(variance, 5)
	variance_cdf = round(variance_cdf, 5)
	quantile_0025 = round(quantile_0025, 5)
	quantile_0975 = round(quantile_0975, 5)
	prior_min = round(prior_min, 5)
	prior_median = round(prior_median, 5)
	prior_max = round(prior_max, 5)
	prior_variance = round(prior_variance, 5)
	prior_quantile_0025 = round(prior_quantile_0025, 5)
	prior_quantile_0975 = round(prior_quantile_0975, 5)

	res_SI = data.frame(param_id, expectation, median, prior_median, variance, variance_cdf, prior_variance, quantile_0025, prior_quantile_0025, quantile_0975, prior_quantile_0975, prior_min, prior_max)

	write.table(res_SI, outfile, append=T, col.names=T, row.names=F, sep='\t', quote=F)
}

###################################
# locus specific model comparison #
###################################
posterior = function(expected, variance, q1, q2, nsims){
        x = rnorm(nsims, expected, sqrt(variance))
        x = x+abs(min(x))
        x = x/max(x)
        x = x/((max(x)-min(x))/(q2-q1))
        x = x+q1
        return(x)
}

i = 0
if( predicted_model_iso_mig$allocation == 'migration' ){ # if the best model contains ongoing migration
	# get the bpfile
	bpfile = read.table(paste('ABC_', nameA, '_', nameB, '/bpfile', sep=''), skip=1, h=F)
	L = ceiling(median(as.numeric(bpfile[1,])))
	nsamA = ceiling(median(as.numeric(bpfile[2,])))
	nsamB = ceiling(median(as.numeric(bpfile[3,])))
	theta = median(as.numeric(bpfile[4,]))
	rho = median(as.numeric(bpfile[5,]))

	# prior distribution from the estimated values
	## model of migration
	i=i+1; N1 = posterior(res_IM$expectation[i], res_IM$variance[i], res_IM$quantile_0025[i], res_IM$quantile_0975[i], nsims_monolocus)
	i=i+1; N2 = posterior(res_IM$expectation[i], res_IM$variance[i], res_IM$quantile_0025[i], res_IM$quantile_0975[i], nsims_monolocus)
	i=i+1; Na = posterior(res_IM$expectation[i], res_IM$variance[i], res_IM$quantile_0025[i], res_IM$quantile_0975[i], nsims_monolocus)
	i=i+1; shape_N_a = posterior(res_IM$expectation[i], res_IM$variance[i], res_IM$quantile_0025[i], res_IM$quantile_0975[i], nsims_monolocus)
	i=i+1; shape_N_b = posterior(res_IM$expectation[i], res_IM$variance[i], res_IM$quantile_0025[i], res_IM$quantile_0975[i], nsims_monolocus)
	i=i+1; prior_Tsplit = posterior(res_IM$expectation[i], res_IM$variance[i], res_IM$quantile_0025[i], res_IM$quantile_0975[i], nsims_monolocus)
	i=i+1; M12 = posterior(res_IM$expectation[i], res_IM$variance[i], res_IM$quantile_0025[i], res_IM$quantile_0975[i], nsims_monolocus)
	i=i+1; shape_M12_a = posterior(res_IM$expectation[i], res_IM$variance[i], res_IM$quantile_0025[i], res_IM$quantile_0975[i], nsims_monolocus)
	i=i+1; shape_M12_b = posterior(res_IM$expectation[i], res_IM$variance[i], res_IM$quantile_0025[i], res_IM$quantile_0975[i], nsims_monolocus)
	i=i+1; M21 = posterior(res_IM$expectation[i], res_IM$variance[i], res_IM$quantile_0025[i], res_IM$quantile_0975[i], nsims_monolocus)
	i=i+1; shape_M21_a = posterior(res_IM$expectation[i], res_IM$variance[i], res_IM$quantile_0025[i], res_IM$quantile_0975[i], nsims_monolocus)
	i=i+1; shape_M21_b = posterior(res_IM$expectation[i], res_IM$variance[i], res_IM$quantile_0025[i], res_IM$quantile_0975[i], nsims_monolocus)
	
	prior_N1 = NULL
	prior_N2 = NULL
	prior_Na = NULL
	prior_M12 = NULL
	prior_M21 = NULL
	for(i in 1:nsims_monolocus){
		scalar_N = rbeta(1, shape_N_a[i], shape_N_b[i])
		prior_N1 = c(prior_N1, N1[i]*scalar_N)
		prior_N2 = c(prior_N2, N2[i]*scalar_N)
		prior_Na = c(prior_Na, Na[i]*scalar_N)
		
		scalar_M12 = rbeta(1, shape_M12_a[i], shape_M12_b[i])
		prior_M12 = c(prior_M12, M12[i]*scalar_M12)
		
		scalar_M21 = rbeta(1, shape_M21_a[i], shape_M21_b[i])
		prior_M21 = c(prior_M21, M21[i]*scalar_M21)
	}
	# mscommand = "msnsam 2*nsam nsims_monolocus -t tbs -r tbs tbs -I 2 nsam nsam 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs"
	prior_IM = cbind(rep(nsamA+nsamB, nsims_monolocus), rep(theta, nsims_monolocus), rep(rho, nsims_monolocus), rep(L, nsims_monolocus), rep(nsamA, nsims_monolocus), rep(nsamB, nsims_monolocus), prior_N1, prior_N2, prior_M12, prior_M21, prior_Tsplit, prior_Tsplit, prior_Na)
	setwd(paste('ABC_', nameA, '_', nameB, sep=''))
	system('mkdir migration')
	setwd('migration')
	write.table(prior_IM, 'prior_IM', col.names=F, row.names=F, sep='\t', quote=F)
	write('# monolocus migration', 'bpfile')
	write(L,'bpfile', append=T)
	write(nsamA, 'bpfile', append=T)
	write(nsamB, 'bpfile', append=T)
	write(theta, 'bpfile', append=T)
	write(rho, 'bpfile', append=T)
	
	commande = paste('cat prior_IM | msnsam tbs ', nsims_monolocus,' -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs | mscalc_2pop.py')
	system(commande)
	setwd("../")

	## model of isolation	
	prior_N1 = NULL
	prior_N2 = NULL
	prior_Na = NULL
	prior_M12 = rep(0, nsims_monolocus)
	prior_M21 = rep(0, nsims_monolocus)
	for(i in 1:nsims_monolocus){
		scalar_N = rbeta(1, shape_N_a[i], shape_N_b[i])
		prior_N1 = c(prior_N1, N1[i]*scalar_N)
		prior_N2 = c(prior_N2, N2[i]*scalar_N)
		prior_Na = c(prior_Na, Na[i]*scalar_N)
	}
	# mscommand = "msnsam 2*nsam nsims_monolocus -t tbs -r tbs tbs -I 2 nsam nsam 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs"
	prior_SI = cbind(rep(nsamA+nsamB, nsims_monolocus), rep(theta, nsims_monolocus), rep(rho, nsims_monolocus), rep(L, nsims_monolocus), rep(nsamA, nsims_monolocus), rep(nsamB, nsims_monolocus), prior_N1, prior_N2, prior_M12, prior_M21, prior_Tsplit, prior_Tsplit, prior_Na)
	system('mkdir isolation')
	setwd('isolation')
	write.table(prior_SI, 'prior_SI', col.names=F, row.names=F, sep='\t', quote=F)
	write('# monolocus isolation', 'bpfile')
	write(L,'bpfile', append=T)
	write(nsamA, 'bpfile', append=T)
	write(nsamB, 'bpfile', append=T)
	write(theta, 'bpfile', append=T)
	write(rho, 'bpfile', append=T)
	
	commande = paste('cat prior_SI | msnsam tbs ', nsims_monolocus,' -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs | mscalc_2pop.py')
	system(commande)
	setwd("../../")

	# locus specific model comparison
	iso = read.table(paste('ABC_', nameA, '_', nameB, '/isolation/ABCstat.txt', sep=''), h=T)
	iso = iso[, -grep('min', colnames(iso))]
	iso = iso[, -grep('max', colnames(iso))]
	mig = read.table(paste('ABC_', nameA, '_', nameB, '/migration/ABCstat.txt', sep=''), h=T)
	mig = mig[, -grep('min', colnames(mig))]
	mig = mig[, -grep('max', colnames(mig))]
	modIndexes = c(rep('isolation', nrow(iso)), rep('migration', nrow(mig)))
	all_models_sim = rbind(iso, mig)
		
	ss_2_remove = c(1)
	for(i in 2:ncol(iso)){
		if( sd(iso[, i])==0 ){
			ss_2_remove = c(ss_2_remove, i)
		}
	}
	for(i in 2:ncol(mig)){
		if( sd(mig[, i])==0 ){
			ss_2_remove = c(ss_2_remove, i)
		}
	}
	ss_2_remove = unique(ss_2_remove)
	ss_obs = read.table(paste('ABC_', nameA, '_', nameB, '/ABCstat_loci.txt', sep=''), h=T)
	obs_saved = ss_obs
	ss_obs = ss_obs[, -grep('min', colnames(ss_obs))]
	ss_obs = ss_obs[, -grep('max', colnames(ss_obs))]
	mod = abcrf(modIndexes~., data = data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
	predicted_model = predict(mod, data.frame(ss_obs), training=data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
	modele = predicted_model$allocation
	posterior_prob = predicted_model$post.prob
	
	res = cbind(obs_saved, modele, posterior_prob)
	write.table(res, paste('ABC_', nameA, '_', nameB, '/results_loci_migration.txt', sep=''), col.names=T, row.names=F, quote=F, sep='\t')
}

# quit the script
rm(list=ls())

