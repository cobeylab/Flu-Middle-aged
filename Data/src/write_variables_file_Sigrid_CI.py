
import read_files
import calc_similarities
import calc_first_exposure
import resample_for_ci
import utils

##################################################
###             Set parameters                 ###
##################################################

#Calculate similarity for birth year from beginY to endY
beginY = 1967
endY = 2016

#people born before 1967 have same similarity with 1966
beginP = 1918
endP = 1966

#oldest age of getting first infection
max_age = 20
mean_atr = 0.328
#number of sequences sampled per year
nSamp = 100

#subjects by birth year or by age
by_byear = 0

####################################################################################
###             Read files                     
####################################################################################


frequencyf_name = "../data/Frequency_from_fluview.csv"
total_Af_name = "../processed_data/consolidated_seasonal_incidence_pandemic.csv"
intensityf_name = "../data/Intensity_from_ILI.csv"

if by_byear == 0:
	scalef_name = "../processed_data/intensity_scalings.csv"
elif by_byear == 1:
	scalef_name = "../processed_data/intensity_scalings_by_byear.csv" 

epif_name = "../data/shih_epitope.csv"
epitopes = read_files.read_epitope_shih(epif_name)

seqf1_name = "../data/aligned_6812_AA.fas"
seqf2_name = "../data/aligned_1217_AA.fas"

#sequences
nSamp = 100
seqByYear = read_files.sequences_byYear(1968, 2011, seqf1_name)
seqByYear += read_files.sequences_byYear(2012, 2017, seqf2_name)
seqByYear = read_files.sample_nSamp(seqByYear, nSamp)

seqBySeason = read_files.sequences_bySeason(1992, 2012, seqf1_name)
seqBySeason += read_files.sequences_bySeason(2013, 2017, seqf2_name)
seqBySeason = read_files.sample_nSamp(seqBySeason, nSamp)


####################################################################################
# read intensity, fraction, scale, total influenza A, and mean ILI

intensity_A = calc_first_exposure.get_intensity_A(intensityf_name)


frac_H3, frac_H2, frac_H1 = \
						calc_first_exposure.get_subtype_fractions(
										frequencyf_name)
									
scale = read_files.read_scales(scalef_name)


total_specimen, total_A, fraction_A = resample_for_ci.read_totalA_file(total_Af_name)


mean_ILI = calc_first_exposure.get_mean_ILI(intensityf_name)


####################################################################################
###  calculate exposure/imprinting probabilities 
####################################################################################

####################################################################################
# calculate imprinting probability


raw_p_first_exp, p_uninf, raw_p_first_exp_h3, raw_p_first_exp_h2, raw_p_first_exp_h1 = \
				calc_first_exposure.calc_raw_p_first_exp(
								mean_atr, intensity_A, scale, 
								frac_H3, frac_H2, frac_H1,
								max_age, 
								p_beginY=1918, p_endY = 2016, 
								raw_beginY = 1918, serum_time = 2017)

p_impH3, p_impH2, p_impH1 = calc_first_exposure.calc_p_impH3N2(
							raw_p_first_exp_h3, raw_p_first_exp_h2, raw_p_first_exp_h1,
							p_beginY=1918, p_endY=2016, raw_beginY=1918)

p_impH3_rmvNaives, p_impH2_rmvNaives, p_impH1_rmvNaives = calc_first_exposure.calc_p_impH3N2_rmvNaive(
							raw_p_first_exp_h3, raw_p_first_exp_h2, raw_p_first_exp_h1,
							p_beginY=1918, p_endY=2016, raw_beginY=1918)

# calculate similarity given exposure to H3

raw_p_first_exp_onlyH3, p_uninf_H3 = \
calc_first_exposure.calc_p_first_exp_H3( 
							mean_atr, intensity_A, scale, 
							frac_H3, max_age, 
							p_beginY=beginY, p_endY=endY, raw_beginY=1918, 
							serum_time = 2017)
								
#For similarity, remove naives			
#p_first_expH3N2 = calc_first_exposure.calc_p_first_expH3N2(
#					raw_p_first_exp_onlyH3)
					
p_first_expH3N2 = calc_first_exposure.calc_p_first_expH3N2_rmvnaive(
					raw_p_first_exp_onlyH3)

######################
					
test_vs = ["FRNT_HA_viruses"]
viruses, vyears = read_files.read_virus(test_vs[0])



def get_similarity(sv, p_first_expH3N2): 
	
	similarity = []
	for y in range(beginY, endY+1) :
		
		similarity_y1 =\
		calc_similarities.calc_Svy_season(y, sv, p_first_expH3N2, max_age, beginY, endY)
		
		similarity.append(similarity_y1)
			
	for y in range(beginP, endP+1):
		similarity.insert(0, similarity[0])
		
	return similarity

	
# virus = 1: 3c2.A
# virus = -1: 3c2.A2
	
sv1 = calc_similarities.calc_sv(viruses[1].split("\n")[1], 
								seqByYear, seqBySeason, [j for i in epitopes for j in i])
sv2 = calc_similarities.calc_sv(viruses[-1].split("\n")[1], 
								seqByYear, seqBySeason, [j for i in epitopes for j in i])
								
similarity_tv1 = get_similarity(sv1, p_first_expH3N2)
similarity_tv2 = get_similarity(sv2, p_first_expH3N2)


####################################################################################
# calculate CI for probability of imprinting to H3, H2, and H1

num_r = 1000

mean_atr_resample = resample_for_ci.get_meanAtr_resample(
					mean_atr_gostic=0.28, se_atr_gostic=0.01, num_r = num_r)
					
total_A_resample, frac_A_resample = resample_for_ci.get_totalA_resample(
									total_specimen, fraction_A, total_A, num_r)
								
frac_H3_resample, frac_H2_resample, frac_H1_resample = \
						resample_for_ci.get_subtype_fractions_resample(
						total_A_resample, frac_H3, frac_H2, frac_H1, num_r = num_r)
						
intensity_A_resample = resample_for_ci.get_intensity_resample(mean_ILI, frac_A_resample, intensity_A, num_r)


p_impH3_r = []	
p_impH3_rmvNaives_r = []	
similarity_tv1_r = []
similarity_tv2_r = []

for i in range(0, num_r):		
	
	# imprinting probability resample
	raw_p_first_exp, p_uninf, raw_p_first_exp_h3, raw_p_first_exp_h2, raw_p_first_exp_h1 = \
				calc_first_exposure.calc_raw_p_first_exp(
								mean_atr_resample[i], intensity_A_resample[i], scale, 
								frac_H3_resample[i], frac_H2_resample[i], frac_H1_resample[i],
								max_age, 
								p_beginY=1918, p_endY = 2016, 
								raw_beginY = 1918, serum_time = 2017)

	p_impH3_r1, p_impH2_r1, p_impH1_r1 = calc_first_exposure.calc_p_impH3N2(
							raw_p_first_exp_h3, raw_p_first_exp_h2, raw_p_first_exp_h1,
							p_beginY=1918, p_endY=2016, raw_beginY=1918)
							
	p_impH3_rmvNaive_r1, p_impH2_rmvNaive_r1, p_impH1_rmvNaive_r1 = calc_first_exposure.calc_p_impH3N2_rmvNaive(
							raw_p_first_exp_h3, raw_p_first_exp_h2, raw_p_first_exp_h1,
							p_beginY=1918, p_endY=2016, raw_beginY=1918)

			
	p_impH3_r.append(p_impH3_r1)
	p_impH3_rmvNaives_r.append(p_impH3_rmvNaive_r1)
	

	# similarity resample
	raw_p_first_exp_onlyH3_r1, p_uninf_H3_r1 = \
	calc_first_exposure.calc_p_first_exp_H3( 
								mean_atr_resample[i], intensity_A_resample[i], scale, 
								frac_H3_resample[i], max_age, 
								p_beginY=beginY, p_endY=endY, raw_beginY=1918, 
								serum_time = 2017)

	#For similarity, remove naives

	p_first_expH3N2_r1 = calc_first_exposure.calc_p_first_expH3N2_rmvnaive(
						raw_p_first_exp_onlyH3_r1)

					
	similarity_tv1_r1 = get_similarity(sv1, p_first_expH3N2_r1)
	similarity_tv2_r1 = get_similarity(sv2, p_first_expH3N2_r1)

	similarity_tv1_r.append(similarity_tv1_r1)
	similarity_tv2_r.append(similarity_tv2_r1)
	
	
#get CI
	
p_impH3_CI1, p_impH3_CI2 = resample_for_ci.get_p_imp_CI(p_impH3_r, num_r)
p_impH3_rmvNaives_CI1, p_impH3_rmvNaives_CI2 = resample_for_ci.get_p_imp_CI(p_impH3_rmvNaives_r, num_r)

similarity_tv1_CI1, similarity_tv1_CI2 = resample_for_ci.get_p_imp_CI(similarity_tv1_r, num_r)
similarity_tv2_CI1, similarity_tv2_CI2 = resample_for_ci.get_p_imp_CI(similarity_tv2_r, num_r)

#############################################################################
# save each resample


def write_resamples(resamplef_name, p_imp):

	resamplef = open(resamplef_name, "w")
	
	resamples = list(range(1, num_r+1))
	resamples = ["resample_" + str(s) for s in resamples]
	head= ( ",".join( resamples ) )
	resamplef.write("yob," + head + ",\n")
	
	for y in range(len(p_imp[0])):
		oneline = str(y+1918) + ","
		for i in range(0, num_r):
			oneline += str(p_imp[i][y]) + ","
		resamplef.write(oneline + "\n")
	
	resamplef.close()
	
####################################################


dataf_name = "../covariate/subtype_imprinting_probabilities_CI.csv"
dataf = open(dataf_name, "w")
dataf.write("yob,\
			imp_H3N2,imp_H3N2_ci1,imp_H3N2_ci2,\
			imp_H3N2_rmvNaives,imp_H3N2_rmvNaives_ci1,imp_H3N2_rmvNaives_ci2,\
			similarity_tv1,similarity_tv1_ci1,similarity_tv1_ci2,\
			similarity_tv2,similarity_tv2_ci1,similarity_tv2_ci2\n")

for y in range(1918-1918, 2016-1918+1):
				
			oneline = str(y+1918) + "," 
			
			oneline += str(p_impH3[y]) + ","
			oneline += str(p_impH3_CI1[y]) + "," + str(p_impH3_CI2[y]) + ","

			oneline += str(p_impH3_rmvNaives[y]) + ","
			oneline += str(p_impH3_rmvNaives_CI1[y]) + "," + str(p_impH3_rmvNaives_CI2[y]) + ","
			
			oneline += str(similarity_tv1[y]) +","
			oneline += str(similarity_tv1_CI1[y]) + "," + str(similarity_tv1_CI2[y]) + ","

			oneline += str(similarity_tv2[y]) +","
			oneline += str(similarity_tv2_CI1[y]) + "," + str(similarity_tv2_CI2[y]) 
			
			dataf.write(oneline + "\n")

dataf.close()




