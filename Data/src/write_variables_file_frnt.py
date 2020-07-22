
import read_files
import calc_similarities
import calc_first_exposure
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
#mean attack rate
mean_atr = 0.328
#number of sequences sampled per year
nSamp = 100

#subjects by birth year or by age (by_byear == 0 means by age)
by_byear = 0

##################################################
###             Read files                     ###
##################################################

#input file
seqf1_name = "../data/aligned_6812_AA.fas"
seqf2_name = "../data/aligned_1217_AA.fas"
epif_name = "../data/shih_epitope.csv"
frequencyf_name = "../data/Frequency_from_fluview.csv"
intensityf_name = "../data/Intensity_from_ILI.csv"
vac_cov_for_1617f_name = "../data/vaccov_subj_1617_pen.csv"


if by_byear == 0:
	scalef_name = "../processed_data/intensity_scalings.csv"
elif by_byear == 1:
	scalef_name = "/processed_data/intensity_scalings_by_byear.csv" 
	
#read epitope files
epitopes = read_files.read_epitope_shih(epif_name)

#sequence samples
nSamp = 100
seqByYear = read_files.sequences_byYear(1968, 2011, seqf1_name)
seqByYear += read_files.sequences_byYear(2012, 2017, seqf2_name)
seqByYear = read_files.sample_nSamp(seqByYear, nSamp)

seqBySeason = read_files.sequences_bySeason(1992, 2012, seqf1_name)
seqBySeason += read_files.sequences_bySeason(2013, 2017, seqf2_name)
seqBySeason = read_files.sample_nSamp(seqBySeason, nSamp)


remove_naives = 0

test_vs = ["FRNT_HA_viruses"]
serum_time = [2017]

#####################################################
###  calculate exposure/imprinting probabilities  ###
#####################################################


##################################################
# read intensity, fraction, scale file

intensity_A = calc_first_exposure.get_intensity_A(intensityf_name)


frac_H3, frac_H2, frac_H1 = \
						calc_first_exposure.get_subtype_fractions(
										frequencyf_name)

scale = read_files.read_scales(scalef_name)

##################################################
# Probability of imprinting to H1, H2, H3

# calculate raw probability of first exposure to influenza A at each year
raw_p_first_exp, p_uninf, raw_p_first_exp_h3, raw_p_first_exp_h2, raw_p_first_exp_h1 = \
			calc_first_exposure.calc_raw_p_first_exp(
							mean_atr, intensity_A, scale, 
							frac_H3, frac_H2, frac_H1,
							max_age, 
							p_beginY=1918, p_endY = 2016, 
							raw_beginY = 1918, serum_time = serum_time[0])


# calculate probability of imprinting to H3, H2, and H1
							
p_impH3, p_impH2, p_impH1 = calc_first_exposure.calc_p_impH3N2(
							raw_p_first_exp_h3, raw_p_first_exp_h2, raw_p_first_exp_h1,
							p_beginY=1918, p_endY=2016, raw_beginY=1918)
		
####################################################
#calculate probability of first exposure to H3, without considering other subtypes
raw_p_first_exp_onlyH3, p_uninf_H3 = \
calc_first_exposure.calc_p_first_exp_H3( 
							mean_atr, intensity_A, scale, 
							frac_H3, max_age, 
							p_beginY=beginY, p_endY=endY, raw_beginY=1918, 
							serum_time = serum_time[0])

#probability of ever being exposed to h3n2
#p_ever_infH3N2 = calc_first_exposure.calc_p_expH3N2(raw_p_first_exp_onlyH3)

p_first_expH3N2 = calc_first_exposure.calc_p_first_expH3N2(
					raw_p_first_exp_onlyH3)


#################################################
#vaccine coverage
vc_dict = read_files.read_vac_cov_for1617(vac_cov_for_1617f_name)


for vdx in range(len(test_vs)):

	viruses, vyears = read_files.read_virus(test_vs[vdx])
	seg_sites, wtAA, WT = utils.get_seg_sites(viruses)
	
	#output file
	if remove_naives == 1:
		dataf_name = "../covariate/imprinting_covariates_remove_naives_" + test_vs[vdx]  +".csv"
	elif remove_naives == 0:
		dataf_name = "../covariate/imprinting_covariates_" + test_vs[vdx] + ".csv"
	
	##################################
	#variables to write

	dataf = open(dataf_name, "w")
	dataf.write("Virus,y,vyear,")
	dataf.write("episim,")	
	dataf.write("vac_cov,")
	dataf.write("imp_h1n1,imp_h2n2,imp_h3n2\n")

	t = serum_time[vdx]

	for v in range(len(viruses)):
		virus = viruses[v].split("\n")[1]
		vname = viruses[v].split("\n")[0]
		vyear = vyears[v]
		
		sv = calc_similarities.calc_sv(virus, seqByYear, seqBySeason, [j for i in epitopes for j in i])
		
		for y in range(beginY, endY+1) :
			age = str(t-y)
					
			#similarity at all epitope using all residues		
			
			similarity_all_epitopes = calc_similarities.calc_Svy_season(y, sv, 
										p_first_expH3N2, max_age, beginY, endY)
													
			vac_cov_for1617 = vc_dict[2017-y]
			
			#subtype first exposure probability
			
			p_imp_h1n1 = p_impH1[y-1918]
			p_imp_h2n2 = p_impH2[y-1918]
			p_imp_h3n2 = p_impH3[y-1918]
			
			
			#write in output file

			oneline = vname + "," + str(y) + "," + vyear + ","			
			oneline += str(similarity_all_epitopes) + ","			
			oneline += str(vac_cov_for1617) + ","		
			oneline += str(p_imp_h1n1) + "," + str(p_imp_h2n2) + "," + str(p_imp_h3n2) 
			
			#people born 1966 or before: 
			
			if y == beginY:
				for py in range(beginP, endP+1) :
					oneline_prev = oneline.split(",")
					oneline_prev[1] = str(py) #birth year
					oneline_prev[4] = str(vc_dict[2017-py])
					
					oneline_prev[5] = str(p_impH1[py-1918]) #probability of imprinting to H1
					oneline_prev[6] = str(p_impH2[py-1918]) #probability of imprinting to H2
					oneline_prev[7] = str(p_impH3[py-1918]) #probability of imprinting to H3
					
					oneline_prev = ",".join(oneline_prev)
					dataf.write(oneline_prev + "\n")

					
			#from beginY to endY
			
			dataf.write(oneline + "\n")

	dataf.close()
