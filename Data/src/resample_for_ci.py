import csv
import math
import numpy
import scipy.stats as stats

def get_meanAtr_resample(mean_atr_gostic, se_atr_gostic, num_r):

	atr_resample = []
	for i in range(0, num_r):
		
		atr = numpy.random.normal(mean_atr_gostic, se_atr_gostic)
		atr_inst = -1.0 * numpy.log(1 - atr)
		atr_resample.append(atr_inst)
	
	return atr_resample

def read_totalA_file(totalAf_name):
	inf = open(totalAf_name, "r")
	
	total_specimen = [0 for i in range(1918,1977)]
	totalA = [0 for i in range(1918,1977)]
	fractionA = [0 for i in range(1918,1977)]
	
	for line in inf:
		if line.find("SEASON") >= 0:
			each = line.split(",")
			col_season = each.index("SEASON END YEAR")
			col_totSpec = each.index("TOTAL SPECIMENS")
			col_totA = each.index("TOTAL A")
			col_fracA = each.index("FRAC A")
			
		else:
			each = line.split(",")
			if each[col_season] == "2009.5":
				continue
			total_specimen.append( int(each[col_totSpec]) )
			totalA.append( int(each[col_totA]) )
			fractionA.append( float(each[col_fracA]) )
			
	return total_specimen, totalA, fractionA

def get_subtype_fractions_resample(totalA, frac_H3, frac_H2, frac_H1, num_r):

	frac_H3_resample = []
	frac_H2_resample = []
	frac_H1_resample = []
	
	for r in range(0, num_r):
		fracs_H3_r1 = []
		fracs_H2_r1 = []
		fracs_H1_r1 = []
		
		for y in range(1918-1918, 2018-1918+1):
			
			if y<1977-1918:
				fracs_H3_r1.append(frac_H3[y])
				fracs_H2_r1.append(frac_H2[y])
				fracs_H1_r1.append(frac_H1[y])			
			else:
				fracs_y1 = numpy.random.multinomial(totalA[r][y], [frac_H3[y], frac_H2[y], frac_H1[y]] )
				fh3 = fracs_y1[0]/sum(fracs_y1)
				fh2 = fracs_y1[1]/sum(fracs_y1)
				fh1 = fracs_y1[2]/sum(fracs_y1)
				fracs_H3_r1.append(fh3)
				fracs_H2_r1.append(fh2)
				fracs_H1_r1.append(fh1)
		
		frac_H3_resample.append(fracs_H3_r1)
		frac_H2_resample.append(fracs_H2_r1)
		frac_H1_resample.append(fracs_H1_r1)
		
	return frac_H3_resample, frac_H2_resample, frac_H1_resample

def get_totalA_resample(total_specimen, fraction_A, total_A, num_r):

	total_A_resample = []
	frac_A_resample = []
	
	for r in range(0, num_r):
		total_A_r1 = []
		frac_A_r1 = []
		
		for y in range(1918-1918, 2018-1918+1):
			
			if y<1977-1918:
				total_A_r1.append(total_A[y])
				frac_A_r1.append(fraction_A[y])
				
			else:
				total_A1 = numpy.random.binomial(total_specimen[y], fraction_A[y], 1 )
				total_A_r1.append(total_A1[0])
				frac_A_r1.append(total_A1[0]/total_specimen[y])
		
		total_A_resample.append(total_A_r1)
		frac_A_resample.append(frac_A_r1)

		
	return total_A_resample, frac_A_resample
	
def get_intensity_resample(mean_ILI, frac_A_resample, intensity_A, num_r):
	intensity_A_resample = []

	for r in range(0, num_r):

		intensity_A_r1 = []
		
		for y in range(1918-1918, 2018-1918):
			
			if y<1977-1918:
				intensity_A_r1.append(intensity_A[y])
			
			else:
				intensity_A1 = mean_ILI[y] * frac_A_resample[r][y]
				intensity_A_r1.append(intensity_A1)	
				
		#add pandemic 2009.5
		intensity_A_r1.append(0.00837)
		
		#normalize
		avg_intensity_A_r1 = sum(intensity_A_r1[1977-1918:])/len(intensity_A_r1[1977-1918:])

		for y in range(1977-1918, 2018-1918+1):
			intensity_A_r1[y] = intensity_A_r1[y]/avg_intensity_A_r1
		
		intensity_A_resample.append(intensity_A_r1)

	return (intensity_A_resample)
	
def get_p_imp_CI(frac_subtype, num_r):	
	p_imp_ci1 = []
	p_imp_ci2 = []
	for y in range(1918-1918, 2017-1918):
		p_imp_y = []
		for i in frac_subtype:
			p_imp_y.append(i[y])
		p_imp_y.sort()
		p_imp_ci1.append(p_imp_y[ int(num_r*0.025)-1 ])
		p_imp_ci2.append(p_imp_y[ int(num_r*0.975)-1 ])

	
	return p_imp_ci1, p_imp_ci2
	
	

	
	
	
	
	
	
	
	
	
	