import csv
import math
import numpy
import scipy.stats as stats



def calc_one_subtype_intensity(fracf_name):

	#Calculate relative H3N2 intensity from sequence fraction

	#read subtype fraction file
	with open(fracf_name, "rU") as fracf:
		fraction_list = list(csv.reader(fracf))

	col_year = fraction_list[0].index("year")
	col_frac = fraction_list[0].index("frac_h3n2")

	fraction_list = numpy.array(fraction_list)

	year = fraction_list[:,col_year][1:]
	year = [int(y) for y in year]

	fraction = fraction_list[:,col_frac][1:]
	fraction = [float(f) for f in fraction]

	#detrend subtype specific sequence collection trend
	#by subtracting slope*y from fraction

	params = stats.linregress(year, fraction)
	slope = params[0]

	fraction_detrend = []
	for y in range(len(year)):
		corrected = fraction[y] - slope*y
		fraction_detrend.append(corrected)

	#normalize intensity by dividing detrended fraction by mean fraction
	mean_frac = sum(fraction_detrend)/len(fraction_detrend)

	intensity = []
	for y in range(len(year)):
		intensity_y = 1.0*fraction_detrend[y]/mean_frac
		intensity.append(intensity_y)

	return intensity
	
def get_intensity_A(intensityf_name):
	#with open(intensityf_name, "rU") as int_frac:
	#	intensities = list(csv.reader(int_frac))	
		
	intensityf = open(intensityf_name, "U")
	intensities = []
	
	for line in intensityf:
		if line.find("season") >= 0:
			continue
		if line.find("2009.5") >= 0:
			continue
		intensity = line.split("\n")[0].split(",")
		intensities.append( float(intensity[6]) )

	
	intensities.append(5.0)
	return intensities
		
def get_mean_ILI(intensityf_name):
	#with open(intensityf_name, "rU") as int_frac:
	#	intensities = list(csv.reader(int_frac))	
		
	intensityf = open(intensityf_name, "U")
	intensities = []
	
	for line in intensityf:
		if line.find("season") >= 0:
			continue
		if line.find("2009.5") >= 0:
			continue
			
		intensity = line.split("\n")[0].split(",")
		if float(intensity[0]) < 1977:
			intensities.append(0)
		else:
			intensities.append( float(intensity[4]) )

	return intensities
	
def get_subtype_fractions(frequencyf_name):
	with open(frequencyf_name, "rU") as int_frac:
		int_frac = list(csv.reader(int_frac))
		
	col_year = int_frac[0].index("year")
	col_frac_h3 = int_frac[0].index("H3N2_fraction")
	col_frac_h2 = int_frac[0].index("H2N2_fraction")
	col_frac_h1 = int_frac[0].index("H1N1_fraction")
	
	int_frac = numpy.array(int_frac)
	
	int_frac = int_frac[int_frac[:, 0].argsort()]
	
	pandemic_rowNum = (numpy.where(int_frac[:, col_year] == "2009.5")[0][0])
	pandemic_row = (int_frac[ pandemic_rowNum ] )
	
	int_frac = numpy.delete(int_frac, pandemic_rowNum, axis=0) #remove 2009.5
	int_frac = numpy.delete(int_frac, -1, axis=0) #remove characters
	int_frac = numpy.delete(int_frac, -1, axis=0) #remove 2017-2018 season 
	
	int_frac = numpy.append(int_frac, [pandemic_row], axis=0)
	
	year = int_frac[:,col_year] 
	year = [float(y) for y in year]	
	
	frac_H3 = int_frac[:,col_frac_h3]
	frac_H3 = [float(f) for f in frac_H3]	
	
	frac_H2 = int_frac[:,col_frac_h2]
	frac_H2 = [float(f) for f in frac_H2]	
	
	frac_H1 = int_frac[:,col_frac_h1]
	frac_H1 = [float(f) for f in frac_H1]		
		

	return (frac_H3, frac_H2, frac_H1)
	

def calc_pyj_nextyear(intensity, mean_atr, max_age, beginY, endY, serum_time):
	p = []
	u = []

	for y in range(beginY-beginY, endY+1-beginY):
		#probability that a cohort gets infeced for the first time this year
		py = [0 for j in range(0, max_age+1)]
		#probability that a cohort has never been infected
		uinf = [0 for j in range(0, max_age+1)]

		#probability of getting infection in the first year after birth (j=1)
		# prob = mean attack rate * intensity
		if beginY+y+1 > serum_time:
			py[1] == 0
		else:
			py1 = mean_atr*intensity[y+1]
			py[1] = py1
			uinf[1] = (1.0-py1)


		#probability of getting infection for year j > 1
		# prob = probability of never been infected * infected this year
		for j in range(2, max_age+1):

			if beginY+y+j > serum_time:
				py[j] = 0
			else:
				pyj = uinf[j-1]*mean_atr*intensity[y+j]
				py[j] = pyj

				uinf[j] = uinf[j-1]*(1.0-mean_atr*intensity[y+j])
		
		if y+max_age <= serum_time-beginY:
			py = [pi/sum(py) for pi in py]
		else:
			py = py

		p.append(py)
		u.append(uinf)

	return p

def calc_raw_p_first_exp(mean_atr, intensity_A, scale, frac_H3, frac_H2, frac_H1,
									 max_age, p_beginY, p_endY, raw_beginY, serum_time):
	# For intensity_A, scale, frac_H3, frac_H2, frac_H1, index -1 is for pandemic.
	p = [] #raw_p_first_exp
	u = [] #p_uninf
	p_h3 = [] #raw_p_first_exp_h3
	p_h2 = [] #raw_p_first_exp_h2
	p_h1 = [] #raw_p_first_exp_h1

	# For this function, you need mean attack rate, intensity A(normalized), 
	# fraction H3, and intensity scaling.
	
	#For each birth year y (y is now index)
	for y in range(p_beginY-raw_beginY, p_endY-raw_beginY+1):
	
		#probability that birth year y get infeced for the first time in each year j
		py = [0 for j in range(0, max_age+1)]
		#probability that birth year y has never been infected in each year j (after season j)
		uinf = [0 for j in range(0, max_age+1)]
		#probability that birth year y get infected to subtype for the first time in each year j
		py_h3 = [0 for j in range(0, max_age+1)]
		py_h2 = [0 for j in range(0, max_age+1)]
		py_h1 = [0 for j in range(0, max_age+1)]

		#probability of getting infection in year j = 0
		# py[0] = mean attack rate * intensity_A * frac_H3 * scale
		if raw_beginY+y > serum_time:
			py[0] = 0
			
			py_h3[0] = 0
			py_h2[0] = 0
			py_h1[0] = 0
			
			uinf[0] = 1.0
			
		else:

			attack_rate = mean_atr * intensity_A[y] * scale[y][0]
			p_infection = 1.0 - math.exp(-attack_rate)		
					
			py[0] = p_infection
			uinf[0] = (1.0-p_infection)
			
			py_h3[0] = py[0] * frac_H3[y]
			py_h2[0] = py[0] * frac_H2[y]
			py_h1[0] = py[0] * frac_H1[y]

		#probability of getting infection for year j > 0
		# prob = probability of never been infected * infected this year
		for j in range(1, max_age+1):

			if raw_beginY+y+j > serum_time:
				py[j] = 0
				
				py_h3[j] = 0
				py_h2[j] = 0
				py_h1[j] = 0
				
				uinf[j] = uinf[j-1]
				
			else:
				if(y+j+raw_beginY) == 2010:
				# If we are dealing with 2010 season ( = 2009-10 season),
				# aggregate pandemic + 2009-10 season
					if (raw_beginY + y) == 2009:
						scale_pandemic = 0.65
					elif (raw_beginY + y) < 2009:
						scale_pandemic = 1.0
					
					
					#attack_rate = mean_atr * intensity_A[-1] * scale_pandemic +\
					#				mean_atr * intensity_A[y+j] * scale[y][j] 
					attack_rate = mean_atr * intensity_A[-1] * scale_pandemic
					p_infection = 1.0 - math.exp(-attack_rate)
					py_pandemic = uinf[j-1]*p_infection
					uinf_pandemic = uinf[j-1]*(1.0 - p_infection)
					
					attack_rate = mean_atr * intensity_A[y+j] * scale[y][j]
					p_infection = 1.0 - math.exp(-attack_rate)
					py[j] = py_pandemic + uinf_pandemic*p_infection
					uinf[j] = uinf_pandemic*(1.0 - p_infection)
					
					#We can just use subtype fraction of 2009-10 season
					#as it is very similar to pandemic
				else:
					attack_rate = mean_atr * intensity_A[y+j]  * scale[y][j]
					
					p_infection = 1.0 - math.exp(-attack_rate)

					py[j] = uinf[j-1]*p_infection
					uinf[j] = uinf[j-1]*(1.0 - p_infection)
					#if raw_beginY + y == 2016:
					#	print(mean_atr, intensity_A[y+j], scale[y][j], py[j])
				
				py_h3[j] = py[j] * frac_H3[y+j]
				py_h2[j] = py[j] * frac_H2[y+j]
				py_h1[j] = py[j] * frac_H1[y+j]

		p.append(py)
		u.append(uinf)
		
		p_h3.append(py_h3)
		p_h2.append(py_h2)
		p_h1.append(py_h1)
		
	return p, u, p_h3, p_h2, p_h1

def calc_p_impH3N2(raw_p_first_exp_h3, raw_p_first_exp_h2, raw_p_first_exp_h1,
			p_beginY, p_endY, raw_beginY):
	
	p_impH3 = []
	p_impH2 = []
	p_impH1 = []
	
	for y in range(p_beginY - raw_beginY, p_endY - raw_beginY + 1):
		
		p_impH3_for_y = sum(raw_p_first_exp_h3[y])
		p_impH2_for_y = sum(raw_p_first_exp_h2[y])
		p_impH1_for_y = sum(raw_p_first_exp_h1[y])
			
		norm = p_impH3_for_y + p_impH2_for_y +p_impH1_for_y
		
		if(y < p_endY - raw_beginY - 20 + 1):
			
			p_impH3_for_y = p_impH3_for_y/norm
			p_impH2_for_y = p_impH2_for_y/norm
			p_impH1_for_y = p_impH1_for_y/norm
		
		p_impH3.append(p_impH3_for_y)
		p_impH2.append(p_impH2_for_y)
		p_impH1.append(p_impH1_for_y)
				
	return p_impH3, p_impH2, p_impH1	
	
def calc_p_impH3N2_rmvNaive(raw_p_first_exp_h3, raw_p_first_exp_h2, raw_p_first_exp_h1,
			p_beginY, p_endY, raw_beginY):
	
	p_impH3 = []
	p_impH2 = []
	p_impH1 = []
	
	for y in range(p_beginY - raw_beginY, p_endY - raw_beginY + 1):
		
		p_impH3_for_y = sum(raw_p_first_exp_h3[y])
		p_impH2_for_y = sum(raw_p_first_exp_h2[y])
		p_impH1_for_y = sum(raw_p_first_exp_h1[y])
		norm = p_impH3_for_y + p_impH2_for_y +p_impH1_for_y
		
		p_impH3_for_y = p_impH3_for_y/norm
		p_impH2_for_y = p_impH2_for_y/norm
		p_impH1_for_y = p_impH1_for_y/norm
		
		p_impH3.append(p_impH3_for_y)
		p_impH2.append(p_impH2_for_y)
		p_impH1.append(p_impH1_for_y)
				
	return p_impH3, p_impH2, p_impH1
	
def calc_p_impG1(p_impH1_1, p_impH2_1):

	p_impG1_1 = []
	for y in range(len(p_impH1_1)):
		p_impG1_1.append(p_impH1_1[y] + p_impH2_1[y])
	
	return (p_impG1_1)	

	
def calc_p_impH1_then_everH2(raw_p_first_exp_h1, p_ever_infH2N2,
			p_beginY, p_endY, raw_beginY, serum_time):
	
	p_impH1_then_everH2 = []

	for y in range(p_beginY - raw_beginY, p_endY - raw_beginY + 1):
		
		py_impH1_then_everH2 = 0
		for j in range(len(raw_p_first_exp_h1[0])):
			
			if raw_beginY+y+j >= serum_time:
				break
				
			py_impH1_then_everH2 += raw_p_first_exp_h1[y][j] * p_ever_infH2N2[y+j]
		
		p_impH1_then_everH2.append(py_impH1_then_everH2)
				
	return p_impH1_then_everH2	
	
def calc_p_first_exp_H3(mean_atr, intensity_A, scale, frac_H3, 
						max_age, p_beginY, p_endY, raw_beginY, serum_time):
	# For intensity_A, scale, frac_H3, frac_H2, frac_H1, index -1 is for pandemic.
	p = [] #raw_p_first_exp
	u = [] #p_uninf

	# For this function, you need mean attack rate, intensity A(normalized), 
	# fraction H3, and intensity scaling.
	
	#For each birth year y (y is now index)
	for y in range(p_beginY - raw_beginY, p_endY - raw_beginY + 1):
	
		#probability that birth year y has never been infected in each year j (after season j)
		uinf = [0 for j in range(0, max_age+1)]
		#probability that birth year y get infected to subtype for the first time in each year j
		py = [0 for j in range(0, max_age+1)]

		#probability of getting infection in year j = 0
		# py[0] = mean attack rate * intensity_A * frac_H3 * scale
		if raw_beginY+y > serum_time:

			py[0] = 0
			uinf[0] = 1
			
		else:

			attack_rate = mean_atr * intensity_A[y] * scale[y][0] * frac_H3[y]
			py0 = 1.0 - math.exp(-attack_rate)	
					
			py[0] = py0
			uinf[0] = (1.0-py0)

		#probability of getting infection for year j > 0
		# prob = probability of never been infected * infected this year
		for j in range(1, max_age+1):

			if raw_beginY+y+j > serum_time:
				py[j] = 0

				uinf[j] = uinf[j-1]
				
			else:

				attack_rate = mean_atr * intensity_A[y+j] * scale[y][j] * frac_H3[y+j]
				
				pyj_given_naive = 1.0 - math.exp(-attack_rate)					
				pyj = uinf[j-1] * pyj_given_naive
				py[j] = pyj

				uinf[j] = uinf[j-1]*(1.0 - pyj_given_naive)
				
		p.append(py)
		u.append(uinf)

		
	return p, u

	
def calc_p_first_exp_H3_old_a(mean_atr, intensity_A, scale, frac_H3, 
						max_age, p_beginY, p_endY, raw_beginY, serum_time):
	# For intensity_A, scale, frac_H3, frac_H2, frac_H1, index -1 is for pandemic.
	p = [] #raw_p_first_exp
	u = [] #p_uninf

	# For this function, you need mean attack rate, intensity A(normalized), 
	# fraction H3, and intensity scaling.
	
	#For each birth year y (y is now index)
	for y in range(p_beginY - raw_beginY, p_endY - raw_beginY + 1):
	
		#probability that birth year y has never been infected in each year j (after season j)
		uinf = [0 for j in range(0, max_age+1)]
		#probability that birth year y get infected to subtype for the first time in each year j
		py = [0 for j in range(0, max_age+1)]

		#probability of getting infection in year j = 0
		# py[0] = mean attack rate * intensity_A * frac_H3 * scale
		if raw_beginY+y > serum_time:

			py[0] = 0
			uinf[0] = 1
			
		else:

			py0 = mean_atr * intensity_A[y] * scale[y][0] * frac_H3[y]
					
			py[0] = py0
			uinf[0] = (1.0-py0)

		#probability of getting infection for year j > 0
		# prob = probability of never been infected * infected this year
		for j in range(1, max_age+1):

			if raw_beginY+y+j > serum_time:
				py[j] = 0

				uinf[j] = uinf[j-1]
				
			else:

				pyj_given_naive = mean_atr * intensity_A[y+j] * scale[y][j] * frac_H3[y+j]
					
				pyj = uinf[j-1] * pyj_given_naive
				py[j] = pyj

				uinf[j] = uinf[j-1]*(1.0 - pyj_given_naive)
				
		p.append(py)
		u.append(uinf)

		
	return p, u

def calc_p_first_expH3N2(raw_p_first_exp_onlyH3, is_NA = 0):

	if is_NA == 0:
		y_before = 1966
	elif is_NA == 1:
		y_before = 1955
		
	
	p_first_expH3N2_rmvnaive = []
	
	for y in range(len(raw_p_first_exp_onlyH3)):
		norm = sum(raw_p_first_exp_onlyH3[y])
		p_first_exp_onlyH3_y = []
		
		for j in range(len(raw_p_first_exp_onlyH3[y])):
			if(y < 2017-y_before-20 and norm != 0):
				p_first_exp_onlyH3_y.append(raw_p_first_exp_onlyH3[y][j]/norm)
			else:
				p_first_exp_onlyH3_y.append(raw_p_first_exp_onlyH3[y][j])
			
		p_first_expH3N2_rmvnaive.append(p_first_exp_onlyH3_y)

	return p_first_expH3N2_rmvnaive
	
def calc_p_first_expH3N2_rmvnaive(raw_p_first_exp_onlyH3):
	p_first_expH3N2_rmvnaive = []

	for y in range(len(raw_p_first_exp_onlyH3)):
		norm = sum(raw_p_first_exp_onlyH3[y])
		p_first_exp_onlyH3_y = []
		
		for j in range(len(raw_p_first_exp_onlyH3[y])):
			if (norm != 0):
				p_first_exp_onlyH3_y.append(raw_p_first_exp_onlyH3[y][j]/norm)
			elif (norm == 0):
				p_first_exp_onlyH3_y.append(raw_p_first_exp_onlyH3[y][j])
			
		p_first_expH3N2_rmvnaive.append(p_first_exp_onlyH3_y)

	return p_first_expH3N2_rmvnaive
	
def calc_p_expH3N2(raw_p_first_exp_onlyH3):
	#The begin year is 1966
	
	p_expH3N2 = []
	
	for y in range(len(raw_p_first_exp_onlyH3)):
				
		p_expH3N2_for_y = sum(raw_p_first_exp_onlyH3[y])
		if (y < 2017 - 1966 - 20):
			#print (y, 2017 - 1966 - 20)
			p_expH3N2_for_y = 1.0
		p_expH3N2.append(p_expH3N2_for_y)
		
	return p_expH3N2
	
def calc_p_ever_expH3N2_v(p_ever_infH3N2, p_ever_vaccinated):
	p_ever_expH3N2 = []
	
	for y in range(len(p_ever_infH3N2)):
		p_ever_expH3N2_y = p_ever_infH3N2[y] + p_ever_vaccinated[y] - p_ever_infH3N2[y]*p_ever_vaccinated[y]
		p_ever_expH3N2.append(p_ever_expH3N2_y)
		
	return p_ever_expH3N2
		
		
def calc_pyj(intensity, mean_atr, max_age, beginY, endY):
	p = []
	u = []

	for y in range(beginY-beginY, endY+1-beginY):
		#probability that a cohort gets infeced for the first time this year
		py = [0 for j in range(0, max_age+1)]
		#probability that a cohort has never been infected
		uinf = [0 for j in range(0, max_age+1)]

		#probability of getting infection in the first year
		# prob = mean attack rate * intensity
		py0 = mean_atr*intensity[y]
		py[0] = py0
		uinf[0] = (1.0-py0)


		#probability of getting infection for year j > y
		# prob = probability of never been infected * infected this year
		for j in range(1, max_age+1):

			if beginY+y+j > endY:
				py[j] = 0
			else:
				pyj = uinf[j-1]*mean_atr*intensity[y+j]
				py[j] = pyj

				uinf[j] = uinf[j-1]*(1.0-mean_atr*intensity[y+j])

		#do not normalize for now

		p.append(py)
		u.append(uinf)

	return p

def calc_pyj_remove_naive(intensity, mean_atr, max_age, beginY, endY):
	p = []
	u = []

	for y in range(beginY-beginY, endY+1-beginY):
		#probability that a cohort gets infeced for the first time this year
		py = [0 for j in range(0, max_age+1)]
		#probability that a cohort has never been infected
		uinf = [0 for j in range(0, max_age+1)]

		#probability of getting infection in the first year
		# prob = mean attack rate * intensity
		py0 = mean_atr*intensity[y]
		py[0] = py0
		uinf[0] = (1.0-py0)


		#probability of getting infection for year j > y
		# prob = probability of never been infected * infected this year
		for j in range(1, max_age+1):

			if beginY+y+j > endY:
				py[j] = 0
			else:
				pyj = uinf[j-1]*mean_atr*intensity[y+j]
				py[j] = pyj

				uinf[j] = uinf[j-1]*(1.0-mean_atr*intensity[y+j])

		#normalize probability of getting first infection
		#by dividing by sum of probabilities for the cohort y

		py = [p/sum(py) for p in py]
		p.append(py)
		u.append(uinf)

	return p

def get_simple_p(max_age, beginY, endY):
	p = []
	for y in range(beginY, endY+1):
		py = [1.0/(max_age+1) for j in range(max_age+1)]
		p.append(py)
	return p
