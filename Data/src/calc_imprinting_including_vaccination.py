import math

def calc_raw_p_first_exp_infection_plus_vaccination(mean_atr, intensity_A, scale, frac_H3, frac_H2, frac_H1,
									 max_age, vac_cov, naive_for_vac,
									 p_beginY, p_endY, raw_beginY, serum_time, is_full):

	# For intensity_A, scale, frac_H3, frac_H2, frac_H1, index -1 is for pandemic.
	p = [] #raw_p_first_exp
	u = [] #p_uninf
	p_h3 = [] #raw_p_first_exp_h3
	p_h2 = [] #raw_p_first_exp_h2
	p_h1 = [] #raw_p_first_exp_h1
	p_vac = []

	# For this function, you need mean attack rate, intensity A(normalized), 
	# fraction H3, and intensity scaling.
	
	#For each birth year y (y is now index)
	for y in range(p_beginY-raw_beginY, p_endY-raw_beginY+1):
	
		#probability that birth year y get vaccinated for the first time in each year j
		py_vac = [0 for j in range(0, max_age+1)]
		#probability that birth year y get exposed for the first time in each year j
		py = [0 for j in range(0, max_age+1)]
		#probability that birth year y has never been exposed in each year j (after season j)
		uinf = [0 for j in range(0, max_age+1)]
		nvac = [0 for j in range(0, max_age+1)]
		#probability that birth year y get infected to subtype for the first time in each year j
		py_h3 = [0 for j in range(0, max_age+1)]
		py_h2 = [0 for j in range(0, max_age+1)]
		py_h1 = [0 for j in range(0, max_age+1)]
		

		#probability of getting infection in year j = 0
		# py[0] = mean attack rate * intensity_A * frac_H3 * scale
		if raw_beginY+y > serum_time:
			py_vac[0] = 0
			py[0] = 0
			
			py_h3[0] = 0
			py_h2[0] = 0
			py_h1[0] = 0
						
			uinf[0] = 1.0
			nvac[0] = 1.0
			
		else:

			p_vaccination = vac_cov[y][0]
			py_vac[0] = 1.0*p_vaccination/1.0
			
			if (is_full==0):
				nvac[0] = 1.0*(1.0 - p_vaccination/1.0)
			elif (is_full==1):
				nvac[0] = 1.0
				
			p_naive_after_vac = 1.0* (1.0 - p_vaccination/1.0)
			
			attack_rate = mean_atr * intensity_A[y] * scale[y][0]
			p_infection = 1.0 - math.exp(-attack_rate)		
					
			py[0] = p_naive_after_vac*p_infection
			uinf[0] = p_naive_after_vac*(1.0 - p_infection)
			
			py_h3[0] = py[0] * frac_H3[y]
			py_h2[0] = py[0] * frac_H2[y]
			py_h1[0] = py[0] * frac_H1[y]

		#probability of getting infection for year j > 0
		# prob = probability of never been infected * infected this year
		for j in range(1, max_age+1):

			if raw_beginY+y+j > serum_time:
				py_vac[j] = 0
				py[j] = 0
				
				py_h3[j] = 0
				py_h2[j] = 0
				py_h1[j] = 0
				
				uinf[j] = uinf[j-1]
				nvac[j] = nvac[j-1]
				
			else:
				if(y+j+raw_beginY) == 2010:
				# If we are dealing with 2010 season ( = 2009-10 season),
				# aggregate pandemic + 2009-10 season
					if (raw_beginY + y) == 2009:
						scale_pandemic = 0.65
					elif (raw_beginY + y) < 2009:
						scale_pandemic = 1.0
					
					p_vaccination = naive_for_vac[j] * vac_cov[y][j]
					if j==1:
						p_vaccination = p_vaccination * 0.8 # only 56 percent are older than 6 months

					py_vac[j] = uinf[j-1] *p_vaccination
					
					if (is_full==0):
						nvac[j] = nvac[j-1]*(1.0 - p_vaccination/nvac[j-1])
					elif (is_full==1):
						nvac[j] = 1.0
						
					p_naive_after_vac = uinf[j-1] * (1.0 - p_vaccination/nvac[j-1])
 
					attack_rate = mean_atr * intensity_A[-1] * scale_pandemic
					p_infection = 1.0 - math.exp(-attack_rate)
					py_pandemic = p_naive_after_vac*p_infection
					uinf_pandemic = p_naive_after_vac*(1.0 - p_infection)
					
					attack_rate = mean_atr * intensity_A[y+j] * scale[y][j]
					p_infection = 1.0 - math.exp(-attack_rate)
					py[j] = py_pandemic + uinf_pandemic*p_infection
					uinf[j] = uinf_pandemic*(1.0 - p_infection)
					
					#We can just use subtype fraction of 2009-10 season
					#as it is very similar to pandemic
				else:
				
					p_vaccination = naive_for_vac[j] * vac_cov[y][j] 
					if j==1:
						p_vaccination = p_vaccination * 0.8 # only 56 percent are older than 6 months

					py_vac[j] = uinf[j-1] * p_vaccination/nvac[j-1]
					
					if (is_full==0):
						nvac[j] = nvac[j-1]*(1.0 - p_vaccination/nvac[j-1])
					elif (is_full==1):
						nvac[j] = 1.0
						

					p_naive_after_vac = uinf[j-1] * (1.0 - p_vaccination/nvac[j-1])
					
					
					attack_rate = mean_atr * intensity_A[y+j] * scale[y][j]					
					p_infection = 1.0 - math.exp(-attack_rate)

					py[j] = p_naive_after_vac*p_infection
					uinf[j] = p_naive_after_vac*(1.0 - p_infection)
						
				py_h3[j] = py[j] * frac_H3[y+j]
				py_h2[j] = py[j] * frac_H2[y+j]
				py_h1[j] = py[j] * frac_H1[y+j]
				

		# Save this birth year's first exposure probabilities
		p_vac.append(py_vac)
		p.append(py)
		u.append(uinf)
		
		p_h3.append(py_h3)
		p_h2.append(py_h2)
		p_h1.append(py_h1)
		
	return p_h3, p_h2, p_h1, p_vac
	

def calc_raw_p_first_expH3_infection_plus_vaccination(mean_atr, intensity_A, scale, frac_H3, 
									 max_age, vac_cov, naive_for_vac,
									 p_beginY, p_endY, raw_beginY, serum_time, is_full):

	# For intensity_A, scale, frac_H3, index -1 is for pandemic.
	p = [] #raw_p_first_exp
	u = [] #p_uninf

	p_vac = []

	# For this function, you need mean attack rate, intensity A(normalized), 
	# fraction H3, and intensity scaling.
	
	#For each birth year y (y is now index)
	for y in range(p_beginY-raw_beginY, p_endY-raw_beginY+1):
		
		#probability that birth year y get vaccinated for the first time in each year j
		py_vac = [0 for j in range(0, max_age+1)]
		#probability that birth year y get exposed to H3 for the first time in each year j
		py = [0 for j in range(0, max_age+1)]
		#probability that birth year y has never been exposed in each year j (after season j)
		uinf = [0 for j in range(0, max_age+1)]
		nvac = [0 for j in range(0, max_age+1)]

		#probability of getting infection in year j = 0
		# py[0] = mean attack rate * intensity_A * frac_H3 * scale
		if raw_beginY+y > serum_time:
			py_vac[0] = 0
			py[0] = 0
						
			uinf[0] = 1.0
			nvac[0] = 1.0
			
		else:

			p_vaccination = vac_cov[y][0] 
			py_vac[0] = 1.0*p_vaccination/1.0
			
			if(is_full==0):
				nvac[0] = 1.0*(1.0 - p_vaccination/1.0)
			elif (is_full == 1):
				nvac[0] = 1.0
				
			p_naive_after_vac = 1.0* (1.0 - p_vaccination/1.0)
			
			attack_rate = mean_atr * intensity_A[y] * scale[y][0] * frac_H3[y]
			p_infection = 1.0 - math.exp(-attack_rate)		
					
			py[0] = p_naive_after_vac*p_infection
			uinf[0] = p_naive_after_vac*(1.0 - p_infection)

				
		#probability of getting infection for year j > 0
		# prob = probability of never been infected * infected this year
		for j in range(1, max_age+1):

			if raw_beginY+y+j > serum_time:
				py_vac[j] = 0
				py[j] = 0
				
				uinf[j] = uinf[j-1]
				nvac[j] = nvac[j-1]
				
			else:
							
				p_vaccination = naive_for_vac[j] * vac_cov[y][j] 
				if j==1:
					p_vaccination = p_vaccination * 0.8 # only 56 percent are older than 6 months

				if (p_vaccination > nvac[j-1]) :
					p_vaccination = nvac[j-1] - 0.000001
					
				py_vac[j] = uinf[j-1] * p_vaccination/nvac[j-1]
				
				if(is_full==0):
					nvac[j] = nvac[j-1]*(1.0 - p_vaccination/nvac[j-1])
				elif (is_full==1):
					nvac[j] = 1.0
					
				p_naive_after_vac = uinf[j-1] * (1.0 - p_vaccination/nvac[j-1])
				
				attack_rate = mean_atr * intensity_A[y+j] * scale[y][j]	* frac_H3[y+j]			
				p_infection = 1.0 - math.exp(-attack_rate)

				py[j] = p_naive_after_vac*p_infection
				uinf[j] = p_naive_after_vac*(1.0 - p_infection)
			
			if y == p_beginY-raw_beginY+44 and is_full==0:
				print (p_vaccination, py_vac[j], nvac[j], p_naive_after_vac)
				print (p_infection, py[j], uinf[j])

		# Save this birth year's first exposure probabilities
		p_vac.append(py_vac)
		p.append(py)
		u.append(uinf)
		
	return p, p_vac
		
def calc_p_impH3N2_and_impVac(raw_p_first_exp_h3, raw_p_first_exp_h2, raw_p_first_exp_h1, raw_p_first_exp_vac,
			p_beginY, p_endY, raw_beginY):
	
	p_impH3 = []
	p_impH2 = []
	p_impH1 = []
	p_impVac = []
	
	for y in range(p_beginY - raw_beginY, p_endY - raw_beginY + 1):
		
		p_impH3_for_y = sum(raw_p_first_exp_h3[y])
		p_impH2_for_y = sum(raw_p_first_exp_h2[y])
		p_impH1_for_y = sum(raw_p_first_exp_h1[y])
		p_impVac_for_y = sum(raw_p_first_exp_vac[y])
		norm = p_impH3_for_y + p_impH2_for_y + p_impH1_for_y + p_impVac_for_y
		

		#if(y < 2017 - p_beginY - 20):
		#	print (p_beginY, 2017 - p_beginY - 20, y)
		#	p_impH3_for_y = p_impH3_for_y/norm
		#	p_impH2_for_y = p_impH2_for_y/norm
		#	p_impH1_for_y = p_impH1_for_y/norm
		#	p_impVac_for_y = p_impVac_for_y/norm
		
		p_impH3.append(p_impH3_for_y)
		p_impH2.append(p_impH2_for_y)
		p_impH1.append(p_impH1_for_y)
		p_impVac.append(p_impVac_for_y)
				
	return p_impH3, p_impH2, p_impH1, p_impVac	
	
	
def calc_p_expH3_or_expVac(raw_p_first_exp_h3_v, raw_p_first_vac_v):
	
	p_expH3_all = []
	p_expVac_all = []
	p_exp_all = []
	
	for y in range(len(raw_p_first_exp_h3_v)):
		
		p_expH3 = sum(raw_p_first_exp_h3_v[y])
		p_expVac = sum(raw_p_first_vac_v[y])
		p_exp = sum(raw_p_first_exp_h3_v[y]) + sum(raw_p_first_vac_v[y])
		
		if y < 2017-1966-20:
			norm = p_exp
		else:
			norm = 1
			
		p_expH3 = p_expH3/norm
		p_expVac = p_expVac/norm
		p_exp = p_exp/norm
				
		p_expH3_all.append(p_expH3)
		p_expVac_all.append(p_expVac)
		p_exp_all.append(p_exp)
	
	return p_exp_all, p_expH3_all, p_expVac_all



def calc_raw_p_first_vaccination(mean_atr, intensity_A, scale, frac_H3, 
									 max_age, vac_cov, naive_for_vac,
									 p_beginY, p_endY, raw_beginY, serum_time):

	# For intensity_A, scale, frac_H3, index -1 is for pandemic.
	p = [] #raw_p_first_exp
	u = [] #p_uninf

	p_vac = []

	# For this function, you need mean attack rate, intensity A(normalized), 
	# fraction H3, and intensity scaling.
	
	#For each birth year y (y is now index)
	for y in range(p_beginY-raw_beginY, p_endY-raw_beginY+1):
	
		#probability that birth year y get vaccinated for the first time in each year j
		py_vac = [0 for j in range(0, max_age+1)]

		#probability that birth year y has never been vaccinated in each year j (after season j)
		nvac = [0 for j in range(0, max_age+1)]


		if raw_beginY+y > serum_time:
			py_vac[0] = 0						
			nvac[0] = 1.0
			
		else:

			p_vaccination = vac_cov[y][0] 
			py_vac[0] = 1.0*p_vaccination/1.0
			
			nvac[0] = 1.0*(1.0 - p_vaccination/1.0)
	
		
		#probability of getting infection for year j > 0
		# prob = probability of never been infected * infected this year
		for j in range(1, max_age+1):

			if raw_beginY+y+j > serum_time:
				py_vac[j] = 0
				nvac[j] = nvac[j-1]
				
			else:
							
				p_vaccination = naive_for_vac[j] * vac_cov[y][j] 
				if j==1:
					p_vaccination = p_vaccination * 0.8 # only 56 percent are older than 6 months

				py_vac[j] = nvac[j-1]*p_vaccination
				nvac[j] = nvac[j-1]*(1.0 - p_vaccination)
				

			#if y == 90:
			#	print (p_vaccination, py_vac[j], nvac[j])
			
		# Save this birth year's first exposure probabilities
		p_vac.append(py_vac)

		
	return  p_vac














	