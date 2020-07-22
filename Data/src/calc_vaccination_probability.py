def p_ever_vaccinated(p_first_vac):
	p_ever_vac = []
	for y in range(len(p_first_vac)):
		p_ever_vac.append(sum(p_first_vac[y]))
		
	return p_ever_vac

def p_first_vaccinated(vac_cov, max_age, p_beginY, p_endY, raw_beginY, serum_time):
	# calculate probability of ever vacciated for each birth year 

	p = []
	u = []

	#For each birth year y (y is now index)
	for y in range(p_beginY-raw_beginY, p_endY-raw_beginY+1):
	
		#probability that birth year y get vaccinated for the first time in each year j
		py = [0 for j in range(0, max_age+1)]
		#probability that birth year y has never been vaccinated in each year j (after season j)
		uinf = [1 for j in range(0, max_age+1)]

		#probability of vaccination for birth year < 2002
		if raw_beginY + y < 2002:
			p.append(py)
			u.append(uinf)
			continue
	
		
		#probability of getting vaccination in year j = 0

		py[0] = 0		
		uinf[0] = 1			

		#probability of getting vaccination in year j = 1
		if raw_beginY+y+1 > serum_time:
			py[1] = 0
			uinf[1] = uinf[j-1]
		else:
			byear = raw_beginY + y
			season = byear + 1
			age = 0
	
			
			for i in range(len(vac_cov)):
				if vac_cov[i][0] == season and vac_cov[i][1] == age:
					vrate = vac_cov[i][2]
					break
				
			pyj_given_naive = vrate * 0.8
			pyj = uinf[0] * pyj_given_naive
			py[1] = pyj

			uinf[1] = uinf[0]*(1.0 - pyj_given_naive)
			
			
		#probability of getting infection for year j >=2
		# prob = probability of never been vaccinated * vaccination rate this season for this birth year
		for j in range(2, max_age+1):

			if raw_beginY+y+j > serum_time:
				py[j] = 0
				uinf[j] = uinf[j-1]
				
			else:
				byear = raw_beginY + y
				season = byear + j
				age1 = season - byear - 1
				age2 = season - byear
				
				for i in range(len(vac_cov)):
					if vac_cov[i][0] == season and vac_cov[i][1] == age1:
						vrate1 = vac_cov[i][2]
						break
				for i in range(len(vac_cov)):
					if vac_cov[i][0] == season and vac_cov[i][1] == age2:
						vrate2 = vac_cov[i][2]
						break
					
					
				pyj_given_naive = vrate1*0.46 + vrate2*0.54
					
				pyj = uinf[j-1] * pyj_given_naive
				py[j] = pyj

				uinf[j] = uinf[j-1]*(1.0 - pyj_given_naive)
				

		p.append(py)
		u.append(uinf)
		
		
	return p, u
	

def p_first_vaccinated_naive_fraction(vac_cov, naive_for_vac, max_age, p_beginY, p_endY, raw_beginY, serum_time):
	# calculate probability of ever vacciated for each birth year 

	p = []
	u = []

	#For each birth year y (y is now index)
	for y in range(p_beginY-raw_beginY, p_endY-raw_beginY+1):
	
		#probability that birth year y get vaccinated for the first time in each year j
		py = [0 for j in range(0, max_age+1)]
		#probability that birth year y has never been vaccinated in each year j (after season j)
		uinf = [1 for j in range(0, max_age+1)]

		#probability of vaccination for birth year < 2002
		if raw_beginY + y < 2002:
			p.append(py)
			u.append(uinf)
			continue
	
		
		#probability of getting vaccination in year j = 0

		py[0] = 0		
		uinf[0] = 1			

		#probability of getting vaccination in year j = 1
		if raw_beginY+y+1 > serum_time:
			py[1] = 0
			uinf[1] = uinf[j-1]
		else:
			byear = raw_beginY + y
			season = byear + 1
			age = 0
	
			
			for i in range(len(vac_cov)):
				if vac_cov[i][0] == season and vac_cov[i][1] == age:
					vrate = vac_cov[i][2]
					break
				
			pyj_given_naive = vrate * 0.8 #fraction of >6 month olds
			pyj = naive_for_vac[1] * pyj_given_naive
			py[1] = uinf[0]*pyj

			uinf[1] = uinf[0]*(1.0 - pyj)
			
			
		#probability of getting infection for year j >=2
		# prob = probability of never been vaccinated * vaccination rate this season for this birth year
		for j in range(2, max_age+1):

			if raw_beginY+y+j > serum_time:
				py[j] = 0
				uinf[j] = uinf[j-1]
				
			else:
				byear = raw_beginY + y
				season = byear + j
				age1 = season - byear - 1
				age2 = season - byear
				
				for i in range(len(vac_cov)):
					if vac_cov[i][0] == season and vac_cov[i][1] == age1:
						vrate1 = vac_cov[i][2]
						break
				for i in range(len(vac_cov)):
					if vac_cov[i][0] == season and vac_cov[i][1] == age2:
						vrate2 = vac_cov[i][2]
						break
					
					
				pyj_given_naive = vrate1*0.46 + vrate2*0.54
					
				pyj = naive_for_vac[j] * pyj_given_naive
				py[j] = uinf[j-1]*pyj

				uinf[j] = uinf[j-1]*(1.0 - pyj)
				

		p.append(py)
		u.append(uinf)
		
		
	return p, u


	
	
	
	
	
	
	