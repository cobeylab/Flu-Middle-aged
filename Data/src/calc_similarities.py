import csv
import numpy
import scipy.stats as stats
#import isGly
import read_files

def calc_p_exposed(y, testv, sequences, epiLoc, p, max_age, beginY, endY):
	exp = 0
	y = y-beginY #y is now used for index

	for j in range(max_age+1):
		if beginY+y+j > endY+1:
			break
		if beginY < 1968 and y == 0 and j == 0:
			continue

		pyj = p[y][j]

		exp += pyj

	return str(exp)
	
def calc_sv(testv, seqByYear, seqBySeason, epiLoc, is_NA = 0):

	if is_NA == 0:
		first_season = 1969
	elif is_NA == 1:
		first_season = 1958
	
	sv = []
	
	for season in range(first_season, 2017+1):
		
		if season < 1992:
			t = season - (first_season - 1)
			
			#year t-1
			nj = len(seqByYear[t-1])
			svj_1 = 0
			for i in range(nj):
				svi = 0
				for s in epiLoc:
					if testv[s-1] == seqByYear[t-1][i][s-1]:
						svi += 1
				svi = 1.0*svi/(len(epiLoc))
				svj_1 += svi
			
			svj_1 = 1.0*svj_1/nj
			
			#year t
			nj = len(seqByYear[t])
			svj_2 = 0
			for i in range(nj):
				svi = 0
				for s in epiLoc:
					if testv[s-1] == seqByYear[t][i][s-1]:
						svi += 1
				svi = 1.0*svi/(len(epiLoc))
				svj_2 += svi
			
			svj_2 = 1.0*svj_2/nj
			
			svj = (1.0/3.5) * svj_1 + (2.5/3.5) * svj_2
			
		elif season >= 1992:
			t = season - 1992 
			nj = len(seqBySeason[t])
			svj = 0
			for i in range(nj):
				svi = 0
				for s in epiLoc:
					if testv[s-1] == seqBySeason[t][i][s-1]:
						svi += 1
				svi = 1.0*svi/(len(epiLoc))
				svj += svi
			
			svj = 1.0*svj/nj	
		
		#print (svj)
		sv.append(svj)
		
	return sv
	
def calc_vv(testv, vacSeqBySeason, epiLoc):

	sv = []
	
	for season in range(1969, 2017+1):
		
		if season < 2003:
			svj = 0
			
		elif season >= 2003:
			t = season - 2003 
			svj = 0

			for s in epiLoc:
				if testv[s-1] == vacSeqBySeason[t][s-1]:
					svj += 1
			svj = svj/(len(epiLoc))
		sv.append(svj)
		
	return sv
	
def calc_Svy_season(y, sv, p_first_expH3N2, max_age, beginY, endY, is_NA = 0):

	Svy = 0
	y = y-beginY #y is now used for index
	
	sequence_endY = 2017
	
	if is_NA == 0:
		sequence_beginY = 1969
		
	elif is_NA == 1:
		sequence_beginY = 1958

		
	for j in range(max_age+1):

		if beginY + y + j > sequence_endY:
			break
		if beginY + y + j < sequence_beginY: 
			continue

		pyj = p_first_expH3N2[y][j] 
			
		Svy += pyj*sv[beginY - sequence_beginY + y+j]

	return str(Svy)

	
def calc_Svy(y, testv, sequences, epiLoc, p, max_age, beginY, endY):

	Svy = 0
	y = y-beginY #y is now used for index
	sequence_beginY = 1968
	sequence_endY = 2017
	for j in range(max_age+1):

		if beginY + y + j > sequence_endY:
			break
		if beginY + y + j < sequence_beginY: 
			continue

		pyj = p[y][j]
		nj = len(sequences[beginY - sequence_beginY + y+j])

		svj = 0
		for i in range(nj):
			svi = 0
			for s in epiLoc:
				if testv[s-1] == sequences[beginY - sequence_beginY +y+j][i][s-1]:
					svi += 1
			svi = 1.0*svi/(len(epiLoc))
			svj += svi
		
		svj = 1.0*svj/nj
		Svy += pyj*svj

	return str(Svy)
	
def calc_Svy_NA(y, testv, sequences, epiLoc, p, max_age, beginY, endY):

	Svy = 0
	y = y-beginY #y is now used for index
	sequence_beginY = 1957
	sequence_endY = 2017
	for j in range(max_age+1):

		if beginY + y + j > sequence_endY:
			break
		if beginY + y + j < sequence_beginY: 
			continue

		pyj = p[y][j]
		nj = len(sequences[beginY - sequence_beginY + y+j])

		svj = 0
		for i in range(nj):
			svi = 0
			for s in epiLoc:
				if testv[s-1] == sequences[beginY - sequence_beginY +y+j][i][s-1]:
					svi += 1
			svi = 1.0*svi/(len(epiLoc))
			svj += svi
		
		svj = 1.0*svj/nj
		Svy += pyj*svj

	return str(Svy)

	
def calc_Ivy(y, allele, sequences, epiLoc, p, max_age, beginY, endY):
	#calculate probability that is imprinted to wild type
	Svy = 0
	y = y-beginY #y is now used for index

	sequence_beginY = 1968
	sequence_endY = 2017
	for j in range(max_age+1):

		if beginY + y + j > sequence_endY:
			break
		if beginY + y + j < sequence_beginY: 
			continue


		pyj = p[y][j]
		nj = len(sequences[beginY - sequence_beginY +y+j])
		svj = 0
		for i in range(nj):
			svi = 0
			for s in epiLoc:
				if allele == sequences[beginY - sequence_beginY +y+j][i][s-1]:
					svi += 1
			svi = 1.0*svi/(len(epiLoc))
			svj += svi

		svj = 1.0*svj/nj
		Svy += pyj*svj

	return str(Svy)

	
def calc_dvySvy(y, testv, sequences, epiLoc, pngsLoc, p, max_age, delta_opt, beginY, endY):
	dvySvy = 0
	y = y-beginY #y is now used for index
	pngs_v = isGly.find_pngs_in_seq(testv, pngsLoc)

	sequence_beginY = 1968
	sequence_endY = 2017
	for j in range(max_age+1):

		if beginY + y + j > sequence_endY:
			break
		if beginY + y + j < sequence_beginY: 
			continue


		pyj = p[y][j]
		nj = len(sequences[beginY - sequence_beginY +y+j])
		dvjsvj = 0
		for i in range(nj):
			svi = 0
			for s in epiLoc:
				if testv[s-1] == sequences[beginY - sequence_beginY +y+j][i][s-1]:
					svi += 1
			svi = 1.0*svi/(len(epiLoc))
			pngs_i = isGly.find_pngs_in_seq(sequences[beginY - sequence_beginY +y+j][i], pngsLoc)
			dvi = isGly.determin_delta(pngs_i, pngs_v, pngsLoc, delta_opt)
			dvjsvj += dvi*svi

		dvjsvj = 1.0*dvjsvj/nj
		dvySvy += pyj*dvjsvj

	return str(dvySvy)

def calc_GEvy(y, testv, sequences, pngsLoc, p, max_age, delta_opt, beginY, endY):
	Gvy = 0
	y = y-beginY #y is now used for index
	pngs_v = isGly.find_pngs_in_seq(testv, pngsLoc)

	sequence_beginY = 1968
	sequence_endY = 2017
	for j in range(max_age+1):

		if beginY + y + j > sequence_endY:
			break
		if beginY + y + j < sequence_beginY: 
			continue


		pyj = p[y][j]
		nj = len(sequences[beginY - sequence_beginY +y+j])

		gvj = 0
		for i in range(nj):
			gvi = 0
			pngs_i = isGly.find_pngs_in_seq(sequences[beginY - sequence_beginY +y+j][i], pngsLoc)
			gvi = isGly.determin_delta(pngs_i, pngs_v, pngsLoc, delta_opt)
			gvj += gvi

		gvj = 1.0*gvj/(nj)
		Gvy += pyj*gvj

	return str(Gvy)

'''
def calc_Gmvy(y, testv, sequences, glyLoc, p, max_age, beginY, endY):
	# What is it doing?
	print ("calc_Gmvy called")
	Gvy = 0
	for j in range(y, y+max_age+1):
		if beginY+y+j > endY+1:
			break
		if beginY < 1968 and y == 0 and j == 0:
			continue

		pyj = p[y][j]
		nj = len(sequences[j])

		gvj = 0
		for i in range(nj):
			gvi = 0
			for s in glyLoc:
				if isGly(sequences[j][i][s]) == 0 and isGly(testv[s]) == 1:
					gvj -= 1
			gvi = 1.0*gvi/len(glyLoc)
			gvj += gvi

		gvj = 1.0*gvj/(nj)
		Gvy += pyj*gvj

	return str(Gvy)
'''